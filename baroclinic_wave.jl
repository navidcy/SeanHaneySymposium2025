using Oceananigans
using Oceananigans.Units
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using Printf

arch = GPU()
resolution = 1 #1//4
Nx = 360 ÷ resolution
Ny = 160 ÷ resolution
Nz = 10

grid = LatitudeLongitudeGrid(arch;
                             size = (Nx, Ny, Nz),
                             halo = (7, 7, 3),
                             latitude = (-80, 80),
                             longitude = (0, 360),
                             z = (-1000, 0))

momentum_advection = WENOVectorInvariant(order=3)
tracer_advection   = WENO(order=5)

buoyancy = SeawaterBuoyancy(equation_of_state=TEOS10EquationOfState())
model = HydrostaticFreeSurfaceModel(; grid,
                                      momentum_advection,
                                      tracer_advection,
                                      coriolis = HydrostaticSphericalCoriolis(),
                                      buoyancy,
                                      tracers=(:T, :S))

Tᵢ(λ, φ, z) = 30 * (1 - tanh((abs(φ) - 40) / 5)) / 2 + rand()
Sᵢ(λ, φ, z) = 28 - 5e-3 * z + rand()
set!(model, T=Tᵢ, S=Sᵢ)

simulation = Simulation(model, Δt=5minutes, stop_time=180days)

function progress(sim)
    T = sim.model.tracers.T
    S = sim.model.tracers.S
    u, v, w = sim.model.velocities

    msg = @sprintf("%d: %s, max|u|: (%.2e, %.2e, %.2e)",
                   iteration(sim), prettytime(sim),
                   maximum(abs, u),
                   maximum(abs, v),
                   maximum(abs, w))

    msg *= @sprintf(", extrema(T): (%.2f, %.2f), extrema(S): (%.2f, %.2f)",
                    minimum(T), maximum(T), minimum(S), maximum(S))

    @info msg

    return nothing
end

add_callback!(simulation, progress, IterationInterval(10))

u, v, w = model.velocities
T = model.tracers.T
S = model.tracers.S
ζ = ∂x(v) - ∂y(u)
fields = (; u, v, w, T, S, ζ)

ow = JLD2OutputWriter(model, fields,
                      filename = "baroclinic_wave.jld2",
                      indices = (:, :, grid.Nz),
                      schedule = TimeInterval(1days),
                      overwrite_existing = true)

simulation.output_writers[:surface] = ow

run!(simulation)

u = FieldTimeSeries("baroclinic_wave.jld2", "u"; backend = OnDisk())
v = FieldTimeSeries("baroclinic_wave.jld2", "v"; backend = OnDisk())
T = FieldTimeSeries("baroclinic_wave.jld2", "T"; backend = OnDisk())
ζ = FieldTimeSeries("baroclinic_wave.jld2", "ζ"; backend = OnDisk())

times = u.times
Nt = length(times)

n = Observable(1)

Tn = @lift begin
    Tn = interior(T[$n])
    view(Tn, :, :, 1)
end

ζn = @lift interior(ζ[$n], :, :, 1)
Tn = @lift (T[$n], :, :, 1)

un = Field{Face, Center, Nothing}(u.grid)
vn = Field{Center, Face, Nothing}(v.grid)
s = Field(sqrt(un^2 + vn^2))

sn = @lift begin
    parent(un) .= parent(u[$n])
    parent(vn) .= parent(v[$n])
    compute!(s)
    sn = interior(s)
    view(sn, :, :, 1)
end

fig = Figure()
axs = Axis(fig[1, 1])
axζ = Axis(fig[2, 1])
axT = Axis(fig[3, 1])

heatmap!(axs, sn)
heatmap!(axζ, ζn)
heatmap!(axT, Tn)

save("snapshot.png", fig)


record(fig, "baroclinic_wave.mp4", 1:Nt, framerate = 8) do nn
    @info nn
    n[] = nn
end
