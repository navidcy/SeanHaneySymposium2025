using Oceananigans
using Oceananigans.Units
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using Printf
using CairoMakie

filename = "baroclinic_wave"

arch = GPU()

resolution = 1
Nx = 360 ÷ resolution
Ny = 160 ÷ resolution
Nz = 10

grid = LatitudeLongitudeGrid(arch;
                             size = (Nx, Ny, Nz),
                             halo = (7, 7, 7),
                             latitude = (-80, 80),
                             longitude = (0, 360),
                             z = (-3000, 0))

momentum_advection = WENOVectorInvariant(order=9)
tracer_advection   = WENO(order=7)
coriolis = HydrostaticSphericalCoriolis()
buoyancy = SeawaterBuoyancy(equation_of_state=TEOS10EquationOfState())

model = HydrostaticFreeSurfaceModel(; grid, coriolis,
                                      buoyancy, tracers=(:T, :S),
                                      momentum_advection, tracer_advection)

Tᵢ(λ, φ, z) = 30 * (1 - tanh((abs(φ) - 45) / 8)) / 2 + rand()
Sᵢ(λ, φ, z) = 28 - 5e-3 * z + rand()
set!(model, T=Tᵢ, S=Sᵢ)

simulation = Simulation(model, Δt=5minutes, stop_time=365days)

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

add_callback!(simulation, progress, IterationInterval(250))

u, v, w = model.velocities
T = model.tracers.T
S = model.tracers.S
ζ = ∂x(v) - ∂y(u)
fields = (; u, v, w, T, S, ζ)

ow = JLD2OutputWriter(model, fields,
                      filename = filename * ".jld2",
                      indices = (:, :, grid.Nz),
                      schedule = TimeInterval(1days),
                      overwrite_existing = true)

simulation.output_writers[:surface] = ow

run!(simulation)

u = FieldTimeSeries(filename * ".jld2", "u"; backend = OnDisk())
v = FieldTimeSeries(filename * ".jld2", "v"; backend = OnDisk())
T = FieldTimeSeries(filename * ".jld2", "T"; backend = OnDisk())
ζ = FieldTimeSeries(filename * ".jld2", "ζ"; backend = OnDisk())

times = u.times
Nt = length(times)

n = Observable(1)

ζn = @lift ζ[$n]
Tn = @lift T[$n]

un = Field{Face, Center, Nothing}(u.grid)
vn = Field{Center, Face, Nothing}(v.grid)
s = Field(sqrt(un^2 + vn^2))

sn = @lift begin
    parent(un) .= parent(u[$n])
    parent(vn) .= parent(v[$n])
    compute!(s)
end

fig = Figure(size=(800, 1200))

kwargs_axis = (xticks = 0:60:360, yticks = -80:40:80)
axs = Axis(fig[2, 1]; kwargs_axis...)
axζ = Axis(fig[3, 1]; kwargs_axis...)
axT = Axis(fig[4, 1]; kwargs_axis...)

hm = heatmap!(axs, sn, colorrange=(0, 2))
Colorbar(fig[2, 2], hm, label = "Surface speed (m s⁻¹)")

hm = heatmap!(axζ, ζn, colormap=:balance, colorrange=(-2e-5, 2e-5))
Colorbar(fig[3, 2], hm, label = "relative vorticity (s⁻¹)")

hm = heatmap!(axT, Tn, colormap=:thermal, colorrange=(0, 32))
Colorbar(fig[4, 2], hm, label = "surface temperature (ᵒC)")

title = @lift @sprintf("%s", prettytime(times[$n]))
fig[1, :] = Label(fig, title, tellwidth=false)

n[] = 1
save(filename * "_initial_snapshot.png", fig)

n[] = Nt
save(filename * "_final_snapshot.png", fig)

record(fig, filename * ".mp4", 1:Nt, framerate = 16) do nn
    @info nn
    n[] = nn
end
