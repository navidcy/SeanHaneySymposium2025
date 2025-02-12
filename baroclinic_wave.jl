using Oceananigans
using Oceananigans.Units
using SeawaterPolynomials: TEOS10EquationOfState
using Printf

arch = GPU()
resolution = 1 #1//4
Nx = 360 ÷ resolution
Ny = 160 ÷ resolution
Nz = 10

grid = LatitudeLongitudeGrid(arch;
                             size = (Nx, Ny, Nz),
                             halo = (7, 7, 7),
                             latitude = (-80, 80),
                             longitude = (0, 360),
                             z = (-1000, 0))

momentum_advection = WENOVectorInvariant()
tracer_advection = WENO(order=7)
buoyancy = SeawaterBuoyancy(equation_of_state=TEOS10EquationOfState())
model = HydrostaticFreeSurfaceModel(; grid, momentum_advection, tracer_advection,
                                    buoyancy, tracers=(:T, :S))

Tatm(λ, φ, z=0) = 30 * cosd(φ)
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

ocean.output_writers[:surface] = ow

run!(simulation)

