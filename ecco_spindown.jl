import ClimaOcean
using Oceananigans
using Oceananigans.Units
using SeawaterPolynomials.TEOS10: TEOS10EquationOfState
using CFTime
using Printf

arch = GPU()
resolution = 1//4
Nx = 360 ÷ resolution
Ny = 160 ÷ resolution
Nz = 10

grid = LatitudeLongitudeGrid(arch;
                             size = (Nx, Ny, Nz),
                             halo = (7, 7, 7),
                             latitude = (-80, 80),
                             longitude = (0, 360),
                             z = (-1000, 0))

bathymetry = ClimaOcean.regrid_bathymetry(grid) # builds gridded bathymetry based on ETOPO1
grid = ImmersedBoundaryGrid(grid, GridFittedBottom(bathymetry))

momentum_advection = WENOVectorInvariant()
tracer_advection = WENO(order=7)
coriolis = HydrostaticSphericalCoriolis()
buoyancy = SeawaterBuoyancy(equation_of_state=TEOS10EquationOfState())
model = HydrostaticFreeSurfaceModel(; grid, momentum_advection, tracer_advection,
                                    coriolis, buoyancy, tracers=(:T, :S))

dates = DateTimeProlepticGregorian(1993, 1, 1)
set!(model, T = ClimaOcean.ECCOMetadata(:temperature; dates),
            S = ClimaOcean.ECCOMetadata(:salinity; dates))

simulation = Simulation(model, Δt=2minutes, stop_time=180days)

function progress(sim)
    T = sim.model.tracers.T
    S = sim.model.tracers.S
    u, v, w = sim.model.velocities

    u = interior(u)
    v = interior(v)
    w = interior(w)

    msg = @sprintf("%d: %s, max|u|: (%.2e, %.2e, %.2e)",
                   iteration(sim), prettytime(sim),
                   maximum(abs, u),
                   maximum(abs, v),
                   maximum(abs, w))

    T = interior(T)
    S = interior(S)

    msg *= @sprintf(", extrema(T): (%.2f, %.2f), extrema(S): (%.2f, %.2f)", 
                    minimum(T), maximum(T), minimum(S), maximum(S))

    @info msg
                        
    return nothing
end

add_callback!(simulation, progress, IterationInterval(100))

u, v, w = model.velocities
T = model.tracers.T
S = model.tracers.S
ζ = ∂x(v) - ∂y(u)
s = sqrt(u^2 + v^2)
fields = (; u, v, w, T, S, s, ζ)

ow = JLD2OutputWriter(model, fields,
                      filename = "ecco_spindown.jld2",
                      indices = (:, :, grid.Nz),
                      schedule = TimeInterval(12hours),
                      overwrite_existing = true)

simulation.output_writers[:surface] = ow

run!(simulation)

