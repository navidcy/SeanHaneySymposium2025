using ClimaOcean
using Oceananigans
using Oceananigans
using Oceananigans.Units
using Printf

arch = GPU()
resolution = 1//4
Nx = 360 ÷ resolution
Ny = 160 ÷ resolution
Nz = 10

grid = LatitudeLongitudeGrid(arch;
                             size = (Nx, Ny, Nz),
                             halo = (7, 7, 3),
                             latitude = (-80, 80),
                             longitude = (0, 360),
                             z = (-1000, 0))

closure = CATKEVerticalDiffusivity()
equation_of_state = TEOS10EquationOfState()
buoyancy = SeawaterBuoyancy(; equation_of_state)
model = HydrostaticFreeSurfaceModel(; grid, buoyancy, closure, tracers=(:T, :S, :e))
ocean = ocean_simulation(grid; closure=vertical_mixing)

Tatm(λ, φ, z=0) = 30 * cosd(φ)
Tᵢ(λ, φ, z) = 30 * (1 - tanh((abs(φ) - 40) / 5)) / 2 + rand()
Sᵢ(λ, φ, z) = 28 - 5e-3 * z + rand()
set!(ocean.model, T=Tᵢ, S=Sᵢ)

atmos_grid = LatitudeLongitudeGrid(arch;
                                   size = (Nx, Ny),
                                   halo = (7, 7),
                                   latitude = (-80, 80),
                                   longitude = (0, 360),
                                   topology = (Periodic, Bounded, Flat))

# Build and run an OceanSeaIceModel (with no sea ice component) forced by JRA55 reanalysis
atmos_times = range(0, 1day, length=24)
atmosphere = PrescribedAtmosphere(atmos_grid, atmos_times)

zonal_wind(λ, φ) = 4 * sind(2φ)^2 - 2 * exp(-(abs(φ) - 12)^2 / 72)
sunlight(λ, φ) = -200 - 600 * cosd(φ)^2

Ta = CenterField(atmos_grid)
ua = CenterField(atmos_grid)
Qsw = CenterField(atmos_grid)

#set!(Ta, (λ, φ) -> 273.15 + 30 * cosd(φ)^2)
set!(Ta, Tatm)
set!(ua, zonal_wind)
set!(Qsw, sunlight)

parent(atmosphere.tracers.T) .= parent(Ta) .+ 273.15
parent(atmosphere.velocities.u) .= parent(ua)
parent(atmosphere.tracers.q) .= 0
parent(atmosphere.downwelling_radiation.shortwave) .= parent(Qsw)

radiation = Radiation(arch)
coupled_model = ClimaOcean.OceanSeaIceModel(ocean; atmosphere, radiation)
simulation = Simulation(coupled_model, Δt=2minutes, stop_time=4 * 360days)

function progress(sim)
    ocean = sim.model.ocean
    T = ocean.model.tracers.T
    S = ocean.model.tracers.S
    u, v, w = ocean.model.velocities

    msg = @sprintf("%d: %s, max|u|: (%.2e, %.2e, %.2e)",
                   iteration(sim), prettytime(sim),
                   maximum(abs, u),
                   maximum(abs, v),
                   maximum(abs, w))

    msg *= @sprintf(", extrema(T): (%.2f, %.2f), extrema(S): (%.2f, %.2f)", 
                    minimum(T), maximum(T), minimum(S), maximum(S))

    ΣQ = simulation.model.fluxes.total.ocean.heat
    Ql = simulation.model.fluxes.turbulent.fields.latent_heat
    Qs = simulation.model.fluxes.turbulent.fields.sensible_heat
    Fv = simulation.model.fluxes.turbulent.fields.water_vapor
    τx = simulation.model.fluxes.turbulent.fields.x_momentum
    τy = simulation.model.fluxes.turbulent.fields.y_momentum

    msg *= @sprintf(", extrema(ΣQ): (%.2e, %.2e), extrema(Ql): (%.2e, %.2e), extrema(Qs): (%.2e, %.2e), max|τ|: (%.2e, %.2e)",
                    minimum(ΣQ),
                    maximum(ΣQ),
                    minimum(Ql),
                    maximum(Ql),
                    minimum(Qs),
                    maximum(Qs),
                    maximum(abs, τx),
                    maximum(abs, τy))
                    
    @info msg
                    
    return nothing
end

add_callback!(simulation, progress, IterationInterval(100))

u, v, w = ocean.model.velocities
T = ocean.model.tracers.T
S = ocean.model.tracers.S
ζ = ∂x(v) - ∂y(u)
fields = (; u, v, w, T, S, ζ)

fluxes = (;
  ΣQ = simulation.model.fluxes.total.ocean.heat,
  Ql = simulation.model.fluxes.turbulent.fields.latent_heat,
  Qs = simulation.model.fluxes.turbulent.fields.sensible_heat,
  Fv = simulation.model.fluxes.turbulent.fields.water_vapor,
  τx = simulation.model.fluxes.turbulent.fields.x_momentum,
  τy = simulation.model.fluxes.turbulent.fields.y_momentum
)

ow = JLD2OutputWriter(ocean.model, merge(fields, fluxes),
                      filename = "ocean2.jld2",
                      indices = (:, :, grid.Nz),
                      schedule = TimeInterval(4days),
                      overwrite_existing = true)

ocean.output_writers[:surface] = ow

run!(simulation)
