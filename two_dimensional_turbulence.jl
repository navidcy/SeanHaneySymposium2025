using Oceananigans, CairoMakie

# The third dimension is "flattened" to reduce the domain from three to two dimensions.
topology = (Periodic, Periodic, Flat)
architecture = CPU() # CPU() works just fine too for this small example.
x = y = (0, 2π)
grid = RectilinearGrid(architecture; size=(256, 256), x, y, topology)

model = NonhydrostaticModel(; grid, advection=WENO(order=5))

ϵ(x, y) = 2rand() - 1 # Uniformly-distributed random numbers between [-1, 1).
set!(model, u=ϵ, v=ϵ)

simulation = Simulation(model; Δt=0.01, stop_time=10)

wizard = TimeStepWizard(cfl=0.8, max_change=1.2, max_Δt=0.5)
simulation.callbacks[:wizard] = Callback(wizard, IterationInterval(10))

using Printf

function progress_message(sim)
    max_abs_u = maximum(abs, sim.model.velocities.u)
    walltime = prettytime(sim.run_wall_time)

    return @info @sprintf("Iteration: %04d, time: %1.3f, Δt: %.2e, max(|u|) = %.1e, wall time: %s\n",
                          iteration(sim), time(sim), sim.Δt, max_abs_u, walltime)
end

add_callback!(simulation, progress_message, IterationInterval(10))

u, v, w = model.velocities
ζ = ∂x(v) - ∂y(u)

filename = "two_dimensional_turbulence"

simulation.output_writers[:fields] = JLD2OutputWriter(model, (; ζ),
                                                      schedule = TimeInterval(0.2),
                                                      filename = filename * ".jld2",
                                                      overwrite_existing = true)


run!(simulation)

heatmap(ζ; colormap = :balance, colorrange = (-3, 3), axis=(; aspect=1))
save(filename * ".png", current_figure())

## Make a movie

ζ_timeseries = FieldTimeSeries(filename * ".jld2", "ζ")
times = ζ_timeseries.times

n = Observable(1)

ζ = @lift ζ_timeseries[$n]

title = @lift "vorticity, t = " * string(round(times[$n], digits=2))

fig = Figure(fontsize=24)

ax = Axis(fig[1, 1]; title, aspect = 1)
heatmap!(ax, ζ; colormap = :balance, colorrange = (-3, 3))

frames = 1:length(times)
record(fig, filename * ".mp4", frames, framerate=18) do i
    n[] = i
end
