using GLMakie, FileIO
using Downloads: download

#earth_img = load(download("https://upload.wikimedia.org/wikipedia/commons/5/56/Blue_Marble_Next_Generation_%2B_topography_%2B_bathymetry.jpg"))
earth_img = load("blue-marble-only-land.png")

#=
using Oceananigans
using Oceananigans.Units
using Oceananigans.ImmersedBoundaries: mask_immersed_field!
using GLMakie
using Printf

function geographic2cartesian(λ, φ, r=1)
    Nλ = length(λ)
    Nφ = length(φ)
    λ = repeat(reshape(λ, Nλ, 1), 1, Nφ)
    φ = repeat(reshape(φ, 1, Nφ), Nλ, 1)

    λ_azimuthal = λ .+ 180  # Convert to λ ∈ [0°, 360°]
    φ_azimuthal = 90 .- φ   # Convert to φ ∈ [0°, 180°] (0° at north pole)

    x = @. r * cosd(λ_azimuthal) * sind(φ_azimuthal)
    y = @. r * sind(λ_azimuthal) * sind(φ_azimuthal)
    z = @. r * cosd(φ_azimuthal)

    return x, y, z
end

filename = "simple_forced_global.jld2"
s = FieldTimeSeries(filename, "s") #; backend = OnDisk())
T = FieldTimeSeries(filename, "T") #; backend = OnDisk())
ζ = FieldTimeSeries(filename, "ζ") #; backend = OnDisk())
=#

land = interior(ζ.grid.immersed_boundary.bottom_height) .== 0
Nt = length(s)
Nx, Ny, Nz = size(s.grid)
for n in 1:Nt
    interior(ζ[n], :, 1:Ny, :)[land] .= NaN
    interior(s[n])[land] .= NaN
    interior(T[n])[land] .= NaN
end

times = s.times
Nt = length(times)

grid = s.grid
λ = λnodes(grid, Center(), Center(), Center())
φ = φnodes(grid, Center(), Center(), Center())
x, y, z = geographic2cartesian(λ, φ)

Nx, Ny, Nz = size(grid)

n = Observable(1)

ζn = @lift interior(ζ[$n], :, :, 1) .* 1e5
sn = @lift interior(s[$n], :, :, 1)
Tn = @lift interior(T[$n], :, :, 1)

fig = Figure(size=(1400, 700))

kw = (elevation= deg2rad(-10), azimuth=deg2rad(90), aspect=:equal)

axs = Axis3(fig[1, 1]; kw...)
axζ = Axis3(fig[1, 2]; kw...)
axT = Axis3(fig[1, 3]; kw...)

N = 1024 ÷ 4 # 2048
θ = LinRange(0, π, N)
φ = LinRange(0, 2π, 2N)
r = 0.98
xbm = [r * cos(φ) * sin(θ) for θ in θ, φ in φ]
ybm = [r * sin(φ) * sin(θ) for θ in θ, φ in φ]
zbm = [r * cos(θ) for θ in θ, φ in φ]

# Create figure
#fig = Figure(size = (800, 600))
#kw = (elevation= deg2rad(-10), azimuth=deg2rad(90), aspect=:equal)

# Plot textured mesh
surface!(axs, xbm, ybm, zbm;
         color = earth_img,
         shading = NoShading,
         backlight = 1.5f0)

surface!(axζ, xbm, ybm, zbm;
         color = earth_img,
         shading = NoShading,
         backlight = 1.5f0)

surface!(axT, xbm, ybm, zbm;
         color = earth_img,
         shading = NoShading,
         backlight = 1.5f0)

sf = surface!(axs, x, y, z, color=sn, colorrange=(0, 0.7), colormap=:viridis) #, nan_color=:lightgray)
Colorbar(fig[2, 1], sf, width=Relative(0.6), vertical=false, label = "Surface speed (m s⁻¹)", labelsize=20)

sf = surface!(axζ, x, y, z, color=ζn, colorrange=(-2, 2), colormap=:balance) #, nan_color=:lightgray)
Colorbar(fig[2, 2], sf, width=Relative(0.6), vertical=false, label = "relative vorticity (10⁻⁵ s⁻¹)", labelsize=20)

sf = surface!(axT, x, y, z, color=Tn, colorrange=(0, 31), colormap=:thermal) #, nan_color=:lightgray)
Colorbar(fig[2, 3], sf, width=Relative(0.6), vertical=false, label = "surface temperature (ᵒC)", labelsize=20)

# surface!(axζ, x, y, z, color=ζn, colorrange=(-2e-5, 2e-5), colormap=:balance, nan_color=:lightgray)
# surface!(axs, x, y, z, color=sn, colorrange=(0, 1), colormap=:viridis, nan_color=:lightgray)
# surface!(axT, x, y, z, color=Tn, colorrange=(0, 1), colormap=:viridis, nan_color=:lightgray)

t = s.times
title = @lift @sprintf("%s days", round((times[$n]+ 1e-6) / days, digits=2))
Label(fig[0, :], title, fontsize=24)

hidedecorations!(axζ)
hidedecorations!(axs)
hidedecorations!(axT)
hidespines!(axζ)
hidespines!(axs)
hidespines!(axT)

rowgap!(fig.layout, 1, Relative(-0.1))
colgap!(fig.layout, 1, Relative(-0.07))
colgap!(fig.layout, 2, Relative(-0.07))

record(fig, "ecco_spinup_bm.mp4", 1:Nt, framerate = 24) do nn
    @info "Drawing frame $nn of $Nt..."
    n[] = nn
end

