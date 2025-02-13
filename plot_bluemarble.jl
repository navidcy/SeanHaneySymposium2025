using GLMakie, FileIO
using Downloads: download

earth_img = load(download("https://upload.wikimedia.org/wikipedia/commons/5/56/Blue_Marble_Next_Generation_%2B_topography_%2B_bathymetry.jpg"))

n = 1024 ÷ 4 # 2048
θ = LinRange(0, π, n)
φ = LinRange(0, 2π, 2 * n)
x = [cos(φ) * sin(θ) for θ in θ, φ in φ]
y = [sin(φ) * sin(θ) for θ in θ, φ in φ]
z = [cos(θ) for θ in θ, φ in φ]

# Create figure
fig = Figure(size = (800, 600))
kw = (elevation= deg2rad(-10), azimuth=deg2rad(90), aspect=:equal)
ax = Axis3(fig[1, 1]; kw...)

# Plot textured mesh
surface!(ax, x, y, z;
    color = earth_texture,
    shading = NoShading,
    backlight = 1.5f0
    )

hidedecorations!(ax)
hidespines!(ax)

fig
