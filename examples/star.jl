include("../generalized_cylinder_inductance.jl")

using Plots, Interpolations

# Example: 5-pointed star (inductance per unit length)
h = 0.02  # height in meters
N = 10   # number of turns
star_points = [
    [0.0, 1.0],
    [0.2245, 0.3090],
    [0.9511, 0.3090],
    [0.3633, -0.1180],
    [0.5878, -0.8090],
    [0.0, -0.3819],
    [-0.5878, -0.8090],
    [-0.3633, -0.1180],
    [-0.9511, 0.3090],
    [-0.2245, 0.3090],
    [0.0, 1.0]
] ./ 100

star_t = LinRange(0, 1, length(star_points))
star_path = linear_interpolation(star_t, star_points)
star_ds = Base.Fix1(derivative, star_path) ∘ Base.Fix2(%, 1)

t = LinRange(0, 1, 101)
ℓs = [generalized_cylinder_inductance_per_length(star_path, h, ti; n = N / h, ds = star_ds)
      for ti in t]

xs = [star_path(ti)[1] for ti in t]
ys = [star_path(ti)[2] for ti in t]
star_plt = plot(
    xs * 100,
    ys * 100,
    line_z = ℓs * 1e6,
    linewidth = 8,
    xlabel = "x (cm)",
    ylabel = "y (cm)",
    legend = false,
    colorbar = true,
    grid = false,
    aspect_ratio = :equal,
    size = (400, 400),
    title = "Inductance Contribution (μH/m)"
)
