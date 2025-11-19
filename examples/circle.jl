include("../generalized_cylinder_inductance.jl")
include("../circular_cylinder_inductance.jl")

using Plots

# Example: right circular cylinder
R = 0.01  # radius in meters
h = 0.02  # height in meters
N = 10   # number of turns
circular_path = (t) -> R .* [cospi(2t), sinpi(2t)]
circular_ds = (t) -> (2R * π) .* [-sinpi(2t), cospi(2t)]

L = generalized_cylinder_inductance(circular_path, h; N, ds = circular_ds)
L_nagoaka = circular_cylinder_inductance(R, h; N)
@assert isapprox(L, L_nagoaka; rtol = 1e-5)

println("Inductance of right circular cylinder: $(1e6 * L) μH")

# Example: right circular cylinder (inductance per unit length)
t = LinRange(0, 1, 101)
ℓs = [generalized_cylinder_inductance_per_length(
          circular_path,
          h,
          ti;
          n = N / h,
          ds = circular_ds
      ) for ti in t]
@assert isapprox(sum(ℓs[1:(end - 1)]) * 2π * R / (length(t) - 1), L; rtol = 1e-5)

xs = [circular_path(ti)[1] for ti in t]
ys = [circular_path(ti)[2] for ti in t]
circle_plt = plot(
    xs * 100,
    ys * 100,
    line_z = ℓs * 1e6,
    linewidth = 8,
    xlabel = "x (cm)",
    ylabel = "y (cm)",
    clim = (ℓs[1] * 1e6 * 0.999, ℓs[1] * 1e6 * 1.001),
    legend = false,
    colorbar = true,
    aspect_ratio = :equal,
    size = (400, 400),
    grid = false,
    title = "Inductance Contribution (μH/m)"
)
