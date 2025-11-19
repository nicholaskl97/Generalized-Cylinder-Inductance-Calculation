include("../generalized_cylinder_inductance.jl")

using Plots

# Example: semi-circular cylinder (inductance per unit length)
R = 0.01  # radius in meters
h = 0.02  # height in meters
N = 10   # number of turns

semi_circular_path = function (t)
    if t < 0.5
        return R .* [cospi(2t), sinpi(2t)]
    else
        return (R * (4t - 3), zero(R))
    end
end
semi_circular_ds = function (t)
    if t < 0.5
        return (2R * π) .* [-sinpi(2t), cospi(2t)]
    else
        return (4R, zero(R))
    end
end

t = LinRange(0, 1, 101)
ℓs = [generalized_cylinder_inductance_per_length(
          semi_circular_path,
          h,
          ti;
          n = N / h,
          ds = semi_circular_ds
      ) for ti in t]

xs = [semi_circular_path(ti)[1] for ti in t]
ys = [semi_circular_path(ti)[2] for ti in t]
semi_circle_plt = plot(
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

# Example: semi-circle with sinusoidal segment (inductance per unit length)
curious_path = function (t)
    if t < 0.5
        return R .* [cospi(2t), sinpi(2t)]
    else
        return [R * (4t - 3), R / 10 * sinpi(10t)]
    end
end

t = LinRange(0, 1, 101)
ℓs = [generalized_cylinder_inductance_per_length(curious_path, h, ti; n = N / h)
      for ti in t]

xs = [curious_path(ti)[1] for ti in t]
ys = [curious_path(ti)[2] for ti in t]
curious_plt = plot(
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

plot(semi_circle_plt, curious_plt, layout = (2, 1))
