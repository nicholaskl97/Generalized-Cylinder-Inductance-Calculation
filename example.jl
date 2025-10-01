include("generalized_cylinder_inductance.jl")
include("circular_cylinder_inductance.jl")

# Example: right circular cylinder
R = 0.01  # radius in meters
h = 0.02  # height in meters
N = 10   # number of turns
circular_path = (t) -> R .* (cospi(2t), sinpi(2t))
circular_ds = (t) -> (2R * π) .* (-sinpi(2t), cospi(2t))

L = generalized_cylinder_inductance(circular_path, h; N, ds = circular_ds)
L_nagoaka = circular_cylinder_inductance(R, h; N)
@assert isapprox(L, L_nagoaka; rtol = 1e-5)

println("Inductance of right circular cylinder: $(1e6 * L) μH")
