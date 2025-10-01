using Elliptic: ellipke

"""
    circular_cylinder_inductance(R, height; N, μ_r)

Calculate the inductance, in H, of a right circular cylindrical conuctor.

# Positional Arguments
  - `R`: The radius of the cylinder, in meters.
  - `height`: The height of the cylinder, in meters.

# Keyword Arguments
  - `N`: The number of turns (default is 1).
  - `μ_r`: The relative permeability of the medium (default is 1).
"""
function circular_cylinder_inductance(R::T, height; N = 1, μ_r = 1) where T<:Real
    return μ_r * T(4e-7) * π * N^2 * (π * R^2) / height * nagaokas_coeff(2R / height)
end

"""
    nagaokas_coeff(u)

Calculate Nagaoka's coefficient for a circular, cylindrical inductor with diameter-to-height
ratio `u`.
"""
function nagaokas_coeff(u::T) where T<:Real
    K, E = T.(ellipke(u^2 / (1 + u^2)))
    return T(4 / 3π) * (sqrt(1 + u^2) * (K - E) / u^2 + sqrt(1 + u^2) * E - u)
end
