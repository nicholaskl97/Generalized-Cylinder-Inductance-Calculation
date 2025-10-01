using HCubature: hcubature
using ForwardDiff: derivative
using LinearAlgebra: norm, ⋅

"""
    generalized_cylinder_inductance(path, height; N, μ_r, ds)

Calculate the inductance, in H, of a generalized cylindrical conductor using formula in:

N. Klugman, J. Vedral and J. Lang, "Self-Inductance of an Extrusion of a Planar Curve," 2020
International Applied Computational Electromagnetics Society Symposium (ACES), Monterey, CA,
USA, 2020, doi: 10.23919/ACES49320.2020.9196194.

# Positional Arguments
  - `s`: A function that takes a single argument `t` in the range `[0, 1]` and returns a 2D
    point (as a tuple or array, in meters) representing the boundary of the bottom face of
    the conductor.
  - `height`: The height of the cylinder, in meters.

  # Keyword Arguments
  - `N`: The number of turns (default is 1).
  - `μ_r`: The relative permeability of the medium (default is 1).
  - `ds`: A function that takes a single argument `t` in the range `[0, 1]` and returns the
    derivative of `s` at `t`. By default, this is computed using automatic differentiation.
  - `norm`, `rtol`, `atol`, `maxevals`, `initdiv`: Additional parameters for numerical
    integration. See `hcubature` documentation for details, as these are passed directly to
    it, except `atol`, which is rescaled so that it is also provided in H.
"""
function generalized_cylinder_inductance(s, height::T; N = 1, μ_r = 1,
        ds = Base.Fix1(derivative, s), norm = norm, rtol = sqrt(eps(T)), atol = 0,
        maxevals = typemax(Int), initdiv = 1) where {T <: Real}
    C = (μ_r * N^2 / height * T(4e-7))

    # Instead of integrating over (t1, t2) ∈ [0, 1]² , we wish to exploit symmetry and
    # integrate over [0, 1] x [0, t1] and multiply by 2, which is more efficient and moves
    # the singularites when t1 = t2 to the boundary (which `hcubature` does not sample).
    # However, `hcubature` does not support variable limits of integration, so we perform
    # the change of variables: u = t2 / t1. This gives us an integral over [0, 1]² with an
    # additional factor of t1 in the integrand.
    return C * hcubature(
        Base.Fix2(inductance_integrand, (s, ds, height)),
        zeros(T, 2),
        ones(T, 2);
        rtol,
        atol = atol / C,
        maxevals,
        initdiv
    )[1]
end

function inductance_integrand(t, p)
    s, ds, h = p

    t1, u = t
    t2 = u * t1

    Δs = norm(s(t1) .- s(t2))

    return (asinh(h / Δs) - h / (Δs + sqrt(h^2 + Δs^2))) * (ds(t1) ⋅ ds(t2)) * t1
end
