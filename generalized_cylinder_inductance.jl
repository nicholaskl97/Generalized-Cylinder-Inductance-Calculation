using HCubature: hcubature, hquadrature
using ForwardDiff: derivative
using LinearAlgebra: norm, ⋅

"""
    generalized_cylinder_inductance(path, height; N, ds, <kwargs>)

Calculate the inductance, in H, of a generalized cylindrical air-core inductor using the
formula in:

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
  - `ds`: A function that takes a single argument `t` in the range `[0, 1]` and returns the
    derivative of `s` at `t`. By default, this is computed using automatic differentiation.
  - `norm`, `rtol`, `atol`, `maxevals`, `initdiv`: Additional parameters for numerical
    integration. See `hcubature` documentation for details, as these are passed directly to
    it, except `atol`, which is rescaled so that it is also provided in H.
"""
function generalized_cylinder_inductance(
        s,
        height::T;
        N = 1,
        ds = Base.Fix1(derivative, s),
        norm = norm,
        rtol = sqrt(eps(T)),
        atol = 0,
        maxevals = typemax(Int),
        initdiv = 1
) where {T <: Real}
    C = (N^2 / height * T(4e-7))

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

"""
    generalized_cylinder_inductance_per_length(path, height, t; n, ds, <kwargs>)

Calculate the inductance per unit length, in H/m, of a generalized cylindrical air-core
inductor using the formula in:

N. Klugman, J. Vedral and J. Lang, "Self-Inductance of an Extrusion of a Planar Curve," 2020
International Applied Computational Electromagnetics Society Symposium (ACES), Monterey, CA,
USA, 2020, doi: 10.23919/ACES49320.2020.9196194.

# Positional Arguments
  - `s`: A function that takes a single argument `t` in the range `[0, 1]` and returns a 2D
    point (as a tuple or array, in meters) representing the boundary of the bottom face of
    the conductor.
  - `height`: The height of the cylinder, in meters.
  - `t`: A parameter in the range `[0, 1]` at which to evaluate the inductance per unit
    length.

# Keyword Arguments
  - `n`: The number of turns per unit length (default is `1 / height`).
  - `ds`: A function that takes a single argument `t` in the range `[0, 1]` and returns the
    derivative of `s` at `t`. By default, this is computed using automatic differentiation.
  - `norm`, `rtol`, `atol`, `maxevals`, `initdiv`: Additional parameters for numerical
    integration. See `hquadrature` documentation for details, as these are passed directly
    to it, except `atol`, which is rescaled so that it is also provided in H/m.
"""
function generalized_cylinder_inductance_per_length(
        s,
        height::T,
        t::T;
        n = one(T) / height,
        ds = Base.Fix1(derivative, s),
        norm = norm,
        rtol = sqrt(eps(T)),
        atol = 0,
        maxevals = typemax(Int),
        initdiv = 1
) where {T <: Real}
    C = (n^2 * height * T(2e-7))

    # Integrate over [0, t] and [t, 1] separately to avoid the singularity at τ = t.
    I1 = if t == zero(T)
        zero(T)
    else
        hquadrature(
            Base.Fix2(inductance_per_length_integrand, (t, s, ds, height)),
            zero(T),
            t;
            rtol,
            atol = atol / C,
            maxevals,
            initdiv
        )[1]
    end

    I2 = if t == zero(T)
        zero(T)
    else
        hquadrature(
            Base.Fix2(inductance_per_length_integrand, (t, s, ds, height)),
            t,
            one(T);
            rtol,
            atol = atol / C,
            maxevals,
            initdiv
        )[1]
    end

    return C * (I1 + I2)
end

function inductance_per_length_integrand(τ, p)
    t, s, ds, h = p

    Δs = norm(s(t) .- s(τ))
    dst = ds(t)
    ŝ_t = dst ./ norm(dst)

    return (asinh(h / Δs) - h / (Δs + sqrt(h^2 + Δs^2))) * (ŝ_t ⋅ ds(τ))
end
