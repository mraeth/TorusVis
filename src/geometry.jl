# ------------------------------------------------------------
# Device geometry abstraction.
#
# Angle/coordinate convention used throughout TorusVis (deliberately different
# from the ad-hoc θ=toroidal/φ=poloidal convention in the original scripts,
# to match ordinary plasma-physics usage):
#   ζ (zeta)  -- toroidal angle, and by construction also the lab-frame
#                azimuthal angle around the z-axis (every concrete geometry's
#                `surface_point` must respect this: x = ... * cos(ζ), y = ... * sin(ζ)).
#   θ (theta) -- poloidal angle within the device's cross-section at a given ζ.
#   ρ (rho)   -- minor-radius-like coordinate labelling nested cross-section
#                contours; ρ = 0 is the magnetic axis. Units are whatever the
#                concrete geometry's own `minor_radius_bounds` says (physical
#                metres for `MillerGeometry`, normalized 0..1 for `FourierGeometry`).
#
# Any device is described purely through this interface -- nothing downstream
# (voxelize.jl, render.jl) contains device-specific code.
# ------------------------------------------------------------

abstract type AbstractGeometry end

"""
    surface_point(g, ζ, ρ, θ) -> (x, y, z)

Lab-frame Cartesian coordinates of the point at toroidal angle `ζ`, minor
radius `ρ`, poloidal angle `θ`. Every implementation must place this point at
lab-frame azimuth `ζ` (i.e. `atan(y, x) == ζ`).
"""
function surface_point(g::AbstractGeometry, ζ, ρ, θ) end

"""
    invert_cross_section(g, ζ, R_cyl, z) -> (ρ, θ)

Inverse of the cross-section shape at toroidal angle `ζ`: given a point at
cylindrical radius `R_cyl = hypot(x, y)` and height `z`, recover the `(ρ, θ)`
that maps to it. Each geometry subtracts its own (possibly ζ-dependent)
magnetic-axis radius internally before inverting -- callers never need to
know where that axis sits.
"""
function invert_cross_section(g::AbstractGeometry, ζ, R_cyl, z) end

"""
    field_pitch(g, ζ) -> Float64

dθ/dζ along a magnetic field line at toroidal angle `ζ` (the reciprocal of the
local safety factor `q`). Kept ζ-dependent in the interface so geometries with
shear can be added later, even though both current geometries return a constant.
"""
function field_pitch(g::AbstractGeometry, ζ) end

"""
    wedge(g) -> Union{Nothing, NTuple{2,Float64}}

`(ζ_A, ζ_B)` toroidal bounds of a wedge cut out of the device body, or
`nothing` if the device is a fully closed torus.
"""
function wedge(g::AbstractGeometry) end

"""
    nfp(g) -> Int

Number of toroidal field periods (discrete rotational symmetry). `1` for an
axisymmetric device.
"""
function nfp(g::AbstractGeometry) end

"""
    minor_radius_bounds(g) -> (ρ_lo, ρ_hi)

Full range of the minor-radius coordinate `ρ`, in whatever units this
geometry uses for it (see module docstring above).
"""
function minor_radius_bounds(g::AbstractGeometry) end

"""
    cut_point(g, ζ, ρ, θ; off=0.0, s=1)

`surface_point(g, ζ, ρ, θ)` nudged by `off*s` along the lab-frame azimuthal
(toroidal) direction at `ζ` -- i.e. off the exact cut plane, to avoid
z-fighting between coincident surfaces. Generic for any `AbstractGeometry`
because `ζ` is always the lab-frame azimuth (see `surface_point`'s contract).
"""
function cut_point(g::AbstractGeometry, ζ, ρ, θ; off=0.0, s=1)
    x, y, z = surface_point(g, ζ, ρ, θ)
    return (x + off * s * (-sin(ζ)), y + off * s * cos(ζ), z)
end

# ------------------------------------------------------------
# Cross-section shape inversion, shared machinery.
#
# Every concrete geometry's cross-section, at a fixed ζ, has the form
# (u, v) = ρ .* (a(θ), b(θ)) for some direction functions a(θ), b(θ) --
# so (u, v) and (a(θ), b(θ)) always share a direction, and θ can be found by
# matching that direction against a precomputed lookup table, then
# ρ = |(u,v)| / |(a(θ),b(θ))|. This is exact and cheap (no root-finding),
# and is reused by both MillerGeometry (one LUT, built once -- its shape
# doesn't depend on ζ) and FourierGeometry (one LUT per cached ζ-bin, since
# its shape does depend on ζ).
# ------------------------------------------------------------
struct ShapeLUT
    θgrid::Vector{Float64}
    angle::Vector{Float64}
    a::Vector{Float64}
    b::Vector{Float64}
end

function build_shape_lut(ab::Function; n=4001)
    θgrid = collect(range(-π, π, length=n))
    a = Vector{Float64}(undef, n)
    b = Vector{Float64}(undef, n)
    for (i, θ) in enumerate(θgrid)
        a[i], b[i] = ab(θ)
    end
    angle = atan.(b, a)
    perm = sortperm(angle)
    return ShapeLUT(θgrid[perm], angle[perm], a[perm], b[perm])
end

function invert_shape(lut::ShapeLUT, u, v)
    α = atan(v, u)
    idx = clamp(searchsortedfirst(lut.angle, α), 2, length(lut.angle))
    α0, α1 = lut.angle[idx-1], lut.angle[idx]
    t = α1 ≈ α0 ? 0.0 : clamp((α - α0) / (α1 - α0), 0.0, 1.0)
    θ = lut.θgrid[idx-1] + t * (lut.θgrid[idx] - lut.θgrid[idx-1])
    a = lut.a[idx-1] + t * (lut.a[idx] - lut.a[idx-1])
    b = lut.b[idx-1] + t * (lut.b[idx] - lut.b[idx-1])
    ρ = hypot(u, v) / hypot(a, b)
    return ρ, θ
end

# ------------------------------------------------------------
# MillerGeometry: axisymmetric, D-shaped (Miller-parametrized) cross-section.
# Reduces to a plain circular torus when κ = 1, δ = 0. TCVGeometry() is a
# preset of this tuned to roughly match the TCV tokamak.
# ------------------------------------------------------------
struct MillerGeometry <: AbstractGeometry
    R0::Float64        # major radius
    a::Float64         # minor radius (edge)
    κ::Float64         # elongation
    δ::Float64         # triangularity
    q::Float64         # safety factor
    wedge_width::Float64
    lut::ShapeLUT      # shape is ζ-independent (axisymmetric), so one LUT suffices
end

function MillerGeometry(; R0::Real, a::Real, κ::Real, δ::Real, q::Real=1.0, wedge_width::Real=π / 2)
    lut = build_shape_lut(θ -> (cos(θ + δ * sin(θ)), κ * sin(θ)))
    return MillerGeometry(Float64(R0), Float64(a), Float64(κ), Float64(δ), Float64(q), Float64(wedge_width), lut)
end

"""
    TCVGeometry(; q=1.0, wedge_width=π/2)

A `MillerGeometry` preset with TCV's approximate parameters
(R0 ≈ 0.88 m, a ≈ 0.25 m, κ ≈ 1.8, δ ≈ 0.4).
"""
TCVGeometry(; q::Real=1.0, wedge_width::Real=π / 2) =
    MillerGeometry(; R0=0.88, a=0.25, κ=1.8, δ=0.4, q=q, wedge_width=wedge_width)

shape_uv(g::MillerGeometry, ρ, θ) = (ρ * cos(θ + g.δ * sin(θ)), g.κ * ρ * sin(θ))

function surface_point(g::MillerGeometry, ζ, ρ, θ)
    u, v = shape_uv(g, ρ, θ)
    R = g.R0 + u
    return (R * cos(ζ), R * sin(ζ), v)
end

invert_cross_section(g::MillerGeometry, ζ, R_cyl, z) = invert_shape(g.lut, R_cyl - g.R0, z)
field_pitch(g::MillerGeometry, ζ) = 1 / g.q
wedge(g::MillerGeometry) = g.wedge_width > 0 ? (g.wedge_width / 2, 2π - g.wedge_width / 2) : nothing
nfp(::MillerGeometry) = 1
minor_radius_bounds(g::MillerGeometry) = (0.0, g.a)

# ------------------------------------------------------------
# FourierGeometry: non-axisymmetric boundary given by a double Fourier
# harmonic series (VMEC-style), for stellarators such as W7-X.
#   R_b(ζ,θ) = Σ Rmn * cos(m θ - n·nfp·ζ)
#   Z_b(ζ,θ) = Σ Zmn * sin(m θ - n·nfp·ζ)
# Interior points (ρ ∈ [0,1]) are a linear scaling from the magnetic axis
# (the m=0 harmonics of Rmn, which are ζ-dependent but θ-independent) out to
# this boundary -- a simple, illustrative interpolation, NOT a true
# nested-flux-surface MHD equilibrium.
# ------------------------------------------------------------
struct FourierGeometry <: AbstractGeometry
    nfp::Int
    Rmn::Vector{NTuple{3,Float64}}   # (m, n, coefficient)
    Zmn::Vector{NTuple{3,Float64}}
    q::Float64
    wedge_width::Float64
    lut_cache::Dict{Int,ShapeLUT}    # lazily built per ζ-bin (shape depends on ζ here)
    nζbins::Int
end

function FourierGeometry(; nfp::Integer, Rmn, Zmn, q::Real=1.0, wedge_width::Real=0.0, nζbins::Integer=64)
    Rmn64 = [(Int(m), Int(n), Float64(c)) for (m, n, c) in Rmn]
    Zmn64 = [(Int(m), Int(n), Float64(c)) for (m, n, c) in Zmn]
    return FourierGeometry(Int(nfp), Rmn64, Zmn64, Float64(q), Float64(wedge_width), Dict{Int,ShapeLUT}(), Int(nζbins))
end

"""
    W7XGeometry(; q=1.0, wedge_width=π/2, nfp=5, nζbins=64)

A `FourierGeometry` with a small, illustrative default harmonic set that
qualitatively resembles W7-X's rotating triangular/bean cross-section over
one field period. This is a schematic stand-in, NOT a real VMEC boundary --
supply your own `Rmn`/`Zmn` via `FourierGeometry(...)` for an accurate shape.

Defaults to a wedge cut (like `TCVGeometry`) so the data cut face is visible;
pass `wedge_width=0.0` for a fully closed device instead.
"""
function W7XGeometry(; q::Real=1.0, wedge_width::Real=π / 2, nfp::Integer=5, nζbins::Integer=64)
    Rmn = [(0, 0, 5.5), (1, 0, 0.55), (1, 1, 0.09), (0, 1, -0.03)]
    Zmn = [(1, 0, 0.55), (1, 1, -0.09), (2, 1, 0.02)]
    return FourierGeometry(; nfp=nfp, Rmn=Rmn, Zmn=Zmn, q=q, wedge_width=wedge_width, nζbins=nζbins)
end

boundary_R(g::FourierGeometry, ζ, θ) = sum(c * cos(m * θ - n * g.nfp * ζ) for (m, n, c) in g.Rmn)
boundary_Z(g::FourierGeometry, ζ, θ) = sum(c * sin(m * θ - n * g.nfp * ζ) for (m, n, c) in g.Zmn)
axis_R(g::FourierGeometry, ζ) = sum(c * cos(-n * g.nfp * ζ) for (m, n, c) in g.Rmn if m == 0)

function surface_point(g::FourierGeometry, ζ, ρ, θ)
    R0 = axis_R(g, ζ)
    R = R0 + ρ * (boundary_R(g, ζ, θ) - R0)
    Z = ρ * boundary_Z(g, ζ, θ)
    return (R * cos(ζ), R * sin(ζ), Z)
end

function _lut_for_zeta(g::FourierGeometry, ζ)
    period = 2π / g.nfp
    ζmod = mod(ζ, period)
    bin = clamp(round(Int, ζmod / period * g.nζbins), 0, g.nζbins - 1)
    return get!(g.lut_cache, bin) do
        ζbin = (bin + 0.5) / g.nζbins * period
        R0 = axis_R(g, ζbin)
        build_shape_lut(θ -> (boundary_R(g, ζbin, θ) - R0, boundary_Z(g, ζbin, θ)))
    end
end

invert_cross_section(g::FourierGeometry, ζ, R_cyl, z) = invert_shape(_lut_for_zeta(g, ζ), R_cyl - axis_R(g, ζ), z)
field_pitch(g::FourierGeometry, ζ) = 1 / g.q
wedge(g::FourierGeometry) = g.wedge_width > 0 ? (g.wedge_width / 2, 2π - g.wedge_width / 2) : nothing
nfp(g::FourierGeometry) = g.nfp
minor_radius_bounds(::FourierGeometry) = (0.0, 1.0)
