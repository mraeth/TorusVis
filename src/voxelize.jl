# ------------------------------------------------------------
# Cartesian voxel grid <-> (ζ, ρ, θ) resampling, driven entirely through the
# `AbstractGeometry` interface (surface_point / invert_cross_section /
# field_pitch) -- no device-specific code below.
#
# Periodicity, deliberately asymmetric between the two directions:
#
# - Toroidally (ζ, `nz`), the box IS periodically filled: rather than tiling
#   physical copies of the small (N_ρ, N_θ, N_ζ) data array, every voxel's
#   toroidal data index is computed via `mod` against the box's own angular
#   width `box_tor_width`, so the same handful of frames repeat, at zero
#   extra memory, as many times as needed to cover the *entire* displayed
#   toroidal section (`ζ_lo` to `ζ_hi`) with no tile seams.
#
# - Poloidally (θ, `ny`), the box is NOT filled/repeated -- only the box's
#   own true angular width `box_pol_width` is populated (rejecting anything
#   outside it), so what's rendered honestly shows how large the flux tube
#   actually is relative to the device, rather than papering over the whole
#   poloidal circle with repeated copies. `field_pitch(g, ζ)` still shifts
#   this single window's centre continuously as ζ sweeps, so it spirals
#   along a magnetic field line as it advances toroidally.
#
# Only the radial (ρ, `nx`) axis is a genuine fixed, non-periodic window --
# that's the real edge-to-core extent of the box.
# ------------------------------------------------------------

"""
    bounding_box(g, ζ_ref, ζ_lo, ζ_hi, ρ1, ρ2, box_pol_width; nζ=60, nρ=8, nθ=16)

Cartesian bounding box swept by the flux tube -- `ζ ∈ [ζ_lo, ζ_hi]`,
`ρ ∈ [ρ1, ρ2]`, and only the box's own poloidal width `box_pol_width`
(centred on the field-line shift at each `ζ`, not the full poloidal circle)
-- found by direct sampling, robust regardless of the geometry's curvature.
`ζ_ref` is the toroidal angle at which the field-line shift is zero (pass the
same value used in `precompute_voxel_map`).
"""
function bounding_box(g::AbstractGeometry, ζ_ref, ζ_lo, ζ_hi, ρ1, ρ2, box_pol_width; nζ=60, nρ=8, nθ=16)
    xmin = ymin = zmin = Inf
    xmax = ymax = zmax = -Inf
    half = box_pol_width / 2
    for ζ in range(ζ_lo, ζ_hi, length=nζ), ρ in range(ρ1, ρ2, length=nρ), θoff in range(-half, half, length=nθ)
        shift = (ζ - ζ_ref) * field_pitch(g, ζ)
        θ = θoff + shift
        x, y, z = surface_point(g, ζ, ρ, θ)
        xmin, xmax = min(xmin, x), max(xmax, x)
        ymin, ymax = min(ymin, y), max(ymax, y)
        zmin, zmax = min(zmin, z), max(zmax, z)
    end
    return (xmin, xmax), (ymin, ymax), (zmin, zmax)
end

"""
    precompute_voxel_map(g, ζ_lo, ζ_hi, ζ_ref, ρ1, ρ2, N_ρ, N_θ, N_ζ,
                          box_pol_width, box_tor_width, xs, ys, zs)

For every voxel centre in the Cartesian grid `(xs, ys, zs)`, finds whether it
falls inside the flux tube -- `ζ ∈ [ζ_lo,ζ_hi]`, `ρ ∈ [ρ1,ρ2]`, and within
`box_pol_width` of the field-line-following poloidal centre at that `ζ`
(this is a real, non-repeated window: points outside it are rejected, so the
rendered volume shows the flux tube's true poloidal size rather than filling
the whole poloidal circle) -- and, if so, which data index of the source
`(N_ρ, N_θ, N_ζ)` array it samples. Toroidally, `kk` IS periodically wrapped
against `box_tor_width` so the same `N_ζ` frames repeat to cover the entire
`[ζ_lo,ζ_hi]` section. `ζ_ref` is the toroidal angle at which the field-line
spiral shift is zero (pass the same value across calls that should share one
continuous field line, and the same value passed to `bounding_box`).

Also marks a one-voxel `border` shell just outside `inside`, carrying the
same data index as its nearest `inside` neighbour but zero alpha: hardware
trilinear filtering (`interpolate=true` in `volume!`) blends every voxel with
its neighbours, and without this, blending real colour into the fully
transparent (0,0,0,0) exterior darkens the whole boundary into a grey/black
halo. Uses the full 26-neighbourhood (face+edge+corner), since trilinear
filtering blends all 8 corners of the voxel cell containing a sample point.
"""
function precompute_voxel_map(g::AbstractGeometry, ζ_lo, ζ_hi, ζ_ref, ρ1, ρ2,
    N_ρ, N_θ, N_ζ, box_pol_width, box_tor_width, xs, ys, zs)
    Nx, Ny, Nz = length(xs), length(ys), length(zs)
    inside = falses(Nx, Ny, Nz)
    ri = zeros(Int, Nx, Ny, Nz)
    tj = zeros(Int, Nx, Ny, Nz)
    kk = zeros(Int, Nx, Ny, Nz)
    half = box_pol_width / 2
    for (k, z) in enumerate(zs), (j, y) in enumerate(ys), (i, x) in enumerate(xs)
        ζ = mod2pi(atan(y, x))
        (ζ_lo <= ζ <= ζ_hi) || continue
        R_cyl = hypot(x, y)
        ρ, θ = invert_cross_section(g, ζ, R_cyl, z)
        (ρ1 <= ρ <= ρ2) || continue

        shift = (ζ - ζ_ref) * field_pitch(g, ζ)
        θrel = mod(θ - shift + π, 2π) - π       # position relative to the field-line-following window centre
        (-half <= θrel <= half) || continue      # outside the flux tube's actual poloidal extent -- reject, don't wrap

        pf = (θrel + half) / box_pol_width * N_θ
        ζf = mod(ζ - ζ_lo, box_tor_width)          # toroidally, DO wrap: fill the whole [ζ_lo,ζ_hi] section
        kf = ζf / box_tor_width * N_ζ
        rf = (ρ - ρ1) / (ρ2 - ρ1) * (N_ρ - 1)

        inside[i, j, k] = true
        ri[i, j, k] = clamp(round(Int, rf) + 1, 1, N_ρ)
        tj[i, j, k] = clamp(round(Int, pf) + 1, 1, N_θ)
        kk[i, j, k] = clamp(round(Int, kf) + 1, 1, N_ζ)
    end

    border = falses(Nx, Ny, Nz)
    neighbours = [(di, dj, dk) for di in -1:1, dj in -1:1, dk in -1:1 if !(di == 0 && dj == 0 && dk == 0)]
    for k in 1:Nz, j in 1:Ny, i in 1:Nx
        inside[i, j, k] && continue
        for (di, dj, dk) in neighbours
            ni, nj, nk = i + di, j + dj, k + dk
            if 1 <= ni <= Nx && 1 <= nj <= Ny && 1 <= nk <= Nz && inside[ni, nj, nk]
                ri[i, j, k] = ri[ni, nj, nk]
                tj[i, j, k] = tj[ni, nj, nk]
                kk[i, j, k] = kk[ni, nj, nk]
                border[i, j, k] = true
                break
            end
        end
    end
    return inside, border, ri, tj, kk
end

"""
    sample_rgba(sl, (lo, hi), inside, border, ri, tj, kk, cmap, alpha_fn) -> Array{RGBAf}

Bakes a `(N_ρ, N_θ, N_ζ)` data slab plus a shared `(lo, hi)` colorrange into a
voxel RGBA volume: colormap-looked-up + `alpha_fn`-shaped opacity inside the
window, fully transparent outside it (`border` voxels still carry real
colour -- see `precompute_voxel_map`). Consumed by `:absorptionrgba`.
"""
function sample_rgba(sl, (lo, hi), inside, border, ri, tj, kk, cmap, alpha_fn)
    n = length(cmap)
    out = Array{RGBAf}(undef, size(inside))
    @inbounds for idx in eachindex(inside)
        if inside[idx] || border[idx]
            v = sl[ri[idx], tj[idx], kk[idx]]
            frac = clamp((v - lo) / (hi - lo + eps(Float32)), 0.0f0, 1.0f0)
            c = cmap[clamp(round(Int, frac * (n - 1)) + 1, 1, n)]
            a = inside[idx] ? Float32(alpha_fn(frac)) : 0.0f0
            out[idx] = RGBAf(red(c), green(c), blue(c), a)
        else
            out[idx] = RGBAf(0.0f0, 0.0f0, 0.0f0, 0.0f0)
        end
    end
    return out
end
