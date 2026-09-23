# ------------------------------------------------------------
# Scene assembly: flat cut-face projection, volumetric field-line-following
# render, and the top-level `visualize` entry point that picks between them.
# Everything here is generic over `AbstractGeometry` -- no device-specific
# code below.
# ------------------------------------------------------------

const BG = RGBf(0.05, 0.05, 0.08)
const ACCENT = RGBf(0.95, 0.73, 0.25)
const TORUS_ALPHA = 0.16f0
const CAP_ALPHA = 0.5f0
const DIMMED_ALPHA = 0.32f0
const SHADING_KW = (; diffuse=Vec3f(0.55), specular=Vec3f(0.35), shininess=24.0f0)

"""
    mount_face!(ax, ax_panel, g, ζ, s, ρ1, ρ2, N_ρ, N_θ, box_pol_width, field, colorrange, cmap_name)

Maps a 2D scalar slice onto an annular window right at cut face `ζ`, drawn on
the torus with a gold outline, and mirrored as a flat heatmap in `ax_panel`.
The window is centred on `θ = 0` (spanning `[-box_pol_width/2, box_pol_width/2]`)
-- this must match `precompute_voxel_map`'s convention exactly (it centres the
flux tube on the field-line shift, which is zero at `ζ = ζ_ref`), or the flat
face and the volume it's supposed to sit flush against end up rotated
relative to each other.

Since the box is poloidally periodic (see `voxelize.jl`'s module docstring),
one dimmed copy is also drawn on either side of the real window, as visual
context showing that the domain repeats -- only the centre copy (the one
mirrored in `ax_panel`) is at full opacity. Skipped if the window is wide
enough that three of it wouldn't fit around the circle without overlapping.

Returns the window's 3D and 2D corner points, used for the connector lines.
"""
function mount_face!(ax, ax_panel, g::AbstractGeometry, ζ, s, ρ1, ρ2, N_ρ, N_θ, box_pol_width,
    field, colorrange, cmap_name)
    half = box_pol_width / 2
    ρs = range(ρ1, ρ2, length=N_ρ)
    θs = range(-half, half, length=N_θ)
    ε = (ρ2 - ρ1) * 0.02

    heatmap!(ax_panel, ρs, θs, field; colormap=cmap_name, colorrange=colorrange, interpolate=true)

    can_tile = 3 * box_pol_width < 2π
    if can_tile
        field3 = lift(v -> hcat(v, v, v), field)
        cmap = to_colormap(cmap_name)
        function tint(v, (lo, hi))
            n, m = size(v)
            out = Matrix{RGBAf}(undef, n, m)
            @inbounds for j in 1:m, i in 1:n
                frac = clamp((v[i, j] - lo) / (hi - lo + eps(Float32)), 0.0f0, 1.0f0)
                c = cmap[clamp(round(Int, frac * (length(cmap) - 1)) + 1, 1, length(cmap))]
                a = (N_θ < j <= 2N_θ) ? 1.0f0 : DIMMED_ALPHA
                out[i, j] = RGBAf(c.r, c.g, c.b, a)
            end
            return out
        end
        field_tinted = lift(tint, field3, colorrange)
        θs_ext = range(-3half, 3half, length=3N_θ)
        P = [cut_point(g, ζ, ρ, θ; off=ε, s=s) for ρ in ρs, θ in θs_ext]
        surface!(ax, getindex.(P, 1), getindex.(P, 2), getindex.(P, 3);
            color=field_tinted, shading=false)
    else
        P = [cut_point(g, ζ, ρ, θ; off=ε, s=s) for ρ in ρs, θ in θs]
        surface!(ax, getindex.(P, 1), getindex.(P, 2), getindex.(P, 3);
            color=field, colormap=cmap_name, colorrange=colorrange, shading=false)
    end

    fine = range(-half, half, length=200)
    inner = [cut_point(g, ζ, ρ1, θ; off=2ε, s=s) for θ in fine]
    outer = [cut_point(g, ζ, ρ2, θ; off=2ε, s=s) for θ in fine]
    edge1 = [cut_point(g, ζ, ρ, -half; off=2ε, s=s) for ρ in (ρ1, ρ2)]
    edge2 = [cut_point(g, ζ, ρ, half; off=2ε, s=s) for ρ in (ρ1, ρ2)]
    for line in (inner, outer, edge1, edge2)
        lines!(ax, Point3f.(line); color=ACCENT, linewidth=2.5)
    end

    corners3 = [cut_point(g, ζ, ρ, θ; off=3ε, s=s)
                for (ρ, θ) in ((ρ1, -half), (ρ2, -half), (ρ2, half), (ρ1, half))]
    corners2 = [Point2f(ρ, θ) for (ρ, θ) in ((ρ1, -half), (ρ2, -half), (ρ2, half), (ρ1, half))]
    return corners3, corners2
end

"""
    mount_volume!(ax, g, ζ_lo, ζ_hi, ζ_ref, ρ1, ρ2, N_ρ, N_θ, N_ζ,
                  box_pol_width, box_tor_width, slab_obs, colorrange, cmap_name, alpha_fn;
                  vox_h=nothing, absorption=150.0f0)

Renders the data as a true 3D volume: a flux tube of the box's own true
poloidal size that extends toroidally over the entire displayed body
(`ζ_lo` to `ζ_hi`), spiralling along a field line of pitch `field_pitch(g, ζ)`.
A single Cartesian voxel grid / single
`volume!` object covers the whole thing -- toroidally, the box's periodicity
(see `voxelize.jl`) is exploited via index wraparound to fill the entire
`ζ_lo..ζ_hi` section with no tile seams; poloidally, only the box's own true
`box_pol_width` is populated (not the whole poloidal circle), so what's
rendered honestly shows how large the flux tube actually is relative to the
device.
"""
function mount_volume!(ax, g::AbstractGeometry, ζ_lo, ζ_hi, ζ_ref, ρ1, ρ2, N_ρ, N_θ, N_ζ,
    box_pol_width, box_tor_width, slab_obs, colorrange, cmap_name, alpha_fn;
    vox_h=nothing, absorption=150.0f0)
    (xmin, xmax), (ymin, ymax), (zmin, zmax) = bounding_box(g, ζ_ref, ζ_lo, ζ_hi, ρ1, ρ2, box_pol_width)
    # Default voxel size is derived from the *physical* (lab-frame) thickness
    # of the ρ-window, not the raw ρ2-ρ1 span -- for MillerGeometry ρ is
    # already physical (so this reduces to the old ρ2-ρ1 behaviour), but for
    # FourierGeometry ρ is normalized 0..1 while the device itself can be
    # physically much larger, so using the raw span there would pick a voxel
    # size many orders of magnitude too fine for the actual bounding box.
    if vox_h === nothing
        ζ_mid = (ζ_lo + ζ_hi) / 2
        p_out = surface_point(g, ζ_mid, ρ2, 0.0)
        p_in = surface_point(g, ζ_mid, ρ1, 0.0)
        radial_thickness = hypot(p_out[1] - p_in[1], p_out[2] - p_in[2], p_out[3] - p_in[3])
        h = radial_thickness / 8
        # Safety net: a very wide toroidal section or poloidal window can
        # still make the bounding box big even though it's just the flux
        # tube's own footprint. Auto-coarsen the default (only the default
        # -- an explicit `vox_h` is respected exactly) to stay under a fixed
        # voxel budget, growing `h` just enough to fit; pass `vox_h`
        # explicitly for finer detail at the cost of memory/time.
        voxel_budget = 8_000_000
        n0 = ((xmax - xmin + 4h) / h) * ((ymax - ymin + 4h) / h) * ((zmax - zmin + 4h) / h)
        if n0 > voxel_budget
            h *= (n0 / voxel_budget)^(1 / 3)
        end
    else
        h = vox_h
    end
    pad = 2h
    Nx = max(2, round(Int, (xmax - xmin + 2pad) / h) + 1)
    Ny = max(2, round(Int, (ymax - ymin + 2pad) / h) + 1)
    Nz = max(2, round(Int, (zmax - zmin + 2pad) / h) + 1)
    xs = range(xmin - pad, xmax + pad, length=Nx)
    ys = range(ymin - pad, ymax + pad, length=Ny)
    zs = range(zmin - pad, zmax + pad, length=Nz)

    inside, border, ri, tj, kk = precompute_voxel_map(g, ζ_lo, ζ_hi, ζ_ref, ρ1, ρ2,
        N_ρ, N_θ, N_ζ, box_pol_width, box_tor_width, xs, ys, zs)
    cmap = to_colormap(cmap_name)
    voxel_rgba = lift((sl, cr) -> sample_rgba(sl, cr, inside, border, ri, tj, kk, cmap, alpha_fn),
        slab_obs, colorrange)

    # interpolate=true gives smooth (trilinearly-filtered) volume rendering,
    # safe only because of the `border` shell computed above. shading=false
    # and transparency=true are both load-bearing: Volume defaults to
    # shading=true (gradient-estimated lighting, which paints a dark rim at
    # our sharp inside/outside boundary) and transparency=false (which,
    # composited against this scene's transparency=true torus shell/caps,
    # also produces a dark rim) -- both previously diagnosed and fixed, and
    # must not regress.
    volume!(ax, (first(xs), last(xs)), (first(ys), last(ys)), (first(zs), last(zs)), voxel_rgba;
        algorithm=:absorptionrgba, absorption=absorption, interpolate=true, shading=false,
        transparency=true)
    return nothing
end

"""
    _default_azimuth(mode, ζ_A, ζ_B)

Camera azimuth default, found empirically -- `Axis3`'s `viewmode=:stretch`
combined with `aspect=:data` badly distorts the naive geometric notion of
"facing" a plane (each axis gets independently rescaled to its own data
range before the camera looks at it), so this can't be derived analytically;
an azimuth sweep was rendered and scored by the projected area of the
mounted window's gold outline to find it.

Same angle for both modes -- `ζ_B`, the far/blank end of the wedge -- so a
`:face` and a `:volume` render of the same geometry line up under the same
camera.
"""
_default_azimuth(mode::Symbol, ζ_A, ζ_B) = ζ_B

function _resolve_ρ_window(g::AbstractGeometry, ρ_window)
    if ρ_window === nothing
        ρlo, ρhi = minor_radius_bounds(g)
        span = ρhi - ρlo
        return (ρhi - span / 8, ρhi)
    else
        return ρ_window
    end
end

"""
    _draw_shell_and_caps!(ax, g, ζ_A, ζ_B, ζbounds)

Draws the semi-transparent torus shell (full poloidal loop, `ζ_A` to `ζ_B`)
and, if the device isn't fully closed (`ζbounds !== nothing`), the two
semi-transparent end caps at `ζ_A`/`ζ_B`. Shared by `visualize` and
`visualize_still` so the two can't drift out of sync.
"""
function _draw_shell_and_caps!(ax, g::AbstractGeometry, ζ_A, ζ_B, ζbounds)
    ρ_edge = minor_radius_bounds(g)[2]
    θt = range(0, 2π, length=80)
    ζt = range(ζ_A, ζ_B, length=200)
    Pshell = [surface_point(g, ζ, ρ_edge, θ) for ζ in ζt, θ in θt]
    shell_base = Makie.to_color(:steelblue3)
    shell_color = RGBAf(red(shell_base), green(shell_base), blue(shell_base), TORUS_ALPHA)
    surface!(ax, getindex.(Pshell, 1), getindex.(Pshell, 2), getindex.(Pshell, 3);
        color=fill(shell_color, size(Pshell)), transparency=true, SHADING_KW...)

    if ζbounds !== nothing
        ρlo, ρhi = minor_radius_bounds(g)
        ρd = range(ρlo, ρhi, length=6)
        αd = range(0, 2π, length=120)
        cap_base = Makie.to_color(:gray40)
        cap_color = RGBAf(red(cap_base), green(cap_base), blue(cap_base), CAP_ALPHA)
        for ζ in (ζ_A, ζ_B)
            P = [surface_point(g, ζ, ρ, α) for ρ in ρd, α in αd]
            surface!(ax, getindex.(P, 1), getindex.(P, 2), getindex.(P, 3);
                color=fill(cap_color, size(P)), transparency=true, SHADING_KW...)
        end
    end
    return nothing
end

"""
    visualize(g, grid, ts, simtime; mode=:volume, kwargs...) -> Figure

Top-level entry point. `mode = :face` draws only the flat annular projection
on the geometry's data cut face (`wedge(g)[1]`, or ζ=0 for a closed device).
`mode = :volume` draws that same face plus a single volumetric render of the
flux tube extending toroidally to fill the *entire* displayed body, spiralling
along a field line as it goes. Poloidally the volume stays the box's own true
size (not stretched to fill the poloidal circle) -- the point is to show how
large the flux tube actually is relative to the device, not to paper over the
whole cross-section.

Keywords:
- `ρ_window`: `(ρ1, ρ2)` minor-radius window to display; defaults to the
  outer 1/8 of `minor_radius_bounds(g)`.
- `box_pol_width`: angular width (rad) the box's `ny` axis physically spans
  -- this is the flux tube's true poloidal size, NOT tiled/repeated.
- `box_tor_width`: angular width (rad) the box's `nz` axis physically spans
  -- this IS tiled/repeated (periodically, via index wraparound) to cover the
  whole toroidal section. Tune both to your simulation's actual
  poloidal/parallel extent; defaults reproduce the original TCV script's scale.
- `colormap`, `alpha_fn`: colour/opacity mapping for the volume render.
- `vox_h`: voxel size for the volume grid; defaults to `span(ρ_window)/8`.
- `fps`, `azimuth_speed`, `record_to`: as before -- `record_to` must end in
  `.mp4` or `.gif`; `nothing` (default) opens an interactive window.
- `azimuth`, `elevation`: camera angles; see `_default_azimuth` for what
  `azimuth=nothing` picks per mode and why.
"""
function visualize(g::AbstractGeometry, grid, ts, simtime;
    mode::Symbol=:volume,
    ρ_window=nothing,
    box_pol_width::Real=2π * 20 / 750,
    box_tor_width::Real=π / 2,
    colormap=:balance,
    alpha_fn=frac -> abs(2frac - 1),
    vox_h=nothing,
    absorption=150.0f0,
    fps=30, azimuth_speed=0.0, record_to=nothing,
    azimuth=nothing, elevation=0.2)

    mode in (:face, :volume) || error("mode must be :face or :volume, got $(repr(mode))")
    if record_to !== nothing && !(endswith(record_to, ".mp4") || endswith(record_to, ".gif"))
        error("record_to must end in \".mp4\" or \".gif\", got \"$record_to\"")
    end

    N_ρ, N_θ, N_ζ = size(ts[1])
    nt = length(ts)
    ρ1, ρ2 = _resolve_ρ_window(g, ρ_window)

    ζbounds = wedge(g)
    ζ_A, ζ_B = ζbounds === nothing ? (0.0, 2π) : ζbounds

    sim = Observable(ts[1])
    clim = Observable(Float32(max(maximum(abs, sim[]), 1.0f-6)))
    colorrange = lift(c -> (-c, c), clim)
    face_density = lift(sl -> sl[:, :, 1], sim)

    frame_label = Observable("t = 1 / $nt")
    fig = Figure(size=(1700, 1000), backgroundcolor=BG, figure_padding=8)
    colgap!(fig.layout, 8)
    rowgap!(fig.layout, 8)
    Label(fig[0, 1:2], frame_label; color=:white, fontsize=20, font=:bold)

    az = azimuth === nothing ? _default_azimuth(mode, ζ_A, ζ_B) : azimuth
    ax = Axis3(fig[1, 1]; aspect=:data, azimuth=az, elevation=elevation, perspectiveness=0.3,
        viewmode=:stretch, protrusions=0, backgroundcolor=BG)
    hidedecorations!(ax)
    hidespines!(ax)
    colsize!(fig.layout, 1, Relative(0.72))

    ax2 = Axis(fig[1, 2]; title="$(ts.label) (flat)", titlecolor=:white,
        backgroundcolor=BG, aspect=AxisAspect(1))
    hidedecorations!(ax2)
    hidespines!(ax2)
    colsize!(fig.layout, 2, Relative(0.28))

    _draw_shell_and_caps!(ax, g, ζ_A, ζ_B, ζbounds)

    cornersA3, cornersA2 = mount_face!(ax, ax2, g, ζ_A, -1, ρ1, ρ2, N_ρ, N_θ, box_pol_width,
        face_density, colorrange, colormap)

    if mode == :volume
        mount_volume!(ax, g, ζ_A, ζ_B, ζ_A, ρ1, ρ2, N_ρ, N_θ, N_ζ, box_pol_width, box_tor_width,
            sim, colorrange, colormap, alpha_fn; vox_h=vox_h, absorption=absorption)
    end

    fig_px(scene, p) = Point2f(Makie.project(scene, p)) .+ Point2f(minimum(viewport(scene)[]))
    function connector_points()
        pts = Point2f[]
        for (c3, c2) in zip(cornersA3, cornersA2)
            push!(pts, fig_px(ax.scene, Point3f(c3)), fig_px(ax2.scene, c2))
        end
        return pts
    end
    connector = Observable(Point2f[])
    linesegments!(fig.scene, connector; space=:pixel, color=(ACCENT, 0.55), linewidth=1.5, linestyle=:dash)

    display(fig)
    sleep(0.1)
    connector[] = connector_points()

    function update!(t)
        f = ts[t]
        sim[] = f
        clim[] = Float32(max(maximum(abs, f), 1.0f-6))
        frame_label[] = "t = $t / $nt"
        ax.azimuth[] += azimuth_speed
        connector[] = connector_points()
    end

    if record_to !== nothing
        record(fig, record_to, 1:nt; framerate=fps) do t
            update!(t)
        end
        return fig
    end

    t = 1
    while events(fig).window_open[]
        update!(t)
        t = mod1(t + 1, nt)
        sleep(1 / fps)
    end
    return fig
end

"""
    visualize_still(g, grid, field::PlasmaCore.ScalarField;
                     mode=:volume, label="field", kwargs...) -> Figure

Renders a single static frame from one `(N_ρ, N_θ, N_ζ)` `ScalarField` --
no `TimeSeries`, no animation, no interactive loop, no `display`. Useful for
a quick preview or a reproducible still image:

```julia
fig = visualize_still(TCVGeometry(), grid, ts[100]; mode=:volume)
save("preview.png", fig)
```

Accepts the same geometry/appearance keywords as `visualize` (`ρ_window`,
`box_pol_width`, `box_tor_width`, `colormap`, `alpha_fn`, `vox_h`,
`absorption`, `azimuth`, `elevation`) minus the ones that only make sense for
a time series (`fps`, `azimuth_speed`, `record_to`). `label` sets the flat
panel's title (in `visualize` this comes from the `TimeSeries`'s own label).
"""
function visualize_still(g::AbstractGeometry, grid, field::PlasmaCore.ScalarField;
    mode::Symbol=:volume,
    ρ_window=nothing,
    box_pol_width::Real=2π * 20 / 750,
    box_tor_width::Real=π / 2,
    colormap=:balance,
    alpha_fn=frac -> abs(2frac - 1),
    vox_h=nothing,
    absorption=150.0f0,
    azimuth=nothing, elevation=0.2,
    label="field")

    mode in (:face, :volume) || error("mode must be :face or :volume, got $(repr(mode))")

    N_ρ, N_θ, N_ζ = size(field)
    ρ1, ρ2 = _resolve_ρ_window(g, ρ_window)

    ζbounds = wedge(g)
    ζ_A, ζ_B = ζbounds === nothing ? (0.0, 2π) : ζbounds

    clim = Float32(max(maximum(abs, field), 1.0f-6))
    colorrange = Observable((-clim, clim))
    slab = Observable(field)
    face_density = lift(sl -> sl[:, :, 1], slab)

    fig = Figure(size=(1700, 1000), backgroundcolor=BG, figure_padding=8)
    colgap!(fig.layout, 8)
    rowgap!(fig.layout, 8)

    az = azimuth === nothing ? _default_azimuth(mode, ζ_A, ζ_B) : azimuth
    ax = Axis3(fig[1, 1]; aspect=:data, azimuth=az, elevation=elevation, perspectiveness=0.3,
        viewmode=:stretch, protrusions=0, backgroundcolor=BG)
    hidedecorations!(ax)
    hidespines!(ax)
    colsize!(fig.layout, 1, Relative(0.72))

    ax2 = Axis(fig[1, 2]; title="$label (flat)", titlecolor=:white,
        backgroundcolor=BG, aspect=AxisAspect(1))
    hidedecorations!(ax2)
    hidespines!(ax2)
    colsize!(fig.layout, 2, Relative(0.28))

    _draw_shell_and_caps!(ax, g, ζ_A, ζ_B, ζbounds)

    cornersA3, cornersA2 = mount_face!(ax, ax2, g, ζ_A, -1, ρ1, ρ2, N_ρ, N_θ, box_pol_width,
        face_density, colorrange, colormap)

    if mode == :volume
        mount_volume!(ax, g, ζ_A, ζ_B, ζ_A, ρ1, ρ2, N_ρ, N_θ, N_ζ, box_pol_width, box_tor_width,
            slab, colorrange, colormap, alpha_fn; vox_h=vox_h, absorption=absorption)
    end

    fig_px(scene, p) = Point2f(Makie.project(scene, p)) .+ Point2f(minimum(viewport(scene)[]))
    connector = Observable(Point2f[])
    linesegments!(fig.scene, connector; space=:pixel, color=(ACCENT, 0.55), linewidth=1.5, linestyle=:dash)

    # As in `visualize`: the layout must settle (sizes/cameras resolved)
    # before `Makie.project` gives correct pixel positions for the
    # connector lines.
    display(fig)
    sleep(0.1)
    pts = Point2f[]
    for (c3, c2) in zip(cornersA3, cornersA2)
        push!(pts, fig_px(ax.scene, Point3f(c3)), fig_px(ax2.scene, c2))
    end
    connector[] = pts

    return fig
end
