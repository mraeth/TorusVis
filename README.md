# TorusVis

A Julia package for visualizing turbulence-simulation data on toroidal
devices. It maps a small, periodic simulation box (as produced by a
flux-tube gyrokinetic/turbulence code) onto a 3D device geometry, either as a
flat 2D projection on a single cut face, or as a true volumetric render of
the flux tube spiralling along a magnetic field line through the device
body.

All units are Larmor radius (length) and cyclotron frequency (time/rate).

## Capabilities

- **Two render modes**, picked via one function: `mode = :face` (a flat
  annular projection at the cut face) or `mode = :volume` (that same face
  plus a true 3D volumetric render).
- **Pluggable device geometry.** Any device is described by six functions
  (`surface_point`, `invert_cross_section`, `field_pitch`, `wedge`, `nfp`,
  `minor_radius_bounds`) — nothing in the data or rendering code is
  device-specific. Two examples ship built in:
  - `TCVGeometry` — an axisymmetric, D-shaped (Miller-parametrized)
    cross-section, tuned to roughly match the TCV tokamak.
  - `W7XGeometry` — a non-axisymmetric cross-section given by a double
    Fourier harmonic series (VMEC-style), illustrative of a W7-X-like
    stellarator. **The default coefficients are schematic, not a real VMEC
    boundary** — supply your own via `FourierGeometry(...)` for an accurate
    shape.
  You can also construct either with arbitrary parameters
  (`MillerGeometry(; R0, a, κ, δ, q, wedge_width)` /
  `FourierGeometry(; nfp, Rmn, Zmn, q, wedge_width)`), or implement the
  six-function interface for an entirely different device.
- **PlasmaCore.jl-based data structures.** `read_series` loads an HDF5
  dataset lazily (one timestep at a time, not eagerly) into a
  `PlasmaCore.TimeSeries` of `PlasmaCore.ScalarField`s, alongside a
  `PlasmaCore.Grid` and a `PlasmaCore.SimulationTime`.
- **Honest flux-tube scale.** The simulation box's poloidal (`ny`) and
  toroidal (`nz`) axes are both periodic, but the renderer treats them
  differently on purpose:
  - **Toroidally**, the box's data *is* tiled (via index wraparound, not
    physical array copies) to fill the *entire* displayed toroidal section —
    the point is to show the flux tube following a field line all the way
    around the device.
  - **Poloidally**, the box's data is *not* tiled — only its own true
    angular width is populated. This is deliberate: it shows how large the
    flux tube actually is relative to the device, rather than smearing it
    across the whole poloidal circle.
  - `field_pitch(g, ζ)` (the reciprocal safety factor `1/q`) continuously
    shifts the poloidal window as ζ advances, so the flux tube spirals along
    a real field line rather than sitting at a fixed poloidal angle.
  - In `mode=:face`, the flat annulus is drawn with one dimmed copy on
    either side of the real (full-opacity) window, as visual context that
    the domain repeats poloidally — only the centre copy is the one
    mirrored in the flat 2D panel.
- **MP4/GIF export** via `record_to`, or an interactive window when it's
  left as `nothing`.
- **Still images from a single `ScalarField`.** `visualize_still` renders one
  frame (no `TimeSeries`, no animation) — pass any single timestep (e.g.
  `ts[100]`) or a hand-built `ScalarField` and `save` the returned `Figure`.

## Package layout

```
TorusVis/
  Project.toml, Manifest.toml
  src/
    TorusVis.jl     -- module: includes + exports
    geometry.jl      -- AbstractGeometry interface, MillerGeometry/TCVGeometry,
                        FourierGeometry/W7XGeometry, ShapeLUT inversion machinery
    data.jl          -- read_series (PlasmaCore-based HDF5 loader), periodic_extend
    voxelize.jl      -- Cartesian voxel grid <-> (ζ,ρ,θ) resampling, geometry-driven
    render.jl        -- mount_face!, mount_volume!, visualize(...)
  scripts/
    tcv_face.jl, tcv_volume.jl, w7x_face.jl, w7x_volume.jl, periodic_extend_demo.jl
```

## Interfaces

### Geometry

```julia
abstract type AbstractGeometry end

surface_point(g, ζ, ρ, θ)        -> (x, y, z)   # lab-frame Cartesian point
invert_cross_section(g, ζ, R_cyl, z) -> (ρ, θ)   # inverse, at cylindrical radius R_cyl & height z
field_pitch(g, ζ)                -> Float64      # dθ/dζ along a field line (1/q)
wedge(g)                          -> (ζ_A, ζ_B) or nothing   # cut bounds, or nothing = closed device
nfp(g)                            -> Int          # toroidal field periods
minor_radius_bounds(g)           -> (ρ_lo, ρ_hi) # full ρ range (metres for Miller, 0..1 for Fourier)
```

Coordinate convention: `ζ` = toroidal angle (always the lab-frame azimuth,
i.e. `atan(y,x) == ζ`), `θ` = poloidal angle, `ρ` = minor-radius-like label
(`ρ=0` at the magnetic axis).

Concrete geometries:

```julia
TCVGeometry(; q=1.0, wedge_width=π/2)
MillerGeometry(; R0, a, κ, δ, q=1.0, wedge_width=π/2)   # general axisymmetric D-shape

W7XGeometry(; q=1.0, wedge_width=π/2, nfp=5, nζbins=64)
FourierGeometry(; nfp, Rmn, Zmn, q=1.0, wedge_width=0.0, nζbins=64)   # general non-axisymmetric
```

`Rmn`/`Zmn` are vectors of `(m, n, coefficient)` triples for
`R(ζ,θ) = Σ Rmn·cos(mθ - n·nfp·ζ)`, `Z(ζ,θ) = Σ Zmn·sin(mθ - n·nfp·ζ)`.
Both `TCVGeometry` and `W7XGeometry` default to a wedge cut (so the data cut
face is visible); pass `wedge_width=0.0` for a fully closed device instead
(no end caps drawn) — `scripts/w7x_volume.jl` demonstrates this case.

### Data

```julia
read_series(path, dsetname; dt=1.0, b0=1.0, Bdir=3, flip_radial=false) ->
    (grid::PlasmaCore.CartGrid, ts::PlasmaCore.TimeSeries, simtime::PlasmaCore.SimulationTime)
```

Reads an HDF5 dataset shaped `(nx, ny, nz, nt)` = (radial, poloidal-in-box,
toroidal-in-box, time). `dt` is the frame spacing in cyclotron units
(`Ω_i⁻¹`); it defaults to `1.0` (unitless frame-index spacing) since these
files don't record a physical Δt — pass the real value once you have it.
`ts[t]` lazily loads and caches one timestep as a `PlasmaCore.ScalarField`;
the whole file is never eagerly materialized.

```julia
periodic_extend(field::PlasmaCore.ScalarField, dim::Int, n::Int) -> PlasmaCore.ScalarField
```

Materializes `n` periodic copies of one frame along `dim` — an explicit,
direct "append arrays to wrap around" utility for scripting/inspection (see
`scripts/periodic_extend_demo.jl`). The renderer itself does not use this;
it wraps indices with `mod` against the original array instead, at zero
extra memory.

### Rendering

```julia
visualize(g::AbstractGeometry, grid, ts, simtime;
    mode=:volume,                    # :face or :volume
    ρ_window=nothing,                # (ρ1,ρ2); defaults to outer 1/8 of minor_radius_bounds(g)
    box_pol_width=2π*20/750,         # true angular width (rad) of the box's ny axis -- NOT tiled
    box_tor_width=π/2,               # angular width (rad) of the box's nz axis -- IS tiled to fill the display
    colormap=:balance,
    alpha_fn=frac -> abs(2frac - 1), # opacity shaping across the colorrange
    vox_h=nothing,                   # voxel size; auto-derived + budget-capped if not given
    absorption=150.0f0,              # :absorptionrgba glow intensity
    fps=30, azimuth_speed=0.0,
    record_to=nothing,               # nothing = interactive window; else path ending in .mp4/.gif
    azimuth=nothing, elevation=0.2,  # camera; azimuth defaults to the same angle for both modes
) -> Figure
```

`box_pol_width`/`box_tor_width` should match your simulation's actual
poloidal/parallel physical extent — the defaults just reproduce the scale of
the original TCV script this package grew out of.

For a single static frame instead of a time series, use `visualize_still`:

```julia
visualize_still(g::AbstractGeometry, grid, field::PlasmaCore.ScalarField;
    mode=:volume, label="field",     # label sets the flat panel's title
    ρ_window=nothing, box_pol_width=2π*20/750, box_tor_width=π/2,
    colormap=:balance, alpha_fn=frac -> abs(2frac - 1),
    vox_h=nothing, absorption=150.0f0,
    azimuth=nothing, elevation=0.2,
) -> Figure
```

It takes the same geometry/appearance keywords as `visualize` minus the ones
that only make sense for a time series (`fps`, `azimuth_speed`,
`record_to`) — no animation, no interactive loop, just one render:

```julia
fig = visualize_still(TCVGeometry(), grid, ts[100]; mode=:volume, label="dn")
save("preview.png", fig)
```

## Installation

TorusVis and its PlasmaCore.jl dependency are distributed through the
[BSLRegistry](https://gitlab.mpcdf.mpg.de/bsl6d/BSLRegistry) (neither is in Julia's General
registry). Add the registry once per machine / cluster account:

```julia
pkg> registry add https://gitlab.mpcdf.mpg.de/bsl6d/BSLRegistry.git
```

then install it like any registered package:

```julia
pkg> add TorusVis
julia> using TorusVis
```

Alternatively, type `using TorusVis` directly in the REPL: if the package isn't installed in the
active environment yet, Julia offers to install it.

Pick up new releases with `pkg> registry up` followed by `pkg> up`. To run the showcase scripts
below, work from a clone of this repository instead (see *First-time setup*).

## Running the showcase scripts

From the `TorusVis/` directory, with `density.h5` one level up (the default;
override with the `TORUSVIS_DATADIR` environment variable):

```sh
cd TorusVis
julia --project=. scripts/tcv_face.jl        # flat annular projection, TCV -- writes a GIF
julia --project=. scripts/w7x_face.jl        # flat annular projection, schematic W7-X -- writes a GIF
julia --project=. scripts/tcv_volume.jl      # volumetric flux tube spiralling around TCV -- writes a GIF
julia --project=. scripts/w7x_volume.jl      # volumetric flux tube, closed stellarator body -- writes a GIF
julia --project=. scripts/periodic_extend_demo.jl   # periodic_extend utility demo (non-interactive)
```

All four `visualize(...)` scripts export a full-length GIF to
`scripts/output/` rather than opening an interactive window — edit their
`record_to=...` argument to change the path or switch to `.mp4`. Any script
can be flipped back to an interactive window by removing `record_to`
(passing `nothing`, the default), e.g.:

```julia
visualize(g, grid, ts, simtime; mode=:volume, record_to="tcv_volume.mp4", fps=24)
```

### First-time setup

If `TorusVis/Manifest.toml` isn't already resolved in your environment:

```sh
cd TorusVis
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

(PlasmaCore.jl comes from the BSLRegistry, so add it first; see
[Installation](#installation). To test against a local PlasmaCore.jl checkout,
`julia --project=. -e 'using Pkg; Pkg.develop(path="/path/to/PlasmaCore.jl")'`
— this only changes the gitignored `Manifest.toml`; `Pkg.free("PlasmaCore")` switches back.)

## Writing your own scene

```julia
using TorusVis

grid, ts, simtime = read_series("density.h5", "dn")

g = TCVGeometry(q=1.0)                 # or W7XGeometry(), or your own MillerGeometry/FourierGeometry
visualize(g, grid, ts, simtime; mode=:volume)
```

Swap in your own device by implementing the six `AbstractGeometry` functions
above for a new struct — everything else (data loading, voxelization,
rendering, both modes) works unchanged.
