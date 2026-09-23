# ------------------------------------------------------------
# HDF5 -> PlasmaCore data structures.
#
# All units in this package are Larmor radius (length) and cyclotron
# frequency (time/rate), per the simulations this is built around.
# ------------------------------------------------------------

"""
    read_series(path, dsetname; dt=1.0, b0=1.0, Bdir=3, flip_radial=false)
        -> (grid::PlasmaCore.CartGrid, ts::PlasmaCore.TimeSeries, simtime::PlasmaCore.SimulationTime)

Read an HDF5 dataset shaped `(nx, ny, nz, nt)` -- radial x poloidal-in-box x
toroidal-in-box x time -- into PlasmaCore data structures:

- `ts` is a `PlasmaCore.TimeSeries` of `PlasmaCore.ScalarField`s that loads one
  `(nx,ny,nz)` timestep from disk on demand and caches it -- the whole file is
  never eagerly materialized.
- `grid` is a `PlasmaCore.Grid` describing the box's own logical
  (ρ, poloidal-window, toroidal-window) axes as a unit cube `[0,1]^3`; the
  physical ρ/θ/ζ mapping is the job of an `AbstractGeometry`, not this grid.
  `Bdir` marks which axis runs along the background field (axis 3, the
  toroidal/nz axis, by default -- consistent with these being flux-tube-like
  boxes elongated along the field).
- `simtime` is a `PlasmaCore.SimulationTime` with frame spacing `dt`
  (ion-cyclotron units, Ω_i⁻¹). `dt=1.0` is a placeholder (unitless
  frame-index spacing) for files without a recorded physical Δt -- pass the
  real value once it's known.

The radial (1st) axis is used as-is by default (`flip_radial=false`): the
simulation's own on-disk ordering already has the hot/core side at the high
radial index, which `mount_face!`/`mount_volume!` place at `ρ2`, the right
side of the flat panel -- matching the simulation's setup. Pass
`flip_radial=true` to reverse it instead, if a dataset's on-disk convention
runs the other way.
"""
function read_series(path::AbstractString, dsetname::AbstractString;
    dt::Real=1.0, b0::Real=1.0, Bdir::Integer=3, flip_radial::Bool=false)
    fid = h5open(path, "r")
    dset = fid[dsetname]
    nx, ny, nz, nt = size(dset)

    loader = t -> begin
        raw = Float32.(dset[:, :, :, t])
        PlasmaCore.ScalarField(flip_radial ? reverse(raw, dims=1) : raw)
    end
    ts = PlasmaCore.TimeSeries(loader, collect(1:nt); label=dsetname)
    finalizer(_ -> close(fid), ts)

    grid = PlasmaCore.Grid([0.0, 0.0, 0.0], [1.0, 1.0, 1.0], [nx, ny, nz], 3,
        Float64(b0), Int(Bdir); type=PlasmaCore.Cart)
    simtime = PlasmaCore.SimulationTime(Float64(dt), Float64(dt) * (nt - 1); nmax=nt)

    return grid, ts, simtime
end

"""
    periodic_extend(field::PlasmaCore.ScalarField, dim::Int, n::Int) -> PlasmaCore.ScalarField

Materialize `n` periodic copies of `field` concatenated along dimension `dim`
-- a direct, explicit "append arrays to wrap around the torus" utility for
scripting/inspection (see `periodic_extend_demo.jl`). The renderer itself does
NOT use this internally: materializing many copies of a `(64,64,16)` frame to
cover a full poloidal+toroidal sweep would be wasteful. Instead `voxelize.jl`
wraps indices with `mod1` against the original, small array -- `dim`'s
periodicity is exploited via index arithmetic, at zero extra memory cost,
regardless of how many times the box needs to repeat to fill the display.
"""
function periodic_extend(field::PlasmaCore.ScalarField, dim::Int, n::Int)
    n >= 1 || throw(ArgumentError("n must be >= 1"))
    data = parent(field)
    return PlasmaCore.ScalarField(cat(ntuple(_ -> data, n)...; dims=dim))
end
