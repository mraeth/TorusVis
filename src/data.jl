# ------------------------------------------------------------
# Deserialized simulation output -> PlasmaCore data structures.
#
# All units in this package are Larmor radius (length) and cyclotron
# frequency (time/rate), per the simulations this is built around.
# ------------------------------------------------------------

"""
    read_series(path) -> (grid::PlasmaCore.Grid, ts::PlasmaCore.TimeSeries, simtime::PlasmaCore.SimulationTime)

Deserializes a `(simtime, grid, ts)` tuple written by the simulation /
post-processing pipeline with `Serialization.serialize`, and returns it
reordered to match `visualize(g, grid, ts, simtime; ...)`'s argument order.

Unlike the old HDF5-based loader this replaces, `grid` here is not a
placeholder unit cube: its configuration-space axes (radial,
poloidal-in-box, toroidal-in-box) carry the box's real physical extent in
Larmor radii, which `visualize`/`visualize_still` use (together with
`rho_larmor`) to size the box against the device.
"""
function read_series(path::AbstractString)
    simtime, grid, ts = deserialize(path)
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
