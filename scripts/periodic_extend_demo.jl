# Demonstrates `periodic_extend`, the explicit "append arrays to wrap around
# the torus" utility, and contrasts it with the mod1-index approach the
# renderer actually uses internally (see voxelize.jl's module docstring).
using TorusVis
using PlasmaCore

datadir = get(ENV, "TORUSVIS_DATADIR", joinpath(@__DIR__, "..", ".."))
grid, ts, simtime = read_series(joinpath(datadir, "density.h5"), "dn")

frame = ts[1]                      # PlasmaCore.ScalarField, size (nx, ny, nz)
println("one frame: ", size(frame))

# Literally materialize 4 periodic copies along the poloidal (2nd) axis --
# useful for inspection/analysis, but note the memory cost: this is 4x the
# data of one frame, and would be 4x *every* frame if done for a whole
# TimeSeries. That's why the renderer instead computes `mod1`-wrapped indices
# against the original small array on the fly (zero extra memory, works for
# an arbitrary number of wraps).
wrapped = periodic_extend(frame, 2, 4)
println("4 periodic copies along dim 2: ", size(wrapped))
@assert size(wrapped) == (size(frame, 1), 4 * size(frame, 2), size(frame, 3))
@assert wrapped[:, 1:size(frame, 2), :] == parent(frame)
println("periodic_extend round-trips the original data in its first block: OK")
