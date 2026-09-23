# Volumetric render on a schematic W7-X-like stellarator, with no wedge cut
# (a fully closed device -- `wedge_width=0`): demonstrates a non-axisymmetric
# geometry whose cross-section shape itself changes with toroidal angle, and
# the closed-device (no end caps) rendering path.
#
# Exports a GIF rather than opening an interactive window -- see
# scripts/output/ for the result.
using TorusVis

datadir = get(ENV, "TORUSVIS_DATADIR", joinpath(@__DIR__, "..", ".."))
grid, ts, simtime = read_series(joinpath(datadir, "density.h5"), "dn")

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)

g = W7XGeometry(q=1.0, wedge_width=0.0)
visualize(g, grid, ts, simtime; mode=:volume, record_to=joinpath(outdir, "w7x_volume.gif"), fps=24)
