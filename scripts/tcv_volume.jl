# Volumetric render on a TCV-like tokamak: the density data extends off the
# cut face into a full 3D volume that spirals along a q=1 field line,
# wrapping all the way around the displayed torus body toroidally (the
# box's nz axis is periodic, so the same handful of frames repeat
# continuously to fill the whole display), while poloidally the flux tube
# keeps its own true (small) size.
#
# Exports a GIF rather than opening an interactive window -- see
# scripts/output/ for the result.
using TorusVis

datadir = get(ENV, "TORUSVIS_DATADIR", joinpath(@__DIR__, "..", ".."))
grid, ts, simtime = read_series(joinpath(datadir, "density.h5"), "dn")

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)

g = TCVGeometry(q=1.0)
visualize(g, grid, ts, simtime; mode=:volume, record_to=joinpath(outdir, "tcv_volume.gif"), fps=24)
