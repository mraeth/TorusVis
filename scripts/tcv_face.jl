# Flat cut-face projection on a TCV-like tokamak: the density fluctuation
# data mapped onto an annular window right at the cut surface, with no
# volumetric extension into the torus body.
#
# Exports a GIF rather than opening an interactive window -- see
# scripts/output/ for the result.
using TorusVis

datadir = get(ENV, "TORUSVIS_DATADIR", joinpath(@__DIR__, "..", ".."))
grid, ts, simtime = read_series(joinpath(datadir, "en.h5"), "en")

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)

g = TCVGeometry()
visualize(g, grid, ts, simtime; mode=:face, record_to=joinpath(outdir, "tcv_face.gif"), fps=24)
