# Flat cut-face projection on a schematic W7-X-like stellarator: same data,
# mounted on a non-axisymmetric (Fourier-harmonic-boundary) device instead of
# TCV's Miller D-shape. The harmonic coefficients here are illustrative, not
# a real VMEC boundary -- see `FourierGeometry`'s docstring.
#
# Exports a GIF rather than opening an interactive window -- see
# scripts/output/ for the result.
using TorusVis

datadir = get(ENV, "TORUSVIS_DATADIR", joinpath(@__DIR__, "..", ".."))
grid, ts, simtime = read_series(joinpath(datadir, "density.h5"), "dn")

outdir = joinpath(@__DIR__, "output")
mkpath(outdir)

g = W7XGeometry()
visualize(g, grid, ts, simtime; mode=:face, record_to=joinpath(outdir, "w7x_face.gif"), fps=24)
