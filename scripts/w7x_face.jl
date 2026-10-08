# Flat cut-face projection on a schematic W7-X-like stellarator: same data,
# mounted on a non-axisymmetric (Fourier-harmonic-boundary) device instead of
# TCV's Miller D-shape. The harmonic coefficients here are illustrative, not
# a real VMEC boundary -- see `FourierGeometry`'s docstring.
#
# Exports a GIF rather than opening an interactive window -- see
# _output/ for the result.
using TorusVis

datadir = get(ENV, "TORUSVIS_DATADIR", joinpath(@__DIR__, "..", ".."))
grid, ts, simtime = read_series(joinpath(datadir, "dn"))

outdir = joinpath(@__DIR__, "..", "_output")
mkpath(outdir)

g = W7XGeometry()
# FourierGeometry's ρ is normalized 0..1, not metres (see minor_radius_bounds'
# docstring), so `rho_larmor`-based auto-sizing doesn't apply here -- give the
# box's angular size explicitly instead (these match the old hardcoded default).
visualize(g, grid, ts, simtime; mode=:face, box_pol_width=2π * 20 / 750, box_tor_width=π / 2,
    record_to=joinpath(outdir, "w7x_face.gif"), fps=24)
