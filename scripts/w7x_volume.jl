# Volumetric render on a schematic W7-X-like stellarator, with no wedge cut
# (a fully closed device -- `wedge_width=0`): demonstrates a non-axisymmetric
# geometry whose cross-section shape itself changes with toroidal angle, and
# the closed-device (no end caps) rendering path.
#
# Exports a GIF rather than opening an interactive window -- see
# _data/ for the result.
using TorusVis

datadir = get(ENV, "TORUSVIS_DATADIR", joinpath(@__DIR__, "..", ".."))
grid, ts, simtime = read_series(joinpath(datadir, "dn"))

outdir = joinpath(@__DIR__, "..", "_data")
mkpath(outdir)

g = W7XGeometry(q=1.0, wedge_width=0.0)
# FourierGeometry's ρ is normalized 0..1, not metres (see minor_radius_bounds'
# docstring), so `rho_larmor`-based auto-sizing doesn't apply here -- give the
# box's angular size explicitly instead (these match the old hardcoded default).
visualize(g, grid, ts, simtime; mode=:volume, box_pol_width=2π * 20 / 750, box_tor_width=π / 2,
    record_to=joinpath(outdir, "w7x_volume.gif"), fps=24)
