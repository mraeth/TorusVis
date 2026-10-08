# Flat cut-face projection on a TCV-like tokamak: the density fluctuation
# data mapped onto an annular window right at the cut surface, with no
# volumetric extension into the torus body.
#
# Exports a GIF rather than opening an interactive window -- see
# _output/ for the result.

using Pkg
Pkg.activate(@__DIR__)
# Pkg.Registry.add("https://gitlab.mpcdf.mpg.de/bsl6d/BSLRegistry")
Pkg.develop(path=joinpath(@__DIR__, ".."))

using TorusVis

using Serialization

datadir = "/home/mraeth@ipp.mpg.de/misc/torussim"
simTime, grid, ts = deserialize(joinpath(datadir, "en.bslLD"))

# grid, ts, simtime = read_series(joinpath(datadir, "en.h5"), "en")

outdir = joinpath(@__DIR__, "..", "_output")
mkpath(outdir)

g = TCVGeometry()
visualize(g, grid, ts, simTime; mode=:face, record_to=joinpath(outdir, "tcv_face.gif"), fps=24)
visualize(g, grid, ts, simTime; mode=:volume, record_to=joinpath(outdir, "tcv_volume.gif"), fps=24)
