module TorusVis

using GLMakie
using HDF5
using PlasmaCore
using Colors: red, green, blue

include("geometry.jl")
include("data.jl")
include("voxelize.jl")
include("render.jl")

export AbstractGeometry, MillerGeometry, TCVGeometry, FourierGeometry, W7XGeometry
export surface_point, invert_cross_section, field_pitch, wedge, nfp, minor_radius_bounds, cut_point
export read_series, periodic_extend
export visualize, visualize_still

end # module TorusVis
