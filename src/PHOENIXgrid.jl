module PHOENIXgrid

using Printf
using CSV # DataFrames
using FITSIO
using Plots
import Statistics: median
import Random: randn
import Interpolations: LinearInterpolation

include("io.jl")
include("fit.jl")
include("mock.jl")

end # module
