module PHOENIXgrid

using Printf
using CSV # DataFrames
using FITSIO
using Plots
import Statistics: median
import Random: randn
import Interpolations: LinearInterpolation

export load_grid_summary, load_one_spectrum, load_grid_of_spectra,
    make_spectrum_mask, mask_spectrum,
    fit_one_spectrum, fit_grid, fill_p, create_prior_grid, marginalise,
    summarise_results, print_results, get_best_fit



include("io.jl")
include("fit.jl")
include("mock.jl")
include("visualise.jl")

end # module
