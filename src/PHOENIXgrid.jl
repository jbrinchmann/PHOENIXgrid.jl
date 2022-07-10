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



abstract type AbstractFitResult end

#
# This fit result type is specific for the fit_grid routine
# Note that I am not storing χ2 because that is recoverable from
# p and min_χ2 through χ2 = min_χ2-ln(p)*2
#
mutable struct GridFitResult <: AbstractFitResult
    p::Array{Float64, 5}
    A::Array{Float64, 5}
    min_χ2::Float64
    i_min::Vector{Int64}
    v::Vector{Float64}
    logg::Vector{Float64}
    Teff::Vector{Float64}
    FeH::Vector{Float64}
    alpha::Vector{Float64}
end


include("io.jl")
include("fit.jl")
include("mock.jl")
include("visualise.jl")

end # module
