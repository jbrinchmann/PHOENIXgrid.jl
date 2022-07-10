"""Get the top-level directory with the data.
This needs to be extended to support environment variables"""
ROOTDIR() = "/data2/jarle/spexxyDATA/PhoenixGrids/"

"""The name of a FITS file with a particular Z, Teff, logg, alpha value"""
function filename(Z, Teff, logg; alpha=0.0)
    TOP = ROOTDIR()*"MuseLsfFITS/PHOENIX-ACES-AGSS-COND-2011/"

    # We need to do a small adjustment when Z == 0
    if (Z==0.0)
        Z = -0.0
    end
    subdir = @sprintf "Z%+3.1f" Z
    suffix = ""
    
    if alpha != 0.0
        suffix = @sprintf ".Alpha=%+4.2f" alpha
    end
    subdir = subdir*suffix
    subdir = subdir*"/"
    
    fname = @sprintf "lte%05d%+4.2f%+3.1f%s.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits" Teff -logg Z suffix

    # TODO: Add a check for testing for the existence

    return TOP*subdir*fname
end

"""Load the CSV file with an overview of the grid"""
load_grid_summary() = CSV.File(ROOTDIR()*"grid.csv")

"""Create a wavelength axis from a header"""
 λ_from_header(h) = [exp((i-h["CRPIX1"])*h["CDELT1"]+h["CRVAL1"]) for i in 1:h["NAXIS1"]]

"""Load a single FITS spectrum"""
function load_one_spectrum(Z, Teff, logg; alpha=0.0, verbose=false)
    fname = filename(Z, Teff, logg; alpha=alpha)

    if isfile(fname)
        fh = FITS(fname)
        header = read_header(fh[1])
        sp = read(fh[1])
        close(fh)
        λ = λ_from_header(header)
        # Note that I store the dlnλ value here
        result =  Dict(:λ => λ, :flux => sp, :dlnλ => header["CDELT1"])
    else
        if verbose
            println("I found nothing in $fname")
        end
        result = nothing
    end

    return result
end

"""Load a grid of all spectra"""
function load_grid_of_spectra(;grid=nothing)
    if grid == nothing
        grid = load_grid_summary()
    end

    # Find the unique values for the various dimensions
    loggs = sort(unique(grid.logg))
    Ng = length(loggs)
    Teffs = sort(unique(grid.Teff))
    NT = length(Teffs)
    FeHs = sort(unique(grid.FeH))
    NFeH = length(FeHs)
    alphas = sort(unique(grid.Alpha))
    Nalpha = length(alphas)

    dims = Dict(:logg => loggs, :Teff => Teffs, :FeH => FeHs,
                :logZ => FeHs, :alpha => alphas)
    
    # The spectra all seem to have the same number of pixels so
    # we can just read the first and use that to set dimensions
    # as well as the wavelength array. However as not all parameter
    # combinations might exist, we do this below in the loop.
    # However we need to define the variables as locals out here

    local spectra
    local λ
    keep = ones(Int8, Ng, NT, NFeH, Nalpha)
    

    # Loop over to read in the necessary files.
    first = true
    for ig=1:Ng
        for iT=1:NT
            for iZ=1:NFeH
                for ia=1:Nalpha
                    tmp = load_one_spectrum(FeHs[iZ], Teffs[iT], loggs[ig]; alpha=alphas[ia])
                    if tmp == nothing
                        keep[ig, iT, iZ, ia] = 0
                    else
                        if first
                            N = length(tmp[:flux])
                            spectra = zeros(N, Ng, NT, NFeH, Nalpha)
                            λ = tmp[:λ]

                            first = false
                        end
                        spectra[:, ig, iT, iZ, ia] = tmp[:flux]
                    end
                end
            end
        end
    end
    dims[:λ] = λ
    
    return λ, spectra, keep, dims
end


"""Create a mask for a spectrum. It is possible to specify regions as well as
lines with widths around to keep.

Examples: 

  regions - a vector of vectors that specifies the regions. E.g:
            regions = [[4750, 5000.0], [6300, 6310], [4370, 4400]]


"""
function make_spectrum_mask(sp; regions=[], lines=[], width=10.0)

    Nλ = length(sp[:λ])
    mask = zeros(Nλ) # To start with, all is masked.

    # We loop over all wavelengths and for each we check whether
    # we are in a region or close to a line to keep.

    for i=1:Nλ

        λ = sp[:λ][i]

        # First check the regions
        for r in regions
            if (λ >= r[1]) & (λ <= r[2])
                mask[i]=1
            end
        end

        # Then check the lines
        for l in lines
            if (λ >= l-width/2.0) & (λ <= l+width/2.0)
                mask[i]=1
            end
        end        
    end

    return mask
end


"""This function gets a mask for the spectrum and then either
    subsets the spectrum or sets the uncertainty in the masked 
    regions to very high values"""
function mask_spectrum(sp; subset=false, kwargs...)
    
    mask = make_spectrum_mask(sp; kwargs...)

    
    if subset
        inds = [i for i in 1:length(mask) if mask[i] == 1]

        spout = Dict(:λ => sp[:λ][inds],
                     :flux => sp[:flux][inds])
        if :dflux in keys(sp)
            spout[:dflux] = sp[:dflux][inds]
        end
    else
        spout = copy(sp)
        
        inds = [i for i in 1:length(mask) if mask[i] == 0]
        spout[:dflux][inds] .= 1e30
    end

    return spout
    
end
