"""Get the top-level directory with the data.
This needs to be extended to support environment variables"""
ROOTDIR() = "/data2/jarle/spexxyDATA/PhoenixGrids/"

"""The name of a FITS file with a particular Z, Teff, logg, alpha value"""
function filename(Z, Teff, logg; alpha=0.0)
    TOP = ROOTDIR()*"MuseLsfFITS/PHOENIX-ACES-AGSS-COND-2011/"
    subdir = @sprintf "Z%+3.1f" Z
    if alpha != 0.0
        suffix = @sprintf ".Alpha=%+4.2f" alpha
        subdir = subdir*suffix
    end

    fname = @sprintf "lte%05d%+4.2f%+3.1F.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits" Teff -logg Z

    # Add a check for testing for the existence

    return fname
end

"""Load the CSV file with an overview of the grid"""
load_grid_summary() = CSV.File(ROOTDIR()*"grid.csv")

"""Create a wavelength axis from a header"""
 λ_from_header(h) = [exp((i-h["CRPIX1"])*h["CDELT1"]+h["CRVAL1"]) for i in 1:h["NAXIS1"]]

"""Load a single FITS spectrum"""
function load_one_spectrum(Z, Teff, logg; alpha=0.0)
    fname = filename(Z, Teff, logg; alpha=alpha)

    fh = FITS(fname)
    header = read_header(fh[1])
    sp = read(fh[1])
    close(fh)
    λ = λ_from_header(header)
    
    return {:λ => λ, :flux => sp}
end

"""Load a grid of all spectra"""
function load_grid_of_spectra(;grid=nothing)
    if grid == nothing
        grid = load_grid_summary()
    end

    # Find the unique values for the various dimensions
    loggs = unique(t.logg)
    teffs = unique(t.Teff)
    FeHs = unique(t.FeH)
    alphas = unique(t.Alpha)


    # The spectra all seem to have the same number of pixels so
    # we can just read the first and use that to set dimensions
    # as well as the wavelength array.

    
end
