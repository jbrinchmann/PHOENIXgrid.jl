#
# Routines for fitting a spectrum against a model. This is Bayesian and gridded.
#

"""Return the scaling from the χ² calculation in the fit_one_spectrum function"""
function scale(f, fM, σ)

    Ω1 = 0.0
    Ω2 = 0.0
    for i=1:length(f)
        s2 = σ[i]^2
        Ω1 = Ω1 + f[i]*fM[i]/s2
        Ω2 = Ω2 + fM[i]^2/s2
    end

    return Ω1/Ω2
end


"""Internal function for interpolation
"""
function _fast_interpolate(λobs, λM, fM, v)

    cspeed = 299792.458
    # Then some useful variables for the interpolation below
    Nobs = length(λobs)
    ΔlogλM = 3e-5 # Fixed for PHOENIX grids
    Δλobs = λobs[2]-λobs[1]
    tmp = (1+v/cspeed)*λM[1]
    fMshifted = zeros(Nobs)

    for iλ=1:Nobs
        # Find the index in fM that is just below or equal to λobs[i]

        z = (1/ΔlogλM)*log(λobs[iλ]/tmp)
        j = floor(Int, z)+1
        frac = (λM[j]-λobs[iλ])/Δλobs
        
        # Calculate the model flux at this observed λ
        fMshifted[iλ] = fM[j]*(1-frac)+fM[j+2]*frac
    end

    return fMshifted
end

"""Fit a model to an observed spectrum
    
As input we expect an observed spectrum with three keys: 
 - :λ is expected to have the wavelength in Å
 - :flux should have the observed flux (in arbitrary units)
 - :dflux should have the uncertainty array (in the same units as :flux

The model spectrum needs to have the same columns except :dflux (it is assumed
to be noise-free). 

We then step through velocities in vmin:vstep:vmax and calculate

        χ² = ∑ (fᵢ-A fMᵢ)²/σᵢ²

where A can be found from dχ²/dA = 0.

We loop over the velocities and shift he model spectrum to that velocity and
interpolate it onto the observed wavelength vector. 

"""
function fit_one_spectrum(sp_obs, sp_M; vmin=-100, vmax=100, vstep=1.0,
                          around=nothing, balmeronly=false)

    # Prepare some variables - some are just for readability below
    velocities = collect(vmin:vstep:vmax)
    λM = sp_M[:λ]
    fM = sp_M[:flux]
    λobs = sp_obs[:λ]
    fobs = sp_obs[:flux]
    dfobs = sp_obs[:dflux]
    NM = length(fM)

    if balmeronly
        around = [4861.325, 6562.8]
    end

    local w
    if around == nothing
        # Weighting is by 1/sigma^2 
        w = 1.0./dfobs
    else
        # In this case weight will be 1/sigma^2 where you are
        # close to around[j] otherwise 0.
        w = zeros(NM)
 #       for j=1:NM
 #           if abs(
 #           
 #       end
    end
    
    
    # I could get this from PhysicalConstants.jl but as it is the only
    # one I expect I need, I won't bother. To be checked.
    cspeed = 299792.458

    χ2 = zeros(length(velocities))
    for i=1:length(velocities)
        
        # We do an inline interpolation to speed things up - note that we
        # have no boundary checks - this could be inserted higher up by
        # testing with max or min velocity

        fMshifted = _fast_interpolate(λobs, λM, fM, velocities[i])
        # Compare the two.
        A = scale(fobs, fMshifted, dfobs)

        χ2[i] = sum(((fobs.-A.*fMshifted).*w).^2)
    end
    
    return Dict(:v=>velocities, :chi2 => χ2)
end


function fit_grid(sp_obs, grid, keep, dims; vmin=-100., vmax=100.0, vstep=1.0,
                  around=nothing, balmeronly=false)

    # The dimensionality of it all
    Nλ, Ng, NT, NZ, Na = size(grid)

    N_spobs = length(sp_obs[:flux])
    
    local χ2, vaxis, p
    min_chi2 = 1e30
    first = true
    for ig=1:Ng
        for iT=1:NT
            for iZ=1:NZ
                for ia=1:Na

                    if keep[ig, iT, iZ, ia] == 0
                        continue
                    end
                    sp_M = Dict(:λ => dims[:λ], :flux => grid[:, ig, iT, iZ, ia])
                    res = fit_one_spectrum(sp_obs, sp_M; vmin=vmin, vmax=vmax,
                                           vstep=vstep, around=around, balmeronly=balmeronly)

                    if first
                        vaxis = res[:v]
                        Nv = length(vaxis)
                        χ2 = fill(-9999.0, Nv, Ng, NT, NZ, Na)
                        p = fill(0.0, Nv, Ng, NT, NZ, Na)
                        first = false
                    end

                    tmp = res[:chi2]
                    χ2[:, ig, iT, iZ, ia] = tmp

                    if minimum(tmp) < min_chi2
                        min_chi2 = minimum(tmp)
                    end
                end
            end
        end
    end


    # I took this out of the loop because we need to subtract off the
    # minimum χ2. 
    p = exp.(.-(χ2.-min_chi2)./2.0)
    
    return vaxis, χ2, p, min_chi2
end


function interp_2d(im_p, im_keep)

    # interpolation along x-axis
    Nx, Ny = size(im_p)

    imout = zeros(Nx, Ny)
    xpos = 1:Nx
    
    for j=1:Ny

        x = [i for i in 1:Nx if (im_keep[i, j] > 0)]
        y = [im_p[i, j] for i in 1:Nx if (im_keep[i, j] > 0)]

        itp = LinearInterpolation(x, y, extrapolation_bc=0)
        imout[:, j] = itp(xpos)
    end

    return imout

end

function fill_p(p, keep)
    # Fill in the missing pieces.

    Nv, Ng, NT, NZ, Na = size(p)
    pfull = zeros(Nv, Ng, NT, NZ, Na)
    for iv=1:Nv
        for iZ=1:NZ
            for ia=1:Na
                im_p = p[iv, :, :, iZ, ia]# [1,:,:,1,1]
                im_keep = keep[:, :, iZ, ia]
                imout = interp_2d(im_p, im_keep)

                pfull[iv, :, :, iZ, ia] = imout
            end
        end
    end

    return pfull
end


function _create_prior(x, settings; norm=false)

    result = zeros(size(x))
    
    k = keys(settings)
    if :max in k
        # Uniform prior
        for i=1:length(x)
            if (x[i] >= settings[:min]) & (x[i] <= settings[:max])
                result[i]=1.0
            else
                result[i] = 0.0
            end
        end
    elseif :σ in k
        μ = settings[:mean]
        σ = settings[:σ]
        z = exp.(-((x.-μ)./σ).^2)
        result = z/(sqrt(2*π)*σ)
    elseif :func in k
        f = settings[:func]
        result = f(x)
    else
        println("No recognised prior!")
    end

    if norm
        # This normalises by summing - this is not proper
        # normalisation!
        result = result/sum(result)
    end

    return result
end
   
        
        


"""create_prior_grid(p; v=, logg=, Teff=, FeH=, alpha=)

This creates a multi-dimensional array with prior information
for this object. The priors should all be provided as Dict's. 

If the dict has keys `:min` and `:max` then a Uniform prior 
is applied between the two values. If :σ, :mean are given then
a Gaussian (N(:mean, :σ^2)) is created, and if it has a key :func,
this is assumed to be a function that can be called to give
a prior value for each value in the dimensions. 

If the prior is set to `nothing`, then nothing is indeed done. 
"""
function create_prior_grid(p, vaxis, dims; v=nothing, logg=nothing,
                           Teff=nothing, FeH=nothing,
                           alpha=nothing)
    
    Nv, Ng, NT, NZ, Na = size(p)

    # The prior needs to be the same size as p
    prior = ones(size(p))

    if v != nothing
        tmp = _create_prior(vaxis, v)
        y = repeat(tmp, 1, Ng, NT, NZ, Na);
        prior = prior.*y
    end

    if logg != nothing
        tmp = _create_prior(dims[:logg], logg)
        y = repeat(tmp, 1, Nv, NT, NZ, Na);
        y = permutedims(y, [2, 1, 3, 4, 5])
        prior = prior.*y
    end

    if Teff != nothing
        tmp = _create_prior(dims[:Teff], Teff)
        y = repeat(tmp, 1, Nv, Ng, NZ, Na);
        y = permutedims(y, [3, 1, 2, 4, 5])
        prior = prior.*y
    end
    
    if FeH != nothing
        tmp = _create_prior(dims[:FeH], FeH)
        y = repeat(tmp, 1, Nv, Ng, NT, Na);
        y = permutedims(y, [4, 1, 2, 3, 5])
        prior = prior.*y
    end

    if alpha != nothing
        tmp = _create_prior(dims[:alpha], alpha)
        y = repeat(tmp, 1, Nv, Ng, NT, NZ);
        y = permutedims(y, [5, 1, 2, 3, 4])
        prior = prior.*y
    end


    return prior
end


sumsqueeze(x; dims=nothing) = dropdims(sum(x, dims=dims), dims=dims)

    
"""Marginalise the output from fit_grid

The routine creates 1D and 2D marginals from the likelihood.
If provided, a prior is applied. 
"""
function marginalise(p, dims, keep; ignore_missing=false, pfull=nothing,
                     prior=nothing)
    #
    # We want to create both 1D and 2D summaries here
    #

    if pfull == nothing
        pfull = fill_p(p, keep)
    end
    # Apply prior if needed
    if prior != nothing
        pfull = pfull*prior
    end
    
    Nv, Ng, NT, NZ, Na = size(pfull)


    #
    # PDFs to create
    #
    P_v = zeros(Nv)
    P_logg = zeros(Ng)
    P_Teff = zeros(NT)
    P_FeH = zeros(NZ)
    P_a = zeros(Na)

    P_vg = zeros(Nv, Ng)
    P_vT = zeros(Nv, NT)
    P_vZ = zeros(Nv, NZ)
    P_va = zeros(Nv, Na)
    P_gT = zeros(Ng, NT)
    P_gZ = zeros(Ng, NZ)
    P_ga = zeros(Ng, Na)
    P_TZ = zeros(NT, NZ)
    P_Ta = zeros(NT, Na)
    P_Za = zeros(NZ, Na)
    
    
    # To speed this up and avoid multiple passes through the arrays, this
    # is done in a single loop. This could of course be shortened if I used
    # a more sophisticated data structure but this is clear and readable so
    # I stick with this.
    for iv=1:Nv
        for ig=1:Ng
            for iT=1:NT
                for iZ=1:NZ
                    for ia=1:Na
                        x = pfull[iv, ig, iT, iZ, ia]
                        if isfinite(x)
                            P_v[iv] = P_v[iv] + x
                            P_logg[ig] = P_logg[ig] + x
                            P_Teff[iT] = P_Teff[iT] + x
                            P_FeH[iZ] = P_FeH[iZ] + x
                            P_a[ia] = P_a[ia] + x
                            
                            P_vg[iv, ig] = P_vg[iv, ig] + x
                            P_vT[iv, iT] = P_vT[iv, iT] + x
                            P_vZ[iv, iZ] = P_vZ[iv, iZ] + x
                            P_va[iv, ia] = P_va[iv, ia] + x
                            
                            P_gT[ig, iT] = P_gT[ig, iT] + x
                            P_gZ[ig, iZ] = P_gZ[ig, iZ] + x
                            P_ga[ig, ia] = P_ga[ig, ia] + x
                            
                            P_TZ[iT, iZ] = P_TZ[iT, iZ] + x
                            P_Ta[iT, ia] = P_Ta[iT, ia] + x
                            
                            P_Za[iZ, ia] = P_Za[iZ, ia] + x
                        end
                    end
                end
            end
        end
    end


    # Finally, insert this into a dict and normalise everything
    res = Dict(:P_v => P_v, :P_logg => P_logg, :P_Teff => P_Teff,
               :P_FeH => P_FeH, :P_Z => P_FeH, :P_a => P_a,
               :P_vg => P_vg, :P_vT => P_vT, :P_vZ => P_vZ,
               :P_va => P_va, :P_gT => P_gT, :P_gZ => P_gZ,
               :P_ga => P_ga, :P_TZ => P_TZ, :P_Ta => P_Ta,
               :P_Za => P_Za)

    for k in keys(res)
        y = res[k]
        res[k] = y/sum(y)
    end

    return res
end


