#
# Routines for fitting a spectrum against a model. This is Bayesian and gridded.
#

global const cspeed = 299792.458



"""Return the scaling from the χ² calculation in the fit_one_spectrum function"""
function scale(f, fM, σ)

    Ω1 = 0.0
    Ω2 = 0.0
    for i=1:length(f)
        if (isfinite(f[i]))
            # Ignore NaNs
            s2 = σ[i]^2
            Ω1 = Ω1 + f[i]*fM[i]/s2
            Ω2 = Ω2 + fM[i]^2/s2
        end
    end
    return Ω1/Ω2
end


"""Internal function for interpolation

This relies on the fact that the PHOENIX grid is equally sampled in log λ. 
There are no assumptions on the sampling of the observed spectrum.

"""
function _fast_interpolate(λobs, λM, fM, v)

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

function ensure_finite_spectrum(sp)

    x = sp[:λ]
    y = sp[:flux]
    dy = sp[:dflux]
    N = length(x)

    xout = [x[i] for i=1:N if (isfinite(y[i])) & (dy[i] > 0)]
    yout = [y[i] for i=1:N if (isfinite(y[i])) & (dy[i] > 0)]
    dyout = [dy[i] for i=1:N if (isfinite(y[i])) & (dy[i] > 0)]
    

    spout = Dict(:λ => xout, :flux => yout, :dflux => dyout)
    return spout
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
                          around=nothing, balmeronly=false, dataisfinite=false)

    # Prepare some variables - some are just for readability below
    velocities = collect(vmin:vstep:vmax)
    λM = sp_M[:λ]
    fM = sp_M[:flux]
    NM = length(fM)

    # The observed spectrum. 
    # Keep only good data unless the caller ensures us that the
    # data are finite.
    if !dataisfinite
        sp_obs = ensure_finite_spectrum(sp_obs)
    end
    λobs = sp_obs[:λ]
    fobs = sp_obs[:flux]
    dfobs = sp_obs[:dflux]
    Nobs = length(fobs)

        

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

    χ2 = zeros(length(velocities))
    Scales = zeros(length(velocities))
    for i=1:length(velocities)
        
        # We do an inline interpolation to speed things up - note that we
        # have no boundary checks - this could be inserted higher up by
        # testing with max or min velocity

        fMshifted = _fast_interpolate(λobs, λM, fM, velocities[i])
        # Compare the two.
        A = scale(fobs, fMshifted, dfobs)
        Scales[i] = A
        χ2[i] = sum(((fobs.-A.*fMshifted).*w).^2)
    end
    
    return Dict(:v=>velocities, :chi2 => χ2, :A => Scales)
end


#
"""fit_grid(sp_obs, grid, keep, dims; kwargs)

Fit a grid of PHOENIX models to the observed spectrum in sp_obs.
The spectrum must be a dict with keys :λ, :flux, :dflux.

The 
"""
function fit_grid(sp_obs, grid, keep, dims; vmin=-100., vmax=100.0, vstep=1.0,
                  around=nothing, balmeronly=false)

    # The dimensionality of it all
    Nλ, Ng, NT, NZ, Na = size(grid)

    # We do not want any NaNs
    sp_obs = ensure_finite_spectrum(sp_obs)
    
    N_spobs = length(sp_obs[:flux])
    
    local χ2, vaxis, p, A
    min_chi2 = 1e30
    # The index of the minimum model
    i_min = [-1, -1, -1, -1, -1]
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
                                           vstep=vstep, around=around, balmeronly=balmeronly,
                                           dataisfinite=true)

                    if first
                        vaxis = res[:v]
                        Nv = length(vaxis)
                        χ2 = fill(-9999.0, Nv, Ng, NT, NZ, Na)
                        p = fill(0.0, Nv, Ng, NT, NZ, Na)
                        A = fill(0.0, Nv, Ng, NT, NZ, Na)
                        first = false
                    end

                    tmp = res[:chi2]
                    χ2[:, ig, iT, iZ, ia] = tmp
                    A[:, ig, iT, iZ, ia] = res[:A]
                    
                    ok = [i for i in 1:length(tmp) if tmp[i] > 0]
                    if length(ok) > 0
                        i_tmp = argmin(tmp[ok])
                        i_min_v = ok[i_tmp]
                        if tmp[i_min_v] < min_chi2
                            min_chi2 = tmp[i_min_v]
                            i_min = [i_min_v, ig, iT, iZ, ia]
                        end
                    end
                end
            end
        end
    end


    # I took this out of the loop because we need to subtract off the
    # minimum χ2. 
    p = exp.(.-(χ2.-min_chi2)./2.0)

    result = GridFitResult(p, A, min_chi2, i_min, vaxis, dims[:logg], dims[:Teff],
                           dims[:FeH], dims[:alpha])

    return result
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
        y = permutedims(y, [2, 3, 1, 4, 5])
        prior = prior.*y
    end
    
    if FeH != nothing
        tmp = _create_prior(dims[:FeH], FeH)
        y = repeat(tmp, 1, Nv, Ng, NT, Na);
        y = permutedims(y, [2, 3, 4, 1, 5])
        prior = prior.*y
    end

    if alpha != nothing
        tmp = _create_prior(dims[:alpha], alpha)
        y = repeat(tmp, 1, Nv, Ng, NT, NZ);
#        println("Prior dimension=",size(prior))
#        println("    y dimension=",size(y))
#        println("  tmp dimension=",size(tmp))
        y = permutedims(y, [2, 3, 4, 5, 1])
#        println("    y dimension=",size(y))

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
        pfull = pfull.*prior
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


"""Convenience function that maps between strings and symbols for the 
    marginalisation results"""
function marginalisation_keymap()
    k = Dict("v" => :P_v, "logg" => :P_logg, "Teff" => :P_Teff,
             "FeH" => :P_FeH, "Z" => :P_FeH, "alpha" => :P_a,
             "vlogg" => :P_vg, "vTeff" => :P_vT, "vFeH" => :P_vZ,
             "vZ" => :P_vZ, "valpha" => :P_va, "loggTeff" => :P_gT,
             "loggZ" => :P_gZ, "loggFeH" => :P_gZ, "loggalpha" => :P_ga,
             "TeffZ" => :P_TZ, "TeffFeH" => :P_TZ, "Teffalpha" => :P_Ta,
             "FeHalpha" => :P_Za, "Zalpha" => :P_Za)
end




"""Given the result from marginalise - calculate quantiles of the PDFs"""
function summarise_results(dims, m)


    q, ss_v = quantiles_epdf(dims[:v], m[:P_v])
    q, ss_g = quantiles_epdf(dims[:logg], m[:P_logg])
    q, ss_T = quantiles_epdf(dims[:Teff], m[:P_Teff])
    q, ss_Z = quantiles_epdf(dims[:FeH], m[:P_FeH])
    q, ss_a = quantiles_epdf(dims[:alpha], m[:P_a])
    
    return Dict(:v => ss_v, :logg => ss_g, :Teff => ss_T,
                :FeH => ss_Z, :alpha => ss_a, :quantiles => q)
end


"""Given the results from the function above - print a summary"""
function print_results(s)

    # This relies on the format used above
    i_p16 = argmin(abs.(s[:quantiles].-0.16))
    i_p50 = argmin(abs.(s[:quantiles].-0.50))
    i_p84 = argmin(abs.(s[:quantiles].-0.84))
    
    @printf "   v = %0.2f [%0.2f, %0.2f]\n" s[:v][i_p50] s[:v][i_p16] s[:v][i_p84]
    @printf "Teff = %0.2f [%0.2f, %0.2f]\n" s[:Teff][i_p50] s[:Teff][i_p16] s[:Teff][i_p84]
    @printf "Fe/H = %0.2f [%0.2f, %0.2f]\n" s[:FeH][i_p50] s[:FeH][i_p16] s[:FeH][i_p84]
    @printf "logg = %0.2f [%0.2f, %0.2f]\n" s[:logg][i_p50] s[:logg][i_p16] s[:logg][i_p84]
    @printf "α/Fe = %0.2f [%0.2f, %0.2f]\n" s[:alpha][i_p50] s[:alpha][i_p16] s[:alpha][i_p84]

end



"""Find the first crossing of the function with target value"""
function lin_interp_solve(x, y, y_wanted; istart=1)
    if length(x) != length(y)
        println("Problem!!!")
        println("x=", x)
        println("y=", y)
        raise
    end
    i_below = -1
    i_above = -1
    previous_diff = y[1]-y_wanted
    for i=istart:length(x)
        # Check first for equality
        if y[i] == y_wanted
            i_below = i
            i_above = i
            break
        else
            diff = y[i]-y_wanted
            if sign(diff) != sign(previous_diff)
                # We have a zero-crossing!
                i_below = i-1
                i_above = i
                break
            end
        end
    end

    if i_below < 1
        solution = nothing
    else
        # Ok we have a solution. If this is exact it is easy
        if (i_below == i_above)
            solution = x[i_below]
        else
            # Linear interpolation needed
            delta_y = y[i_above]-y[i_below]
            delta_x = x[i_above]-x[i_below]

            solution = x[i_below] + (y_wanted-y[i_below])*delta_x/delta_y
        end
    end

    return solution
end

"""Find quantiles for a given probability distribution"""
function quantiles_epdf(x_pdf, pdf; quantiles=[0.025, 0.16, 0.5, 0.84, 0.975])

    N_q = length(quantiles)
    cpdf = cumsum(pdf, dims=1)
    cpdf = cpdf/maximum(cpdf)
    values = zeros(N_q)


    # Check for finiteness
    if !all(isfinite.(cpdf))
        println("PDF with NaNs!")
        values = [NaN for i in 1:N_q]
    else
        for i=1:N_q
            tmp = lin_interp_solve(x_pdf, cpdf, quantiles[i])
            if tmp == nothing
                values[i] = -99999.0
            else
                values[i] = tmp
            end
            
        end
    end
    
    return (quantiles, values)
end

function stringify_q(q)
    tmp = @sprintf "%.1f" 100*q
    return replace(tmp, "." => "p")
end
    

#
# Functions to get back spectra
#

"""get_best_fit

Create a best-fit spectrum using the model. This function just gets the best fit
model."""
function get_best_fit(dims, g, A, i_min; λobs=4750:1.25:9300)


    i_v, i_g, i_T, i_Z, i_a = i_min

    # The best-fit velocity
    best_v = dims[:v][i_v]
    
    # Get the spectrum at velocity=0.0
    fM = g[:, i_g, i_T, i_Z, i_a]

    # Interpolate it in velocity
    fMobs = _fast_interpolate(λobs, dims[:λ], fM, best_v)
    
    # Scale it by the best-fit scale
    fMobs = fMobs .* A[i_v, i_g, i_T, i_Z, i_a]

    return Dict(:λ => λobs, :flux => fMobs, :scl=>A[i_v,i_g,i_T, i_Z, i_a],
                :modelflux => fM)
end



"""get_fitted_spectrum_pdf

Create a best-fit spectrum using full likelihood. This loops over all models
and constructs a median spectrum and the quantiles around this.  However it does
this per velocity and then combines again over velocity. 
"""
function get_best_fit(dims, g, A, p; λobs=4750:1.25:9300)


    


    
    return Dict(:λ => λobs, :flux => fMobs, :scl=>A[i_v,i_g,i_T, i_Z, i_a],
                :modelflux => fM)
end


