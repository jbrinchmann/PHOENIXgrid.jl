# Functions to create mock observations.

"""Create a fake observed spectrum

If median S/N is given, then a flat uncertainty array is created. 

"""
function mock_observed_spectrum(sp; v_los=0.0, median_sn=10.0,
                                σ=nothing, λobs=nothing)

    if λobs == nothing
        λobs = collect(4750:1.25:9300) # Canonical MUSE wavelength array
    end
    
    # First shift in velocity
    fMshifted = _fast_interpolate(λobs, sp[:λ], sp[:flux], v_los)

    # Then add noise
    if σ==nothing
        med_signal = median(fMshifted)
        σ = med_signal/median_sn
    end

    N = length(λobs)
    fobs = zeros(N)
    dfobs = zeros(N)
    ϵ = zeros(N)
    for i=1:N
        if length(σ) == 1
            dfobs[i] = σ
        else
            dfobs[i] = σ[i]
        end
        ϵ[i] = randn(Float64)*dfobs[i]
        fobs[i] = fMshifted[i] + ϵ[i]
    end            
        

    return Dict(:λ => λobs, :flux=>fobs, :dflux=>dfobs, :ϵ=>ϵ)
    
end
