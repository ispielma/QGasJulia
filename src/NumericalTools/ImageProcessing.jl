"""
ImageProcessing

This module is as grab-bag of tools ported from python for image analysis, typically of
cold atom data
"""
module ImageProcessing

    import LinearAlgebra as LA
    import FFTW
    import Distributions as Dist
    import Statistics as Stat

    import QGas.NumericalTools.ArrayDimensions as AD

    """
    ImageInfo
    
    contains the basic information requed to interpert cold atom images
    """ 
    
    _filter_width = 1.0 # a default width
    struct ImageInfo{dim}
        npnts::Tuple{Vararg{Int64, dim}}
        dxs::Tuple{Vararg{Float64, dim}}
        units::Tuple{Vararg{String, dim}}
        filter_width::Float64 # width of gaussian for data filtering
        filter_function::Function # Function to use for filtering
        filter::Array{Float64, dim} # filter for use in low-pass for recovery of bad data
        ranges::Tuple{Vararg{AbstractVector, dim}} # Intended to be 1D ranges
        ndrange::AD.NDRange{dim} # Mesh-grid like object
    end 
    function ImageInfo(
            npnts::Tuple{Vararg{Int64, dim}}, 
            dxs::Tuple{Vararg{Float64, dim}}, 
            units::Tuple{Vararg{String, dim}}; 
            filter_width=_filter_width) where dim
        ranges = Tuple(range(-dx*npnt/2, dx*npnt/2, length=npnt) for (dx, npnt) in zip(dxs, npnts))
        ndrange = AD.NDRange(ranges)

        # Future option for more filter choices
        filter = FFTW.fftshift([exp( -0.5.*sum( (xy./filter_width).^2 )) for xy in ndrange]) # Build filter here
        
        return ImageInfo{dim}(npnts, dxs, units, filter_width, gauss_window, filter, ranges, ndrange)
    end
    function ImageInfo(
            npnts::Tuple{Vararg{Int64, dim}},
            dxs::Tuple{Vararg{Float64, dim}}; 
            filter_width=_filter_width) where dim
        units = Tuple("" for i in 1:dim)
        return ImageInfo(npnts, dxs, units); filter_width=filter_width
    end
    function ImageInfo(npnts::Tuple{Vararg{Int64, dim}}; filter_width=_filter_width) where dim
        dsx = Tuple(1.0 for i in 1:dim)
        return ImageInfo(npnts, dxs; gauss_width=gauss_width)
    end

    """
        tukey_window(f, a)

    f : defines the quantity on which the window is evaluated where the window cuts off
        for |f| ≥ 1/2.  Note that often the tukey is defined instead from 0 to 1 rather
        than -1/2 to 1/2

    a : tukey window parameter with 0 being a box window and 1 being a cos window

    scale : x axis scale factor.  I use this for radial tukeys that make a window starting 
        at zero and naturally ending at 1.  In this case one would set scale = 0.5

        min : y minimum of the window (usually 0.0)

        max : y maximum of the window (usually 1.0)

    """
    function tukey_window(f::Number, a; scale=1.0, min=0.0, max=1.0)
        g = abs(f*scale)
        if g > 0.5
            window = 0
        elseif g < 0.5*(1 - a) 
            window = 1
        else
            window = ( 1-cos(pi*(1-2*g)/a) ) / 2
        end
    
        # Now scale to fit in the mix / max range
        window = (max-min)*window + min

        return window
    end
    function tukey_window(f::AbstractArray, a; scale=1.0, min=0.0, max=1.0)
    
        window = similar(f, Float64) 
    
        for i in CartesianIndices(f)
            window[i] = tukey_window(f[i], a; scale=scale, min=min, max=max)
        end
    
        return window
    end

    """
        gauss_window(w, a)

    f : defines the quantity on which the window is evaluated.

    σ : gaussian width

    min : y minimum of the window (usually 0.0)

    max : y maximum of the window (usually 1.0)

    """
    function gauss_window(f::Number, σ; min=0.0, max=1.0)

        window = exp( -0.5*sum( (f/σ)^2 ))
    
        # Now scale to fit in the mix / max range
        window = (max-min)*window + min

        return window
    end
    function gauss_window(f::AbstractArray, σ; min=0.0, max=1.0)
    
        window = similar(f, Float64) 
    
        for i in CartesianIndices(f)
            window[i] = gauss_window(f[i], σ; min=min, max=max)
        end
    
        return window
    end

    """
        conv_with_weights(data, w, f)
    
    computes the convilution of a dataset with a function including the weight factors w.  The norm of data will be 
        conserved.

    data: array to be convolved
    w : weight factor (equal to 1/σ, with uncertainty σ)
    f : filter array (equal in shape to data), does not need to be normalized 
    """
    function conv_with_weights(data::AbstractArray, w::AbstractArray, f::AbstractArray)

        f_fft = FFTW.rfft(f)


        # Compute numerator
        d = size(data)[1] # needed for irfft
        ans = FFTW.irfft( FFTW.rfft(data .* w.^2) .* f_fft, d)

        # Compute denominator
        ans ./= FFTW.irfft( FFTW.rfft(w.^2) .* f_fft, d)

        return ans
    end

    """
    preprocess_probe(probe, dark, filter)

    preprocesses a probe beam

    probe : probe image
    dark : dark image
    filter : the fancy version convolves the data with filter to make a smoothed version of the data (weighted by the uncertanties)
        and fills in the invalid (negative) data with the filtered data.
    """
    function preprocess_probe(probe::Array{Float64}, dark::Array{Float64})
        probe -= dark

        return probe
    end
    # These methods convert to Float64 only if needed
    preprocess_probe(probe::Array{Float64}, dark::AbstractArray) = preprocess_probe(probe, convert(Array{Float64}, dark))
    preprocess_probe(probe::AbstractArray, dark::Array{Float64}) = preprocess_probe(convert(Array{Float64}, probe), dark)
    preprocess_probe(probe::AbstractArray, dark::AbstractArray) = preprocess_probe(convert(Array{Float64}, probe), convert(Array{Float64}, dark))

    function filter_probe(probe, filter)

        weights = [p <= 0.0 ? 0.0 : sqrt(p) for p in probe]
        probe_smooth = conv_with_weights(probe, weights, filter)
        mask = weights .== 0
        probe[mask] .= probe_smooth[mask]

        return probe
    end

    preprocess_probe(probe, dark, filter::AbstractArray) = filter_probe(preprocess_probe(probe, dark), filter)
    preprocess_probe(probe, dark, imginfo::ImageInfo) = preprocess_probe(probe, dark, imginfo.filter)

    """
        average_dark(darks)

    return the average of the dark frames as well as the average noise
    """
    function average_dark(darks::Array{Float64})
        mean = Stat.mean(darks; dims=1)
        std = Stat.std(darks; mean=mean, dims=1)

        # For some reason mean and std do not reduce the dimension of the array and have
        # a 1 sized dimension on the first index
        reshape(mean, size(mean)[2:end] )
        reshape(std, size(std)[2:end] )

        return mean[1,:,:],  std[1,:,:]
    end
    average_dark(darks::AbstractArray) = average_dark(convert(Array{Float64}, darks))

    """
        α(I, I_bar, ΔI2)

    Generates a single realization of a field consistant with:

    I   The observed field
    I_bar   The known mean field
    ΔI2     The pixel-by-pixel variance variance

    """
    function α_realization(I, I_bar, ΔI2; addnoise=false)
        κ_min = 0.01

        # Field variance
        Δα = sqrt.(2 .* I_bar .* ((1 .+ ΔI2 ./ (4 .* I_bar.^2)).^0.5 .- 1))

        α = sqrt.(I_bar) .+ 0im
        α .+= (I - I_bar) ./ (2 .* α)
        
        # Now add noise  
        if addnoise
            α .+= 1.0im * Dist.rand.(Dist.Normal.(0, Δα./2))
        end

        return α
    end
end