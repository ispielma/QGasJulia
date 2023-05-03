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

    #=
     ######  ######## ########  ##     ##  ######  ########  ######
    ##    ##    ##    ##     ## ##     ## ##    ##    ##    ##    ##
    ##          ##    ##     ## ##     ## ##          ##    ##
     ######     ##    ########  ##     ## ##          ##     ######
          ##    ##    ##   ##   ##     ## ##          ##          ##
    ##    ##    ##    ##    ##  ##     ## ##    ##    ##    ##    ##
     ######     ##    ##     ##  #######   ######     ##     ######
    =#

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
        filter = gauss_window(ndrange, filter_width)
        filter = FFTW.fftshift(filter / sum(filter))
        
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

    #=
       ##      ## #### ##    ## ########   #######  ##      ##  ######
       ##  ##  ##  ##  ###   ## ##     ## ##     ## ##  ##  ## ##    ##
       ##  ##  ##  ##  ####  ## ##     ## ##     ## ##  ##  ## ##
       ##  ##  ##  ##  ## ## ## ##     ## ##     ## ##  ##  ##  ######
       ##  ##  ##  ##  ##  #### ##     ## ##     ## ##  ##  ##       ##
       ##  ##  ##  ##  ##   ### ##     ## ##     ## ##  ##  ## ##    ##
        ###  ###  #### ##    ## ########   #######   ###  ###   ######
    =#

    # Here I have a series of distributions and windows.  The key distinction between the 
    # two functions are that windows have a peak value of 1 (i.e. inside the window), while
    # distributions are normalized to 1

    # Plan: for each item the "base" function will be the distribution as it will be the most
    # simple to express, and the normed version will be a different function that calls that 
    # and then appends the norm.
    # I will be careful to minimize the number of times the norm is computed in the case of vector
    # operations (i.e don't comput it inside the scalar method if it is called from the vector method)

    """
        shift_window(f, min, max)

    shifts a function that has a domain [0,1] to be in the domain [min,max]
    """
    shift_window(f, min, max) = (max-min).*f .+ min

    """
        tukey_window(f, a)

    f : defines the quantity on which the window is evaluated where the window cuts off
        for |f| ≥ 1/2.  Note that often the tukey is defined instead from 0 to 1 rather
        than -1/2 to 1/2

    a : tukey window parameter with 0 being a box window and 1 being a cos window

    w : x axis width.  I use this for radial tukeys that make a window starting 
        at zero and naturally ending at 1.  In this case one would set width = 2

    """
    function tukey_window(f::Number, a)
        f = abs(f)
        if f > 0.5
            window = 0
        elseif f < 0.5*(1 - a) 
            window = 1
        else
            window = ( 1-cos(pi*(1-2*f)/a) ) / 2
        end

        return window
    end 

    tukey_window(f, a, w) = tukey_window(f./w, a)
    tukey_window(f, a, w, cen) = tukey_window(f-cen, w, a)
    
    # Multi dimension methods
    tukey_window(fs::Tuple{Vararg{T, dim}}, a) where {T, dim} = tukey_window(sqrt(sum(fs.^2)), a)
    tukey_window(fs::Tuple{Vararg{T, dim}}, a, w::Number) where {T, dim} = tukey_window(fs ./ w, a)
    tukey_window(fs::Tuple{Vararg{T, dim}}, a, ws) where {T, dim} = tukey_window(Tuple( f / w for (f, w) in zip(fs, ws)), a)
    tukey_window(fs::Tuple{Vararg{T, dim}}, a, ws, cens) where {T, dim} = tukey_window(Tuple( f - c for (f, c) in zip(fs, cens)), a, ws)

    # version for many points
    tukey_window(f::AbstractArray, args...) = [tukey_window(f[i], args...) for i in CartesianIndices(f) ]


    """
        gauss_window(w, a)

    f : defines the quantity on which the window is evaluated.

    w : gaussian width

    """

    gauss_window(f::Number) = exp(-0.5*f^2 )

    gauss_window(f, w) = gauss_window(f./w)
    gauss_window(f, w, cen) = gauss_window(f-cen, w)
    
    # Multi dimension methods
    gauss_window(fs::Tuple{Vararg{T, dim}}) where {T, dim} = gauss_window(sqrt(sum(fs.^2)))
    gauss_window(fs::Tuple{Vararg{T, dim}}, w::Number) where {T, dim} = gauss_window(fs ./ w)
    gauss_window(fs::Tuple{Vararg{T, dim}}, ws) where {T, dim} = gauss_window(Tuple( f / w for (f, w) in zip(fs, ws)))
    gauss_window(fs::Tuple{Vararg{T, dim}}, ws, cens) where {T, dim} = gauss_window(Tuple( f - c for (f, c) in zip(fs, cens)), ws)

    # version for many points
    gauss_window(f::AbstractArray, args...) = [gauss_window(f[i], args...) for i in CartesianIndices(f) ]


    
    #=
       ########  #######   #######  ##        ######
          ##    ##     ## ##     ## ##       ##    ##
          ##    ##     ## ##     ## ##       ##
          ##    ##     ## ##     ## ##        ######
          ##    ##     ## ##     ## ##             ##
          ##    ##     ## ##     ## ##       ##    ##
          ##     #######   #######  ########  ######
    =#

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

    function PSD(Field; shift=false)
        psd = abs.(FFTW.rfft(Field)).^2
    
        if shift
            psd = FFTW.fftshift(psd)
        end
    
        return psd
    end

    #=
        ######   #######  ##       ########        ###    ########  #######  ##     ##  ######
       ##    ## ##     ## ##       ##     ##      ## ##      ##    ##     ## ###   ### ##    ##
       ##       ##     ## ##       ##     ##     ##   ##     ##    ##     ## #### #### ##
       ##       ##     ## ##       ##     ##    ##     ##    ##    ##     ## ## ### ##  ######
       ##       ##     ## ##       ##     ##    #########    ##    ##     ## ##     ##       ##
       ##    ## ##     ## ##       ##     ##    ##     ##    ##    ##     ## ##     ## ##    ##
        ######   #######  ######## ########     ##     ##    ##     #######  ##     ##  ######
    =#

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
    filter_probe(probe, imginfo::ImageInfo) = filter_probe(probe, imginfo.filter)


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