"""
    FieldTools

This module is as grab-bag of tools ported from python that have been useful for working 
with real and complex fields
"""
module FieldTools

    import QGas.NumericalTools.ArrayDimensions as AD
    import LinearAlgebra as LA

    """
        field_from_magnitude(magnitude)

    Return a field given a magnitude.
    """
    function field_from_magnitude(mag::Number) # Method for a single number
        return real(mag) < 0 ? 0 : sqrt(real(mag))
    end
    function field_from_magnitude(mag::Array{<:Number})
        field = zeros(ComplexF64, size(mag))

        for i in CartesianIndices(mag)
            field[i] = field_from_magnitude(mag[i])
        end

        return field
    end

    """
        normalize_field!(ϕ::Array{<:Number}; 
                            norm=1.0, 
                            adims::Union{AD.Dimensions, Nothing}=nothing, 
                            mode::String="field")

    Returns the field and the normalization coefficient.

    field : Field to be normalized

    norm : value to be normalized to

    dims = ArrayDimensions.Dimensions object associated with Field, used obtain 
        the integration measure

    Mode :
        "field" assume object is a complex field and norm to ∑ |ϕ|^2
        "magnitude" assume object is a distribution and use ∑ ϕ
    """
    function normalize_field!(ϕ::Array{<:Number}; 
                            norm=1.0, 
                            adims::Union{AD.Dimensions, Nothing}=nothing, 
                            mode::String="field")

        dx = (isnothing(adims) ? 1.0 :  AD.DeltasProd(adims)) / norm

        if mode == "magnitude"
            total = sum(ϕ) * dx
        elseif mode == "field"
            total = LA.norm(ϕ) * sqrt(dx) # LA.norm has sqrt built in
        else
            msg = "invalid mode = $(mode)"
            throw( ArgumentError(msg) )
        end

        ϕ ./= total

        return ϕ, total
    end

end
