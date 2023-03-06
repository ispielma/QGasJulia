module FieldSimulation
    import QGas.NumericalTools.ArrayDimensions as AD

    # Should all of the field names be lower case to be standard?
    Base.@kwdef mutable struct Field
        # Config items
        BasePosition::Float64 = 0.0
        Position::Float64 = 0.0
        norm::Float64 = 1.0
        MaxSpatialFrequency::Float64 = NaN
        RandomSeed::Int64 = 0 # interperted as no seed

        # Data and support
        adims::AD.Dimensions = AD.Dimensions()
        BaseField::Array{Float64} = []
        CurrentField::Array{Float64} = []

        # Functions to be evaluated after each update
        UpdateFunction::Function = identity
        UpdateFunction_args = tuple()
        UpdateFunction_kwargs = Dict{String, Any}()

        # Postprocessing Functions
        PostProcessFunction::Function = identity
        PostProcessFunction_args = tuple()
        PostProcessFunction_kwargs = Dict{String, Any}()
        PostProcessFunction_results = tuple()
    end



end 
