module QGas
    export FileIO
    export AtomicPhysics    
    export NumericalTools    
    export Simulation    

    include("AtomicPhysics.jl")

    include("FileIO.jl")

    include("NumericalTools/NumericalTools.jl")

    include("Simulation/Simulation.jl")

    # Check to see if this is being run from the command line
    if abspath(PROGRAM_FILE) == @__FILE__
        print("Running from command line")
    end
end
