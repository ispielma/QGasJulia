module QGas
export FileIO

import Formatting

include("FileIO.jl")

# Check to see if this is being run from the command line
if abspath(PROGRAM_FILE) == @__FILE__
    print("Running from command line")
end

end # module BasicImageProcessing_2022_11_16
