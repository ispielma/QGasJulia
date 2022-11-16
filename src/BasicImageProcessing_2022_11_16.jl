module BasicImageProcessing_2022_11_16

import HDF5

greet() = print("Hello World!")

function OpenData(filename::String)
    HDF5.h5open(filename, "r") do fid
        image_grp = fid["images"]
        MOT_x_grp = image_grp["MOT_x/fluorescence"] # ["MOT", "MOT_dark", "cMOT", "cMOT_dark"]
        MOT_y_grp = image_grp["MOT_y/absorption"] # ["MOT_TOF", "MOT_TOF_dark", "MOT_TOF_probe"]
        MOT_z_grp = image_grp["MOT_z/fluorescence"] # ["MOT", "MOT_dark", "cMOT", "cMOT_dark"]
        global MOT_x = HDF5.read(MOT_x_grp["MOT"])
        global MOT_dark_x = HDF5.read(MOT_x_grp["MOT_dark"])
    end

    MOT_fl_x = convert(Array{Int16}, MOT_x) - convert(Array{Int16}, MOT_dark_x)
    
    return MOT_fl_x

end 

end # module BasicImageProcessing_2022_11_16
