module FileIO

import HDF5
import Format

#=
 ######   ######## ##    ## ######## ########     ###    ##
##    ##  ##       ###   ## ##       ##     ##   ## ##   ##
##        ##       ####  ## ##       ##     ##  ##   ##  ##
##   #### ######   ## ## ## ######   ########  ##     ## ##
##    ##  ##       ##  #### ##       ##   ##   ######### ##
##    ##  ##       ##   ### ##       ##    ##  ##     ## ##
 ######   ######## ##    ## ######## ##     ## ##     ## ########
=#

"""
dict_to_h5

recursivly iterates through d and saves contents into the h5 reference
"""
function dict_to_h5(h5ref::HDF5.H5DataStore, d::Dict{String, Any})
    for (key, item) in d

        # Clear out old data
        if key in keys(h5ref)
            HDF5.delete_object(h5ref, key)
        end

        # recurse if needed
        if typeof(item) == Dict{String, Any}
            HDF5.create_group(h5ref, key)

            dict_to_h5(h5ref[key], item)
        else
            HDF5.write_dataset(h5ref, key, item)
        end
        
    end
end
dict_to_h5(filename::AbstractString, d::Dict{String, Any}) = HDF5.h5open( f -> dict_to_h5(f, d), filename, "cw")

"""
h5_to_dict

Returns a dictionary containg the data from the h5 file
"""
h5_to_dict(filename) = HDF5.h5open(h5_to_dict, filename, "r")
h5_to_dict(h5ref::HDF5.H5DataStore) = read(h5ref)

#=
##          ###    ########   ######   ######  ########  #### ########  ########
##         ## ##   ##     ## ##    ## ##    ## ##     ##  ##  ##     ##    ##
##        ##   ##  ##     ## ##       ##       ##     ##  ##  ##     ##    ##
##       ##     ## ########   ######  ##       ########   ##  ########     ##
##       ######### ##     ##       ## ##       ##   ##    ##  ##           ##
##       ##     ## ##     ## ##    ## ##    ## ##    ##   ##  ##           ##
######## ##     ## ########   ######   ######  ##     ## #### ##           ##
=#

struct LabscriptConfig
    experiment_shot_storage::String
    output_folder_format::String
    filename_prefix_format::String
    extension::String
    file_fstring::String
end 
LabscriptConfig(experiment_shot_storage, output_folder_format, filename_prefix_format; extension=".h5") = LabscriptConfig(
    experiment_shot_storage, 
    output_folder_format, 
    filename_prefix_format,
    extension, 
    labscript_file_fstring(experiment_shot_storage, output_folder_format, filename_prefix_format; extension=extension)
)

"""
    LabscriptSequence

Basic information used to define an indificual sequence consisting of a collection of shots
"""
struct LabscriptSequence
    script_basename::String
    year::Integer
    month::Integer
    day::Integer
    index::Integer
    shots::Vector{Integer}
end 
# LabscriptSequence(script_basename, year, month, day, index, shots) = LabscriptSequence(script_basename, year, month, day, index, Vector(shots))

"""
    ImageGroup
    
Information needed to locate a specific image
"""
struct ImageGroup
    orientation::String # Orientation label for saved image.
    label::String # Label of saved image (ignore if empty)
    images::Vector{String} # Identifiers of saved images.
end

raw"""
    labscript_file

Creates a typical labscript style prefix for a file including the path

An example file might be:
    "2024-03-22_0010_BEC_elongatedTrap_FK_PEP_fieldMonitor_066.h5"

and for the path:
    "~/Documents/data/RbChip/Labscript/BEC_elongatedTrap_FK_PEP_fieldMonitor/2024/03/22/0010"

Notice the partly redundent information where in this example we had:
    experiment_shot_storage = "~/Documents/data/RbChip/Labscript"
    script_basename = "BEC_elongatedTrap_FK_PEP_fieldMonitor"
    year = 2024
    month = 03
    day = 22
    sequence_index = 10
    shot = 66

It is important to use format strings such as "0>4d" that will make strings like 0111 or 0001 to correspond
to the configuration information usually in the labscript config such as:

    output_folder_format = %%Y\%%m\%%d\{sequence_index:04d}
    filename_prefix_format = %%Y-%%m-%%d_{sequence_index:04d}_{script_basename}

Here we will use a slightly different format based on the python f-string.  For the example above we would have

    experiment_shot_storage = "~/Documents/data/RbChip/Labscript"
    output_folder_format = "{year:04d}/{month:02d}/{day:02d}/{sequence_index:04d}"
    filename_prefix_format = "{year:04d}-{month:02d}-{day:02d}_{sequence_index:04d}_{script_basename}_{shot:04d}"
    script_basename = "BEC_elongatedTrap_FK_PEP_fieldMonitor"
    year = 2024
    month = 03
    day = 22
    sequence_index = 10
    shot=66
"""
function labscript_file_fstring(
    experiment_shot_storage::String,
    output_folder_format::String,
    filename_prefix_format::String;
    extension::String=".h5"
    )

    ks = (
        "script_basename",
        "year",
        "month",
        "day",
        "sequence_index",
        "shot"
    )

    experiment_shot_storage = _fstring_key_to_index(experiment_shot_storage, 0, ks)

    output_folder_format = _fstring_key_to_index(output_folder_format, 0, ks)

    filename_prefix_format = _fstring_key_to_index(filename_prefix_format, 0, ks)

    return "$(experiment_shot_storage)/$(output_folder_format)/$(filename_prefix_format)$(extension)"
end

"""
fstring_key_to_index

Takes a string formatted like a python f-string using the kwarg format and replaces with numerical indices like the args format.
"""
fstring_key_to_index(s::String, args...; kwargs...) =  _fstring_key_to_index(s, length(args),  keys(kwargs))

"""
    _fstring_key_to_index

Helper function

    n_args : Number of arguments
    ks : tuple or vector of strings
"""
function _fstring_key_to_index(s::String, n_args::Integer, ks)
    n_ks = length(ks)
    rules = Tuple("{"*string(ks[idx]) => "{"*string(idx+n_args) for idx in 1:n_ks)
    replace(s, rules...)
end

"""
file_name_array : this function returns an array of files that should be evaulated
    Takes as parameters a list of DataFile structs
    with keyword parameters
    extension: the filename
    formatstring: a python style format string directing how to format the index
        "0>4d" will make strings like 0111 or 0001
"""
function file_name_array(sequences::Vector{LabscriptSequence}, labconfig::LabscriptConfig)

    # Now build the array of file names
    filelist = Array{String}(undef, 0)

    for sequence in sequences	
        
        # And here is the list for this element of files.  
        newfiles = [Format.format(
            labconfig.file_fstring, 
            sequence.script_basename, 
            sequence.year,
            sequence.month,
            sequence.day,
            sequence.index,
            i) for i in sequence.shots]
        append!(filelist, newfiles)
    end
    
    return filelist
end
file_name_array(file::LabscriptSequence, labconfig::LabscriptConfig) = file_name_array([file], labconfig)

"""Get previously saved image from the h5 file.

h5_file : h5 reference or a file name

Args:
    orientation (str): Orientation label for saved image.
    label (str): Label of saved image.
    image (str): Identifier of saved image.

KWargs:
    cast = type to convert to if not nothing

Raises:
    Exception: If the image or paths do not exist.

Returns:
    2-D image array.
"""
function get_image(h5_file::HDF5.File, orientation::String, label::String, image::String; cast=Float64)

    h5path = "images"
    if ~haskey(h5_file, h5path)
        throw( KeyError("File does not contain any images") )
    end

    if ~haskey(h5_file[h5path], orientation)
        throw( KeyError("File does not contain any images with orientation $(orientation)") )
    end
    h5path = h5path*"/"*orientation

    if ~haskey(h5_file[h5path], label) # no label is OK
        throw( GenericException("File does not contain any images with label $(label)") )
    end
    h5path = h5path*"/"*label

    if ~haskey(h5_file[h5path], image)
        throw( GenericException("Image $(image) not found in file") )
    end

    data = HDF5.read(h5_file[h5path][image]) 

    data = ( cast == nothing ? data : convert(Array{cast}, data) )

    return data
end
function get_image(filename::String, orientation::String, label::String, image::String)
    image = HDF5.h5open(filename, "r") do h5_file
        get_image(h5_file, orientation, label, image)
    end    
    return image
end

"""Get previously saved image from the h5 file.

h5_file : h5 file name

Args:
    orientation (str): Orientation label for saved image.
    label (str): Label of saved image.
    image (str): Identifier of saved image.

Returns:
    Dict of 2-D image arrays.
"""
function get_images(h5_file::HDF5.File, imagegroups::Vector{ImageGroup})
    local imagedict = Dict{String, Dict}()
    for imagegroup in imagegroups
        if ~(imagegroup.orientation in keys(imagedict))
            imagedict[imagegroup.orientation] = Dict{String, Dict}()
        end

        if ~(imagegroup.label in keys(imagedict[imagegroup.orientation]))
            imagedict[imagegroup.orientation][imagegroup.label] = Dict{String, Any}()
        end

        for image in imagegroup.images
            imagedict[imagegroup.orientation][imagegroup.label][image] = get_image(h5_file, imagegroup.orientation, imagegroup.label, image)
        end
    end
    return imagedict
end
function get_images(filename::String, imagegroups::Vector{ImageGroup})

    imagedict = HDF5.h5open(filename, "r") do h5_file
        imagedict = get_images(h5_file, imagegroups)
    end

    return imagedict
end
get_images(filename::Union{HDF5.File, String}, imagegroup::ImageGroup) = get_images(filename, [imagegroup])

end # module FileIO