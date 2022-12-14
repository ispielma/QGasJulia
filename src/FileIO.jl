module FileIO

import HDF5
import Formatting

struct DataFiles
    prefix::String # Beginning of filename 
    initialindex::Integer # First index 
    finalindex::Integer # Last index
    extension::String # End of filename (like h5 or txt)
end 
DataFiles(prefix::String, initialindex::Integer, finalindex::Integer; extension::String=".h5") = DataFiles(prefix, initialindex, finalindex, extension)

struct ImageGroup
    orientation::String # Orientation label for saved image.
    label::String # Label of saved image (ignore if empty)
    images::Vector{String} # Identifiers of saved images.
end

"""
FileNameArray : this function returns an array of files that should be evaulated
    Takes as parameters a list of DataFile structs
    with keyword parameters
    extension: the filename
    formatstring: a python style format string directing how to format the index
        "0>4d" will make strings like 0111 or 0001
"""
function FileNameArray(files::Vector{DataFiles}; formatstring::String="d")


    # Now build the tuple of file names
    filelist = Array{String}(undef, 0)

    indexformat = Formatting.format("{{1}}{{2:{}}}{{3}}", formatstring)

    for file in files	
    
        if (file.finalindex < file.initialindex)
            throw( KeyError("FileNameArray $(file.filename): initial index $(file.initialindex) must be smaller than Final index $(file.finalindex)") )
        end
        
        # And here is the list for this element of files.  
        newfiles = [Formatting.format(indexformat, file.prefix, i, file.extension) for i in file.initialindex:file.finalindex]
        append!(filelist, newfiles)
    end
    
    return filelist
end
FileNameArray(file::DataFiles; kwargs...) = FileNameArray([file]; kwargs...)

"""Get previously saved image from the h5 file.

h5_file : h5 reference or a file name

Args:
    orientation (str): Orientation label for saved image.
    label (str): Label of saved image.
    image (str): Identifier of saved image.

Raises:
    Exception: If the image or paths do not exist.

Returns:
    2-D image array.
"""
function getimage(h5_file::HDF5.File, orientation::String, label::String, image::String)

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

    return HDF5.read(h5_file[h5path][image]) 
end
function getimage(filename::String, orientation::String, label::String, image::String)
    (image) = HDF5.h5open(filename, "r") do h5_file
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
function getimages(h5_file::HDF5.File, imagegroups::Vector{ImageGroup})
    local imagedict = Dict{String, Dict}()
    for imagegroup in imagegroups
        if ~(imagegroup.orientation in keys(imagedict))
            imagedict[imagegroup.orientation] = Dict{String, Dict}()
        end

        if ~(imagegroup.label in keys(imagedict[imagegroup.orientation]))
            imagedict[imagegroup.orientation][imagegroup.label] = Dict{String, Any}()
        end

        for image in imagegroup.images
            imagedict[imagegroup.orientation][imagegroup.label][image] = getimage(h5_file, imagegroup.orientation, imagegroup.label, image)
        end
    end
    return imagedict
end
function getimages(filename::String, imagegroups::Vector{ImageGroup})

    imagedict = HDF5.h5open(filename, "r") do h5_file
        imagedict = getimages(h5_file, imagegroups)
    end

    return imagedict
end
getimages(filename::Union{HDF5.File, String}, imagegroup::ImageGroup) = getimages(filename, [imagegroup])

end # module FileIO