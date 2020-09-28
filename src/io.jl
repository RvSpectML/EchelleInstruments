"""
Shared code for file io.  Currently, a convenience wrapper for FITSIO.

Author: Eric Ford
Created: August 2020
"""

using DataFrames, FITSIO

""" `read_header( filename )`
Read header from FITS file and return Dict with contents.
Optional inputs:
- hdu: Specifies which HDU to read from the FITS file.  (Default: 1)
"""
function read_header(fn::String; header_idx::Integer = 1)
    #println("# Reading: ",fn, " hdu= ",header_idx)
    f = FITS(fn)
    @assert 1<=header_idx<=length(f)
    #@assert read_key(f[header_idx],"SKY-OBJ")[1] == "Solar"
    hdr = FITSIO.read_header(f[header_idx])
    metadata = Dict(zip(map(k->Symbol(k),hdr.keys),hdr.values))
end

""" `read_fits_header( filename )`
Read header from FITS file and return Dict with contents.
Optional inputs:
- hdu: Specifies which HDU to read from the FITS file.  (Default: 1)
"""
function read_fits_header(fn::String; header_idx::Integer = 1)
    read_header(fn,header_idx=header_idx)
end

""" `read_metradata_from_fits( filename, fields )`
Read metadata in FITS header for specified keys and return data as a Dict.  `fields` can be an array of symbols or strings.
Optional inputs:
- Passing both fields (an array of symbols) and fields_str (an array of strings) as named parameters allows for differences in string used in FITS file and Symbol used in resulting Dict.
- hdu: Specifies which HDU to read from the FITS file.  (Default: 1)
"""
function read_metradata_from_fits
end

function read_metadata_from_fits(fn::String, fields::Array{Symbol,1} ; hdu::Integer = 1)
    fields_str=string.(fields)
    read_metadata_from_fits(fn,hdu=hdu,fields=fields,fields_str=fields_str)
end

function read_metadata_from_fits(fn::String, fields_str::AbstractArray{AS,1} ; hdu::Integer = 1)  where { AS<:AbstractString }
    fields = map(f->Symbol(f),fields_str)
    read_metadata_from_fits(fn,hdu=hdu,fields=fields,fields_str=fields_str)
end

function read_metadata_from_fits(fn::String; fields::Array{Symbol,1}, fields_str::AbstractArray{AS,1}, hdu::Integer = 1 )  where { AS<:AbstractString }
    @assert length(fields) == length(fields_str)
    @assert 1 <= hdu <= 3
    f = FITS(fn)
    hdr = FITSIO.read_header(f[hdu])
    # Check that header has all expected fields
    #println(fn)
    for field in fields_str
        #println("  ",field,": ",typeof(hdr[field]))
        @assert findfirst(isequal(field),keys(hdr)) != nothing
    end
    values = map(s->hdr[s],fields_str)

    df = Dict{Symbol,Any}(zip(fields,values))
    return df
end
