"""
   IO functions for the HARPS spectrograph
    https://www.eso.org/sci/facilities/lasilla/instruments/harps.html
Author: Alex Wise and collaborators
Created: April 2021
"""

"""Create Dataframe containing filenames and key data for all files YYYY-MM-DDThh:mm:ss.mss (year-month-day T hour:minute:second:milisecond in directory"""
function make_manifest(data_path::String ; max_spectra_to_use::Int = 1000 )
    df_filenames = EchelleInstruments.make_manifest(data_path,r"^2010-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}.\d{3}$")

    df_files = DataFrame(read_metadata(df_filenames.Filename[1]*"/e2ds.fits"))
    keys = propertynames(df_files)
    allowmissing!(df_files, keys[map(k->k∉[:Filename, :bjd, :target, :airmass],keys)] )

    if length(df_filenames.Filename) >= 2
        map(fn->add_metadata_from_fits!(df_files,fn*"/e2ds.fits"),df_filenames.Filename[2:end])
    end

    df_files

end

"""Create Dict containing filename and default metadata from file."""
function read_metadata(fn::Union{String,FITS})
    fields_to_save = metadata_symbols_default(HARPS2D())
    fields_str_to_save = metadata_strings_default(HARPS2D())
    dict = read_metadata_from_fits(fn,fields=fields_to_save,fields_str=fields_str_to_save)
    check_metadata_fields_expected_present(dict,fields_to_save)
    dict[:Filename] = fn
    if haskey(dict,:airmass) && typeof(dict[:airmass]) == String
       dict[:airmass] = parse(Float64,dict[:airmass])
    end
    return dict
end

function add_metadata_from_fits!(df::DataFrame, fn::String)
    println("# Reading metadata from ", fn)
    metadata_keep = read_metadata(fn)
    #reject file from df if it is a calibration (only keep spectra)
    #if metadata_keep[:target] != "Solar" return df end
    #calibration_target_substrings = ["Etalon","LFC","ThArS","LDLS","FlatBB"]
    #if any(map(ss->contains( metadata_keep[:target],ss),calibration_target_substrings)) 
    #   return df
    #end

    for k in propertynames(df)
        if typeof(metadata_keep[k]) == Missing continue end
        if typeof(metadata_keep[k]) == Nothing
            metadata_keep[k] = missing
            continue
        end
        if (eltype(df[!,k]) <: Union{Real,Union{<:Real,Missing}}) && !(typeof(metadata_keep[k]) <: Union{Real,Union{<:Real,Missing}})
            metadata_keep[k] = missing
        end
    end
    #println("metadata_keep = ", metadata_keep)
    if haskey(metadata_keep,:airmass) && typeof(metadata_keep[:airmass]) == String
       metadata_keep[:airmass] = parse(Float64,metadata_keep[:airmass])
    end
    push!(df, metadata_keep)
    return df
end




""" Read HARPS wavelength file.
    Note:HARPS wavelengths are stored in Float-32 format, which reduces wavelength precision.
    This function restores that precision using a smooth polynomial fit """
function read_λ(fn::String)
    global harps_data_path
    λ_fits = FITS(harps_data_path*"/waves/"*fn)
    λ_original = convert(Array{Float64,2},FITSIO.read(λ_f[1]))
    nx, norder = size(λ_original)
    λ_new = zeros(nx,norder)
    for i in 1:norder
        waves_fit = Polynomials.fit(1:nx,λ_original[:,i],2)
        λ_new[:,i] = waves_fit.(1:nx)
    end
    return λ_new
end

function read_λ(fn::String, orders_to_read::AR ) where  { AR<:AbstractRange }
    global harps_data_path
    λ_fits = FITS(harps_data_path*"/waves/"*fn)
    λ_original = convert(Array{Float64,2},FITSIO.read(λ_f[1],:,orders_to_read))
    nx, norder = size(λ_original)
    λ_new = zeros(nx,norder)
    for i in 1:norder
        waves_fit = Polynomials.fit(1:nx,λ_original[:,i],2)
        λ_new[:,i] = waves_fit.(1:nx)
    end
    return λ_new
end

""" Read HARPS data from FITS file, and return in a Spectra2DBasic object."""
function read_data   end

function read_data(f::FITS, metadata::Dict{Symbol,Any} )
    @assert length(f) >= 2
    @assert all(map(i->typeof(f[i]) <: FITSIO.ImageHDU,2:length(f)))
    λ = read_λ(metadata[:wfile])
    flux = FITSIO.read(f[1])
    var = NaNMath.sqrt.(flux)
    metadata[:normalization] = :raw
    spectrum = Spectra2DBasic(λ, flux, var, HARPS2D(), metadata=metadata)
    apply_doppler_boost!(spectrum,metadata)
    return spectrum
end


function read_data(f::FITS, metadata::Dict{Symbol,Any}, orders_to_read::AR ) where  { AR<:AbstractRange }
    @assert length(f) >= 2
    @assert all(map(i->typeof(f[i]) <: FITSIO.ImageHDU,2:length(f)))
    λ = read_λ(metadata[:wfile], orders_to_read)
    flux = FITSIO.read(f[1],:,orders_to_read)
    var = NaNMath.sqrt.(flux)
    metadata[:normalization] = :raw
    spectrum = Spectra2DBasic(λ, flux, var, HARPS2D(), metadata=metadata)
    apply_doppler_boost!(spectrum,metadata)
    return spectrum
end


function read_data(fn::String, metadata::Dict{Symbol,Any} )
    f = FITS(fn)
    read_data(f, metadata)
end

function read_data(fn::String)
    f = FITS(fn)
    hdr = FITSIO.read_header(f[1])
    metadata = Dict(zip(map(k->Symbol(k),hdr.keys),hdr.values))
    read_data(f, metadata)
end

function read_data(dfr::DataFrameRow{DataFrame,DataFrames.Index})
    fn = dfr.Filename
    metadata = Dict(zip(keys(dfr),values(dfr)))
    read_data(fn,metadata)
end

function read_data(fn::String, orders_to_read::AR ) where  { AR<:AbstractRange }
    f = FITS(fn)
    metadata = read_metadata(f)
    read_data(f,metadata, orders_to_read)
end

function read_data(dfr::DataFrameRow{DataFrame,DataFrames.Index}, orders_to_read::AR ) where  { AR<:AbstractRange }
    fn = dfr.Filename
    metadata = Dict(zip(keys(dfr),values(dfr)))
    read_data(fn,metadata, orders_to_read)
end

function apply_barycentric_correction!(λ::AbstractArray{T,1}, z::Real) where { T<:Real }

end

""" Read normalization function for approximate continuum+blaze correction
TODO: Figure out mapping between these orders and what's stored in the FITS file before using. """
function read_normalization(fn::String = joinpath(pkgdir(EchelleInstruments),"data","neid","HARPS_normalization_function.csv") )
    CSV.read(fn, DataFrame)
end
