"""
   IO functions for the HARPS-N spectrograph
    https://plone.unige.ch/HARPS-N/
Author: Alex Wise and collaborators
Created: April 2021
"""

"""Create Dataframe containing filenames and key data for all files YYYY-MM-DDThh:mm:ss.mss (year-month-day T hour:minute:second:milisecond in directory"""
function make_manifest(data_path::String ; max_spectra_to_use::Int = 1000 )
    df_filenames = EchelleInstruments.make_manifest(data_path,r"^20\d{2}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}.\d{3}$")

    df_files = DataFrame(read_metadata(joinpath(df_filenames.Filename[1],"s2d.fits")))
    keys = propertynames(df_files)
    allowmissing!(df_files, keys[map(k->k∉[:Filename, :bjd, :target, :airmass],keys)] )

    if length(df_filenames.Filename) >= 2
        map(fn->add_metadata_from_fits!(df_files,joinpath(fn,"s2d.fits")),df_filenames.Filename[2:end])
    end

    RV_data = CSV.read(joinpath(data_path,"harpn_sun_release_timeseries_2015-2018.csv"), DataFrame)

    # berv = RV_data[:, [:berv]]
    berv_bary_to_helio = RV_data[:, [:berv_bary_to_helio]] .* -1.0
    rv_diff_extinction = RV_data[:, [:rv_diff_extinction]] .* -1.0

    df_files = hcat(df_files, berv_bary_to_helio, rv_diff_extinction)
    append!(keys,[:ssb_rv, :diff_ext_rv])
    rename!(df_files,keys)

    df_files

end

"""Create Dict containing filename and default metadata from file."""
function read_metadata(fn::Union{String,FITS})
    fields_to_save = metadata_symbols_default(HARPSN2D())
    fields_str_to_save = metadata_strings_default(HARPSN2D())
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


""" Read HARPS-N data from FITS file, and return in a Spectra2DBasic object."""
function read_data   end

function read_data(f::FITS, metadata::Dict{Symbol,Any} )
    @assert length(f) >= 2
    @assert all(map(i->typeof(f[i]) <: FITSIO.ImageHDU,2:length(f)))
    img_idx = Dict(map(i->first(read_key(f[i],"EXTNAME"))=>i,2:length(f)))
    λ, flux, var  = FITSIO.read(f[img_idx["WAVEDATA_AIR_BARY"]]), FITSIO.read(f[img_idx["SCIDATA"]]), FITSIO.read(f[img_idx["ERRDATA"]])
    metadata[:normalization] = :raw
    spectrum = Spectra2DBasic(λ, flux, var, HARPSN2D(), metadata=metadata)
    apply_doppler_boost!(spectrum,metadata)
    return spectrum
end


function read_data(f::FITS, metadata::Dict{Symbol,Any}, orders_to_read::AR ) where  { AR<:AbstractRange }
    @assert length(f) >= 2
    @assert all(map(i->typeof(f[i]) <: FITSIO.ImageHDU,2:length(f)))
    img_idx = Dict(map(i->first(read_key(f[i],"EXTNAME"))=>i,2:length(f)))
    λ, flux, var  = FITSIO.read(f[img_idx["WAVEDATA_AIR_BARY"]]), FITSIO.read(f[img_idx["SCIDATA"]]), FITSIO.read(f[img_idx["ERRDATA"]])
    metadata[:normalization] = :raw
    spectrum = Spectra2DBasic(λ, flux, var, HARPSN2D(), metadata=metadata)
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
function read_normalization(fn::String = joinpath(pkgdir(EchelleInstruments),"data","neid","HARPSN_normalization_function.csv") )
    CSV.read(fn, DataFrame)
end
