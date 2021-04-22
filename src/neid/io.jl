"""
   IO functions for the NEID spectrograph
   https://neid.psu.edu/
Author: Eric Ford and collaborators
Created: August 2020
"""

"""Create Dataframe containing filenames and key data for all files neid*.fits in directory"""
function make_manifest(data_path::String ; max_spectra_to_use::Int = 1000 )
    df_filenames = EchelleInstruments.make_manifest(data_path,r"^neidL1_\d+[T\.]\d+\.fits$")
    #=
    dir_filelist = readdir(data_path,join=true)
    idx_spectra = map(fn->occursin(r"^neidL1_\d+[T\.]\d+\.fits$", last(split(fn,'/')) ),dir_filelist)
    spectra_filelist = dir_filelist[idx_spectra]
    =#
    #=
    df_files = DataFrame(Filename = String[], target = String[], bjd = Float64[], ssbz=Float64[] )
    map(fn->add_metadata_from_fits!(df_files,fn),spectra_filelist)
    df_files
    =#
    #@assert length(spectra_filelist) >= 1
    #df_files = DataFrame(read_metadata(spectra_filelist[1]))
    df_files = DataFrame(read_metadata(df_filenames.Filename[1]))
    keys = propertynames(df_files)
    allowmissing!(df_files, keys[map(k->k∉[:Filename, :bjd, :target],keys)] )

    #if length(spectra_filelist) >= 2
    if length(df_filenames.Filename) >= 2
        map(fn->add_metadata_from_fits!(df_files,fn),df_filenames.Filename[2:end])
    end
    #=
    for i in 2:length(spectra_filelist)
        try
            add_metadata_from_fits!(df_files,spectra_filelist[i])

        catch
            println("# Problem reading ", spectra_filelist[i])
            fields_to_save = metadata_symbols_default(NEID2D())
            fields_str_to_save = metadata_strings_default(NEID2D())

            metadata_tmp = read_metadata_from_fits(spectra_filelist[i],fields=fields_to_save,fields_str=fields_str_to_save)
            println("Before: ", metadata_tmp)
            #println("After : ", metadata_tmp)
            #push!(df_files,metadata_tmp)
        end
    end
    =#
    df_files

end

"""Create Dict containing filename and default metadata from file."""
function read_metadata(fn::Union{String,FITS})
    fields_to_save = metadata_symbols_default(NEID2D())
    fields_str_to_save = metadata_strings_default(NEID2D())
    dict = read_metadata_from_fits(fn,fields=fields_to_save,fields_str=fields_str_to_save)
    check_metadata_fields_expected_present(dict,fields_to_save)
    dict[:Filename] = fn
    return dict
end

function add_metadata_from_fits!(df::DataFrame, fn::String)
    metadata_keep = read_metadata(fn)
    #if metadata_keep[:target] != "Solar" return df end
    calibration_target_substrings = ["Etalon","LFC","ThArS","LDLS","FlatBB"]

   if any(map(ss->contains( metadata_keep[:target],ss),calibration_target_substrings)) 
      return df
   end


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
    push!(df, metadata_keep)
    return df
end


""" Read NEID data from FITS file, and return in a Spectra2DBasic object."""
function read_data   end

function read_data(f::FITS, metadata::Dict{Symbol,Any} )
    @assert length(f) >= 2
    @assert all(map(i->typeof(f[i]) <: FITSIO.ImageHDU,2:length(f)))
    img_idx = Dict(map(i->first(read_key(f[i],"EXTNAME"))=>i,2:length(f)))
    λ, flux, var  = FITSIO.read(f[img_idx["SCIWAVE"]]), FITSIO.read(f[img_idx["SCIFLUX"]]), FITSIO.read(f[img_idx["SCIVAR"]])
    metadata[:normalization] = :raw
    spectrum = Spectra2DBasic(λ, flux, var, NEID2D(), metadata=metadata)
    apply_doppler_boost!(spectrum,metadata)
    return spectrum
end


function read_data(f::FITS, metadata::Dict{Symbol,Any}, orders_to_read::AR ) where  { AR<:AbstractRange }
    @assert length(f) >= 2
    @assert all(map(i->typeof(f[i]) <: FITSIO.ImageHDU,2:length(f)))
    img_idx = Dict(map(i->first(read_key(f[i],"EXTNAME"))=>i,2:length(f)))
    λ, flux, var  = FITSIO.read(f[img_idx["SCIWAVE"]],:,orders_to_read), FITSIO.read(f[img_idx["SCIFLUX"]],:,orders_to_read), FITSIO.read(f[img_idx["SCIVAR"]],:,orders_to_read)
    metadata[:normalization] = :raw
    spectrum = Spectra2DBasic(λ, flux, var, NEID2D(), metadata=metadata)
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

""" Read space delimited file with differential extinction corrections, interpolate to bjd's in df and insert into df[:,diff_ext_rv]. """
function read_differential_extinctions!(fn::String, df::DataFrame, df_time_col::Symbol = :bjd)
    df_diff_ext = CSV.read(fn, DataFrame, delim='\t')
    @assert any(isequal(:JD),propertynames(df_diff_ext))
    @assert any(isequal(:delta_vr),propertynames(df_diff_ext))
    diff_ext = LinearInterpolation(df_diff_ext.JD, df_diff_ext.delta_vr)
    df[!,:diff_ext_rv] = -diff_ext(df[!,df_time_col])/100 # cm/s -> m/s
    return df
end

""" Read normalization function for approximate solar continuum+blaze correction
TODO: Figure out mapping between these orders and what's stored in the FITS file before using. """
function read_solar_normalization(fn::String = joinpath(pkgdir(EchelleInstruments),"data","neid","NEID_normalization_function.csv") )
    CSV.read(fn, DataFrame)
end

""" Read calibrator data """
function read_cal_data(f::FITS, metadata::Dict{Symbol,Any}, orders_to_read::AR )  where { AR<:AbstractRange  }
    @assert length(f) >= 2
    @assert all(map(i->typeof(f[i]) <: FITSIO.ImageHDU,2:length(f)))
    img_idx = Dict(map(i->first(read_key(f[i],"EXTNAME"))=>i,2:length(f)))
    λ, flux, var  = FITSIO.read(f[img_idx["CALWAVE"]],:,orders_to_read), FITSIO.read(f[img_idx["CALFLUX"]],:,orders_to_read), FITSIO.read(f[img_idx["CALVAR"]],:,orders_to_read)
    metadata[:normalization] = :raw
    spectrum = Spectra2DBasic(λ, flux, var, NEID2D(), metadata=metadata)
    apply_doppler_boost!(spectrum,metadata)
    return spectrum
end

function read_cal_data(fn::String, metadata::Dict{Symbol,Any} )
    f = FITS(fn)
    read_cal_data(f, metadata)
end

function read_cal_data(fn::String)
    f = FITS(fn)
    hdr = FITSIO.read_header(f[1])
    metadata = Dict(zip(map(k->Symbol(k),hdr.keys),hdr.values))
    read_cal_data(f, metadata)
end
