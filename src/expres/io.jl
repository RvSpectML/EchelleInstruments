"""
   IO functions for the EXPRES spectrograph
      http://exoplanets.astro.yale.edu/expresBlog/
      https://ui.adsabs.harvard.edu/abs/2016SPIE.9908E..6TJ/abstract
Author: Eric Ford
Created: August 2020
"""

"""Create Dataframe containing filenames and key data for all files neid*.fits in directory"""
function make_manifest(data_path::String)
    dir_filelist = readdir(data_path,join=true)
    idx_spectra = map(fn->occursin(r"^\d+_\d+\.\d+\.fits$", last(split(fn,'/')) ),dir_filelist)
    spectra_filelist = dir_filelist[idx_spectra]
    @assert length(spectra_filelist) >= 1
    df_files = DataFrame(read_metadata(spectra_filelist[1]))
    if length(spectra_filelist) >= 2
        map(fn->add_metadata_from_fits!(df_files,fn),spectra_filelist[2:end])
    end
    # Date transition Based on EXPRES webpage
    df_files[!,:expres_epoch] = map(b-> b ? 5 : 4, df_files[!,:bjd] .> EXPRES.jd2mjd(datetime2julian(DateTime(2019,8,4))))
    df_files
end

"Return modified Julian date based on input Julian date"
function jd2mjd(jd::Real)
    @assert jd > 2400000.5  # There's no EPRV data from that long ago!
    mjd = jd - 2400000.5
    return mjd
end

"""Create Dict containing filename and default metadata from FITS file headers."""
function read_metadata(f::FITS, filename::String)
    fields_to_save = metadata_symbols_default(EXPRES2D())
    fields_str_to_save = metadata_strings_default(EXPRES2D())
    dict1 = read_metadata_from_fits(f,hdu=1,fields=fields_to_save,fields_str=fields_str_to_save)
    dict1[:Filename] = filename
    dict2 = read_metadata_from_fits(f,hdu=2,fields=fields_to_save,fields_str=fields_str_to_save)
    dict3 = read_metadata_from_fits(f,hdu=3,fields=fields_to_save,fields_str=fields_str_to_save)
    dict = merge(dict1,dict2,dict3)
    check_metadata_fields_expected_present(dict,fields_to_save)
    return dict
end

"""Create Dict containing filename and default metadata from FITS file headers."""
function read_metadata(fn::String)
    f = FITS(fn)
    read_metadata(f, fn)
end

""" `add_metadata_from_fits!(df, filename)`
Return updated dataframe after adding metadata from FITS file header."""
function add_metadata_from_fits!(df::DataFrame, fn::String)
    metadata_keep = read_metadata(fn)
    push!(df, metadata_keep)
    return df
end

#=
import NaNMath

function running_mean(x::AbstractArray{T1,1}; smooth_length::Integer ) where { T1 <: Real }
    n = length(x)
    smooth = zeros(T1,n)
    for i in 1:n
        v = view(x,max(1,i-smooth_length):min(i+smooth_length,n))
        smooth[i] = NaNMath.sum(v) / (2*smooth_length+1)
        #=
        num = zero(T1)
        denom = zero(T1)
        for y in v
            if y != NaN
                num += y
                denom += 1
            end
        end
        smooth[i] = num/denom
        =#
    end
    return smooth
end

=#

default_smooth_blaze_degree = 8
default_smooth_blaze_length = 60

function smooth_blaze(x::AbstractArray{T1,1}; degree::Integer = default_smooth_blaze_degree, smooth_length::Integer = default_smooth_blaze_length) where { T1 <: Real }
    poly = Polynomials.fit(1:length(x),x,degree,weights=max.(zero(x),x) )
    poly.(collect(range(1.0,length=length(x))))
end

function smooth_blaze(x::AbstractArray{T1,2}; degree::Integer = default_smooth_blaze_degree, smooth_length::Integer = default_smooth_blaze_length) where { T1 <: Real }
    return x  # WARNING bypassing this for now
    smooth_2d = fill(NaN,size(x))
    for i in 1:size(x,2)
        num_pix_nonnan = length(minmax_col_nonnan[i])
        deg = num_pix_nonnan > degree-2 ? degree : num_pix_nonnan-2
        if num_pix_nonnan > 1
            #println("# Smoothing order idx ", i, " pixels ", minmax_col_nonnan[i])
            smooth_2d[minmax_col_nonnan[i],i] .= smooth_blaze(view(x,minmax_col_nonnan[i],size(x,2)), degree=deg, smooth_length=smooth_length)
        end
    end
    smooth_2d
end


""" Read EXPRES data from FITS file, and return in a Spectra2DBasic object."""
function read_data   end

function read_data(f::FITS, metadata::Dict{Symbol,Any}; normalization::Symbol = :raw, use_excalibur::Bool = true )
    λ, spectrum, uncertainty = FITSIO.read(f["optimal"],"bary_wavelength"), FITSIO.read(f["optimal"],"spectrum"), FITSIO.read(f["optimal"],"uncertainty")
    if use_excalibur
        # For pixels where a presumably more accurate wavelength is avaliable, overwrite it.
        excalibur_mask = haskey(metadata,:excalibur_mask) ? metadata[:excalibur_mask] : FITSIO.read(f["optimal"],"excalibur_mask")
        λ_excalibur = FITSIO.read(f["optimal"],"bary_excalibur")
        λ[excalibur_mask] .= λ_excalibur[excalibur_mask]
    end
    if normalization == :blaze
        flux = spectrum
        var = uncertainty.^2
        metadata[:normalization] = :blaze
    elseif normalization == :raw
        # Restore fluxes to include the blaze function and scale uncertainties appropriately
        blaze = haskey(metadata,:blaze) ? metadata[:blaze] : FITSIO.read(f["optimal"],"blaze")
        blaze_smoothed = smooth_blaze(blaze)
        flux = spectrum.*blaze_smoothed
        # Since EXPRES pipeline returns standard deviation rather than variance
        var = (uncertainty.*blaze_smoothed).^2
        metadata[:normalization] = :raw
    elseif normalization == :continuum
        #blaze = haskey(metadata,:blaze) ? metadata[:blaze] : FITSIO.read(f["optimal"],"blaze")
        #blaze_smoothed = smooth_blaze(blaze)
        continuum = haskey(metadata,:continuum) ? metadata[:continuum] : FITSIO.read(f["optimal"],"continuum")
        flux = spectrum./continuum
        var = (uncertainty./continuum).^2
        metadata[:normalization] = :continuum
    else
        @error "# Reading data directly with normalization " * string(normalization) * " is not implemented."
    end
    Spectra2DBasic(λ, flux, var, EXPRES2D(), metadata=metadata)
end

function read_data(fn::String, metadata = Dict{Symbol,Any}();
                        normalization::Symbol = :raw, use_excalibur::Bool = true,
                        store_min_data::Bool = false, store_blaze::Bool = !store_min_data,
                        store_continuum::Bool = !store_min_data, store_tellurics::Bool = !store_min_data,
                        store_pixel_mask::Bool = !store_min_data, store_excalibur_mask::Bool = use_excalibur,
                        store_all_metadata::Bool = !store_min_data )
    # println("# Reading data from ", fn)
    f = FITS(fn)
    if store_all_metadata || isempty(metadata)
        read_metadata(f,fn)
        hdr1 = FITSIO.read_header(f[1])
        metadata1 = Dict(zip(map(k->Symbol(k),hdr1.keys),hdr1.values))
        hdr2 = FITSIO.read_header(f[2])
        metadata2 = Dict(zip(map(k->Symbol(k),hdr2.keys),hdr2.values))
        metadata_combo = merge(metadata,metadata1,metadata2)
    else
        metadata_combo = copy(metadata)
    end

    if store_blaze
        metadata_combo[:blaze] = FITSIO.read(f["optimal"],"blaze")
    end
    if store_continuum
        metadata_combo[:continuum] = FITSIO.read(f["optimal"],"continuum")
    end
    if store_tellurics
        metadata_combo[:tellurics] = FITSIO.read(f["optimal"],"tellurics")
    end
    if store_pixel_mask
        metadata_combo[:pixel_mask] = FITSIO.read(f["optimal"],"pixel_mask")
    end
    if store_excalibur_mask
        metadata_combo[:excalibur_mask] = FITSIO.read(f["optimal"],"excalibur_mask")
    end
    read_data(f,metadata_combo, normalization=normalization, use_excalibur=use_excalibur)
end

""" Read EXPRES data from FITS a file using Filename from a DataFrameRow.
    Adds other columns from DataFrameRow as metadata. """
function read_data(dfr::DataFrameRow{DataFrame,DataFrames.Index};
                            normalization::Symbol = :raw, use_excalibur::Bool = true,
                            store_min_data::Bool = false, store_blaze::Bool = !store_min_data,
                            store_continuum::Bool = !store_min_data, store_tellurics::Bool = !store_min_data,
                            store_pixel_mask::Bool = !store_min_data, store_excalibur_mask::Bool = use_excalibur,
                            store_all_metadata::Bool = !store_min_data)
    fn = dfr.Filename
    metadata = Dict(zip(keys(dfr),values(dfr)))
    read_data(fn,metadata,  normalization=normalization, use_excalibur=use_excalibur,
                            # store_min_data=store_min_data,
                            store_blaze=store_blaze, store_continuum=store_continuum, store_tellurics=store_tellurics,
                            store_pixel_mask=store_pixel_mask, store_excalibur_mask=store_excalibur_mask,
                            store_all_metadata=store_all_metadata )
end

""" Read only EXPRES data from FITS file, and leave metadata empty."""
function read_data_only(f::FITS)
    metadata = Dict{Symbol,Any}()
    read_data(f, metadata)
end

""" Read only EXPRES data from FITS file, and leave metadata empty."""
function read_data_only(fn::String)
    f = FITS(fn)
    read_data_only(f)
end
