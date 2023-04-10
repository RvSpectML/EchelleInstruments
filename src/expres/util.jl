
""" Normalize spectrum based on blaze model (also tries to include black body function) from FITS file. """
function blaze_normalize_spectrum!(spectrum::ST) where { ST<:AbstractSpectra }
    @assert haskey(spectrum.metadata,:continuum)
    if spectrum.metadata[:normalization] == :raw
        blaze_smoothed = smooth_blaze(blaze)
        spectrum.flux ./= blaze_smoothed
        spectrum.var ./= (blaze_smoothed ) .^2
        spectrum.metadata[:normalization] = :blaze
    elseif spectrum.metadata[:normalization] == :continuum
        spectrum.flux .*= spectrum.metadata[:continuum]
        spectrum.var .*= (spectrum.metadata[:continuum] ) .^2
        spectrum.metadata[:normalization] = :blaze
    elseif spectrum.metadata[:normalization] == :blaze
        # do nothing
    else
        @error "Normalizing from " * string(spectrum.metadata.normalization) * " to blaze isn't implemented yet"
    end
    return spectrum
end

""" Normalize each spectrum based on blaze model (also tries to include black body function) from FITS files. """
function blaze_normalize_spectra!(spectra::AS) where { ST<:AbstractSpectra, AS<:AbstractArray{ST} }
    for spectrum in spectra
        blaze_normalize_spectrum!(spectrum)
    end
    return spectra
end

""" Normalize spectrum based on continuum model from FITS file. """
function continuum_normalize_spectrum!(spectrum::ST) where { ST<:AbstractSpectra }
    @assert haskey(spectrum.metadata,:continuum)
    if spectrum.metadata[:normalization] == :raw
        blaze_smoothed = smooth_blaze(blaze)
        spectrum.flux ./= spectrum.metadata[:continuum] .* blaze_smoothed
        spectrum.var ./= (spectrum.metadata[:continuum].* blaze_smoothed ) .^2
        spectrum.metadata[:normalization] = :continuum
    elseif spectrum.metadata[:normalization] == :blaze
        spectrum.flux ./= spectrum.metadata[:continuum]
        spectrum.var ./= (spectrum.metadata[:continuum] ) .^2
        spectrum.metadata[:normalization] = :continuum
    elseif spectrum.metadata[:normalization] == :continuum
        # do nothing
    else
        @error "Normalizing from " * string(spectrum.metadata.normalization) * " to continuum isn't implemented yet"
    end
    return spectrum
end

""" Normalize each spectrum based on continuum model from FITS files. """
function continuum_normalize_spectra!(spectra::AS) where { ST<:AbstractSpectra, AS<:AbstractArray{ST} }
    for spectrum in spectra
        continuum_normalize_spectrum!(spectrum)
    end
    return spectra
end

""" `filter_line_list( linelist, inst ; λmin, λmax )`
Return DataFrame based on linelist, but filtered to only include lines between λmin and λmax.
Defaults for λmin and λmax specified are in traits for EXPRES.
"""
function filter_line_list(df::DataFrame, inst::IT ;
                            λmin::Real = default_filter_line_list_λmin, λmax::Real = default_filter_line_list_λmax ) where {
                                   IT<:EXPRES.AnyEXPRES }
   df |> @filter(λmin <= _.lambda <= λmax) |>
    #    @filter( _.lambda < 6000.0 ) |>                       # Avoid tellurics at redder wavelengths
    #    @filter( _.lambda >6157 || _.lambda < 6155  ) |>   # Avoid "line" w/ large variability
    DataFrame
end

""" `find_worst_telluric_in_each_chunk( chunk_list_timeseries, spectra )`
Return an array of the worst (smallest) telluric for each chunk at any time for the provided ChunkListTimeseries.
"""
function find_worst_telluric_in_each_chunk( clt::AbstractChunkListTimeseries, data::AbstractArray{AS,1} )  where { T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,2}, A2<:AbstractArray{T2,2}, A3<:AbstractArray{T3,2}, IT<:AnyEXPRES, AS<:Spectra2DBasic{T1,T2,T3,A1,A2,A3,IT}  }
   num_lines = num_chunks(clt)
   num_obs = length(clt)
   min_telluric_model_one_obs = ones(num_lines, num_obs )
   min_telluric_model_all_obs  = ones(num_lines)
   for ch_idx in 1:num_chunks(clt)
       for t_idx in 1:num_obs
          view_indices = clt.chunk_list[t_idx].data[ch_idx].λ.indices
          cols = view_indices[1]
          order = view_indices[2]
          min_telluric_model_one_obs[ch_idx, t_idx] = minimum(view(data[t_idx].metadata[:tellurics], cols, order))
      end # times
      min_telluric_model_all_obs[ch_idx] = minimum(min_telluric_model_one_obs[ch_idx, :])
  end # lines
  return min_telluric_model_all_obs
end


"""
   `make_clean_line_list_from_tellurics(line_list, expres_data; Δv_to_avoid_tellurics )`
Returns a new line list that excludes lines with telluric contamination.
Inputs:
- `line_list`:  Dataframe containing field lambda
- `expres_data`:  Array of spectra
- `Δv_to_avoid_tellurics`:  in m/s
Outputs:
- `line_list_without_tellurics`:   DataFrame with fields: `lambda`, `weight`, `lambda_lo`, and `lambda_hi`.
Warning: Currently, assumes a tellurics value in metadata for each spectra, such as is provided by EXPRES.
"""
function make_clean_line_list_from_tellurics(line_list::DataFrame, expres_data::DT; Δv_to_avoid_tellurics::Real = default_Δv_to_avoid_tellurics,
            v_center_to_avoid_tellurics::Real = 0.0, telluric_threshold::Real = 1
               ) where { T1<:Real, A1<:AbstractArray{T1}, T2<:Real, A2<:AbstractArray{T2}, T3<:Real, A3<:AbstractArray{T3}, IT<:EXPRES.AnyEXPRES, ST<:Spectra2DBasic{T1,T2,T3,A1,A2,A3,IT}, DT<:AbstractArray{ST,1} }

   @assert 0.5*RvSpectMLBase.max_bc_earth_rotation <= Δv_to_avoid_tellurics <= 4*RvSpectMLBase.max_bc
   line_list_to_search_for_tellurics = copy(line_list)
   # TODO: Check sign convention for v_center_to_avoid_tellurics
   line_list_to_search_for_tellurics.lambda_lo = line_list_to_search_for_tellurics.lambda./calc_doppler_factor(Δv_to_avoid_tellurics).*calc_doppler_factor(v_center_to_avoid_tellurics)
   line_list_to_search_for_tellurics.lambda_hi = line_list_to_search_for_tellurics.lambda.*calc_doppler_factor(Δv_to_avoid_tellurics).*calc_doppler_factor(v_center_to_avoid_tellurics)
   chunk_list_timeseries = RvSpectMLBase.make_chunk_list_timeseries_around_lines(expres_data,line_list_to_search_for_tellurics)
   line_list_to_search_for_tellurics.min_telluric_model_all_obs = find_worst_telluric_in_each_chunk( chunk_list_timeseries, expres_data)
   line_list_no_tellurics_df = line_list_to_search_for_tellurics |> @filter(_.min_telluric_model_all_obs >= telluric_threshold ) |> DataFrame
end

default_min_Δv_clean = 8000.0

function find_ranges_with_tellurics_in_order(spectrum::ST, order::Integer; telluric_threshold::Real = 1, min_Δv_clean::Real = default_min_Δv_clean, verbose::Bool = false) where { T1<:Real, A1<:AbstractArray{T1}, T2<:Real, A2<:AbstractArray{T2}, T3<:Real, A3<:AbstractArray{T3}, IT<:EXPRES.AnyEXPRES, ST<:Spectra2DBasic{T1,T2,T3,A1,A2,A3,IT}  }

    @assert haskey(spectrum.metadata,:tellurics)
    telluric_ranges = Vector{Tuple{Int64,Int64}}(undef,0)

    c = RvSpectMLBase.speed_of_light_mps
    all_wavelengths = range(first(spectrum.λ[:,order]),stop=last(spectrum.λ[:,order]),length=size(spectrum.λ,1) )
    tellurics = spectrum.metadata[:tellurics]
    start = findfirst(x->!isnan(x),view(tellurics,:,order) )
    stop = findlast(x->!isnan(x),view(tellurics,:,order) )
    if isnothing(start) || isnothing(stop)
        if verbose   println("# The entire order ", order, " was NaNs!?!")   end
        return Vector{Tuple{eltype(spectrum.λ), eltype(spectrum.λ)} }(undef,0)
    end
    @assert stop >= start +1
    in_telluric = tellurics[start,order] < telluric_threshold
    idx_start_this_telluric = in_telluric ? start : 0
    idx_stop_this_telluric = 0
    for i in start+1:stop
        if in_telluric
            if !(tellurics[i,order] < telluric_threshold)  # exitted telluric
                idx_stop_this_telluric = i-1
                #println("# Order = " , order, "  Adding pixels = ", idx_start_this_telluric, " - ", idx_stop_this_telluric, " λ = ", spectrum.λ[idx_start_this_telluric,order], " - ", spectrum.λ[idx_stop_this_telluric,order] )
                push!(telluric_ranges, (idx_start_this_telluric,idx_stop_this_telluric) )
                #push!(telluric_ranges,(spectrum.λ[idx_start_this_telluric,order], spectrum.λ[idx_stop_this_telluric,order]) )
                in_telluric = false
            end
        else
            if (tellurics[i,order] < telluric_threshold)  # entered telluric
                idx_start_this_telluric = i
                if length(telluric_ranges) >= 1
                    idx_last_start = last(telluric_ranges)[1]
                    idx_last_stop = last(telluric_ranges)[2]
                    if spectrum.λ[idx_start_this_telluric,order] - spectrum.λ[idx_last_stop,order] <= min_Δv_clean/c * spectrum.λ[idx_start_this_telluric,order]
                        idx_start_this_telluric = first(last(telluric_ranges))
                        pop!(telluric_ranges)
                    end
                end
                in_telluric = true
            end
        end
    end
    if in_telluric
        idx_stop_this_telluric = stop
        #println("# Order = " , order, "  Adding pixels = ", idx_start_this_telluric, " - ", idx_stop_this_telluric, " λ = ", spectrum.λ[idx_start_this_telluric,order], " - ", spectrum.λ[idx_stop_this_telluric,order] )
        push!(telluric_ranges,(idx_start_this_telluric,idx_stop_this_telluric) )
    end
    #lambda_range = map(r->(spectrum.λ[first(r),order], spectrum.λ[last(r),order] ), telluric_ranges)
    DataFrame(:lambda_lo=>map(r->spectrum.λ[first(r),order], telluric_ranges), :lambda_hi=>map(r->spectrum.λ[last(r),order], telluric_ranges) )
end

import EchelleInstruments.merge_sorted_wavelength_ranges

function find_ranges_with_tellurics(spectrum::ST; min_order::Integer = 1, max_order::Integer = size(spectrum.λ,2), min_Δv_clean::Real = default_min_Δv_clean, telluric_threshold::Real = 1 ) where {
            T1<:Real, A1<:AbstractArray{T1}, T2<:Real, A2<:AbstractArray{T2}, T3<:Real, A3<:AbstractArray{T3}, IT<:EXPRES.AnyEXPRES, ST<:Spectra2DBasic{T1,T2,T3,A1,A2,A3,IT}  }
    @assert haskey(spectrum.metadata,:tellurics)
    c = RvSpectMLBase.speed_of_light_mps

    #list_of_ranges = Vector{Tuple{eltype(spectrum.λ),eltype(spectrum.λ)}}(undef,0)
    list_of_ranges = DataFrame()
    for ord in min_order:max_order
        ranges_for_order = find_ranges_with_tellurics_in_order(spectrum,ord, min_Δv_clean=min_Δv_clean, telluric_threshold=telluric_threshold)
        if size(ranges_for_order,1) >= 1
            append!(list_of_ranges,ranges_for_order)
        end
    end

    if size(list_of_ranges,1)<=1   return list_of_ranges   end
    sort!(list_of_ranges, :lambda_lo )
    non_overlapping_ranges = merge_sorted_wavelength_ranges(list_of_ranges,min_Δv_clean=min_Δv_clean)
    return non_overlapping_ranges
end

#import EchelleInstruments.choose_obs_idx_for_init_guess

function choose_obs_idx_for_init_guess(df::DataFrame, inst::IT ) where { IT<:EXPRES.AnyEXPRES }
   @assert size(df,1) >= 1
   @assert hasproperty(df,:snr_prelim)
   idx = findmax(df.snr_prelim)[2]
   @assert 1 <= idx <= size(df,1)
   return idx
end
