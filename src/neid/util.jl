function filter_line_list(df::DataFrame, inst::IT ; λmin::Real = default_λmin, λmax::Real = default_λmax ) where { IT<:NEID.AnyNEID }
   df |> @filter(λmin <= _.lambda <= λmax) |>
    #    @filter( _.lambda < 6000.0 ) |>                       # Avoid tellurics at redder wavelengths
    #    @filter( _.lambda >6157 || _.lambda < 6155  ) |>   # Avoid "line" w/ large variability
    DataFrame
end

"""  read_telluric_ranges( filename )
Return DataFrame (keys lambda_lo and lambda_hi) with wavelength ranges to be avoided as containing tellurics
based on provided CSV file.  Assumes path is included in filename.
"""
function read_telluric_ranges(fn::String)
    @assert occursin(r"\.csv$",fn)
    #println("Trying to read from >",fn,"<.")
    @assert isfile(fn)
    df = CSV.read(fn, DataFrame)
    @assert size(df,1) >= 1
    @assert hasproperty(df,:lambda_lo)
    @assert hasproperty(df,:lambda_hi)
    return df
end


"""
   `make_clean_line_list_from_tellurics(line_list, neid_data; Δv_to_avoid_tellurics )`
Returns a new line list that excludes lines with telluric contamination.
Inputs:
- `line_list`:  Dataframe containing field lambda
- `neid_data`:  Array of spectra
- `Δv_to_avoid_tellurics`:  in m/s
Outputs:
- `line_list_without_tellurics`:   DataFrame with fields: `lambda`, `weight`, `lambda_lo`, and `lambda_hi`.
Warning: Currently, reads in a list of wavelength ranges with tellurics and this file was created based on tellurics metadata in EXPRES data for HR 101501.
"""
function make_clean_line_list_from_tellurics(line_list::DataFrame, neid_data::DT ;
    Δv_to_avoid_tellurics::Real = 0.0, v_center_to_avoid_tellurics::Real = 0.0, tellurics_filename::String = joinpath(pkgdir(EchelleInstruments),"data/neid/telluric_ranges.csv")
               ) where { T1<:Real, A1<:AbstractArray{T1}, T2<:Real, A2<:AbstractArray{T2}, T3<:Real, A3<:AbstractArray{T3}, IT<:NEID.AnyNEID, ST<:Spectra2DBasic{T1,T2,T3,A1,A2,A3,IT}, DT<:AbstractArray{ST,1} }
   @assert 0.5*RvSpectMLBase.max_bc_earth_rotation <= Δv_to_avoid_tellurics <= 4*RvSpectMLBase.max_bc
   line_list_to_search_for_tellurics = copy(line_list)
   line_list_to_search_for_tellurics.lambda_lo = line_list_to_search_for_tellurics.lambda./calc_doppler_factor(Δv_to_avoid_tellurics).*calc_doppler_factor(v_center_to_avoid_tellurics)
   line_list_to_search_for_tellurics.lambda_hi = line_list_to_search_for_tellurics.lambda.*calc_doppler_factor(Δv_to_avoid_tellurics).*calc_doppler_factor(v_center_to_avoid_tellurics)
   println("Hello")
   #chunk_list_timeseries = RvSpectMLBase.make_chunk_list_timeseries_around_lines(neid_data,line_list_to_search_for_tellurics)

   telluric_ranges = read_telluric_ranges(tellurics_filename)
   # TODO: Check sign convention for v_center_to_avoid_tellurics
   #telluric_ranges.lambda_lo .*= calc_doppler_factor(v_center_to_avoid_tellurics) / calc_doppler_factor(Δv_to_avoid_tellurics)
   #telluric_ranges.lambda_hi .*= calc_doppler_factor(v_center_to_avoid_tellurics) * calc_doppler_factor(Δv_to_avoid_tellurics)

   line_list_no_tellurics_df = line_list_to_search_for_tellurics |> # @filter(λmin <= _.lambda <= λmax) |>
      @filter( !RvSpectMLBase.is_in_wavelength_range_list(_.lambda, list=telluric_ranges ) ) |>
      DataFrame
end

#= expres
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
=#



function choose_obs_idx_for_init_guess(df::DataFrame, inst::IT )  where { IT<:NEID.AnyNEID }
   @assert size(df,1) >= 1
   return 1
   # TODO: Replace with something that estimates signal to noise
   #=
   @assert hasproperty(df,:snr_prelim)
   idx = findmax(df.snr_prelim)[2]
   @assert 1 <= idx <= size(df,1)
   return idx
   =#
end
