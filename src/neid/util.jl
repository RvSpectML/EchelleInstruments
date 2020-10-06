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
   @assert 100 <= Δv_to_avoid_tellurics <= 50000
   telluric_ranges = read_telluric_ranges(tellurics_filename)
   # TODO: Check sign convention for v_center_to_avoid_tellurics
   telluric_ranges.lambda_lo .*= calc_doppler_factor(v_center_to_avoid_tellurics) / calc_doppler_factor(Δv_to_avoid_tellurics)
   telluric_ranges.lambda_hi .*= calc_doppler_factor(v_center_to_avoid_tellurics) * calc_doppler_factor(Δv_to_avoid_tellurics)

   line_list |> # @filter(λmin <= _.lambda <= λmax) |>
      @filter( !RvSpectMLBase.is_in_wavelength_range_list(_.lambda, list=telluric_ranges ) ) |>
      DataFrame
end


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
