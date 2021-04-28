function filter_line_list(df::DataFrame, inst::IT ; λmin::Real = default_λmin, λmax::Real = default_λmax ) where { IT<:HARPS.AnyHARPS }
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

#=
#old neid functions

function make_λ_list_for_bad_columns(line_list::DataFrame, neid_data::DT )  where {
            T1<:Real, A1<:AbstractArray{T1}, T2<:Real, A2<:AbstractArray{T2}, T3<:Real, A3<:AbstractArray{T3},
             IT<:NEID.AnyNEID, ST<:Spectra2DBasic{T1,T2,T3,A1,A2,A3,IT}, DT<:AbstractArray{ST,1} }
    order_list = 1:size(first(neid_data).λ,2)
    df_bad_col_λs = DataFrame(:order=>Int[], :lambda_lo=>Float64[], :lambda_hi=>Float64[])
    for order in order_list
        λ_extrema = NaNMath.extrema(first(neid_data).λ[:,order])
        if isnan(first(λ_extrema)) ||  isnan(last(λ_extrema))   continue   end
        for bcr in bad_col_ranges(NEID2D(), order)
            pixlo = first(bcr)
            pixhi = pixlo + 1
            if pixhi > size(first(neid_data).λ,1)
                pixhi -= 1
                pixlo -= 1
            end
            doppler_factors = map(obsid-> haskey(neid_data[obsid].metadata,:doppler_factor) ? neid_data[obsid].metadata[:doppler_factor] : 1 , 1:length(neid_data))
            Δλ_pixel = (first(neid_data).λ[pixhi,order] - first(neid_data).λ[pixlo,order]) * doppler_factors[1]
            (λ_lo, λ_hi) = mapreduce(obsid->extrema(neid_data[obsid].λ[bcr,order]) .* doppler_factors[obsid],
                               (a,b) -> (min(a[1],b[1]), max(a[2],b[2])),   1:length(neid_data) )
            λ_lo -= Δλ_pixel/2
            λ_hi += Δλ_pixel/2
            push!(df_bad_col_λs, Dict(:order=>order, :lambda_lo=>λ_lo, :lambda_hi=>λ_hi) )
        end
    end
    return df_bad_col_λs
end

function make_good_orders_pixels_df(neid_data::DT ; orders::A4 = orders_to_use_default(first(neid_data).inst) ) where {
                T1<:Real, A1<:AbstractArray{T1}, T2<:Real, A2<:AbstractArray{T2}, T3<:Real, A3<:AbstractArray{T3},
                 IT<:NEID.AnyNEID, ST<:Spectra2DBasic{T1,T2,T3,A1,A2,A3,IT}, DT<:AbstractArray{ST,1},
                 T4 <: Integer, A4<:AbstractArray{T4} }
    @warn "Deprecated use make_good_orders_pixels_df(inst) instead."
    df = DataFrame(order=Int[],pixels=UnitRange[])
    for ord in orders
        pixel_array = EchelleInstruments.calc_complement_index_ranges(get_pixel_range(NEID2D(),ord),EchelleInstruments.NEID.bad_col_ranges(NEID2D(),ord))
        order_array = repeat([ord],length(pixel_array))
        df_tmp = DataFrame(:order=>order_array,:pixels=>pixel_array)
        append!(df,df_tmp)
    end
    return df
end

function make_good_orders_pixels_df(inst::IT ; orders::A4 = orders_to_use_default(inst) ) where {
                 IT<:NEID.AnyNEID, T4 <: Integer, A4<:AbstractArray{T4} }
    df = DataFrame(order=Int[],pixels=UnitRange[])
    for ord in orders
        pixel_array = EchelleInstruments.calc_complement_index_ranges(get_pixel_range(NEID2D(),ord),EchelleInstruments.NEID.bad_col_ranges(NEID2D(),ord))
        order_array = repeat([ord],length(pixel_array))
        df_tmp = DataFrame(:order=>order_array,:pixels=>pixel_array)
        append!(df,df_tmp)
    end
    return df
end
=#

#=
function make_clean_line_list_from_bad_columns(line_list::DataFrame, neid_data::DT ;
                Δv_to_avoid_tellurics::Real = 0.0, v_center_to_avoid_tellurics::Real = 0.0
                )  where { T1<:Real, A1<:AbstractArray{T1}, T2<:Real, A2<:AbstractArray{T2}, T3<:Real, A3<:AbstractArray{T3}, IT<:NEID.AnyNEID, ST<:Spectra2DBasic{T1,T2,T3,A1,A2,A3,IT}, DT<:AbstractArray{ST,1} }

                if size(df_bad_col_λs,1) < 1 return line_list end
                line_list_to_search_for_tellurics = copy(line_list)
                line_list_to_search_for_tellurics.lambda_lo = line_list_to_search_for_tellurics.lambda./calc_doppler_factor(Δv_to_avoid_tellurics).*calc_doppler_factor(v_center_to_avoid_tellurics)
                line_list_to_search_for_tellurics.lambda_hi = line_list_to_search_for_tellurics.lambda.*calc_doppler_factor(Δv_to_avoid_tellurics).*calc_doppler_factor(v_center_to_avoid_tellurics)

                line_list_no_tellurics_df = line_list_to_search_for_tellurics |> # @filter(λmin <= _.lambda <= λmax) |>
                   @filter( !RvSpectMLBase.is_in_wavelength_range_list(_.lambda, order=_.order, list=df_bad_col_λs ) ) |>
                   DataFrame
end
=#

#=
#old neid functions


"""
   `make_clean_line_list_from_tellurics(line_list, neid_data; Δv_to_avoid_tellurics )`
Returns a new line list that excludes lines with telluric contamination.
Inputs:
- `line_list`:  Dataframe containing field lambda
# - `neid_data`:  Array of spectra
- `Δv_to_avoid_tellurics`:  in m/s
- `v_center_to_avoid_tellurics` : in m/s
- `tellurics_filename`:  filename with wavelength ranges affected by tellurics
Outputs:
- `line_list_without_tellurics`:   DataFrame with fields: `lambda`, `weight`, `lambda_lo`, and `lambda_hi`.
Warning: Currently, reads in a list of wavelength ranges with tellurics and this file was created based on tellurics metadata in EXPRES data for HR 101501.
"""
function make_clean_line_list_from_tellurics(line_list::DataFrame, neid_data::DT ;
    Δv_to_avoid_tellurics::Real = 0.0, v_center_to_avoid_tellurics::Real = 0.0, tellurics_filename::String = joinpath(pkgdir(EchelleInstruments),"data/neid/telluric_ranges.csv")
               ) where { T1<:Real, A1<:AbstractArray{T1}, T2<:Real, A2<:AbstractArray{T2}, T3<:Real, A3<:AbstractArray{T3}, IT<:NEID.AnyNEID, ST<:Spectra2DBasic{T1,T2,T3,A1,A2,A3,IT}, DT<:AbstractArray{ST,1} }
   @assert 0.5*RvSpectMLBase.max_bc_earth_rotation <= Δv_to_avoid_tellurics <= 4*RvSpectMLBase.max_bc

   if hasproperty(line_list,:lambda_lo) &&  hasproperty(line_list,:lambda_hi)
      line_list_to_search_for_tellurics = copy(line_list)
   else
      line_list_to_search_for_tellurics = RvSpectMLBase.add_line_boundaries_to_line_list(line_list, Δv_to_avoid_tellurics=Δv_to_avoid_tellurics, v_center_to_avoid_tellurics=v_center_to_avoid_tellurics )
   end
   #chunk_list_timeseries = RvSpectMLBase.make_chunk_list_timeseries_around_lines(neid_data,line_list_to_search_for_tellurics)

   if !hasproperty(line_list, :order)
       order_info = get_order_info(neid_data) # , orders_to_use=orders_to_use)
       #println("# order_info contains: ", names(order_info))
       #println("# line_list_df contains: ", names(line_list_df))
       line_list_to_search_for_tellurics = assign_lines_to_orders(line_list_to_search_for_tellurics, order_info, Δv_to_avoid_tellurics=Δv_to_avoid_tellurics, v_center=v_center_to_avoid_tellurics )
   end

   telluric_ranges = read_telluric_ranges(tellurics_filename)   # already includes range due to BCs
   # TODO: Check sign convention for v_center_to_avoid_tellurics
   # telluric_ranges.lambda_lo .*= calc_doppler_factor(v_center_to_avoid_tellurics) # / calc_doppler_factor(Δv_to_avoid_tellurics)
   # telluric_ranges.lambda_hi .*= calc_doppler_factor(v_center_to_avoid_tellurics) # * calc_doppler_factor(Δv_to_avoid_tellurics)

   #println("# make_clean_line_list_from_tellurics: λ_list_for_bad_cols_df /= calc_doppler_factor ")
   λ_list_for_bad_cols_df = make_λ_list_for_bad_columns(line_list_to_search_for_tellurics, neid_data)
   λ_list_for_bad_cols_df.lambda_lo .*= calc_doppler_factor(v_center_to_avoid_tellurics) / calc_doppler_factor(Δv_to_avoid_tellurics)
   λ_list_for_bad_cols_df.lambda_hi .*= calc_doppler_factor(v_center_to_avoid_tellurics) * calc_doppler_factor(Δv_to_avoid_tellurics)

   line_list_no_tellurics_df = line_list_to_search_for_tellurics |> # @filter(λmin <= _.lambda <= λmax) |>
      @filter( !RvSpectMLBase.is_in_wavelength_range_list(_.lambda, order=_.order, list=λ_list_for_bad_cols_df ) ) |>
      @filter( !RvSpectMLBase.is_in_wavelength_range_list(_.lambda, list=telluric_ranges ) ) |>
      DataFrame

end
=#

#=
function add_line_boundaries_to_line_list(line_list::DataFrame; Δv_to_avoid_tellurics::Real = 0.0, v_center_to_avoid_tellurics::Real = 0.0 )
    @assert 0.5*RvSpectMLBase.max_bc_earth_rotation <= Δv_to_avoid_tellurics <= 4*RvSpectMLBase.max_bc
    line_list_to_search_for_tellurics = copy(line_list)
    line_list_to_search_for_tellurics.lambda_lo = line_list_to_search_for_tellurics.lambda./calc_doppler_factor(Δv_to_avoid_tellurics).*calc_doppler_factor(v_center_to_avoid_tellurics)
    line_list_to_search_for_tellurics.lambda_hi = line_list_to_search_for_tellurics.lambda.*calc_doppler_factor(Δv_to_avoid_tellurics).*calc_doppler_factor(v_center_to_avoid_tellurics)
    return line_list_to_search_for_tellurics
end
=#

function choose_obs_idx_for_init_guess(df::DataFrame, inst::IT )  where { IT<:HARPS.AnyHARPS }
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
