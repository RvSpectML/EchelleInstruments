function calc_mean_Δv(λ::AbstractVector{T1}, λ_min::Real, λ_max::Real ) where { T1<:Real }
    idx_min = searchsortedfirst(λ,λ_min)
    idx_max = searchsortedlast(λ,λ_max)
    @assert 1 <= idx_min <= length(λ)
    @assert 1 <= idx_max <= length(λ)
    Δv = log(λ[idx_max]/λ[idx_min])/(idx_max-idx_min) * RvSpectMLBase.speed_of_light_mps
end


""" `find_ranges_with_tellurics( spectra )`
Return Dataframe of non-overlapping, sorted wavelength ranges that were marked as having tellurics in any of provided spectra.
Calls instrument-specific `find_ranges_with_tellurics` on each spectrum.
"""
function find_ranges_with_tellurics(spectra::AS; min_order::Integer = 1, max_order::Integer = size(first(spectra).λ,2), min_Δv_clean::Real = default_min_Δv_clean, telluric_threshold::Real = 1 ) where { ST<:AbstractSpectra, AS<:AbstractArray{ST} }
    inst_mod = get_inst_module(first(spectra).inst)
    list_of_ranges = mapreduce(s->inst_mod.find_ranges_with_tellurics(s, min_order=min_order, max_order=max_order, min_Δv_clean=min_Δv_clean, telluric_threshold=telluric_threshold), append!, spectra)
    #return list_of_ranges
    if size(list_of_ranges,1)<=1   return list_of_ranges   end
    #sort!(list_of_ranges, by=x->x[1] )
    sort!(list_of_ranges, :lambda_lo )
    non_overlapping_ranges = merge_sorted_wavelength_ranges(list_of_ranges, min_Δv_clean=min_Δv_clean)
    return non_overlapping_ranges
end


""" `merge_sorted_wavelength_ranges( df; min_Δv_clean )`
Return Dataframe of non-overlapping, sorted wavelength ranges (keys lambda_lo and lambda_hi) that is (approximately) the union of all the wavelength ranges in the input dataframe.
min_Δv_clean limits results included in output.
Input dataframe are assumed to be sorted by :lambda_lo
"""
function merge_sorted_wavelength_ranges(df::DataFrame; min_Δv_clean::Real = default_min_Δv_clean)
    @assert hasproperty(df, :lambda_lo)
    @assert hasproperty(df, :lambda_hi)
    @assert issorted(df[!,:lambda_lo])
    if size(df,1) <= 1    return df    end
    c_mps = RvSpectMLBase.speed_of_light_mps
    non_overlapping_ranges = DataFrame(:lambda_lo=>[], :lambda_hi=>[])
    last_range_lo = df[1,:lambda_lo]
    last_range_hi = df[1,:lambda_hi]
    for i in 2:size(df,1)
        this_range = df[i,:]
        if this_range.lambda_lo - last_range_hi > min_Δv_clean/c_mps*last_range_hi
            push!(non_overlapping_ranges, Dict(:lambda_lo=>last_range_lo, :lambda_hi=>last_range_hi))
            last_range_lo = this_range.lambda_lo
            last_range_hi = this_range.lambda_hi
        else
            last_range_hi = this_range.lambda_hi
        end
    end
    push!(non_overlapping_ranges, Dict(:lambda_lo=>last_range_lo, :lambda_hi=>last_range_hi))
    return non_overlapping_ranges
end

""" `calc_complement_wavelength_ranges( df_in ; lambda_start, lambda_stop )`
Return DataFrame with the complement of wavelength ranges in input DataFrame.
Inputs:
- `df_in`: DataFrame containing `lambda_lo` and `lambda_hi`
Optional Inputs:
- `λ_start`:   Output range should start at `λ_start` rather than first entry in df
- `λ_stop`:   Output range should end at `λ_stop` rather than first entry in df
Returns:
- `df_out`: DataFrame containing `lambda_lo` and `lambda_hi`
"""
function calc_complement_wavelength_ranges(df_in::DataFrame; λ_start::Real = df_in[1,:lambda_lo], λ_stop::Real=df_in[end,:lambda_hi] )
    @assert hasproperty(df_in,:lambda_lo)
    @assert hasproperty(df_in,:lambda_hi)
    df_out = λ_start < df_in[1,:lambda_lo] ? DataFrame( :lambda_lo=>[λ_start], :lambda_hi=>[df_in[1,:lambda_lo]] ) :
                                                     DataFrame( :lambda_lo=>[], :lambda_hi=>[] )
    append!(df_out, DataFrame( :lambda_lo=>df_in[1:end-1,:lambda_hi],  :lambda_hi=>df_in[2:end,:lambda_lo] ) )
    if  df_in[end,:lambda_hi] < λ_stop
        append!(df_out, Dict( :lambda_lo=>df_in[end,:lambda_hi], :lambda_hi=>λ_stop ) )
    end
    df_out
end

#global already_printed_stuff = 0

""" `make_ranges_without_tellurics( telluric_list ; lambda_start, lambda_stop, min_Δv, max_Δv )`
Return DataFrame with the complement of wavelength ranges in input DataFrame.
Inputs:
- `telluric_list`: DataFrame containing `lambda_lo` and `lambda_hi` for each wavelength range to be excluded due to tellurics
Optional Inputs:
- `min_Δv`:  Don't include chunks that are smaller than `min_Δv` (2*max barycentric correction)
- `max_Δv`:  Split up any chunks larger than `max_Δv` (Inf)
- `λ_start`:   Output range should start at `λ_start` rather than first entry in telluric_list
- `λ_stop`:   Output range should end at `λ_stop` rather than last entry in telluric_list
"""
function make_ranges_without_tellurics(telluric_list::DataFrame; min_Δv::Real = 2*RvSpectMLBase.max_bc, max_Δv::Real =Inf, λ_start::Real = telluric_list[1,:lambda_lo], λ_stop::Real=telluric_list[end,:lambda_hi] )
    @assert issorted(telluric_list.lambda_lo)
    c = RvSpectMLBase.speed_of_light_mps
    df_full = calc_complement_wavelength_ranges(telluric_list,λ_start=λ_start,λ_stop=λ_stop)
    df_large_chunks = df_full |> @filter((_.lambda_hi-_.lambda_lo)*2/(_.lambda_hi+_.lambda_lo) >= min_Δv/c) |> DataFrame
    df_nice_sized_chunks =  df_large_chunks |> @filter((_.lambda_hi-_.lambda_lo)*2/(_.lambda_hi+_.lambda_lo) <= max_Δv/c) |> DataFrame
    if max_Δv < Inf
        df_too_large_chunks = df_large_chunks |> @filter((_.lambda_hi-_.lambda_lo)*2/(_.lambda_hi+_.lambda_lo) > max_Δv/c) |> DataFrame
        for row in eachrow(df_too_large_chunks)
            #global already_printed_stuff
            #if already_printed_stuff > 100   continue   end
            #  TODO: Pay attention to where spectrum has natural breaks (order boundaries, no lines, etc.)
            Δv_to_be_divided = (row.lambda_hi-row.lambda_lo)*2/(row.lambda_hi+row.lambda_lo) * c
            println("# row: λ= ", row.lambda_lo, " - ", row.lambda_hi, " Δv_to_be_divided= ",Δv_to_be_divided)
            Δv_to_be_divided = log(row.lambda_hi/row.lambda_lo)  * c
            num_chunks = ceil(Int, Δv_to_be_divided/max_Δv )
            Δv_new_chunks = Δv_to_be_divided/num_chunks
            println("# Δv_to_be_divided = ", Δv_to_be_divided, " num_chunks = ", num_chunks, " Δv_new_chunks= ",Δv_new_chunks)
            new_λ_lo = exp.(log(row.lambda_lo) .+ Δv_new_chunks./c .* (0:(num_chunks-1)))
            new_λ_hi = exp.(log(row.lambda_lo) .+ Δv_new_chunks./c .* (1:num_chunks))
            append!(df_nice_sized_chunks,Dict( :lambda_lo=>new_λ_lo, :lambda_hi=>new_λ_hi ) )
            #already_printed_stuff += 1
        end
        sort!(df_nice_sized_chunks, :lambda_lo )
    end
    return df_nice_sized_chunks
    #return df_large_chunks
end

function calc_complement_index_ranges(idx_all::AR1, idx_bad::AA1 ) where { AR1 <: AbstractRange, AR2 <: AbstractRange, AA1 <: AbstractArray{AR2} }
    function is_in_any_range(x::eltype(AR1), r::eltype(AA1) )
        any(map(b->b.start<=x<=b.stop, idx_bad))
    end
    is_not_in_any_range(x::eltype(AR1), r::eltype(AA1) ) = !is_in_any_range(x,r)
    idx_out = Array{typeof(idx_all),1}()
    start = 0
    stop = 0
    while stop < last(idx_all)
        start = findnext(x->is_not_in_any_range(x,idx_bad), idx_all,stop+1)
        if start == nothing
            break
        end
        stop = findnext(x->is_in_any_range(x,idx_bad), idx_all, start)
        if stop == nothing
            stop = length(idx_all)
        else
            stop -= 1
        end
        push!(idx_out,idx_all[start]:idx_all[stop])
    end
    return idx_out
end


function make_bad_pixel_list_into_ranges(x::V) where { T<:Integer, V<:AbstractVector{T} }
    if length(x) == 0
        return []
    elseif length(x) == 1
        return [x[1]:x[1]]
    else
        tmp = findall(x[2:end].-x[1:end-1].!=1)
        pix_start = vcat( x[1], x[tmp.+1] )
        pix_stop = vcat( x[tmp], x[length(x)] )
        return map(p->p[1]:p[2], zip(pix_start,pix_stop) )
    end
end
