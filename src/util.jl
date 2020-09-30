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
function merge_sorted_wavelength_ranges(df::DataFrame; min_Δv_clean::Real = default_min_Δv_clean) where { T<:Real }
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

""" `calc_complement_wavelength_ranges( df ; lambda_start, lambda_stop )`
Return DataFrame with the compliment of wavelength ranges in input DataFrame.
Optionally, specify that output ranges should start/stop beyond wavelength range of input DataFrame.
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

""" `calc_complement_wavelength_ranges( df ; lambda_start, lambda_stop )`
Return DataFrame with the compliment of wavelength ranges in input DataFrame.
Optionally, specify that output ranges should start/stop beyond wavelength range of input DataFrame.
"""
function make_ranges_without_tellurics(telluric_list::DataFrame; min_Δv::Real = 60000, λ_start::Real = telluric_list[1,:lambda_lo], λ_stop::Real=telluric_list[end,:lambda_hi] )
    df_full = calc_complement_wavelength_ranges(telluric_list,λ_start=λ_start,λ_stop=λ_stop)
    df_large_chunks = df_full |> @filter((_.lambda_hi-_.lambda_lo)*2/(_.lambda_hi+_.lambda_lo) >= min_Δv/RvSpectMLBase.speed_of_light_mps) |> DataFrame
    return df_large_chunks
end

""" `is_in_wavelength_range_list(λ; list )`
Return true if λ is between lambda_lo and lambda_hi for any row in list
"""
function is_in_wavelength_range_list(λ::Real; list::DataFrame  )
    @assert hasproperty(list, :lambda_lo)
    @assert hasproperty(list, :lambda_hi)
    idx =  searchsortedfirst(list[:,:lambda_hi], λ)
    return idx>size(list,1) || !(list[idx,:lambda_lo]<=λ<=list[idx,:lambda_hi]) ?  false : true
end


""" `is_in_wavelength_range_list(λ_lo, λ_hi; list )`
Return true if there is overlap between (λ_lo, λ_hi) and lambda_lo and lambda_hi for any row in list
# TODO: test
"""
function is_in_wavelength_range_list(λ_lo::Real, λ_hi::Real; list::DataFrame  )
    @assert λ_lo < λ_hi
    @assert hasproperty(list, :lambda_lo)
    @assert hasproperty(list, :lambda_hi)
    idx =  searchsortedfirst(list[:,:lambda_hi], λ_lo)
    if idx>size(list,1)    return false  end
    if λ_lo<=list[idx,:lambda_hi] &&  λ_hi>=list[idx,:lambda_lo]
        return true
    else
        return false
    end
end

""" `break_chunk_into_chunks_without_tellurics( chunk, df_telluric_ranges )`
WIP
"""
function break_chunk_into_chunks_without_tellurics(chunk::AbstractChunkOfSpectrum, tellurics::DataFrame; min_chunk_length::Integer = 6)
    telluric_mask = is_in_wavelength_range_list.(chunk.λ, list=tellurics)
    new_chunks = Vector{typeof(chunk)}(undef,0)
    idx_start = findfirst(.!telluric_mask)
    idx_stop = findlast(.!telluric_mask)
    if isnothing(idx_start) || isnothing(idx_stop)
        return new_chunks
    end
    idx_lo = idx_start
    idx_hi = idx_start
    while idx_hi < idx_stop
        idx_lo = findfirst( view(telluric_mask,idx_hi:idx_stop) )
        if isnothing(idx_lo)
            break
        end
        idx_lo += idx_lo-1
        idx_hi = findfirst( view(telluric_mask,idx_lo:idx_stop) )
        if isnothing(idx_hi)
            idx_hi = idx_stop
        else
            idx_hi += idx_lo-2
        end
        if idx_hi-idx_lo > min_chunk_length # TODO Make min_chunk_length
            @warn "Not implemented yet"
            push!(new_chunks,ChunkOfSpectrum())
        end
    end
    return new_chunks
end
