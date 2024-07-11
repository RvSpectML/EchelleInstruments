"""
   Traits for the NEID spectrograph.
   https://neid.psu.edu/
Author: Eric Ford and collaborators
Created: August 2020
Updated: March 2021 for DRP v0.7
         May 2021 for DRP v1.0
         November 2021 for DRP v1.1
"""

""" Delegates loading of code specifying types essential to the package.  """

import RvSpectMLBase: min_order, max_order, min_pixel_in_order, max_pixel_in_order, min_pixel, max_pixel

#= Pre-ship
min_order(::NEID2D) = 1
max_order(::NEID2D) = 90
=#
# At KPNO
#min_order(::NEID2D) = 1
min_order(::NEID2D) = 4  # DRP 0.7
#max_order(::NEID2D) = 126 # DRS 0.5
#max_order(::NEID2D) = 122 # DRS 0.6
#max_order(::NEID2D) = 118 # DRS 0.7
max_order(::NEID2D) = 122 # DRS 1.0

min_pixel_in_order(inst::NEID2D) = 1
max_pixel_in_order(inst::NEID2D) = 9216

min_pixel(::NEID1D) = 1
max_pixel(::NEID1D) = (max_order(NEID1D())-min_order(NEID1D())+1)*9216 # TODO: Update once know size of NEID's 1d extracted spectra


import RvSpectMLBase: bad_col_ranges
#bad_col_ranges(inst::NEID2D, ord::Int) = [439:449, 1934:1943, 6714:6714] # DRP v0.6
#=
function bad_col_ranges(inst::NEID2D, ord::Int)   # DRP v0.7
    if ord == 1      return [1:1712, 1938:1938]
    elseif ord == 2  return [1:866, 1938:1938]
    elseif ord == 3  return [1:152, 438:450, 1938:1938]
    elseif 4 <= ord <= 44   return [438:450, 1938:1938]
    elseif 45 <= ord <= 46   return [438:450, 1938:1938, 6714:6714]
    elseif 47 <= ord <= 122   return [438:450, 1932:1946, 6714:6714]
    end
end
=#
#=
function bad_col_ranges(inst::NEID2D, ord::Int)   # DRP v1.1
    if ord == 1      return [1:1718, 1936:1940, 7603:9216]
    elseif ord == 2  return [1:871, 1936:1940]
    elseif ord == 3  return [1:156, 437:451, 1936:1940, 7791:9216]
    elseif ord == 4  return [1:230, 437:451, 1936:1940, 7895:9216]
    elseif ord == 5  return [1:676, 1936:1940, 8499:9216]
    elseif ord == 6  return [1:2396, 8330:9216]
    elseif ord == 7  return [1:3157, 8052:9216]
    elseif ord == 8  return [437:451, 1936:1940, 8452:9216]
    elseif ord == 9  return [437:451, 1936:1940, 5449:8302]
    elseif 10 <= ord <= 34   return [437:451, 1936:1940]
    elseif 35 <= ord <= 36   return [437:451, 1936:1940, 4584:4584]
    elseif 37 <= ord <= 44   return [437:451, 1936:1940, 4584:4584, 5404:5404]
    elseif 45 <= ord <= 45   return [437:451, 1936:1940, 6714:6714]
    elseif ord == 46         return [437:451, 1936:1940, 3199:3199, 6714:6714]
    elseif 47 <= ord <= 60   return [437:451, 1936:1948, 3199:3199, 6714:6714]
    elseif 61 <= ord <= 120  return [437:451, 1930:1948, 6714:6714]
    elseif ord == 121        return [437:451, 986:4213, 6714:6714]
    elseif ord == 122        return [1:9216]
    end
end
=#

#=
function bad_col_ranges(inst::NEID2D, ord::Int)   # DRP v1.2
    if ord == 1      return [1:1718, 1936:1940, 7603:9216]
    elseif ord == 2  return [1:871, 1936:1940]
    elseif ord == 3  return [1:156, 437:451, 1936:1940, 7791:9216]
    elseif ord == 4  return [1:230, 437:451, 1936:1940, 7895:9216]
    elseif ord == 5  return [1:676, 1936:1940, 8499:9216]
    elseif ord == 6  return [1:2396, 8330:9216]
    elseif ord == 7  return [1:3157, 8052:9216]
    elseif ord == 8  return [437:451, 1936:1940, 8452:9216]
    elseif ord == 9  return [437:451, 1936:1940, 5449:8302]
    elseif 10 <= ord <= 34   return [437:451, 1936:1940]
    elseif 35 <= ord <= 36   return [437:451, 1936:1940, 4584:4584]
    elseif 37 <= ord <= 44   return [437:451, 1936:1940, 4584:4584, 5404:5404]
    elseif 45 <= ord <= 45   return [437:451, 1936:1940, 6714:6714]
    elseif 46 <= ord <= 60   return [437:451, 1936:1948, 3199:3199, 6714:6714]
    elseif 61 <= ord <= 120  return [437:451, 1930:1948, 6714:6714]
    elseif ord == 121   return UnitRange[437:451, 1209:4876, 6714:6714]
    elseif ord == 122   return UnitRange[1:9216]
    end
end
=#

function bad_col_ranges(inst::NEID2D, ord::Int)   # DRP v1.3
    if ord == 1   return UnitRange[1:1735, 1936:1940, 7825:9216]
    elseif ord == 2   return UnitRange[1:911, 914:914, 916:916, 1936:1940, 8263:9216]
    elseif ord == 3   return UnitRange[1:223, 225:228, 437:451, 1936:1940, 7817:9216]
    elseif ord == 4   return UnitRange[437:451, 1936:1940, 7810:9216]
    elseif ord == 5   return UnitRange[1:756, 1936:1940, 9179:9216]
    elseif ord == 6   return UnitRange[1:3111, 8037:9216]
    elseif ord == 7   return UnitRange[1:3295, 8122:9216]
    elseif ord == 8   return UnitRange[437:451, 1936:1940, 5636:7896]
    elseif ord == 9   return UnitRange[437:451, 1936:1940, 5068:8143]
    elseif 10 <= ord <= 34   return [437:451, 1936:1940]
    elseif 35 <= ord <= 36   return [437:451, 1936:1940, 4584:4584]
    elseif 37 <= ord <= 44   return [437:451, 1936:1940, 4584:4584, 5404:5404]
    elseif 45 <= ord <= 45   return [437:451, 1936:1940, 6714:6714]
    elseif 46 <= ord <= 60   return [437:451, 1936:1948, 3199:3199, 6714:6714]
    elseif 61 <= ord <= 108  return [437:451, 1930:1948, 6714:6714]
    elseif ord == 109        return UnitRange[437:451, 1930:1948, 3962:3962, 6714:6714]
    elseif 110 <= ord <= 120  return [437:451, 1930:1948, 6714:6714]
    elseif ord == 121   return UnitRange[1:1503, 1930:1948, 6129:9216]
    elseif ord == 122   return UnitRange[1:9216]
    end
end


import RvSpectMLBase: orders_to_use_default, min_col_default, max_col_default

#= Pre-ship
#orders_to_use_default(inst::NEID2D) = 1:52   # Avoiding redder orders due to tellurics
orders_to_use_default(inst::NEID2D) = 1:71   # Avoiding 72 because of NaNs in solar data
min_col_default(::NEID2D, ord::Integer) = 451              # Avoiding smaller columns due to NaNs
#min_col_default(::NEID2D) = 2000              # Avoiding smaller columns due to lower flux and distortions
max_col_default(::NEID2D, ord::Integer) = 9216 - (min_col_default(NEID2D(),ord)-1)   # Avoiding larger columns for symmetry
=#
# At KPNO
#orders_to_use_default(::NEID2D) = 55:113 # everything plaussibly usable as of DRS 0.5
#orders_to_use_default(::NEID2D) = 56:111 # everything plaussibly usable as of DRS 0.6
#orders_to_use_default(::NEID2D) = 60:90  # relatively safe
orders_to_use_default(::NEID2D) = 56:98   # avoid worst of tellurics
#=
function min_col_default(::NEID2D, ord::Integer)
#    return 1445 +500 # DRS 0.6, avoiding NaN in cols 1435-1444
    if ord == 55
        return 787
    else
        return 500
    end
end
=#
#max_col_default(::NEID2D, ord::Integer) = 8429  # DRS 0.5
#max_col_default(::NEID2D, ord::Integer) = 6214  # DRS 0.6, avoiding NaN in col 6215
#max_col_default(::NEID2D, ord::Integer) = 9215  # DRS 0.6, 0.7
min_col_default(::NEID2D, ord::Integer) = 1500  # Avoid where continuum normalization effected by edges/scattered light
max_col_default(::NEID2D, ord::Integer) = 8000  # Avoid where continuum normalization effected by edges/scattered light

import RvSpectMLBase: get_pixel_range
function get_pixel_range(inst::NEID2D, ord::Integer)
    minc = max(min_col_default(inst, ord)) #, min_col_excalibur(inst,order), min_col_nonnan(inst,order))
    maxc = min(max_col_default(inst, ord)) #, max_col_excalibur(inst,order), max_col_nonnan(inst,order))
    return minc:maxc
end

import RvSpectMLBase: metadata_symbols_default, metadata_strings_default
#= Pre-ship headers
#metadata_symbols_default(::AnyNEID) = Symbol[:bjd, :target, :ssbz]
#metadata_strings_default(::AnyNEID) = String["OBSJD", "SKY-OBJ", "SSBZ000"]
=#
# Headers at KPNO
metadata_symbols_default(::AnyNEID) = Symbol[:bjd, :target, :exptime, :airmass, :ssbz,
                            :wavecal, :driftfun, :drpversion, :drpextsnr,
                            :envotmp, :envohum, :envwinds, :envwindd, :dqlevel0, :dqlevel1, :dqlevel2 ]
metadata_strings_default(::AnyNEID) = String["OBSJD", "OBJECT", "EXPTIME", "AIRMASS",  "SSBZ100",
                            "WAVECAL", "DRIFTFUN", "E_VER", "EXTSNR", 
                            "ENVOTMP", "ENVOHUM", "ENVWINDS", "ENVWINDD", "DQLEVEL0", "DQLEVEL1", "DQLEVEL2" ]

import RvSpectMLBase: default_ccf_mask_v_width
default_ccf_mask_v_width(::AnyNEID) = 620.953

import RvSpectMLBase: get_inst_module
get_inst_module(::AnyNEID) = NEID

import RvSpectMLBase: get_λ_range
function get_λ_range(data::CLT) where { T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,2}, A2<:AbstractArray{T2,2}, A3<:AbstractArray{T3,2},
                                       IT<:AnyNEID, CLT<:Spectra2DBasic{T1,T2,T3,A1,A2,A3,IT} }
   (λmin, λmax) = NaNMath.extrema(data.λ)
   return (min=λmin, max=λmax)
end

#=
# Pre-ship
default_λmin = 3950.0  # Based on HD solar data from PSU, should generalize
default_λmax = 9500.0  #
=#
# At KPNO
#default_λmin = 5040.04045429369  # Based on one FITS file, should generalize
#default_λmax = 9810.539033567005 # Based on one FITS file, should generalize
# default_λmin = 5145.0 # DRP v0.6
# default_λmax = 9822.0 # DRP v0.6
default_λmin = 3571.0 # DRP v0.7 rounded down
default_λmax = 11048.0 # DRP v0.7
