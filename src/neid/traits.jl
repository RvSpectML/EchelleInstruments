"""
   Traits for the NEID spectrograph
   https://neid.psu.edu/
Author: Eric Ford and collaborators
Created: August 2020
"""

""" Delegates loading of code specifying types essential to the package.  """

import RvSpectMLBase: min_order, max_order, min_pixel_in_order, max_pixel_in_order, min_pixel, max_pixel

min_order(::NEID2D) = 1
max_order(::NEID2D) = 90
min_pixel_in_order(inst::NEID2D) = 1
max_pixel_in_order(inst::NEID2D) = 9216

min_pixel(::NEID1D) = 1
max_pixel(::NEID1D) = 90*9216 # TODO: Update once know size of NEID's 1d extracted spectra

import RvSpectMLBase: orders_to_use_default, min_col_default, max_col_default
#orders_to_use_default(inst::NEID2D) = 1:52   # Avoiding redder orders due to tellurics
orders_to_use_default(inst::NEID2D) = 1:71   # Avoiding 72 because of NaNs in solar data
min_col_default(::NEID2D, ord::Integer) = 451              # Avoiding smaller columns due to NaNs
#min_col_default(::NEID2D) = 2000              # Avoiding smaller columns due to lower flux and distortions
max_col_default(::NEID2D, ord::Integer) = 9216 - (min_col_default(NEID2D(),ord)-1)   # Avoiding larger columns for symmetry

import RvSpectMLBase: get_pixel_range
function get_pixel_range(inst::NEID2D, ord::Integer)
    minc = max(min_col_default(inst, ord)) #, min_col_excalibur(inst,order), min_col_nonnan(inst,order))
    maxc = min(max_col_default(inst, ord)) #, max_col_excalibur(inst,order), max_col_nonnan(inst,order))
    return minc:maxc
end

import RvSpectMLBase: metadata_symbols_default, metadata_strings_default
metadata_symbols_default(::AnyNEID) = Symbol[:bjd, :target, :ssbz]
metadata_strings_default(::AnyNEID) = String["OBSJD", "SKY-OBJ", "SSBZ000"]

import RvSpectMLBase: default_ccf_mask_v_width
default_ccf_mask_v_width(::AnyNEID) = 620.953

import RvSpectMLBase: get_inst_module
get_inst_module(::AnyNEID) = NEID

import RvSpectMLBase: get_λ_range
function get_λ_range(data::CLT) where { T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,2}, A2<:AbstractArray{T2,2}, A3<:AbstractArray{T3,2},
                                       IT<:AnyNEID, CLT<:Spectra2DBasic{T1,T2,T3,A1,A2,A3,IT} }
   (λmin, λmax) = extrema(data.λ)
   return (min=λmin, max=λmax)
end

default_λmin = 3950.0  # Based on HD solar data from PSU, should generalize
default_λmax = 9500.0  #
