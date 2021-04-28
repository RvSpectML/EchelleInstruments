"""
   Traits for the HARPS spectrograph.
    https://www.eso.org/sci/facilities/lasilla/instruments/harps.html
Author: Alex Wise and collaborators
Created: April 2021
"""

""" Delegates loading of code specifying types essential to the package.  """

import RvSpectMLBase: min_order, max_order, min_pixel_in_order, max_pixel_in_order, min_pixel, max_pixel


min_order(::HARPS2D) = 1
max_order(::HARPS2D) = 72

min_pixel_in_order(inst::HARPS2D) = 1
max_pixel_in_order(inst::HARPS2D) = 4096

min_pixel(::HARPS1D) = 1
max_pixel(::HARPS1D) = 300000 # TODO: Update once know size of HARPS's 1d extracted spectra


import RvSpectMLBase: bad_col_ranges

function bad_col_ranges(inst::HARPS2D, ord::Int)
    []#placeholder function to be updated if/when we decide certain pixels should be excluded from analysis
end


import RvSpectMLBase: orders_to_use_default, min_col_default, max_col_default


orders_to_use_default(::HARPS2D) = 7:71   # avoid last order due to oxygen tellurics, first 7 due to poor SNR

min_col_default(::HARPS2D, ord::Integer) = 500  # Avoid where continuum normalization effected by edges/scattered light
max_col_default(::HARPS2D, ord::Integer) = 3500  # Avoid where continuum normalization effected by edges/scattered light

import RvSpectMLBase: get_pixel_range
function get_pixel_range(inst::HARPS2D, ord::Integer)
    minc = max(min_col_default(inst, ord))
    maxc = min(max_col_default(inst, ord))
    return minc:maxc
end

import RvSpectMLBase: metadata_symbols_default, metadata_strings_default
# HARPS headers
metadata_symbols_default(::AnyHARPS) = Symbol[:bjd, :target, :exptime, :airmass, :ssb_rv_kmps, :wfile]
metadata_strings_default(::AnyHARPS) = String["MJD-OBS", "OBJECT", "EXPTIME", "AIRMASS", "ESO DRS BERV", "ESO DRS CAL TH FILE"]

import RvSpectMLBase: default_ccf_mask_v_width
default_ccf_mask_v_width(::AnyHARPS) = 820.0 #approx one HARPS pixel

import RvSpectMLBase: get_inst_module
get_inst_module(::AnyHARPS) = HARPS

#not sure what this does right now
import RvSpectMLBase: get_λ_range
function get_λ_range(data::CLT) where { T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,2}, A2<:AbstractArray{T2,2}, A3<:AbstractArray{T3,2},
                                       IT<:AnyHARPS, CLT<:Spectra2DBasic{T1,T2,T3,A1,A2,A3,IT} }
   (λmin, λmax) = NaNMath.extrema(data.λ)
   return (min=λmin, max=λmax)
end

default_λmin = 3800.0
default_λmax = 6900.0
