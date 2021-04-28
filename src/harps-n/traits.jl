"""
   Traits for the HARPS-N spectrograph.
    https://plone.unige.ch/HARPS-N/
Author: Alex Wise and collaborators
Created: April 2021
"""

""" Delegates loading of code specifying types essential to the package.  """

import RvSpectMLBase: min_order, max_order, min_pixel_in_order, max_pixel_in_order, min_pixel, max_pixel


min_order(::HARPSN2D) = 1
max_order(::HARPSN2D) = 69

min_pixel_in_order(inst::HARPSN2D) = 1
max_pixel_in_order(inst::HARPSN2D) = 4096

min_pixel(::HARPSN1D) = 1
max_pixel(::HARPSN1D) = 300000 # TODO: Update once know size of HARPS-N's 1d extracted spectra


import RvSpectMLBase: bad_col_ranges

function bad_col_ranges(inst::HARPSN2D, ord::Int)
    [] #placeholder function to be updated if/when we decide certain pixels should be excluded from analysis
end


import RvSpectMLBase: orders_to_use_default, min_col_default, max_col_default


orders_to_use_default(::HARPSN2D) = 1:69 

min_col_default(::HARPSN2D, ord::Integer) = 500  # Avoid where continuum normalization effected by edges/scattered light
max_col_default(::HARPSN2D, ord::Integer) = 3500  # Avoid where continuum normalization effected by edges/scattered light

import RvSpectMLBase: get_pixel_range
function get_pixel_range(inst::HARPSN2D, ord::Integer)
    minc = max(min_col_default(inst, ord))
    maxc = min(max_col_default(inst, ord))
    return minc:maxc
end

import RvSpectMLBase: metadata_symbols_default, metadata_strings_default
# HARPS-N headers
metadata_symbols_default(::AnyHARPSN) = Symbol[:bjd, :exptime, :airmass]
metadata_strings_default(::AnyHARPSN) = String["MJD-OBS", "EXPTIME", "AIRMASS"]

import RvSpectMLBase: default_ccf_mask_v_width
default_ccf_mask_v_width(::AnyHARPSN) = 820.0 #approx one HARPS-N pixel

import RvSpectMLBase: get_inst_module
get_inst_module(::AnyHARPSN) = HARPSN

#not sure what this does right now
import RvSpectMLBase: get_λ_range
function get_λ_range(data::CLT) where { T1<:Real, T2<:Real, T3<:Real, A1<:AbstractArray{T1,2}, A2<:AbstractArray{T2,2}, A3<:AbstractArray{T3,2},
                                       IT<:AnyHARPSN, CLT<:Spectra2DBasic{T1,T2,T3,A1,A2,A3,IT} }
   (λmin, λmax) = NaNMath.extrema(data.λ)
   return (min=λmin, max=λmax)
end

default_λmin = 3800.0
default_λmax = 6900.0
