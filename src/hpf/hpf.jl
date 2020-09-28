"""
   Delegates loading functions & traits for the HPF spectrograph
Author: Eric Ford
Created: Sept 2020
"""

"""
Module providing types and traits and customized functions for the HPF Spectrograph.
"""
module HPF
using RvSpectMLBase
using DataFrames, FITSIO

#type EXPRES <: AbstractInstrument end
struct HPF1D <: AbstractInstrument1D end
struct HPF2D <: AbstractInstrument2D end
const AnyHPF = Union{HPF1D,HPF2D}
export HPF1D, HPF2D, AnyHPF

#include("traits.jl")
export min_order, max_order, min_pixel_in_order, max_pixel_in_order
export orders_to_use_default, min_col_default, max_col_default
export metadata_symbols_default, metadata_strings_default

#include("io.jl")
#export read_data, read_solar_data
#export make_manifest,
# read_header not exported to avoid conflict with FITSIO.read_header

#import RvSpectMLBase: get_inst_module
get_inst_module(::AnyHPF) = HPF

end
