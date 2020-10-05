"""
   Delegates loading functions & traits for the NEID spectrograph
   https://neid.psu.edu/
Author: Eric Ford and collaborators
Created: August 2020
"""

"""
Module providing types and traits and customized functions for the NEID Spectrograph.
"""
module NEID
using RvSpectMLBase
import RvSpectMLBase: AbstractInstrument, AbstractInstrument1D, AbstractInstrument2D
import ..EchelleInstruments
import ..EchelleInstruments: read_manifest, read_metadata_from_fits, check_metadata_fields_expected_present
#import ..EchelleInstruments: read_header
import ..EchelleInstruments: default_Δv_to_avoid_tellurics

using CSV, DataFrames, Query
using FITSIO
using Interpolations

""" Trait for 1D Extracted spectra from NEID """
struct NEID1D <: AbstractInstrument1D end

""" Trait for 2D Extracted spectra from NEID """
struct NEID2D <: AbstractInstrument2D end

""" Trait to specify any 1D or 2D Extracted spectra from NEID """
const AnyNEID = Union{NEID1D,NEID2D}

export NEID, NEID1D, NEID2D, AnyNEID

# traits.jl imports from RvSpectMBase on its own
include("traits.jl")
export min_order, max_order, min_pixel_in_order, max_pixel_in_order
export orders_to_use_default, min_col_default, max_col_default
export metadata_symbols_default, metadata_strings_default
export default_ccf_mask_v_width
export get_inst_module

include("io.jl")
# not exported, since don't have parameter that specializes to instrument

#import RvSpectMLBase: continuum_normalize_spectrum, continuum_normalize_spectra  # Not yet implemented
#import RvSpectMLBase: filter_line_list
#import RvSpectMLBase: find_worst_telluric_in_each_chunk # Not yet implemented
import RvSpectMLBase: make_clean_line_list_from_tellurics
include("util.jl") # TODO
#export continuum_normalize_spectrum!, continuum_normalize_spectra! # Not yet implemented
export filter_line_list
export make_clean_line_list_from_tellurics

end
