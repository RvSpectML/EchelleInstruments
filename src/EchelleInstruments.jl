"""
Delegates loading of code with functions and parameters specific to different instruments.
Subdirectories of src/instruments include provide functions specialized for each instrument,
typically the file I/O and pre-processing, so data ends up in a common format.
src/instruments/common.jl provides routines that can be shared by instruments.
"""
module EchelleInstruments

using RvSpectMLBase
using DataFrames, CSV, Query
using FITSIO

import RvSpectMLBase: default_ccf_mask_v_width  # So this function can be overriden for each instrument

include("common.jl")   # Mostly trait functions to be specialized by instruments
#function get_inst_module end

include("files.jl")
export make_manifest, read_manifest, write_manifest

include("io.jl")
export  read_metadata_from_fits, read_fits_header
# don't export read_header due to conflict with FITSIO

include("util.jl")
# TODO: Decide what to export from util.jl
#export add_line_boundaries_to_line_list, expand_line_boundaries_in_line_list  # moved to RvSpectMLBase

include("expres/expres.jl")
import .EXPRES: EXPRES1D, EXPRES2D, AnyEXPRES #, get_inst_module
#import .EXPRES: get_inst_module, filter_line_list, find_worst_telluric_in_each_chunk
export EXPRES, EXPRES1D, EXPRES2D, AnyEXPRES

include("harps-n/harps-n.jl")
export HARPSN, HARPSN1D, HARPSN2D, AnyHARPSN
import .HARPSN: HARPSN1D, HARPSN2D, AnyHARPSN #, get_inst_module
#=
#import .HARPSN: get_inst_module, filter_line_list, find_worst_telluric_in_each_chunk
=#

include("harps/harps.jl")
export HARPS, HARPS1D, HARPS2D, AnyHARPS
import .HARPS: HARPS1D, HARPS2D, AnyHARPS #, get_inst_module
#=
#import .HARPSN: get_inst_module, filter_line_list, find_worst_telluric_in_each_chunk
=#

include("hpf/hpf.jl")
export HPF, HPF1D, HPF2D, AnyHPF
import .HPF: HPF1D, HPF2D, AnyHPF #, get_inst_module
#=
#import .HPF: get_inst_module, filter_line_list, find_worst_telluric_in_each_chunk
=#

include("neid/neid.jl")
import .NEID: NEID1D, NEID2D, AnyNEID
#import .NEID: get_inst_module, filter_line_list #, find_worst_telluric_in_each_chunk
export NEID, NEID1D, NEID2D, AnyNEID


# TODO: Add more instruments:  HARPS-N, HPF, etc.

end
