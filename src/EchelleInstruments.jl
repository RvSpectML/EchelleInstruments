"""
Delegates loading of code with functions and parameters specific to different instruments.
Subdirectories of src/instruments include provide functions specialized for each instrument,
typically the file I/O and pre-processing, so data ends up in a common format.
src/instruments/common.jl provides routines that can be shared by instruments.
"""
module EchelleInstruments

using RvSpectMLBase

include("common.jl")   # Mostly trait functions to be specialized by instruments
#function get_inst_module end

include("files.jl")
export make_manifest

include("io.jl")
export read_manifest, read_metadata_from_fits
# don't export read_header due to conflict with FITSIO


include("neid/neid.jl")
import .NEID: NEID1D, NEID2D, AnyNEID
#import .NEID: get_inst_module, filter_line_list #, find_worst_telluric_in_each_chunk
export NEID, NEID1D, NEID2D, AnyNEID

include("expres/expres.jl")
import .EXPRES: EXPRES1D, EXPRES2D, AnyEXPRES #, get_inst_module
#import .EXPRES: get_inst_module, filter_line_list, find_worst_telluric_in_each_chunk
export EXPRES, EXPRES1D, EXPRES2D, AnyEXPRES

# TODO: Add more instruments:  HARPS-N, HPF, etc.
#=
include("harps-n/harps-n.jl")
export HARPSN, HARPSN1D, HARPSN2D, AnyHARPSN
import .HARPSN: HARPSN1D, HARPSN2D, AnyHARPSN #, get_inst_module
#import .HARPSN: get_inst_module, filter_line_list, find_worst_telluric_in_each_chunk
=#

end
