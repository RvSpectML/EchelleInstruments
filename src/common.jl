default_Δv_to_avoid_tellurics = RvSpectMLBase.max_bc # 29785.8905  # ~ 2pi AU/year in m/s

max_num_spectra = 10000  # arbitrary limit

#import RvSpectMLBase: make_clean_line_list_from_tellurics
#export filter_line_list
#export make_clean_line_list_from_tellurics
""" For generic functions on spectra that are to be overloaded for each instrument, see [RvSpectMLBase](https://rvspectml.github.io/RvSpectMLBase.jl/stable/)'s [src/instrumenets/common.jl](https://github.com/RvSpectML/RvSpectMLBase.jl/blob/main/src/instruments/common.jl)."""

#=
# This is now in RvSpectmLBase
# Declare functions that should be specialized for each instrument here, so they can be imported into their instrument's module.
"Returns minimum allowable order for the specified instrument"
function min_order end

"Returns maximum allowable order for the specified instrument"
function max_order end

"Returns minimum allowable pixel index for the specified instrument and order index"
function min_pixel_in_order end

"Returns maximum allowable pixel index for the specified instrument and order index"
function max_pixel_in_order end

"Returns minimum allowable pixel index for the specified instrument and any order"
function min_pixel end

"Returns maximum allowable pixel index for the specified instrument and any order"
function max_pixel end

"""Returns range of order indices to be used by default for the specified instrument"""
function orders_to_use_default end

"Returns minimum pixel index to be used by default for the specified instrument"
function min_col_default end

"Returns maximum pixel index to be read from header into metadata by default for the specified instrument"
function max_col_default end

"Returns array of symbols to be stored from metadata in header by default for the specified instrument"
function metadata_symbols_default end

"Returns array of strings to be read from header into metadata by default for the specified instrument"
function metadata_strings_default end

"Returns the maximum range of wavelengths (in solar system barycenter frame) to be used for given spectrum"
function get_λ_range end

"""" `filter_line_list`( df, inst )
Returns a dataframe filtered to only include lines appropriate for the specified instrument.
Inputs:
- df: Dataframe containing a column with key `:lambda` for wavelengths to be considered.
- inst: instrument whose properties will be used in choosing which lines to keep.
Optional Inputs:
- λmin: Overide the instrument's default minimum wavelength
- λmax: Overide the instrument's default maximum wavelength
"""
function filter_line_list end

"Returns the module containing code specific to the specified instrument"
function get_inst_module end

"Returns the default velocity width (in m/s) of a tophat CCF mask for the specified instrument"
function default_ccf_mask_v_width end

# Trait-based functions that provide defaults (can be overwritten by instrument-specific versions)
" Returns range of all allowable orders for the specified instrument"
function orders_all end

orders_all(inst::AbstractInstrument2D) = min_order(inst):max_order(inst)

" Returns range of all allowable pixels for the specified instrument"
function pixels_all end

pixels_all(inst::AbstractInstrument2D) = min_pixels_in_order(inst):max_pixel_in_order(inst)
pixels_all(inst::AbstractInstrument1D) = min_pixel(inst):max_pixel(inst)

" Returns maximum number of pixels in a spectrum for the specified instrument"
function max_pixels_in_spectra end

max_pixels_in_spectra(inst::AbstractInstrument1D) = length(pixels_all(inst))
max_pixels_in_spectra(inst::AbstractInstrument2D) = (max_order(inst)-min_order(inst)+1) * (max_pixel_in_order(inst)-min_pixel_in_order(inst)+1)

" Returns maximum number of pixels in a chunk from the specified instrument"
min_pixels_in_chunk(inst::AbstractInstrument) = 6

=#
