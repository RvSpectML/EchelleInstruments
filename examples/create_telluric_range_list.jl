include("read_expres_data_101501.jl")

#df_for_one_order_of_one_spectra = EXPRES.find_ranges_with_tellurics_in_order(all_spectra[5], 81)
#df_for_one_spectra = EXPRES.find_ranges_with_tellurics(all_spectra[5])

# Overwrite barycentric wavelengths with observatory-frame wavelengths
using FITSIO
map(obsid->all_spectra[obsid].λ .= FITSIO.read(FITS(df_files_use.Filename[obsid])["optimal"],"wavelength"), 1:length(all_spectra) )

wavelength_ranges_with_tellurics = EchelleInstruments.find_ranges_with_tellurics(all_spectra, min_Δv_clean=EXPRES.default_min_Δv_clean)

#using CSV
#CSV.write(joinpath(output_dir,"telluric_ranges.csv"),wavelength_ranges_with_tellurics)

min_λ_min = minimum(map(i->minimum(all_spectra[i].λ[.!isnan.(all_spectra[i].metadata[:tellurics]) .& (all_spectra[i].metadata[:tellurics] .== 1) ]), 1:length(all_spectra) ))
max_λ_max = maximum(map(i->maximum(all_spectra[i].λ[.!isnan.(all_spectra[i].metadata[:tellurics]) .& (all_spectra[i].metadata[:tellurics] .== 1) ]), 1:length(all_spectra) ))

wavelength_ranges_without_tellurics = EchelleInstruments.make_ranges_without_tellurics(wavelength_ranges_with_tellurics,
      λ_start = min_λ_min, λ_stop = max_λ_max, min_Δv= 1000)

wavelength_ranges_without_tellurics2 = EchelleInstruments.make_ranges_without_tellurics(wavelength_ranges_with_tellurics,
            λ_start = min_λ_min, λ_stop = max_λ_max, min_Δv= 1000, max_Δv = 1000000)

#CSV.write(joinpath(output_dir,"telluric_free_ranges_min_Δv=1000.csv"),wavelength_ranges_without_tellurics)

# Useful to compare how much extra information we'd get by reducing Δv
mapreduce(r->(r.lambda_hi-r.lambda_lo)*2/(r.lambda_hi+r.lambda_lo),+,eachrow(wavelength_ranges_without_tellurics))
