verbose = true
 if verbose && !isdefined(Main,:RvSpectMLBase)   println("# Loading RvSpecMLBase")    end
 using RvSpectMLBase
 if verbose && !isdefined(Main,:EchelleInstruments)   println("# Loading EchelleInstruments")    end
 using EchelleInstruments, EchelleInstruments.EXPRES
 using DataFrames, Query


# USER: You must create a data_paths.jl file in one of the default_paths_to_search listed below. It need only contain one line:
# expres_data_path = "/path/to/EXPRES/data/not/including/target_subdir"
target_subdir = "10700"   # USER: Replace with directory of your choice
 fits_target_str = "10700"
 output_dir = "examples/output/"
 paths_to_search_for_param = [pwd(),"examples",joinpath(pkgdir(RvSpectMLBase),"..","RvSpectML","examples"), "/gpfs/group/ebf11/default/ebf11/expres"]
 # NOTE: make_manifest does not update its paths_to_search when default_paths_to_search is defined here, so if you change the line above, you must also include "paths_to_search=default_paths_to_search" in the make_manifest() function call below
 pipeline_plan = PipelinePlan()
 dont_make_plot!(pipeline_plan, :movie)

reset_all_needs!(pipeline_plan)
if need_to(pipeline_plan,:read_spectra)
   if verbose println("# Finding what data files are avaliable.")  end
   eval(read_data_paths(paths_to_search=paths_to_search_for_param))
   @assert isdefined(Main,:expres_data_path)
   df_files = make_manifest(expres_data_path, target_subdir, EXPRES )

   if verbose println("# Reading in customized parameters from param.jl.")  end
   eval(code_to_include_param_jl(paths_to_search=paths_to_search_for_param))

   if verbose println("# Reading in ", size(df_files_use,1), " FITS files.")  end
   @time all_spectra = map(row->EXPRES.read_data(row,store_min_data=true, store_tellurics=true, normalization=:continuum),eachrow(df_files_use))
   #@time all_spectra = map(row->EXPRES.read_data(row,store_min_data=true, store_tellurics=true, normalization=:raw),eachrow(df_files_use))
   #@time all_spectra = map(row->EXPRES.read_data(row,store_min_data=true, store_tellurics=true, normalization=:blaze),eachrow(df_files_use))
   #RvSpectML.discard_blaze(all_spectra)
   #RvSpectML.discard_continuum(all_spectra)  # Seeing if affects CCFs
   GC.gc()
   dont_need_to!(pipeline_plan,:read_spectra)
   all_spectra
end
