""" `make_manifest(data_path::String, target_subdir::String, Inst::Module; [opts] )`
Returns a dataframe containing a list of files to be read and some metadata (e.g., observation times)

# Optional arguements
- verbose = true
"""
function make_manifest(data_path::String, target_subdir::String, Inst::Module; verbose::Bool = true)
   #=
   if Inst == EXPRES
      #@assert isdefined(Main,:expres_data_path)
      data_path = joinpath(expres_data_path,target_subdir)
   elseif Inst == NEID
      #  @assert isdefined(Main,:solar_data_path)
      data_path = joinpath(solar_data_path,target_subdir)
   else
      @error "Specified instrument didn't match a module name, ", string(Inst), "."
   end
   =#
   if !isdir(data_path)
      @error "Can't access data directory ", data_path, "."
   end
   if !isdir(joinpath(data_path,target_subdir))
      @error "Can't access target data subdirectory ", joinpath(data_path,target_subdir), "."
   end

   if verbose  println("# Creating manifest of files to process.")    end
   df_files = Inst.make_manifest(joinpath(data_path,target_subdir))
   if size(df_files,1) < 1
      @error("Did not find any files in ", joinpath(data_path,target_subdir) ,".")
   end
   #=
   df_files = df_files |> @filter( _.target == target_fits ) |> DataFrame
   if size(df_files,1) < 1
      @error("Did not find any files with target field matching ", target_fits,".")
   end
   if size(df_files,1) > max_spectra_to_use
      df_files = df_files |> @take(max_spectra_to_use) |> DataFrame
   end
   =#
   return df_files
end

"""Read manifest containing `filename`, `bjd`, `target`, and optionally additional metadata from CSV file. """
function read_manifest(fn::String)
    df = CSV.read(fn,DataFrame,threaded=false)
    @assert hasproperty(df,:filename)
    @assert hasproperty(df,:bjd)
    @assert hasproperty(df,:target)
    @assert size(df,1) >= 1
    return df
end

"""Write manifest containing `filename`, `bjd`, `target`, and optionally additional metadata from CSV file. """
function write_manifest(fn::String, df)
    @assert hasproperty(df,:filename)
    @assert hasproperty(df,:bjd)
    @assert hasproperty(df,:target)
    @assert size(df,1) >= 1
    CSV.write(fn,df)
end
