""" Functions that were used with NEID data products prior to shipping to KPNO """
module PreShip
using DataFrames, FITSIO

""" Read NEID (non-solar) data from FITS file, and return in a Spectra2DBasic object."""
function read_data   end

function read_data(fn::String, metadata::Dict{Symbol,Any} )
    f = FITS(fn)
    @assert read_key(f[1],"SKY-OBJ")[1] != "Solar"
    λ, flux, var  = FITSIO.read(f["SCIWAVE"]), FITSIO.read(f["Sci Flux"]), FITSIO.read(f["Sci Variance"])
    metadata[:normalization] = :raw
    Spectra2DBasic(λ, flux, var, NEID2D(), metadata=metadata)
end

function read_data(fn::String)
    f = FITS(fn)
    @assert read_key(f[1],"SKY-OBJ")[1] != "Solar"
    hdr = FITSIO.read_header(f[1])
    metadata = Dict(zip(map(k->Symbol(k),hdr.keys),hdr.values))
    λ, flux, var  = FITSIO.read(f["SCIWAVE"]), read(f["Sci Flux"]), read(f["Sci Variance"])
    metadata[:normalization] = :raw
    Spectra2DBasic(λ, flux, var, NEID2D(), metadata=metadata)
end

function read_data(dfr::DataFrameRow{DataFrame,DataFrames.Index})
    fn = dfr.Filename
    metadata = Dict(zip(keys(dfr),values(dfr)))
    read_data(fn,metadata)
end

""" Read NEID Solar data from FITS file, and return in a Spectra2DBasic object."""
function read_solar_data
end

function read_solar_data(fn::String, metadata::Dict{Symbol,Any} )
    f = FITS(fn)
    @assert read_key(f[1],"SKY-OBJ")[1] == "Solar"
    λ, flux, var  = read(f["SKYWAVE"]), read(f["Sky Flux"]), read(f["Sky Variance"])
    Spectra2DBasic(λ, flux, var, NEID2D(), metadata=metadata)
end

function read_solar_data(fn::String)
    f = FITS(fn)
    @assert read_key(f[1],"SKY-OBJ")[1] == "Solar"
    hdr = FITSIO.read_header(f[1])
    metadata = Dict(zip(map(k->Symbol(k),hdr.keys),hdr.values))
    λ, flux, var  = FITSIO.read(f["SKYWAVE"]), FITSIO.read(f["Sky Flux"]), FITSIO.ead(f["Sky Variance"])
    Spectra2DBasic(λ, flux, var, NEID2D(), metadata=metadata)
end

function read_solar_data(dfr::DataFrameRow{DataFrame,DataFrames.Index})
    fn = dfr.Filename
    metadata = Dict(zip(keys(dfr),values(dfr)))
    read_solar_data(fn,metadata)
end

""" Read CSV of NEID drift corrections, interpolate to bjd's in df and insert into df[:,drift]. """
function read_drift_corrections!(fn::String, df::DataFrame, df_time_col::Symbol = :bjd)
    drift_corrections = CSV.read(fn, DataFrame, header=["bjd", "sci_drift", "cal_drift"]);
    @assert any(isequal(:bjd),propertynames(drift_corrections))
    @assert any(isequal(:cal_drift),propertynames(drift_corrections))
    drift_interp = LinearInterpolation(drift_corrections[!,:bjd],drift_corrections[!,:cal_drift])
    df[!,:drift] = -drift_interp.(df[!,df_time_col])
    return df
end

""" Read CSV of NEID barycentric corrections, interpolate to bjd's in df and insert into df[:,ssb_rv]. """
function read_barycentric_corrections!(fn::String, df::DataFrame, df_time_col::Symbol = :bjd)
    ssb_corrections = CSV.read(fn, DataFrame, header=["bjd","rv_ssb"], select=[1,2], types=types=[Float64,Float64], datarow=2, silencewarnings=true);
    @assert any(isequal(:bjd),propertynames(ssb_corrections))
    @assert any(isequal(:rv_ssb),propertynames(ssb_corrections))
    ssb_interp = LinearInterpolation(ssb_corrections.bjd, ssb_corrections.rv_ssb)
    df[!,:ssb_rv] = ssb_interp.(df[!,df_time_col])
    return df
end

end  # module PreShip
