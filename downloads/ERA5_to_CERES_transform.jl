# This file converts the downloaded ERA5 data (which are stored in 1 year per file)
# into the format of CERES data, which is a single file for all data,
# latitude coordinate goes from -90 to 90, latitude spacing is between cell edges,
# and latitude and longitude coordinates have attributes
# that define them as CF-compatible standards.
# Notice that to use CDO one final transformation needs to be done which permutes the
# dimensions so that time is the final dim

using DrWatson
@quickactivate ""
using ClimateBase, NCDatasets
include(srcdir("xarray_loading.jl"))
Time = ClimateBase.Time
EBAF_TOA = datadir("CERES", "EBAF_TOA.nc")
R = ncread(EBAF_TOA, "toa_sw_all_mon")
CERES_LATS = dims(R, Lat).val
CERES_META = dims(R, Lat).metadata
dimensionality = 2
ncdetails(desktop("ERA5_original", "ERA5_monthly_$(dimensionality)D_reanalysis_2001.nc"))

# This needs to install xarray and dask from Conda and loads the .nc file as xarray.
using PyCall
xr = pyimport("xarray")
ERA5_files = desktop("ERA5_original", "ERA5_monthly_$(dimensionality)D_reanalysis_*.nc")
xa = xr.open_mfdataset(ERA5_files)
all_vars = collect(xa.data_vars)
W = climarray_from_xarray(xa, all_vars[1], "test")

function to_CERES_latitude(X, lats = CERES_LATS, meta = CERES_META)
    W = reverse(X; dims = Lat) # ERA5 have different ordering from CERES
    newW = W[Lat(1:length(lats))] # drop 1 latitude index
    for i ∈ 1:length(lats)
        newW[Lat(i)] .= 0.5(W[Lat(i)] .+ W[Lat(i+1)])
    end
    newdims = Any[dims(X)...]
    newlat = Lat(CERES_LATS; metadata = meta)
    newdims[dimnum(X, Lat)] = newlat
    # Simple renaming because `reverse` does not conserve name and attributes
    return ClimArray(newW, (newdims...,); name = X.name, attrib = X.attrib)
end

# %% Start creating the dataset
global_attributes = xa.attrs

function add_dims_to_ncfile!(ds::NCDatasets.AbstractDataset, dimensions::Tuple)
    dnames = [ClimateBase.dim_to_commonname(d) for d in dimensions]
    for (i, d) ∈ enumerate(dnames)
        haskey(ds, d) && continue
        println("writing dimension $(d)...")
        v = dimensions[i].val
        # this conversion to DateTime is necessary because of bugs in NCDatasets.jl
        eltype(v) == Date && (v = DateTime.(v))
        l = length(v)
        defDim(ds, d, l) # add dimension entry
        # write dimension values as a variable as well (mandatory)
        # TODO: add check here that attribs has the "units" quantity (mandatory)
        meta = dimensions[i].metadata
        meta isa Dict || (meta = Dict())
        if isempty(meta) && haskey(ClimateBase.DEFAULT_ATTRIBS, d)
            @warn "Dimension $d has no attributes, adding default attributes (mandatory)."
            meta = ClimateBase.DEFAULT_ATTRIBS[d]
        end
        defVar(ds, d, v, (d, ); attrib = meta)
    end
end

NCDataset(desktop("ERA5_original", "ERA5_$(dimensionality)D.nc"), "c"; attrib = global_attributes) do ds
    for (i, fieldname) in enumerate(all_vars)
        println("processing variable $fieldname...")
        W = climarray_from_xarray(xa, fieldname)
        println("converting to CERES format...")
        X = to_CERES_latitude(W)
        # X = W
        if hasdim(X, Dim{:expver})
            # Stupid Copernicus interface cannot drop `expver` dimension on download
            # So we drop it manually, because it is completely useless for us...
            X = X[Dim{:expver}(1)]
        end
        add_dims_to_ncfile!(ds, dims(X))
        println("writing the CF-variable...")
        attrib = X.attrib
        dnames = [ClimateBase.dim_to_commonname(d) for d in dims(X)]
        data = Array(X)
        @show (fieldname, summary(data), dnames)
        defVar(ds, fieldname, data, (dnames...,); attrib)
    end
end

# %% Test that it worked
X = ncread(desktop("ERA5_original", "ERA5_$(dimensionality)D.nc"), all_vars[1])
