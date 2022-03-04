# This file simply combines the downloaded ERA5 .nc data (which are stored in 1 year per file)
# into a single file with all the years. It is more convenient to work with data this way.

# WARNING!!! because of a bug in NCDatasets.jl https://github.com/Alexander-Barth/NCDatasets.jl/issues/95
# The merging of different .nc files that have different offset and scale factors
# does not work. ERA5 always exports each individual file with different offset and scale factors.
# Real_Value = (Stored_Value X scale_factor) + add_offset

# So, the solutions are to either load each individual file and then merge,
# or use Python's `xarray` package to load, which doesn't have this bug.

using DrWatson
@quickactivate "MinimallyFittingCloudiness"
era5dir = desktop("ERA5_original")

# WRONG WAY, with NCDatasets:
using ClimateBase, NCDatasets, Glob
cd(era5dir)

for f in ("ERA5_monthly_3D", "ERA5_monthly_2D")
    era5 = glob(f*"*")
    alldata = NCDataset(era5; aggdim = "time")
    write(f*".nc", alldata)
end

# CORRECT WAY, with NCDatasets:

# CORRECT WAY, with `xarray`:
using PyCall
xr = pyimport("xarray")
ERA5_files = joinpath(era5dir, "ERA5_monthly_3D_*")
xa = xr.open_mfdataset(ERA5_files)
xa.to_netcdf(joinpath(era5dir, "ERA5_3D.nc"))
