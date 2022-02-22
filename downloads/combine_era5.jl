# This file simply combines the downloaded ERA5 .nc data (which are stored in 1 year per file)
# into a single file with all the years. It is more convenient to work with data this way.

using DrWatson
@quickactivate "MinimallyFittingCloudiness"
using ClimateBase, NCDatasets, Glob

era5dir = desktop("ERA5_original")
cd(era5dir)

for f in ("ERA5_monthly_3D", "ERA5_monthly_2D")
    era5 = glob(f*"*")
    alldata = NCDataset(era5; aggdim = "time")
    write(f*".nc", alldata)
end