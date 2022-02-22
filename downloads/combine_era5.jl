# This file simply combines the downloaded ERA5 .nc data (which are stored in 1 year per file)
# into a single file with all the years. It is more convenient to work with data this way.

using DrWatson
@quickactivate ""
using ClimateBase, NCDatasets, Glob

era5dir = desktop("ERA5_original")
cd(era5dir)

era5 = glob("ERA5_monthly_2D_reanalysis*")
alldata = NCDataset(era5; aggdim = "time")
write("ERA5_2D.nc", alldata)