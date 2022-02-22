# See `how_to_download_ERA5.md`.
# This script is specifically for data with 3 spatial dimensions, as in here:
# https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels-monthly-means?tab=form

using DrWatson
@quickactivate ""

YEARS = string.(2001:2020)
main = "reanalysis-era5-pressure-levels"
times = [lpad(x, 2, '0')*":00" for x in 0:1:23]
producttype = "reanalysis"

using PyCall
cdsapi = pyimport("cdsapi")

c = cdsapi.Client()
datapath = desktop("ERA5_original")

for year in YEARS
savename = joinpath(datapath, "ERA5_hourly_W_$year.nc")
mkpath(dirname(savename))
isfile(savename) && continue

config =
Dict(
    "product_type" => producttype,
    "variable" => ["vertical_velocity"],
    "pressure_level" => ["500"],
    "year" => year,
    "month" => [lpad(x, 2, '0') for x in 1:12],
    "day" => [lpad(x, 2, '0') for x in 1:31],
    "time" => times,
    "expver" => "1",
    "format" => "netcdf",
    "grid" => string.([1.0, 1.0]),
)

try # ensure that remaining downloads will be triggered if one goes wrong.
    c.retrieve(main, config, savename)
catch err
    println("For year $year :")
    showerror(stdout, err)
    println()
end

end
