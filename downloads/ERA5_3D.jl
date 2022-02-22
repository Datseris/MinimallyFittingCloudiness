# See `how_to_download_ERA5.md`.
# This script is specifically for data with 3 spatial dimensions, as in here:
# https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels-monthly-means?tab=form

using DrWatson
@quickactivate ""
YEARS = string.(2001:2020)


main = "reanalysis-era5-pressure-levels-monthly-means"
times = "00:00"
producttype = "monthly_averaged_reanalysis"

using PyCall
cdsapi = pyimport("cdsapi")

c = cdsapi.Client()
datapath = desktop("ERA5_original")

for year in YEARS
savename = joinpath(datapath, "ERA5_3D_$year.nc")
mkpath(dirname(savename))
isfile(savename) && continue

config =
Dict(
    "product_type" => producttype,
    "variable" => [
        # "fraction_of_cloud_cover",
        "relative_humidity",
        "specific_humidity",
        "temperature",
        "vertical_velocity",
    ],
    # "pressure_level" => [ # height dependence version
    #     "1", "5", "20",
    #     "70", "150", "225",
    #     "350", "400", "500",
    #     "650", "700", "775",
    #     "850", "925", "1000",c
    # ],
    "pressure_level" => ["500", "700", "1000"],
    "year" => year,
    "month" => [lpad(x, 2, '0') for x in 1:12],
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
