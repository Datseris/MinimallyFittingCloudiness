# See `how_to_download_ERA5.md`.
# This script is specifically for data with 2 spatial dimensions (lon lat), as in here:
# https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-monthly-means?tab=form

using DrWatson
@quickactivate ""
HOURLY = false # obtain hourly data or monthly
YEARS = 2001:2020

if HOURLY
    main = "reanalysis-era5-single-levels"
    times = [lpad(x, 2, '0')*":00" for x in 0:1:23]
    producttype = "reanalysis"
else
    main = "reanalysis-era5-single-levels-monthly-means"
    times = "00:00"
    producttype = "monthly_averaged_reanalysis"
end

using PyCall
cdsapi = pyimport("cdsapi")

c = cdsapi.Client()
datapath = desktop("ERA5_original")

for year in YEARS
sampling = HOURLY ? "hourly" : "monthly"
savename = joinpath(datapath, "ERA5_$(sampling)_2D_reanalysis_$year.nc")
mkpath(dirname(savename))

config =
Dict(
    "product_type" => producttype,
    "variable" => [
        # Thermodynamics:
        "2m_temperature", "surface_pressure",
        # Clouds:
        "low_cloud_cover", "total_cloud_cover",
        # "total_column_cloud_ice_water", "total_column_cloud_liquid_water",
        # Radiation, surface:
        "surface_net_solar_radiation", "surface_net_solar_radiation_clear_sky", 
        "surface_net_thermal_radiation", "surface_net_thermal_radiation_clear_sky", 
        "surface_solar_radiation_downward_clear_sky", "surface_solar_radiation_downwards",
        "surface_thermal_radiation_downward_clear_sky", "surface_thermal_radiation_downwards",
        # Radiation, TOA:
        "toa_incident_solar_radiation","top_net_solar_radiation", "top_net_solar_radiation_clear_sky",
        "top_net_thermal_radiation", "top_net_thermal_radiation_clear_sky",
        # Water:
        "total_column_water_vapour", "total_precipitation", "evaporation",
        # Vegetation:
        # "high_vegetation_cover", "low_vegetation_cover",
        # Wind:
        "10m_wind_speed",
    ],
    "year" => year,
    "month" => [lpad(x, 2, '0') for x in 1:12],
    "day" => [lpad(x, 2, '0') for x in 1:31],
    "time" => times,
    "format" => "netcdf",
    "expver" => 1,
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


# Merge the different files.