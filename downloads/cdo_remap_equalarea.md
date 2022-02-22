To remap data using CDO on mistral, first I check if the file structure is compatible with CDO:
```
cdo showformat IN.nc
```
If not (typically because time needs to be first dimension), then I transform them:


```
module load nco
ncpdq -a time,lon,lat,level IN.nc OUT.nc
```

Now I can do the remapping (hopefully) to an equal area grid:
```
cdo remapbil,gea250 -setgridtype,lonlat IN.nc OUT.nc
```
