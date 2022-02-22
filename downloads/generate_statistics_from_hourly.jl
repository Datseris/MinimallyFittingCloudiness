# In this file downloaded HOURLY data from ERA5 are converted to higher
# moments or other statistics. These are all grouped into a single file.

using DrWatson
@quickactivate "MinimallyFittingCloudiness"

using NCDatasets, ClimateBase, StatsBase
include(srcdir("statistics.jl"))
Time = ClimateBase.Ti

path = raw"C:\Users\m300808\Desktop\ERA5_original"
filereg = "ERA5_hourly_W"
files = readdir(path; join=true)
files = filter(file -> occursin(filereg, file), files)
name = "w"
outname = "W"

# Create new time dimension
ds = NCDataset(files[1])
thourly = ds["time"] |> collect |> ClimateBase.vector2range
tnew = range(Date(yearmonth(thourly[1])..., 15); step = Month(1), length = 12*length(files))
# Create the remaining dimensions
w = ncread(files[1], name, (:, :, 1:2))
odims = otherdims(w, (Time,))
newdims = (odims..., Time(tnew))

# Create the arrays to be populated with various monents/statistics
nf(x) = count(<(0), x)/length(x) # negative fraction
reducing_fs = (mean, std, skewness, nf)
fields = map(
    f -> ClimArray(zeros(size(newdims)), newdims; name = "$(outname)_"*string(f)), 
    reducing_fs
)

for (j, file) in enumerate(files)
    @show j/length(files)
    # Because the hourly data of one full year do not fit in memory,
    # we need to load only parts of them...
    s = ncsize(file, name)
    v = 2 # selected every `v`th hour
    Ω = ncread(file, name, (:, :, 1:v:s[3]))

    # Do monthly aggregation using ClimateBase
    for (f, F) in zip(reducing_fs, fields)
        if f == mean || f == std
            F_mon = monthlyagg(Ω, f; mday=15)
        else
            F_mon = monthlyagg_custom(Ω, f; mday=15)
        end
        k = 1 + 12*(j-1)
        F[Time(k:k+11)] .= F_mon
    end
end

@assert !any(iszero, fields[2])

# Alrighty, now add them to the file
outfile = joinpath(path, "ERA5_$(outname)_statistics.nc")
a = Dict(
    "Conventions" => "CF-1.6", 
    "Note" => "Monthly statistics of $(name) from hourly data"
)
@tag!(a; storepatch = false)
ncwrite(outfile, fields; globalattr = a)
