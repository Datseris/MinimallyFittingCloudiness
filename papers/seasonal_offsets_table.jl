# This file is expected to run after the `best_fit_` files have run.
# It utilizes the offset containers to make a table of offsets.
using DrWatson
@quickactivate "MinimallyFittingCloudiness"
using DataFrames, PrettyTables
# Make dataframe First
df = DataFrame(type = String[], hemisp = String[], loc = String[], mean = Float64[])

for (name, val) in Iterators.flatten((seasonal_offsets_albedo, seasonal_offsets_lcre))
    type, hemisp, loc = split(name, ' ')
    push!(df, (; type, hemisp, loc, mean = val))
end

latex_out = pretty_table(df, backend = Val(:latex))