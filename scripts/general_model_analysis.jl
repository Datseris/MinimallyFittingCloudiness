#=
This script loads the central dataframe(s) generated with `general_model_fit.jl`
and then performs some standardized analysis like identifying minima, etc.
=#
using DrWatson
@quickactivate "AlbedoBounds"
include(scriptsdir("fields_definition.jl"));
using DataFrames


# %%
F_name = :CRElw
filename = datadir("modelfits", "general", "$(F_name)_all_fits.jld2")
df = wload(filename)["df"]
F = field_dictionary[F_name]
column_to_use = :timezonal_fit_error
df2 = sort(df, column_to_use)

# zonal plot of best 3
for i in 1:3
    row = df2[i, :]
    Ps = map(p -> getindex(field_dictionary,p), row.predictors)
    eval_model_equations(row.expression, length(row.predictors))
    OCEAN_MASK = O .> row.ocean_mask_perc
    Ftz = timezonalmean(F, OCEAN_MASK, row.maxdeg)
    Pstz = map(P -> timezonalmean(normfield(P), OCEAN_MASK, row.maxdeg), Ps)
    Mtz = model(row.timezonal_fit_params, Pstz...)
    figure()
    lats = sind.(gnv(dims(Ftz, Lat)))
    plot(lats, gnv(Ftz); label = F_name)
    plot(lats, gnv(Mtz); label = "M, nrmse = $(round(row.timezonal_fit_error;sigdigits=3))", ls = "-.")
    gca().set_title(
        string(row.predictors)*"\n"*string(row.expression)
    )
    tight_layout(); legend();
end