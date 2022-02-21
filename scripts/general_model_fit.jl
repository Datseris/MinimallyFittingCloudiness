# First we load the necessary packages, files, and arrays in memory
using DrWatson
@quickactivate "AlbedoBounds"
include(scriptsdir("predictors", "cloudiness_predictors_definition.jl"));

###################################################################### #src
# # Inputs/definitions of predictors (user input is here)
###################################################################### #src

# %% #src
# ## Predictors
# Define the predictors as a tuple of `Symbol`s.
predictors = (:Ω_nf, :ECTEI, :WS10)

# ## Field to be predicted
# Here you can use any of `C, CREsw, CRElw, CF`, or whatever else.
Φ = C

# ## Model definition
# Here we express the model's inner code as a String, and later
# we'll use Julia's metaprogramming to actually make it a runnable function
model_expression = "p[1]*x1 + p[2]*x2*(1-x1) + p[3]*x3"

# ## Fit constraints
# We want to do two limitations:
# * Fit only over ocean
# * Fit within a range of latitudes only
# To this end, we define:
MAXDEG = 80 # fit only within ± MAXDEG
ocean_mask_perc = 50 # points with % ≥ than this are considered "ocean"

######################################################################## #src
# # Model fit code (do not touch anything here!)
######################################################################## #src
close("all")
Ps = map(p -> getindex(field_dictionary,p), predictors)

OCEAN_MASK = O .≥ ocean_mask_perc

# Here we use Julia's metaprogramming capabilities to actually instantiate
# the model function. The code also guarantees that the function definition
# will be overwritten for the given amount of predictors we have.
# I am extremely prooud of the way I've set this up.

function eval_model_equations(model_expression, N)
    parsed_expression = Meta.parse(model_expression)
    _vars_expand = [Symbol("x", "$i") for i in 1:N]
    eval(:(model(p, $(_vars_expand...), args...) = @. $(parsed_expression)))
end
eval_model_equations(model_expression, length(predictors))

# To make the parameters of the resulting fit comperable
# with each other, we will normalize the predictors to approximately [0,1] 
# based on spatial variability
function normfield(x)
    mi, ma = extrema((timemean(x)))
    return (x .- mi) ./ (ma - mi)
end
# Ps = map(x -> ClimArray(normfield(x); name = x.name), Ps)

# Zonal means weighted the ocean fraction
oceany_zonal(X, OCEAN_MASK) = to_non_nan_lats(timemean(zonalmean(X, OCEAN_MASK)))
plot_zonal_averages([oceany_zonal(X, OCEAN_MASK) for X in (Φ, Ps...)])

# In the following we will examine the quality of the fit defined by the
# function `model` (see below). For each subsection, we will
# peform the fit of the parameters of `model` to a given form
# of data: full, zonally averaged, zonally+temporally averaged.

# We obtain the necessary parameters from the fit and always then generate a
# **full** spatiotemporal field that is the model fit, by using the
# spatiotemporal fields of the predictors.

include(srcdir("modelling", "general.jl"))

# model(p, fields) = model(p, fields...)
NP = 8 # allow up to `NP` parameters
p0 = ones(NP) # initial parameters don't matter, 10 random seeds are taken
pl = -maximum(Φ)*10ones(NP) # lower bounds
pu = maximum(Φ)*10ones(NP) # upper bounds

function ocean_masked(X, OCEAN_MASK, MAXDEG)
    latitude_bounding = Coord(Lat(Between(-MAXDEG, MAXDEG)))
    Xbounded = X[latitude_bounding]
    bounded_mask_selection = OCEAN_MASK[latitude_bounding]
    return Xbounded[bounded_mask_selection]
end
ocean_masked(X) = ocean_masked(X, OCEAN_MASK, MAXDEG)

# ## Fit full data
X = ocean_masked(Φ)
Ys = [ocean_masked(P) for P in Ps]
@time M, err, p, = perform_fit(model, X, Ys, p0; lower = pl, upper = pu)
full_fit_error = err
full_fit_params = p
println("   Original NRMSE for FULL fit: ", err)
println("   MAXDEG=$(MAXDEG), OCEAN_FRAC=$(ocean_mask_perc)")
println("   fit parameters: ", p, '\n')
timemean_error, mean_correlation, timeseries_errors = 
    plot_total_model_comparison(Φ, Ps, model, p, OCEAN_MASK; MAXDEG)


    # ## Fit zonal+time mean data
function timezonalmean(Φ, OCEAN_MASK, MAXDEG)
    to_non_nan_lats(
        timemean(
            zonalmean(Φ[Coord(Lat((-MAXDEG)..(MAXDEG)))], OCEAN_MASK[Coord(Lat((-MAXDEG)..(MAXDEG)))])
        )
    )
end
Ftz = timezonalmean(Φ, OCEAN_MASK, MAXDEG)
Pstz = map(P -> timezonalmean(P, OCEAN_MASK, MAXDEG), Ps)
Mtz, errtz, ptz, = perform_fit(model, Ftz, Pstz, p0; lower = pl, upper = pu, n = 100)
timezonal_fit_error = errtz
timezonal_fit_params = ptz
println("   -----------------------------------\n")
println("   Original NRMSE for TIME+ZONAL+MEAN fit: ", errtz)
println("   with parameters: ", ptz, '\n')
println("   -----------------------------------")
println("   -----------------------------------")
# Plot resulting fit
fig = figure(); ax = gca();
ax.plot(sind.(gnv(dims(Ftz, Lat))), Ftz.data; label = Φ.name)
ax.plot(sind.(gnv(dims(Mtz, Lat))), Mtz.data; label = "M, nrmse = $(round(errtz;sigdigits=3))", ls = "-.")
(Ftz, Pstz, model, ptz)
fig.suptitle("fit using only zonal+time mean data\nMAXDEG=$(MAXDEG), OCEAN_FRAC=$(ocean_mask_perc)")
ax.legend()
fig.tight_layout(pad=0.3)


# dd results of the fit into the general model fit dataframe for the given field
using DataFrames

F_name = Φ.name
filename = datadir("modelfits", "general", "$(F_name)_all_fits.jld2")
df = isfile(filename) ? load(filename)["df"] : DataFrame(
    ## Column names of an empty dataframe to initialize
    predictors = [], expression = String[], ocean_mask_perc = Float64[],
    maxdeg = Float64[], full_fit_error = Float64[], full_fit_params = Vector{Float64}[], 
    timemean_error = Float64[], mean_correlation = Float64[], timeseries_errors = Vector{Float64}[],
    timezonal_fit_error = Float64[], timezonal_fit_params = Vector{Float64}[],
)

# New row to save to the dataframe
# (order of elements must match order declared in DataFrame initialization)
newrow = (
    predictors, model_expression, Float64(ocean_mask_perc), Float64(MAXDEG),
    full_fit_error, full_fit_params, timemean_error, mean_correlation, timeseries_errors,
    timezonal_fit_error, timezonal_fit_params,
)



# Add current entry only if it doesn't already exist
tocheck = newrow[1:4]
if !(any(x -> isequal(Tuple(x), tocheck), Tables.namedtupleiterator(df[!, 1:4])))
    push!(df, newrow)
    wsave(filename, @strdict(df))
end

# df2 = df[:, [:predictors, :expression, :timezonal_fit_error]]
df2 = select(df, :predictors, :expression, :timezonal_fit_error, 
:timeseries_errors => ByRow(median) => :seasonal_error)

sort!(df2, :timezonal_fit_error)


# %% Contribution of each factors
f1 = @. p[1] + p[2]*Ps[1] # Ω-
f2 = @. p[3]*Ps[2] # eis/ectei/temperature
f3 = @. p[4]*Ps[3] # windspid

for f in (f1,f2,f3)
    plot_spatial_map_with_mask_lat(f, OCEAN_MASK, MAXDEG;
    # vmin = -0.05, vmax = 0.3
    )
end


# %% Make it a notebook #src
using Literate #src
Literate.notebook(@__FILE__, projectdir("notebooks"); execute=false, credit=false) #src
