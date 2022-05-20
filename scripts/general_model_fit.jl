# # Introduction
# This file performs the standardized nonlinear least squares fit on a prescribed target
# field, a given set of predictors, and an expression for how to combine them. It calculates
# various error measures and stores them in a dataframe.

###################################################################### #src
# # Inputs/definitions of predictors (user input is here)
###################################################################### #src
# First we load the necessary packages, files, and arrays in memory
using DrWatson
@quickactivate "MinimallyFittingCloudiness"
include(scriptsdir("fields_definitions.jl"));
# %% #src
# ## Predictors
# Define the predictors as a tuple of `Symbol`s, which can then access
# the dictionary defined in the file `fields_definition.jl`
predictors = (:Ω_std, :Ω_mean, :Tsfc)
predictors = (:Ω_mean, :ECTEI)

# ## Field to be predicted
# Symbol containing the name of the field.
predicted = :L
predicted = :C

# ## Model definition
# Here we express the model's inner code as a String, and later
# we'll use Julia's metaprogramming to actually make it a runnable function.
# Predictors are always expressed as `x1, x2, ...`, and the order corresponds
# to the order of the `predictors` variable.
model_expression = "p[1]*x1 + p[2]*x2 + p[3]*x3"
model_expression = "p[1]*50*(tanh(p[2]*x1 + p[3]*x2) + 1)"

# ## Fit constraints
# We want to do two limitations:
# * Fit only over ocean
# * Fit within a range of latitudes only
# To this end, we define:
MAXDEG = 70  # fit only up to MAXDEG
MINDEG = -70 # fit only down to MINDEG
ocean_mask_perc = 50 # points with % ≥ than this are considered "ocean"

# User input stops here.

######################################################################## #src
# # Model fit code
######################################################################## #src
Φ = field_dictionary[predicted]
Ps = map(p -> getindex(field_dictionary, p), predictors)
OCEAN_MASK = field_dictionary[:O] .≥ ocean_mask_perc
include(srcdir("fitting", "general.jl"))
include(srcdir("fitting", "masking.jl"))

close("all") #src

# Here we use Julia's metaprogramming capabilities to actually instantiate
# the model function. The code also guarantees that the function definition
# will be overwritten for the given amount of predictors we have.
# (the details of this are in `src/fitting/general.jl`)
eval_model_equations(model_expression, predictors)

# This has now created a `model` function, `model(p, x1, x2, ...)` that
# evaluates exactly the string of the `model_expression` variable
model

# Zonal means weighted with the ocean fraction
plot_zonal_averages([oceany_zonalmean(X, OCEAN_MASK) for X in (Φ, Ps...)])

# We obtain the necessary parameters from the fit and always then generate a
# **full** spatiotemporal field that is the model fit, by using the
# spatiotemporal fields of the predictors.
NP = 4 # allow up to `NP` parameters
p0 = ones(NP) # initial parameters don't matter, 10 random seeds are taken
pl = -1000ones(NP) # lower bounds
pu = 1000ones(NP) # upper bounds

# ## Fit full data
println("\n\n---Starting fit")
Y = ocean_masked(Φ, OCEAN_MASK, MAXDEG, MINDEG)
Xs = [ocean_masked(P, OCEAN_MASK, MAXDEG, MINDEG) for P in Ps]
@time M, err, p, = perform_fit(model, Y, Xs, p0; lower = pl, upper = pu, show_info=false, n = 2)
full_fit_error = err
full_fit_params = p
@show predictors
@show model_expression
println("   Original NRMSE for FULL fit: ", err)
# println("   Mean-normalized error: ", μrmse(M, Y))
println("   fit parameters: ", p, '\n')
timemean_error, mean_correlation, timeseries_errors, timezonal_full_error =
    plot_total_model_comparison(Φ, Ps, model, p, OCEAN_MASK; MAXDEG, MINDEG)


# ## Fit zonal+time mean data only
Ftz = maskedtimezonalmean(Φ, OCEAN_MASK, MAXDEG, MINDEG)
Pstz = map(P -> maskedtimezonalmean(P, OCEAN_MASK, MAXDEG, MINDEG), Ps)
w = cosd.(gnv(dims(Ftz, Lat))) # weights
Mtz, errtz, ptz, = perform_fit(model, Ftz, Pstz, p0;
lower = pl, upper = pu, n = 100, w)
timezonal_fit_error = errtz
timezonal_fit_params = ptz
println("   \n")
println("   Original NRMSE for TIME+ZONAL+MEAN fit: ", errtz)
# println("   Mean-normalized error for TIMEZONAL mean fit: ", μrmse(Ftz, Mtz))
println("   with parameters: ", ptz, '\n')
println("   --------------------------------------------------------------\n")

# Plot resulting fit
fig = figure(); ax = gca();
ax.plot(sind.(gnv(dims(Ftz, Lat))), Ftz.data; label = Φ.name)
ax.plot(sind.(gnv(dims(Mtz, Lat))), Mtz.data; label = "M, nrmse = $(round(errtz;sigdigits=3))", ls = "-.")
fig.suptitle("fit using only zonal+time mean data\nMAXDEG=$(MAXDEG), OCEAN_FRAC=$(ocean_mask_perc)")
ax.legend()
fig.tight_layout(pad=0.3)
fig

# ## Add fit resuls to database
using DataFrames

filename = datadir("modelfits", "general", "$(Φ.name)_all_fits.jld2")
df = isfile(filename) ? load(filename)["df"] : DataFrame(
    ## Column names of an empty dataframe to initialize
    predictors = [], expression = String[], ocean_mask_perc = Float64[],
    maxdeg = Float64[], full_fit_error = Float64[], full_fit_params = Vector{Float64}[],
    timemean_error = Float64[], mean_correlation = Float64[],
    timeseries_errors = Vector{Float64}[], timezonal_full_error = Float64[],
    timezonal_fit_error = Float64[], timezonal_fit_params = Vector{Float64}[],
)

# New row to save to the dataframe
# (order of elements must match order declared in DataFrame initialization)
newrow = (
    predictors, model_expression, Float64(ocean_mask_perc), Float64(MAXDEG),
    full_fit_error, full_fit_params, timemean_error, mean_correlation,
    timeseries_errors, timezonal_full_error,
    timezonal_fit_error, timezonal_fit_params,
)

# Add current entry only if it doesn't already exist
tocheck = newrow[1:4]
if !(any(x -> isequal(Tuple(x), tocheck), Tables.namedtupleiterator(df[!, 1:4])))
    push!(df, newrow)
    wsave(filename, @strdict(df))
end

# Subselect dataframe for presentation
df2 = select(df, :predictors, :expression, :timezonal_full_error,
:timeseries_errors => ByRow(median) => :seasonal_error)

sort!(df2, :timezonal_full_error)


# %% Make it a notebook #src
using Literate #src
Literate.notebook(@__FILE__, projectdir("notebooks"); execute=false, credit=false) #src
