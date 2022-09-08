# This script makes the plot for the parameter sensitivity
# for the Supplement, using the random subsampling approach.
# Notice that it requires the `scripts/distribution_of_paramfit_values.jl`
# file to have run and produced all the data. This script only loads data, doesn't produce.
using DrWatson
@quickactivate "MinimallyFittingCloudiness"
using ClimateBase
include(srcdir("plotting", "recipes.jl"))
include(scriptsdir("fields_definitions.jl"))
include(srcdir("fitting", "general.jl"))
include(srcdir("fitting", "fit_errors.jl"))

space_sampling_percentage = 0.25
time_sampling_percentage = 1.0
random_samples = 2000
MAXDEG = 70  # fit only up to MAXDEG
MINDEG = -70 # fit only down to MINDEG
ocean_mask_perc = 50 # points with % ≥ than this are considered "ocean"
WHICH_ERROR_MEASURE = 1 # full error, space error, etc...
# How we will obtain error for bottom plots
# see source code of `main_error_measures` for possibilities
error_obtainer = (x, y, m) -> spatial_nrmse(x, y, m)
error_obtainer = (x, y, m) -> nrmse(x[gnv(m)], y[gnv(m)])

# best fits for each field to compare random sampling errors against
best_fit_expression = (
    ((p, x1, x2) -> @. p[1]*50*(tanh(p[2]*x1 + p[3]*x2) + 1)),
    ((p, x1, x2, x3) -> @. p[1]*x1 + p[2]*x2 + p[3]*x3),
)
best_fit_params = (
    [0.40704233122044947, 6.870173231618394, 0.0803937826587933],
    [42.67918994446545, 208.88564817976484, 0.06557681258808708],
)
best_fit_predictors = (
    (:Ω_mean, :ECTEI),
    (:Ω_std, :Ω_mean, :Tsfc)
)

# %% Start plotting and generating distributions
PyPlot.ioff()
fig = figure(constrained_layout=true; figsize = (figx, 2figy))
subfigs = fig.subfigures(2, 1; wspace=0.07, height_ratios=[1.5, 1.0])
axs_params = subfigs[1].subplots(2, 3)
axs_nrmses = subfigs[2].subplots(1, 2)


OCEAN_MASK = (field_dictionary[:O] .> ocean_mask_perc)[Coord(Lat((-MAXDEG)..(MAXDEG)))]

for (i, predicted) in enumerate((:C, :L)) # iterate over predictants
    input_data = @strdict(
        predicted, space_sampling_percentage,
        random_samples, time_sampling_percentage
    )
    data = load(datadir("parameter_variability", savename(input_data, "jld2")))
    @unpack p_distributions = data

    for (j, pdist) in enumerate(p_distributions)
        ax = axs_params[i,j]
        ax.hist(pdist, 20)
        dispersion = std(pdist)/mean(pdist)
        ax.text(0.05, 0.9, "d = $(round(dispersion; sigdigits = 3))",
            transform = ax.transAxes, zorder = 99, fontsize = 18,
        )
        ax.xaxis.set_major_locator(plt.MaxNLocator(5))
        ax.tick_params(axis="x", labelsize=18)
        ax.set_yticks([])
        if i == 2
            ax.set_xlabel("\$p_$(j)\$")
        end
    end

    # NRMSEs of random parameter samplings
    e = zeros(random_samples) # errors of predicted field
    Ps = map(
        X -> getindex(field_dictionary, X)[Coord(Lat((-MAXDEG)..(MAXDEG)))],
        best_fit_predictors[i]
    )
    expression = best_fit_expression[i]
    M0 = get_model_instance(expression, Ps, best_fit_params[i])
    # Now compute errors versus this "best model version" by randomly sampling parameters
    for k in 1:random_samples
        randomp = [p[k] for p in p_distributions]
        Mp = get_model_instance(expression, Ps, randomp)
        # Can change this error to whatever we want
        err = error_obtainer(M0, Mp, OCEAN_MASK)
        e[k] = err
    end
    # Plot distribution
    ax = axs_nrmses[i]
    ax.hist(e, 50)
    ax.set_xlabel("NRMSE(\$M_b\$, \$M_r\$), \$$(predicted)\$")
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax.tick_params(axis="x", labelsize=18)
    ax.set_yticks([])
end


# subfigs[1].subplots_adjust(wspace = 0.1, left = 0.05, right = 0.95)
wsave(papersdir("plots", "fit_stability"), fig)

PyPlot.ion()