# This file performs the standardized nonlinear least squares fit on a prescribed target
# field, a given set of predictors, and an expression for how to combine them.
# In complement to the main fit script `general_model_fit.jl`, this file
# estimates the variability of parameters when the fit is applied in sub-samples
# of space or time. See `general_model_fit.jl` for comments explaining in
# detail how such fitting scripts work.
using DrWatson
@quickactivate "MinimallyFittingCloudiness"
using ClimateBase
include(srcdir("fitting", "general.jl"))
include(srcdir("masking.jl"))
include(scriptsdir("fields_definitions.jl"));

# Model setup
predicted = :L

if predicted == :C
    predictors = (:Ω_mean, :ECTEI)
    model_expression = "p[1]*50*(tanh(p[2]*x1 + p[3]*x2) + 1)"
elseif predicted == :L
    predictors = (:Ω_std, :Ω_mean, :Tsfc)
    model_expression = "p[1]*x1 + p[2]*x2 + p[3]*x3"
end

# Fit constrains
MAXDEG = 70  # fit only up to MAXDEG
MINDEG = -70 # fit only down to MINDEG
ocean_mask_perc = 50 # points with % ≥ than this are considered "ocean"

# Subsampling in time or space
# (percentage of data sub-sampled)
space_sampling_percentage = 0.25
time_sampling_percentage = 1.0
# If false, the sampling is done sequantially in time. Otherwise randomly
# (sampling in space is always done randomly)
randomized_sampling = false
random_samples = 2000 # how many random samples to do

# %% Run code
Φ = field_dictionary[predicted]
Ps = map(p -> getindex(field_dictionary, p), predictors)
# For masking we do space mean as sub sampling happens orthogonally
# in space and time
OCEAN_MASK = timemean(field_dictionary[:O]) .≥ ocean_mask_perc
bounded_mask_selection = findall(!iszero, gnv(OCEAN_MASK[Coord(Lat(Between(MINDEG, MAXDEG)))]))
total_masked_coords = length(bounded_mask_selection)


eval_model_equations(model_expression, predictors)
eval_model_equations(model_expression, predictors)
NP = 3 # allow up to `NP` parameters
p0 = ones(NP) # initial parameters don't matter to my knowledge

# initial space mask (for accelerating subsampling)
function subsample_indices(Y)
    # First sample in time:
    n = round(Int, size(Y, Tim)*time_sampling_percentage)
    if randomized_sampling
        time_idxs = sample(1:t, n; replace = false)
    else
        time_idxs = 1:n
    end
    # then in space (always here, as this does ocean masking as well)
    # (notice that this always shuffles space, but we don't care)
    n = round(Int, space_sampling_percentage*total_masked_coords)
    coord_idxs = sample(bounded_mask_selection, n; replace = false)
    return time_idxs, coord_idxs
end

function subsample(fields)
    time_idxs, coord_idxs = subsample_indices(fields[1])
    return map(f -> f[Time(time_idxs), Coord(coord_idxs)], fields)
end

# Loop over randomized sampling
p_distributions = [zeros(random_samples) for _ in 1:NP]
errors = zeros(random_samples)

using ProgressMeter
@showprogress 0.2 "Random sampling..." for i in 1:random_samples
    Y, Xs... = subsample((Φ, Ps...))
    M, err, p, = perform_fit(model, Y, Xs, p0; show_info=false, n = 2)
    errors[i] = err
    for (j, v) in enumerate(p)
        p_distributions[j][i] = v
    end
end

input_data = @strdict(
    predicted, space_sampling_percentage,
    random_samples, time_sampling_percentage
)
output_data = merge(input_data, @strdict(p_distributions, errors))
wsave(datadir("parameter_variability", savename(input_data, "jld2")), output_data)

# %% produce the plot
# We need something to compare
fig, axs = subplots(1,3)
for (j, pdist) in enumerate(p_distributions)
    ax = axs[j]
    ax.hist(pdist, 50)
    ax.set_xlabel("\$p_$(j)\$")
    ax.set_yticks([])
end
fig.suptitle("\$$(predicted)\$, Space subs. = $(space_sampling_percentage), Time subs. = $(time_sampling_percentage), samples = $(random_samples)")
fig.tight_layout()


wsave(plotsdir("parameter_variability", savename(input_data, "png")), fig)
