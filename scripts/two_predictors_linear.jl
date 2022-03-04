# This file performs a zero-intercept linear fit of 2 predictors by going through
# all possible 2-predictor combination from a given pool of predictors.
# The model saves various error measures into a dataframe.

using DrWatson
@quickactivate "MinimallyFittingCloudiness"
include(scriptsdir("fields_definitions.jl"));
include(srcdir("fitting", "masking.jl"))
import LsqFit, GLM

predictors = (:Ω_mean, :Ω_std, :Ω_nf, :WS10, :Tsfc, :EIS, :ECTEI, :V, :q700)
predicted = :C
MAXDEG = 80 # fit only within ± MAXDEG
ocean_mask_perc = 50 # points with % ≥ than this are considered "ocean"
# %% #src
allow_intercept = false

fitter_used = GLM # doesn't make a difference in the end.
n = length(predictors)
OCEAN_MASK = field_dictionary[:O] .≥ ocean_mask_perc
Y = field_dictionary[predicted]
Ymasked = ocean_masked(Y, OCEAN_MASK, MAXDEG)
Ytz = maskedtimezonalmean(Y, OCEAN_MASK, MAXDEG)
if allow_intercept
    Xmasked = hcat(Ymasked, Ymasked, ones(eltype(Ymasked), length(Ymasked)))
    linearmodel = (x, p) -> @. p[1]*x[:, 1] + p[2]*x[:, 2] + p[3]
else
    Xmasked = hcat(Ymasked, Ymasked)
    linearmodel = (x, p) -> @. p[1]*x[:, 1] + p[2]*x[:, 2]
end
full_error, timezonal_error, seasonal_error, mean_pearson = [ones(n, n) for _ in 1:4]


for i in 1:n
    X1 = field_dictionary[predictors[i]]
    X1masked = ocean_masked(X1, OCEAN_MASK, MAXDEG)
    X1tz = maskedtimezonalmean(X1, OCEAN_MASK, MAXDEG)
    Xmasked[:, 1] .= X1masked

    for j in (i+1):n

        X2 = field_dictionary[predictors[j]]
        X2masked = ocean_masked(X2, OCEAN_MASK, MAXDEG)
        X2tz = maskedtimezonalmean(X2, OCEAN_MASK, MAXDEG)
        Xmasked[:, 2] .= X2masked

        # Fit using LsqFit
        if fitter_used == LsqFit
            modelfit = LsqFit.curve_fit(linearmodel, Xmasked, Ymasked, [0.0, 0.0, 0.0])
            p1, p2 = modelfit.param
            p3 = allow_intercept ? modelfit.param[3] : 0
        elseif fitter_used == GLM
            modelfit = GLM.lm(Xmasked, Ymasked)
            p1, p2 = GLM.coef(modelfit)
            p3 = allow_intercept ? GLM.coef(modelfit)[3] : 0
        end
        @show (i, j, p1, p2)

        # Fit using GLM
        # Extract error measures
        M = @. p1*X1 + p2*X2 + p3
        Mmasked = @. p1*X1masked + p2*X2masked + p3
        Mtz = @. p1*X1tz + p2*X2tz + p3
        full_error[i,j] = nrmse(Mmasked, Ymasked)
        timezonal_error[i, j] = nrmse(Mtz, Ytz)
        seasonal_error[i, j] = seasonal_nrmse(Y, M, OCEAN_MASK)
        mean_pearson[i, j] = 1 - correlationmap(Y, M, OCEAN_MASK)[2]
    end
end

# Save output
wsave(datadir("modelfits","linear_2_predictors", "$(predicted)_intercept=$(allow_intercept).jld2"),
    @strdict(
        predictors, predicted, full_error,
        timezonal_error, seasonal_error, mean_pearson,
        MAXDEG, ocean_mask_perc,
    )
)