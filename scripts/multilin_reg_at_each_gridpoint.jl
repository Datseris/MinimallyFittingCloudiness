#=
This file performs a multilinear regression of a given field `F`
over given predictor fields `Ps`. The regression is done individually at each
spatial point. It then plots the coefficients of each predictor at each 
spatial point. The coefficients are normalized versus the average timeseries std
of the predictor.
=#

###################################################################### #src
# # Inputs/definitions of predictors
###################################################################### #src
using DrWatson
@quickactivate "AlbedoBounds"
include(scriptsdir("predictors", "cloudiness_predictors_definition.jl"));

# %% #src
# ## Predictors
Ps = (Ω_mean, EIS)

# ## Field to be predicted
F = C
# F = CREsw
# F = CRElw

MAXDEG = 70 # fit only within ± MAXDEG
ocean_mask_perc = 50 # points with % ≥ than this are considered "ocean"

###################################################################### #src
# # Linear fit code
###################################################################### #src
OCEAN_MASK_MEAN = timemean(O) .> ocean_mask_perc

import GLM

# function linear_regression_versus_time(F::ClimArray, Ps::Tuple)
# Prepare output array
coord = dims(F, Coord)
pnames = [String(P.name) for P in Ps]
coefs = Dim{:coef}(pnames)

attrib = Dict("Description" => 
"Contains coeffients of linear fit for a given predictor and for a given space point.
Coefficeints are normalized by temporal std of predictor, averaged over all grid cells.")

output = ClimArray(fill(NaN, length(coord), length(coefs)), (coord, coefs); attrib)

# Calculate std of predictors for coefficient normalization
stds = [spacemean(timeagg(std, P)) for P in Ps]
# stds[1] = 1

xdata = zeros(size(F, Time), length(Ps))
ydata = zeros(size(F, Time))

# loop over coordinates for coefficients
for (cj, c) in enumerate(coord)
    c[2] > MAXDEG && continue
    OCEAN_MASK_MEAN[Coord(cj)] || continue
    for (i, P) in enumerate(Ps)
        xdata[:, i] .= P[Coord(cj)] ./ stds[i]
    end
    ydata .= F[Coord(cj)]
    out = GLM.lm(xdata, ydata)
    output[(Coord(cj))] .= GLM.coef(out)
    # if i want confidence intervals:
    # ci05, ci95 = GLM.confint(out)[i, :] for the i-th predictor
end

# Okay let's plot first coefficient
for (i, pn) in enumerate(pnames)
    fig, ax = earthscatter(output[coef = At(pn)]; cmap = :PRGn)
    ax.set_title(pn*", linear coefficient (normalized)")
end
