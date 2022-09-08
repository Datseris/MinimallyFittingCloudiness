# This file defines various error measures, as well as convenience functions to compute
# error measures between fields. Typically one of the two fields is a "fit" of the
# other, hence the file is named `fit_errors.jl`.
using StatsBase, Statistics

###########################################################################################
# # Standardized errors
###########################################################################################
function mse(x, y)
    m = length(x)
    @assert m == length(y)
    if any(isnan, x)
        j = findall(!isnan, x)
        @inbounds mse = sum(abs2(x[i] - y[i]) for i in j) / length(j)
    else
        @inbounds mse = sum(abs2(x[i] - y[i]) for i in 1:m) / m
    end
    return mse
end
mse(x, y::Real) = mse(x, fill(y, length(x)))

"""
    nrmse(x, y) → e
Return the normalized root mean square error of the "fit" `y` into data `x`.
This number is the relative error of `y` to `x` versus `mean(x)` to `x`, i.e.
if `e < 1` the fit `y` is better than using `mean(x)` as a fit.
"""
nrmse(x,y) = nrmse(vec(x), vec(y))
function nrmse(x::AbstractVector, y::AbstractVector)
    m = length(x)
    if any(isnan, x)
        j = findall(!isnan, x)
        mean_out = mean(x[j])
        _mse     = mse(x[j], y[j])
        @inbounds msemean = sum(abs2(x[i] - mean_out) for i in j) / length(j)
    else
        mean_out = mean(x)
        _mse     = mse(x, y)
        @inbounds msemean = sum(abs2(x[i] - mean_out) for i in 1:m) / m
    end
    nrmse   = sqrt(_mse / msemean)
    return nrmse
end

μrmse(x, y) = sqrt(mse(x, y))/mean(dropnan(x))

###########################################################################################
# # Related with time dimension
###########################################################################################
"""
    seasonal_nrmse(C, M, MASK = tues(C))
Return the median of the nrmses of the seasonal cycles of the two spatiotemporal fields.
Seasonal cycles are computed over the given keyword `latzones`.
The function removes the mean of each seasonal cycle, hence only caring about similarity
in the seasonal variability, not its mean.
"""
function seasonal_nrmse(C, M, MASK = trues(C); latzones = (-90, -30, 0, 30, 90))
    timeseries_errors = Float64[]
    for j in 1:length(latzones)-1
        l1, l2 = latzones[j], latzones[j+1]
        Csel = C[Coord(Lat(l1..l2))]
        Wsel = MASK[Coord(Lat(l1..l2))]
        Msel = M[Coord(Lat(l1..l2))]
        Csel = spacemean(Csel, Wsel)
        Msel = spacemean(Msel, Wsel)
        demeaned = []
        for (n, out) in enumerate((Csel, Msel))
            dates, vals = seasonality(out)
            m = mean.(vals)
            m = m .- mean(m)
            push!(demeaned, m)
        end
        push!(timeseries_errors, nrmse(demeaned[1], demeaned[2]))
    end
    return median(timeseries_errors)
end

"""
    correlationmap(C, M, MASK = trues(C)) → X, mean_cor
Calculate the correlation map of the two spatiotemporal fields based on `Statistics.cor`.
Return the correlation map and the mean correlation while dropping `NaN`.
"""
function correlationmap(F, M, MASK = trues(F))
    FF = copy(F); MM = copy(M)
    # remove NaNs
    X = timemean(MM)
    for i in spatialidxs(X)
        X[i] = Statistics.cor(FF[i...], MM[i...])
    end
    ZZZ = timemean(F, MASK)
    i = findall(isnan, ZZZ)
    X[i] .= NaN
    aaa = spacemean(dropnan(X))
    return X, aaa
end

###########################################################################################
# # Related with space and overall errors
###########################################################################################
"""
    spatial_nrmse(C, M, MASK = trues(C))
Return the spatial error (i.e., error between time-averaged version of the spatiotemporal
fields `C, M`).
"""
function spatial_nrmse(C, M, MASK = trues(C))
    Cmap = timemean(C, MASK)
    Mmap = timemean(M, MASK)
    return nrmse(Cmap, Mmap)
end

"""
    main_error_measures(C, M, MASK = trues(C))
Return the main error measures between the spatiotemporal fields `C, M` which are
used in the study. Specifically:
1. Full NRMSE
1. Spatial NRMSE
1. Median of seasonal NRMSEs
1. Mean of correlation map
In all cases, the given mask is properly taken into account.
All three fields must have identical dimensional layout.
"""
function main_error_measures(C, M, MASK = trues(C))
    fullerr = nrmse(C[gnv(MASK)], M[gnv(MASK)])
    spaterr = spatial_nrmse(C, M, MASK)
    timeerr = seasonal_nrmse(C, M, MASK)
    correrr = correlationmap(C, M, MASK)[2]
    return (fullerr, spaterr, timeerr, correrr)
end
