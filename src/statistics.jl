using Statistics, Dates, StatsBase, Random
export linreg, regularize, nrmse

"""
    linreg(x, y) -> s, o, fit, x
Return best linear fit so that y ≈ s*x + o.
If `x` is a date vector, `1:length(x)` is used instead.
"""
linreg(x::AbstractArray{<:Dates.AbstractDateTime}, y) = linreg(1:length(x), y)
function linreg(x, y)
    X = ones(length(x), 2)
    X[:, 1] .= x
    c = X\y
    s, o = c
    f = s .* x .+ o
    return s, o, f, x
end

regularize(x) = regularize!(copy(x))
regularize!(x) = (x .= (x .- mean(x))./std(x))


export nrmse, rmse
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

using Statistics: mean

"""
    rmse(x, y) → e
Return the root mean square error `e` of the "fit" `y` into data `x`.
"""
rmse(x, y) = sqrt(mse(x, y))

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



"Convenience function to compute seasonal nrmse as the median of four zones."
function seasonal_nrmse(C, M, OCEAN_MASK = ones(C); latzones = (-90, -30, 0, 30, 90),
    subtract_mean = true)
    timeseries_errors = Float64[]
    for j in 1:length(latzones)-1
        l1, l2 = latzones[j], latzones[j+1]
        Csel = C[Coord(Lat(l1..l2))]
        Wsel = OCEAN_MASK[Coord(Lat(l1..l2))]
        Msel = M[Coord(Lat(l1..l2))]
        Csel = spacemean(Csel, Wsel)
        Msel = spacemean(Msel, Wsel)
        if subtract_mean
            Csel = Csel .- mean(Csel)
            Msel = Msel .- mean(Msel)
        end
        push!(timeseries_errors, nrmse(Csel, Msel))
    end
    return median(timeseries_errors)
end

function correlationmap(F, M, OCEAN_MASK)
    FF = copy(F); MM = copy(M)
    FF[.!OCEAN_MASK] .= NaN
    MM[.!OCEAN_MASK] .= NaN
    X = timemean(MM)
    for i in spatialidxs(X)
        X[i] = Statistics.cor(FF[i...], MM[i...])
    end
    aaa = spacemean(dropnan(X))
    return X, aaa
end


using StatsBase
# Oh damn, `monthlyagg(Ω, skewness)` does not work, both because DimData.jl doesn't extend
# it, but also because StatsBase.jl doesn't have a `dims` method for it.
# So I wrote this custom function instead
function monthlyagg_custom(A, f=skewness; mday = 15)
    t0 = dims(A, Time).val
    startdate = Date(year(t0[1]), month(t0[1]), mday)
    finaldate = Date(year(t0[end]), month(t0[end]), mday+1)
    t = startdate:Month(1):finaldate
    tranges = temporalrange(t0, Dates.month)
    other = otherdims(A, Time)
    n = A.name == Symbol("") ? A.name : Symbol(A.name, ", $(string(f))")
    B = ClimArray(zeros(eltype(A), length.(other)..., length(t)), (other..., Time(t)), n)
    for i in 1:length(tranges)
        for j in otheridxs(A, Time)
            B[Time(i), j...] = f(Array(view(A, Time(tranges[i]), j...)))
        end
    end
    return B
end

function to_zonalmean(X)
	if hasdim(X, Lat) && length(dims(X)) == 1
		return X
	elseif hasdim(X, Coord)
		Z = collapse(mean, X, Coord)
		return zonalmean(Z)
	elseif hasdim(X, Lat)
		Z = collapse(mean, X, Lat)
		return Z
	end
end
