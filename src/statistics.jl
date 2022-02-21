using Statistics, Dates, Entropies, StatsBase, Random
export linreg, regularize, findnearest

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


export periodic_correlation
"""
	periodic_correlation(x, y, lags, demean=true)
Calculate the correlation function of periodic signals both of which have period
their length (e.g. functions of longitude).
"""
function periodic_correlation(x,y,τ, demean=true)
	if length(unique(x)) == length(unique(y)) == 1
		# Dividing by `std` if `x` is constant leads to `Inf` results.
		return ones(eltype(x), length(τ))
	end
    if demean
        x = x .- mean(x)
        y = y .- mean(y)
    end
    N = length(x)
    r = [sum(x[n]*y[mod1(n+m, N)] for n in 1:N)/N for m in τ]
    r ./ (std(x)*std(y)) # follow normalization of `StatsBase.crosscor`
end

"""
	relativeMI(x, y, ε, N = 1000)
Return the "relative" mutual information between `x, y` with binning `ε`.
The MI is measured with respect to the null hypothesis that `x, y` are uncorrelated, and
the number returned in fact is ``|m-\\mu|/2\\sigma`` with μ, σ the mean and standard
deviation of the null distribution of shuffled `x, y`. If the returned value is greater
than 1, then you can claim significant correlation.
"""
function relativeMI(x, y, ε, N = 1000)
	x, y = copy(Array(x)), copy(Array(y))
	Hx = genentropy(Dataset(x), ε)
    Hy = genentropy(Dataset(y), ε)
    Hxy = genentropy(Dataset(x,y), ε)
	mi = Hx + Hy - Hxy
	null = zeros(N)
	for i in 1:N
		shuffle!(x); shuffle!(y)
		Hxy = genentropy(Dataset(x,y), ε)
		null[i] = Hx + Hy - Hxy
	end
	μ, σ = mean(null), std(null)
	return abs(mi - μ)/(2σ)
end


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


using StatsBase
# Oh damn, `monthlyagg(Ω, skewness)` does not work, both because DimData.jl doesn't extend
# it, but also because StatsBase.jl doesn't have a `dims` method for it
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

"""
	to_non_nan_lats(A)
Given `A` that must have a latitude dimension, return `B` which has the latitudes of
`A` that *do not* contain `NaN` values. Notice that if `A` has a time dimension,
a time-average is done first to return latitudes that would be non-NaN even after
time averaging.
"""
function to_non_nan_lats(A)
	B = hasdim(A, Ti) ? timemean(A) : A
	is = findall(!isnan, B)
	return A[Lat(is)]
end

function dropnan(x)
	is = findall(!isnan, x)
	x[is]
end


using LinearAlgebra
"""
    findnearest(val, A)
Return the index of `A` which has value nearest to `val`.
"""
function findnearest(val, A)
    i = 1
    d = norm(val .- A[i])
    @inbounds for j in 1:length(A)
        dd = norm(val .- A[j])
        if dd < d
            i = j
            d = dd
        end
    end
    return i
end

"""
	hemispheric_means_lats(A::ClimArray, W = ones(A); l1 = 0, l2 = 90) → nh, sh
Same as `hemispheric_means` for a `Coord` array, but the averaging is further
limited to latitudes from `l1` to `l2`.
"""
function hemispheric_means_lats(A::ClimArray, W = ones(A); l1 = 0, l2 = 90)
	@assert size(A) == size(W)
    nhi, shi = ClimateBase.hemisphere_indices(A)
	i = findall(x -> l1 ≤ x[2] ≤ l2, dims(A, Coord).val[nhi])
	nhi = nhi[i]
	j = findall(x -> -l2 ≤ x[2] ≤ -l1, dims(A, Coord).val[shi])
	shi = shi[j]
    nh = dropagg(mean, A[Coord(nhi)], Coord, W[Coord(nhi)])
    sh = dropagg(mean, A[Coord(shi)], Coord, W[Coord(shi)])
    return nh, sh
end
