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
