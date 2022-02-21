import LsqFit
# Notice that the tolerance related keywords for LsqFit can be found here:
# https://github.com/JuliaNLSolvers/LsqFit.jl/blob/master/src/levenberg_marquardt.jl#L16-L28

# Dispatch helpers
function multimodel(x::Matrix, p, model)
    return model(p, eachcol(x)...)
end

function multimodel(x::Vector{<:AbstractArray}, p, model)
    return model(p, x...)
end

# TODO: There should be only one call to `curve_fit`. The `perform_fit`
# should dispatch accordingly and not be duplicated.

function perform_fit(model, c::ClimArray, xs, p0;
        lower = fill(-Inf, length(p0)),
        upper = fill(Inf,  length(p0)),
        n = 10,
    )
    fitter = (x, p) -> multimodel(x, p, model)
    M = ClimArray(copy(c); name = "Model fit", attrib = copy(p0))
    xdata = hcat([vec(x) for x in xs]...)
    ydata = vec(c)
    err = Inf
    bestp = p0
    for m in 1:n
        p = m == 1 ? p0 : random_params(p0, lower, upper)
        modelfit = LsqFit.curve_fit(fitter, xdata, ydata, p; upper, lower)
        pfit = modelfit.param
        cerr = LsqFit.mse(modelfit)
        if cerr < err
            err = cerr
            bestp = pfit
        end
        # if m > 1
        #     @show bestp
        #     @show err
        # end
    end
    M.attrib .= bestp
    M .= model(bestp, xs...)
    M, nrmse(c, M), bestp
end

function perform_fit(model, c::Vector, xs, p0;
        lower = fill(-Inf, length(p0)),
        upper = fill(Inf,  length(p0)),
        n = 10,
    )
    fitter = (x, p) -> multimodel(x, p, model)
    M = copy(c)
    xdata = hcat(xs...)
    ydata = vec(c)
    err = Inf
    bestp = p0
    for m in 1:n
        p = m == 1 ? p0 : random_params(p0, lower, upper)
        modelfit = LsqFit.curve_fit(fitter, xdata, ydata, p; upper, lower, maxIter = 500)
        pfit = modelfit.param
        cerr = LsqFit.mse(modelfit)
        if cerr < err
            err = cerr
            bestp = pfit
        end
        # if m > 1
        #     @show p
        #     @show bestp
        #     @show err
        # end
    end
    M .= model(bestp, xs...)
    M, nrmse(c, M), bestp
end

function random_params(p0, lower, upper)
    p = copy(p0)
    for (i, (l, u)) in enumerate(zip(lower, upper))
        if !isinf(lower[i]) && !isinf(upper[i])
            p[i] = rand()*(u-l) + l
        end
    end
    return p
end


# here we fit two models sequentially; the second at the residuals
function perform_fit(model1, model2, c, xs, p0;
        ub1 = fill(Inf,  length(p0)),
        ub2 = fill(Inf,  length(p0)),
        lb1 = fill(-Inf, length(p0)),
        lb2 = fill(-Inf, length(p0)),
        p2 = copy(p0),
    )
    M1, err1, p1 = perform_fit(model1, c, xs, p0; lower=lb1, upper=ub1)
    R = c .- M1
    M2, err2, p2 = perform_fit(model2, R, Ps, p2; lower=lb2, upper=ub2)
    M = M1 .+ M2
    return M, nrmse(c,M), p1, p2
end
