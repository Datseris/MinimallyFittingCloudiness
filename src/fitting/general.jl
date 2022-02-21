import LsqFit
# The tolerance related keywords for LsqFit can be found here:
# https://github.com/JuliaNLSolvers/LsqFit.jl/blob/master/src/levenberg_marquardt.jl#L16-L28

##########################################################################################
# Model generation
##########################################################################################
eval_model_equations(model_expression, p) = eval_model_equations(model_expression, length(p)

function eval_model_equations(model_expression, N::Int)
    parsed_expression = Meta.parse(model_expression)
    _vars_expand = [Symbol("x", "$i") for i in 1:N]
    eval(:(model(p, $(_vars_expand...), args...) = @. $(parsed_expression)))
end

function multimodel(x::Matrix, p, model)
    return model(p, eachcol(x)...)
end

function multimodel(x::Vector{<:AbstractArray}, p, model)
    return model(p, x...)
end

##########################################################################################
# Actual fitting
##########################################################################################
"""
    perform_fit(model, F, Xs, p0; kwargs...) → M, e, p
Fit the given `model` function, which takes as an input predictors `Xs`, against
data `F`, using nonlinear least squares optimization from LsqFit.jl.

Return the fitted model, the NRMSE, and the best parameter values.

The fit is the best ouf of `n` different and random initial parameter configurations
within the space allowed by `lower` and `upper`.
"""
function perform_fit(model, c, xs, p0;
        lower = fill(-Inf, length(p0)),
        upper = fill(Inf,  length(p0)),
        n = 10, show_info = false,
    )
    fitter = (x, p) -> multimodel(x, p, model)
    M = init_empty_model(c)
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
        if show_info
            @show (m, err)
            @show bestp
        end
    end
    M .= model(bestp, xs...)
    (M isa ClimArray) && (M.attrib .= bestp)
    return M, nrmse(c, M), bestp
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

init_empty_model(c::ClimArray) = ClimArray(copy(c); name = "Model fit", attrib = copy(p0))
init_empty_model(c) = copy(c)