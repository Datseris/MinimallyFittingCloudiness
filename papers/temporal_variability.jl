using DrWatson
@quickactivate "AlbedoBounds"
include(scriptsdir("predictors", "cloudiness_predictors_definition.jl"));


albedo_expression = "p[1]*x1 + p[2]*x2*(1 - x1)"
albedo_predictors = (:Ω_nf, :ECTEI)
albedo_params = [49.3111, 1.94444]
lcre_expression = "p[1]*x1 + p[2]*x2 + p[3]*x3"
lcre_predictors = (:Ω_mean, :Ω_std, :q700)
lcre_params = [80.7788, 113.346953, 3.1855]
fieldnames = ("\$C\$", "\$L\$")
ocean_mask_perc = 50 # points with % ≥ than this are considered "ocean"
latzones = (-90, -30, 0, 30, 90)

# %%
OCEAN_MASK = O .≥ ocean_mask_perc
function eval_model_equations(model_expression, N)
    parsed_expression = Meta.parse(model_expression)
    _vars_expand = [Symbol("x", "$i") for i in 1:N]
    eval(:(model(p, $(_vars_expand...), args...) = @. $(parsed_expression)))
end

function get_model_instance(expression, predictors, p)
    eval_model_equations(expression, length(predictors))
    Ps = map(p -> getindex(field_dictionary,p), predictors)
    return model(p, Ps...)
end


fig = figure()

for (i, F) in enumerate((C, L))

# Create model
expression = (albedo_expression, lcre_expression)[i]
predictors = (albedo_predictors, lcre_predictors)[i]
params = (albedo_params, lcre_params)[i]
M = get_model_instance(expression, predictors, params)

# Correlation map plot
ax = subplot2grid((8, 2), (1 - 1, i - 1); rowspan = 4, projection=DEFPROJ)
# spatial correlation plot
FF = copy(F); MM = copy(M)
FF[.!OCEAN_MASK] .= NaN
MM[.!OCEAN_MASK] .= NaN
X = timemean(MM)
for i in spatialidxs(X)
    X[i] = Statistics.cor(FF[i...], MM[i...])
end
aaa = spacemean(dropnan(X))
earthscatter!(ax, X; s = 30, vmin = -1, vmax = 1, cmap = "PRGn", add_colorbar = false)
ax.set_title("$(fieldnames[i]), avg. corr = $(rdspl(aaa, 3))")

# Hemispheric timeseries (new)
latzones = (-90, -30, 0, 30, 90)

for j in 1:length(latzones)-1

    ax = subplot2grid((8, 2), (4 - 1 + 5 - j, i - 1))

    l1, l2 = latzones[j], latzones[j+1]
    Csel = F[Coord(Lat(l1..l2))]
    Wsel = OCEAN_MASK[Coord(Lat(l1..l2))]
    Msel = M[Coord(Lat(l1..l2))]
    Csel = spacemean(Csel, Wsel)
    Msel = spacemean(Msel, Wsel)

    for (n, out) in enumerate((Csel, Msel))    
        dates, vals = seasonality(out)
        m = mean.(vals)
        # m = m .- mean(m)
        push!(m, m[1])
        v = std.(vals)
        push!(v, v[1])
        ax.plot(1:13, m; color = "C$(n-1)", lw = 2, label = n == 1 ? "F" : "M",
        ls = n == 1 ? "-" : "--")
        ax.fill_between(1:13, m.-v, m.+v; color = "C$(n-1)", alpha = 0.5)
    end
    ax.set_xticks(1:3:13)
    if j ≠ 1
        setp(ax.get_xticklabels(), visible=false)
    end
end

end


# wsave(papersdir("plots", basename(@__FILE__)[1:end-3]), fig)
