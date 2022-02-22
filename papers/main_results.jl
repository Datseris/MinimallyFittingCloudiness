using DrWatson
@quickactivate "AlbedoBounds"
include(scriptsdir("predictors", "cloudiness_predictors_definition.jl"));
include(srcdir("modelling", "general.jl"))

# %%
input_fields = (C, CRElw)
input_predictors = ((:Ω_nf, :ECTEI), (:Ω_mean, :Ω_std, :q700))
input_expressions = ("p[1]*x1 + p[2]*x2*(1-x1)", "p[1]*x1 + p[2]*x2 + p[3]*x3")
fieldnames = ["\$C\$",  "\$L\$"]
input_params = [[49.3111, 1.94444], [80.7788, 113.346953, 3.1855]]
vmins = [0, 0]
vmaxs = [40, 70]
MAXDEG = 80 # fit only within ± MAXDEG
ocean_mask_perc = 50 # points with % ≥ than this are considered "ocean"
PROJ = ccrs.Mollweide()
levels = 21
latzones = (-90, -30, 0, 30, 90)

# MAIN #########################
OCEAN_MASK = O .≥ ocean_mask_perc
function eval_model_equations(model_expression, N)
    parsed_expression = Meta.parse(model_expression)
    _vars_expand = [Symbol("x", "$i") for i in 1:N]
    eval(:(model(p, $(_vars_expand...), args...) = @. $(parsed_expression)))
end

function get_model_instance(expression, predictors, p)
    eval_model_equations(expression, length(predictors))
    Ps = map(P -> getindex(field_dictionary, P), predictors)
    return model(p, Ps...)
end

function timezonalmean(Φ, OCEAN_MASK, MAXDEG)
    to_non_nan_lats(
        timemean(
            zonalmean(Φ[Coord(Lat((-MAXDEG)..(MAXDEG)))], OCEAN_MASK[Coord(Lat((-MAXDEG)..(MAXDEG)))])
        )
    )
end

function correlationmap(F, M)
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
    
close("all")
# Initialize figure and axis for temporal results
rows = 32 # colobar is always 1 row length. rows must be divisible by 2(length(latzones)-1)
figtemp = figure(figsize = (18, 10))
axCtemp = subplot2grid((rows, 2), (0, 0), rowspan=rows÷2-1, fig = figtemp)
axLtemp = subplot2grid((rows, 2), (0, 1), rowspan=rows÷2-1, fig = figtemp)
ax_cbar_temp = subplot2grid((rows, 2), (rows÷2, 0), colspan=2, fig = figtemp)
axCtemp_zones = []
axLtemp_zones = []
remrows = rows÷2
rowhei = remrows÷(length(latzones)-1)
for j in 1:length(latzones)-1
    for (k, cont) in enumerate((axCtemp_zones, axLtemp_zones))
        ax = subplot2grid((rows, 2), (rows - j*rowhei + 1, k-1); rowspan=rowhei, fig=figtemp)
        # ax.text(0.5, 0.5, "latzone=$(latzones[j])-$(latzones[j+1])")
        j != 1 && setp(ax.get_xticklabels(), visible=false)
        push!(cont, ax)
    end
end
figtemp.tight_layout()
axCtemp = subplot2grid((rows, 2), (0, 0), rowspan=rows÷2-1, fig = figtemp, projection = PROJ)
axLtemp = subplot2grid((rows, 2), (0, 1), rowspan=rows÷2-1, fig = figtemp, projection = PROJ)
coraxis = (axCtemp, axLtemp)
seaaxis = (axCtemp_zones, axLtemp_zones)




for (i, Φ) in enumerate(input_fields)
predictors = input_predictors[i]
expression = input_expressions[i]
Φname = fieldnames[i]
params = input_params[i]
M = get_model_instance(expression, predictors, params)
Ps = map(P -> getindex(field_dictionary, P), predictors)

# Initialize figure and axis for spatial results
rows = 20
offset = 2 # for latitudes
fig = figure(figsize = (18, 5))
axΦ = subplot2grid((rows, 3), (0, 0), rowspan=rows-1)
axM = subplot2grid((rows, 3), (0, 1), rowspan=rows-1)
axZ = subplot2grid((rows, 3), (0, 2), rowspan=rows-offset)
ax_cbar = subplot2grid((rows, 3), (rows-1, 0), colspan=2)
fig.tight_layout()
axΦ = subplot2grid((rows, 3), (0, 0), rowspan=rows-1, projection=PROJ)
axM = subplot2grid((rows, 3), (0, 1), rowspan=rows-1, projection=PROJ)


# Actual plotting
# Spatial maps
latticks = [-70, -30, -0, 30, 70]

coords = gnv(dims(Φ, Coord))
lon = [l[1] for l in coords]
lat = [l[2] for l in coords]

vmin = vmins[i]
vmax = vmaxs[i]
lvls = range(vmin, vmax, length = levels)
cmap = matplotlib.cm.get_cmap(:YlGnBu_r, length(lvls)-1)

Φmap = timemean(Φ, OCEAN_MASK)
Mmap = timemean(M, OCEAN_MASK)

sckwargs = (
    transform = LONLAT, cmap, vmin, vmax, s = 6,
)
axΦ.scatter(lon, lat; c = gnv(Φmap), sckwargs...)
axM.scatter(lon, lat; c = gnv(Mmap), sckwargs...)
axΦ.set_title(Φname*", CERES")
axM.set_title(Φname*", MODEL")
axΦ.coastlines()
axM.coastlines()
for ax in (axΦ, axM)
    gl = ax.gridlines(alpha = 1)
    gl.ylocator = matplotlib.ticker.FixedLocator(latticks)
    gl.xlocator = matplotlib.ticker.FixedLocator([-90, 0, 90])
    # gl.xlines = false
end

# Colorbar
norm = matplotlib.colors.Normalize(vmin=lvls[1], vmax=lvls[end])
cb = matplotlib.colorbar.ColorbarBase(ax_cbar, cmap, norm, extend="both", orientation="horizontal")
cb.set_ticks(lvls[1:max(1, levels÷5):end])

# zonal plot
Ftz = timezonalmean(Φ, OCEAN_MASK, MAXDEG)
Pstz = map(P -> timezonalmean(P, OCEAN_MASK, MAXDEG), Ps)

pl = -maximum(Φ)*10ones(8) # lower bounds
pu = maximum(Φ)*10ones(8) # upper bounds
p0 = zero(pl)
Mtz, errtz, ptz, = perform_fit(model, Ftz, Pstz, p0; lower = pl, upper = pu, n = 100)
Mtzfull = timezonalmean(M, OCEAN_MASK, MAXDEG)


Φlats = sind.(gnv(dims(Ftz, Lat)))

axZ.plot(gnv(Ftz), Φlats; label = Φname, lw = 3)
axZ.plot(gnv(Mtz), Φlats; label = "\$M_z\$", ls = "-.", alpha = 0.7, lw = 2.5) 
axZ.plot(gnv(Mtzfull), Φlats; label = "\$M\$", ls = "--", alpha = 0.7, lw = 2.5)

axZ.set_yticks(sind.(latticks))
axZ.set_yticklabels(string.(latticks) .* "ᵒ")
axZ.legend(fontsize = 18)


fig.subplots_adjust(left = 0.01, right = 0.99, top = 0.99)

add_identifiers!(fig, (axΦ, axM, axZ))
# wsave(papersdir("plots", "results_$(i)"), fig)

# Temporal correlation map
cormap, aaa = correlationmap(Φ, M)
earthscatter!(coraxis[i], cormap; s = 30, vmin = -1, vmax = 1, cmap = "PRGn", add_colorbar = false, levels = 11)
coraxis[i].set_title("$(fieldnames[i]), avg. corr = $(rdspl(aaa, 3))")
if i == 1 # Colorbar
    lvls = range(-1, 1; length = 11)
    cmap = matplotlib.cm.get_cmap("PRGn", length(lvls)-1)
    norm = matplotlib.colors.Normalize(vmin=-1, vmax=1)
    cb = matplotlib.colorbar.ColorbarBase(ax_cbar_temp, cmap, norm, ticklocation="top",extend="both", orientation="horizontal")
    cb.set_ticks(lvls)
end
# Zonal mean plots
for j in 1:length(latzones)-1

    ax = seaaxis[i][j]

    l1, l2 = latzones[j], latzones[j+1]
    Csel = Φ[Coord(Lat(l1..l2))]
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
    ax.set_xlim(1,13)
end

end

# Spacing of temporal plot
figtemp.subplots_adjust(left = 0.06, right = 0.98, top = 0.98, bottom = 0.05, wspace = 0.1, hspace = 0.5)
wsave(papersdir("plots", "results_temporal"), figtemp)

