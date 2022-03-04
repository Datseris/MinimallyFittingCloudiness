using DrWatson
@quickactivate "MinimallyFittingCloudiness"
include(scriptsdir("fields_definitions.jl"));
include(srcdir("fitting", "masking.jl"))
include(srcdir("fitting", "general.jl"))

# %%
input_fields = (C, L)
input_predictors = ((:Ω_nf, :ECTEI), (:Ω_std, :Ω_mean, :q700))
input_expressions = ("p[1]*x1 + p[2]*x2*(1-x1)", "p[1]*x1 + p[2]*x2 + p[3]*x3")
fieldnames = ["\$C\$",  "\$L\$"]
input_params = [[48.88837, 2.08326], [90.63965, 0.05912, 2.69408]]
vmins = [0, 0]
vmaxs = [40, 70]
MAXDEG = 80 # fit only within ± MAXDEG
ocean_mask_perc = 50 # points with % ≥ than this are considered "ocean"
PROJ = ccrs.Mollweide()
levels = 21
latzones = (-90, -30, 0, 30, 90)

# MAIN #########################
OCEAN_MASK = O .≥ ocean_mask_perc

close("all")
# Initialize figure and axis for temporal results
lzl = length(latzones) - 1
map_height = 12 # colorbar is always 1 height
zonal_height = 5 # we have four of these
rows = map_height + 1 + lzl*zonal_height
figtemp = figure(figsize = (18, 12))
axCtemp = subplot2grid((rows, 2), (0, 0), rowspan=map_height, fig = figtemp, projection = PROJ)
axLtemp = subplot2grid((rows, 2), (0, 1), rowspan=map_height, fig = figtemp, projection = PROJ)
ax_cbar_temp = subplot2grid((rows, 2), (map_height+1, 0), colspan=2, fig = figtemp)
axCtemp_zones = []
axLtemp_zones = []
figtemp.subplots_adjust(left = 0.06, right = 0.98, top = 0.98, bottom = 0.05, wspace = 0.1, hspace = 0.5)

for j in 1:length(latzones)-1
    for (k, cont) in enumerate((axCtemp_zones, axLtemp_zones))
        ax = subplot2grid((rows, 2), (rows - j*zonal_height + 1, k-1);
        rowspan=zonal_height, fig=figtemp)
        # ax.text(0.5, 0.5, "latzone=$(latzones[j])-$(latzones[j+1])")
        j != 1 && setp(ax.get_xticklabels(), visible=false)
        push!(cont, ax)
    end
end
coraxis = (axCtemp, axLtemp)
seaaxis = (axCtemp_zones, axLtemp_zones)


for (i, Φ) in enumerate(input_fields)
predictors = input_predictors[i]
expression = input_expressions[i]
Φname = fieldnames[i]
params = input_params[i]
Ps = map(P -> getindex(field_dictionary, P), predictors)
M = get_model_instance(expression, Ps, params)

# Initialize figure and axis for spatial results
rows = 20
fig = figure(figsize = (18, 5))
axΦ = subplot2grid((rows, 3), (0, 0), rowspan=rows-1, projection=PROJ)
axM = subplot2grid((rows, 3), (0, 1), rowspan=rows-1, projection=PROJ)
axD = subplot2grid((rows, 3), (0, 2), rowspan=rows-1, projection=PROJ)
ax_cbar = subplot2grid((rows, 3), (rows-1, 0), colspan=2)
ax_cbar_diff = subplot2grid((rows, 3), (rows-1, 2), colspan=2)

# Actual plotting
# Spatial maps
latticks = [-70, -30, -0, 30, 70]

vmin = vmins[i]
vmax = vmaxs[i]
lvls = range(vmin, vmax, length = levels)
cmap = matplotlib.cm.get_cmap(:YlGnBu_r, length(lvls)-1)

Φmap = timemean(Φ, OCEAN_MASK)
Mmap = timemean(M, OCEAN_MASK)
coords = gnv(dims(Φmap, Coord))
lon = [l[1] for l in coords]
lat = [l[2] for l in coords]


sckwargs = (
    transform = LONLAT, cmap, vmin, vmax, s = 7,
)
axΦ.scatter(lon, lat; c = gnv(Φmap), sckwargs...)
axM.scatter(lon, lat; c = gnv(Mmap), sckwargs...)
axΦ.set_title(Φname*", CERES")
axM.set_title(Φname*", MODEL")
# Colorbar
norm = matplotlib.colors.Normalize(vmin=lvls[1], vmax=lvls[end])
cb = matplotlib.colorbar.ColorbarBase(ax_cbar, cmap, norm, extend="both", orientation="horizontal")
cb.set_ticks(lvls[1:max(1, levels÷5):end])

# Difference map
vdiff = (vmax - vmin)/4
vmin = -vdiff; vmax = vdiff
lvls = range(vmin, vmax, length = levels)
cmap = matplotlib.cm.get_cmap(:YlGnBu_r, length(lvls)-1)
cmap = matplotlib.cm.get_cmap(:PRGn, length(lvls)-1)
sckwargs = (
    transform = LONLAT, cmap, vmin, vmax, s = 7,
)
axD.set_title(Φname*", CERES - MODEL")
norm = matplotlib.colors.Normalize(vmin=lvls[1], vmax=lvls[end])
cb = matplotlib.colorbar.ColorbarBase(ax_cbar_diff, cmap, norm, extend="both", orientation="horizontal")
cb.set_ticks(lvls[1:max(1, levels÷5):end])
axD.scatter(lon, lat; c = gnv(Φmap .- Mmap), sckwargs...)

for ax in (axΦ, axM, axD)
    ax.coastlines()
    gl = ax.gridlines(alpha = 1)
    gl.ylocator = matplotlib.ticker.FixedLocator(latticks)
    gl.xlocator = matplotlib.ticker.FixedLocator([-90, 0, 90])
    # gl.xlines = false
end

fig.subplots_adjust(left = 0.01, right = 0.99, top = 0.99, wspace = 0.05)

add_identifiers!(fig, (axΦ, axM, axD))
wsave(papersdir("plots", "results_$(i)"), fig)

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
    j == 1 && ax.set_xticklabels(["JAN", "APR", "JUL", "OCT", "JAN"])
end

end

# Spacing of temporal plot
figtemp.subplots_adjust(left = 0.06, right = 0.98, top = 0.97, bottom = 0.05, wspace = 0.1, hspace = 0.5)
wsave(papersdir("plots", "results_temporal"), figtemp)

