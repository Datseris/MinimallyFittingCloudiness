using DrWatson
@quickactivate "MinimallyFittingCloudiness"
include(scriptsdir("fields_definitions.jl"));
include(srcdir("fitting", "masking.jl"))
include(srcdir("fitting", "general.jl"))

predictors = (:Ω_nf, :ECTEI)
predicted = :C
Φname = "\$C\$"
expression = "p[1]*x1 + p[2]*x2*(1 - x1)"
params = [48.88837, 2.08326]
MAXDEG = 80 # fit only within ± MAXDEG
ocean_mask_perc = 50 # points with % ≥ than this are considered "ocean"
limits = (0, 40)
function generate_contributions(p, Ps)
    c1 = p[1]*Ps[1]
    z = 1 .- Ps[1]
    c2 = p[2]*Ps[2]
    c3 = p[2] .* Ps[2] .* Ps[1]
    # return c1, c3, Ps[2]
    return c1, c2, c3
end
contribution_titles = ("\$p_1 U\$", "\$p_2 I\$", "\$- p_2 I U\$")

# %% The rest of this will be made on script
PROJ = ccrs.Mollweide()
close("all")
segment = 20
font_size = 18
title_size = 22
segment_gap = 4
latzones = (-90, -30, 0, 30, 90)

function add_colorbar!(ax_cbar, cmap, lvls = nothing; set_ticks = true)
    norm = matplotlib.colors.Normalize(vmin=lvls[1], vmax=lvls[end])
    cb = matplotlib.colorbar.ColorbarBase(ax_cbar, cmap, norm, extend="both", orientation="horizontal")
    if set_ticks
        cb.set_ticks(lvls[1:max(1, levels÷5):end])
    end
    cb.ax.tick_params(labelsize=font_size)
end

fig = figure(figsize = (18, 12))
# Segment 1: maps of ceres, model, diff
axY = subplot2grid((3segment + 2segment_gap, 3), (0, 0), rowspan=segment-1, projection=PROJ)
axM = subplot2grid((3segment + 2segment_gap, 3), (0, 1), rowspan=segment-1, projection=PROJ)
axD = subplot2grid((3segment + 2segment_gap, 3), (0, 2), rowspan=segment-1, projection=PROJ)
ax_cbar = subplot2grid((3segment + 2segment_gap, 3), (segment-1, 0), colspan=2)
ax_cbar_diff = subplot2grid((3segment + 2segment_gap, 3), (segment-1, 2))
# Segment 2: maps of contributions
axc1 = subplot2grid((3segment + 2segment_gap, 3), (segment+segment_gap, 0), rowspan=segment-1, projection=PROJ)
axc1_cbar = subplot2grid((3segment + 2segment_gap, 3), (2segment + segment_gap-1, 0))
axc2 = subplot2grid((3segment + 2segment_gap, 3), (segment+segment_gap, 1), rowspan=segment-1, projection=PROJ)
axc2_cbar = subplot2grid((3segment + 2segment_gap, 3), (2segment + segment_gap-1, 1))
axc3 = subplot2grid((3segment + 2segment_gap, 3), (segment+segment_gap, 2), rowspan=segment-1, projection=PROJ)
axc3_cbar = subplot2grid((3segment + 2segment_gap, 3), (2segment + segment_gap-1, 2))
# Segment 3: time variability
ax_tropics = subplot2grid((3segment + 2segment_gap, 3), (2segment + 2segment_gap, 0), rowspan=segment)
ax_extra = subplot2grid((3segment + 2segment_gap, 3), (2segment + 2segment_gap, 1), rowspan=segment)
axcorr = subplot2grid((3segment + 2segment_gap, 3), (2segment +2segment_gap, 2), rowspan=segment-1, projection=PROJ)
axcorr_cbar = subplot2grid((3segment + 2segment_gap, 3), (3segment + 2segment_gap -1, 2))

fig.subplots_adjust(
    left = 0.01, right = 0.99, top = 0.99, wspace = 0.05, hspace = 0.2, bottom = 0.05
)

OCEAN_MASK = field_dictionary[:O] .> ocean_mask_perc
Y = field_dictionary[predicted]
Ps = map(P -> getindex(field_dictionary, P), predictors)
M = get_model_instance(expression, Ps, params)
c1, c2, c3 = generate_contributions(p, Ps)

latticks = [-70, -30, -0, 30, 70]
levels = 21
vmin, vmax = limits
lvls = range(vmin, vmax; length = levels)
cmap = matplotlib.cm.get_cmap(:YlGnBu_r, length(lvls)-1)

Ymap = timemean(Y, OCEAN_MASK)
Mmap = timemean(M, OCEAN_MASK)
coords = gnv(dims(Ymap, Coord))
lon = [l[1] for l in coords]
lat = [l[2] for l in coords]

sckwargs = (
    transform = LONLAT, cmap, vmin, vmax, s = 7,
)
axY.scatter(lon, lat; c = gnv(Ymap), sckwargs...)
axY.set_title(Φname*", CERES", size = title_size)
axM.scatter(lon, lat; c = gnv(Mmap), sckwargs...)
axM.set_title(Φname*", MODEL", size = title_size)
add_colorbar!(ax_cbar, cmap, lvls)
# Difference map
vdiff = (vmax - vmin)/4
vmin = -vdiff; vmax = vdiff
lvls = range(vmin, vmax, length = levels)
cmap = matplotlib.cm.get_cmap(:PRGn, length(lvls)-1)
sckwargs = (
    transform = LONLAT, cmap, vmin, vmax, s = 7,
)
axD.set_title(Φname*", CERES - MODEL", size = title_size)
axD.scatter(lon, lat; c = gnv(Ymap .- Mmap), sckwargs...)
add_colorbar!(ax_cbar_diff, cmap, lvls)

### Segment 2
for i in 1:3
    c = (c1, c2, c3)[i]
    conmap = timemean(c, OCEAN_MASK)
    axc = (axc1, axc2, axc3)[i]
    axc_cbar = (axc1_cbar, axc2_cbar, axc3_cbar)[i]
    lvls = range(minimum(dropnan(conmap)), maximum(dropnan(conmap)); length = levels)
    cmap = matplotlib.cm.get_cmap(:viridis, length(lvls)-1)
    sckwargs = (
        transform = LONLAT, cmap, s = 7,
    )
    axc.scatter(lon, lat; c = gnv(conmap), sckwargs...)
    axc.set_title(contribution_titles[i], size = title_size)
    add_colorbar!(axc_cbar, cmap, lvls; set_ticks = false)
end

### Segment 3
# Temporal correlation map
cormap, aaa = correlationmap(Y, M, OCEAN_MASK)
earthscatter!(axcorr, cormap; s = 30, vmin = -1, vmax = 1, cmap = "RdBu", add_colorbar = false, levels)
axcorr.set_title("Pearson corr., avg.=$(rdspl(aaa, 3))", size = title_size)
cmap = matplotlib.cm.get_cmap(:RdBu, levels-1)
add_colorbar!(axcorr_cbar, cmap, range(-1, 1; length = levels); set_ticks = true)
# Seasonal timeseries

for j in 1:length(latzones)-1

    ax = j ∈ (1, 4) ? ...

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


# for all map plots
for ax in (axY, axM, axD, axc1, axc2, axc3, axcorr)
    ax.coastlines()
    gl = ax.gridlines(alpha = 0.5)
    gl.ylocator = matplotlib.ticker.FixedLocator(latticks)
    gl.xlocator = matplotlib.ticker.FixedLocator([-90, 0, 90])
end
