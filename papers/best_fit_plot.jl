PyPlot.ioff()
PROJ = ccrs.Mollweide()
segment = 20
font_size = 22
title_size = 32
segment_gap = 8
latzones = (-90, -30, 0, 30, 90)

function add_colorbar!(ax_cbar, cmap, lvls = nothing; set_ticks = true)
    norm = matplotlib.colors.Normalize(vmin=lvls[1], vmax=lvls[end])
    cb = matplotlib.colorbar.ColorbarBase(ax_cbar, cmap, norm, extend="both", orientation="horizontal")
    if set_ticks
        cb.set_ticks(lvls[1:max(1, levels÷5):end])
    end
    cb.ax.tick_params(labelsize=font_size)
end

fig = figure(figsize = (18, 14))
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

OCEAN_MASK = (field_dictionary[:O] .> ocean_mask_perc)[Coord(Lat((-MAXDEG)..(MAXDEG)))]
Y = field_dictionary[predicted][Coord(Lat((-MAXDEG)..(MAXDEG)))]
Xs = map(X -> getindex(field_dictionary, X)[Coord(Lat((-MAXDEG)..(MAXDEG)))], predictors)
M = get_model_instance(expression, Xs, params)
c1, c2, c3 = generate_contributions(params, Xs)

latticks = [-70, -30, -0, 30, 70]
levels = 21
vmin, vmax = limits
lvls = range(vmin, vmax; length = levels)
cmap = matplotlib.cm.get_cmap(:viridis, length(lvls)-1)

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
    cvmin = minimum(dropnan(conmap))
    cvmax = maximum(dropnan(conmap))
    ccmap = contrib_cmaps[i]
    if ccmap == :BrBG
        cvmax = max(cvmax, -cvmin)
        cvmin = -cvmax
    end
    lvls = range(cvmin, cvmax; length = levels)
    cmap = matplotlib.cm.get_cmap(ccmap, length(lvls)-1)
    sckwargs = (
        transform = LONLAT, cmap, s = 7, vmin = minimum(lvls), vmax = maximum(lvls),
    )
    axc.scatter(lon, lat; c = gnv(conmap), sckwargs...)
    axc.set_title(contribution_titles[i], size = title_size)
    add_colorbar!(axc_cbar, cmap, lvls; set_ticks = false)
end

### Segment 3
# Temporal correlation map
cormap, aaa = correlationmap(Y, M, OCEAN_MASK)
earthscatter!(axcorr, cormap; s = 7, vmin = -1, vmax = 1, cmap = "RdBu", add_colorbar = false, levels)
axcorr.set_title("Pearson corr., avg.=$(rdspl(aaa, 3))", size = title_size)
cmap = matplotlib.cm.get_cmap(:RdBu, levels-1)
add_colorbar!(axcorr_cbar, cmap, range(-1, 1; length = levels); set_ticks = true)

# Seasonal timeseries
for j in 1:length(latzones)-1

    ax = j ∈ (1, 4) ? ax_extra : ax_tropics
    ax.set_title(j ∈ (1, 4) ? "extratropics" : "tropics"; size = title_size)
    cna = isodd(j) ? 0 : 2 # color n add

    l1, l2 = latzones[j], latzones[j+1]
    Ysel = Y[Coord(Lat(l1..l2))]
    Wsel = OCEAN_MASK[Coord(Lat(l1..l2))]
    Msel = M[Coord(Lat(l1..l2))]
    Ysel = spacemean(Ysel, Wsel)
    Msel = spacemean(Msel, Wsel)

    for (n, out) in enumerate((Ysel, Msel))
        label = (n == 2 ? "MODEL" : "CERES") * " " * (j < 3 ? "SH" : "NH")
        dates, vals = seasonality(out)
        m = mean.(vals)
        m = m .- mean(m)
        isodd(j) && (m .+= seasonal_offset)
        push!(m, m[1])
        v = std.(vals)
        push!(v, v[1])
        ax.plot(1:13, m; color = "C$(n-1 + cna)", lw = 2, label,
        ls = n == 1 ? "-" : "--")
        ax.fill_between(1:13, m.-v, m.+v; color = "C$(n-1 + cna)", alpha = 0.5)
    end
    ax.legend(ncol = 2, fontsize = 22,
    handlelength = 1, handletextpad = 0.5, )
    ax.set_xlim(1,13)
    ax.set_ylim(seasonal_limits)
    ax.set_xticks(2:3:13)
    # ax.set_xticklabels(["JAN", "APR", "JUL", "OCT", "JAN"], size = font_size)
    ax.set_xticklabels(["FEB", "MAY", "JUL", "NOV"], size = font_size)
    # ax.locator_params(axis='y', nbins=6)
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    ax.tick_params(axis="y"; labelsize=font_size)
    j ∈ (1, 4) && ax.axes.yaxis.set_ticklabels([])
end

# for all map plots
for ax in (axY, axM, axD, axc1, axc2, axc3, axcorr)
    ax.coastlines()
    gl = ax.gridlines(alpha = 0.5)
    gl.ylocator = matplotlib.ticker.FixedLocator(latticks)
    gl.xlocator = matplotlib.ticker.FixedLocator([-90, 0, 90])
end

fig.subplots_adjust(
    left = 0.03, right = 0.99, top = 0.97, wspace = 0.05, hspace = 0.5, bottom = 0.05
)

wsave(papersdir("plots", "results_$(predicted).pdf"), fig; transparent = false, dpi = 100)
PyPlot.ion()