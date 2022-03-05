# This script does something rather simple:
# compares the NRMSE of our model versus that of ERA5
# For our model we have the error already in `general_model_fit.jl`

using DrWatson
@quickactivate "MinimallyFittingCloudiness"
include(scriptsdir("fields_definitions.jl"));
include(srcdir("fitting", "masking.jl"))
include(srcdir("fitting", "general.jl"))
include(srcdir("cloudiness.jl"))
MAXDEG = 80 # fit only within ± MAXDEG
ocean_mask_perc = 50 # points with % ≥ than this are considered "ocean"
OCEAN_MASK = field_dictionary[:O] .≥ ocean_mask_perc

# Some plotting helper functions
function plot_zonal_comparison!(CERES_tz, ERA5_tz, MODEL_tz, MODEL_tz2=nothing; ax = gca())
    alpha = 0.8
    latticks = [-70, -30, -0, 30, 70]
    Φlats = sind.(gnv(dims(CERES_tz, Lat)))
    e1 = nrmse(CERES_tz, ERA5_tz)
    e2 = nrmse(CERES_tz, MODEL_tz)
    ax.plot(Φlats, gnv(CERES_tz); label = "CERES", lw = 3)
    ax.plot(Φlats, gnv(ERA5_tz);  label = "ERA5, \$e\$=$(rdspl(e1, 3))", ls = "--", alpha, lw = 2.5)
    ax.plot(Φlats, gnv(MODEL_tz); label = "FIT, \$e\$=$(rdspl(e2, 3))", ls = "-.", alpha, lw = 2.5)
    if MODEL_tz2 !== nothing
        e3 = nrmse(CERES_tz, MODEL_tz2)
        ax.plot(Φlats, gnv(MODEL_tz2); label = "FTZ, \$e\$=$(rdspl(e3, 3))",
        c = "C4", ls = ":", alpha, lw = 2.5)
    end
    ax.set_xticks(sind.(latticks))
    ax.set_xticklabels(string.(latticks) .* "ᵒ")
end

function compare_field_maps(CERES_FIELDS, ERA5_FIELDS)
    PROJ = ccrs.Mollweide()

    coords = gnv(dims(CERES_FIELDS[1], Coord))
    lon = [l[1] for l in coords]
    lat = [l[2] for l in coords]
    rows = 20
    levels = 21

    for (C, E) in zip(CERES_FIELDS, ERA5_FIELDS)
        fig = figure(figsize = (14, 5))
        axC = subplot2grid((rows, 3), (0, 0), rowspan=rows-1, projection=PROJ)
        axE = subplot2grid((rows, 3), (0, 1), rowspan=rows-1, projection=PROJ)
        axD = subplot2grid((rows, 3), (0, 2), rowspan=rows-1, projection=PROJ)
        ax_cbar = subplot2grid((rows, 3), (rows-1, 0), colspan=2)
        ax_cbar_diff = subplot2grid((rows, 3), (rows-1, 2), colspan=2)
        axC.set_title(C.name)

        tmean(X) = hasdim(X, Time) ? timemean(X) : X
        Cmap, Emap = tmean.((C, E))
        Dmap = Cmap .- Emap
        s = 7
        vmin = minimum(Cmap); vmax = maximum(Cmap)
        lvls = range(vmin, vmax, length = levels)
        cmap = matplotlib.cm.get_cmap(:YlGnBu_r, length(lvls)-1)
        axC.scatter(lon, lat; c = gnv(Cmap), s, cmap, vmin, vmax, transform = LONLAT)
        axE.scatter(lon, lat; c = gnv(Emap), s, cmap, vmin, vmax, transform = LONLAT)
        norm = matplotlib.colors.Normalize(vmin=lvls[1], vmax=lvls[end])
        cb = matplotlib.colorbar.ColorbarBase(ax_cbar, cmap, norm,
        extend="both", orientation="horizontal")
        cb.set_ticks(lvls[1:max(1, levels÷5):end])
        cb.ax.tick_params(labelsize=18)

        vdiff = (vmax - vmin)/4
        vmin = -vdiff; vmax = vdiff
        lvls = range(vmin, vmax, length = levels)
        cmap = matplotlib.cm.get_cmap(:PRGn, length(lvls)-1)
        axD.scatter(lon, lat; c = gnv(Dmap), s, cmap, vmin, vmax, transform = LONLAT)
        norm = matplotlib.colors.Normalize(vmin=lvls[1], vmax=lvls[end])
        cb = matplotlib.colorbar.ColorbarBase(ax_cbar_diff, cmap, norm,
        extend="both", orientation="horizontal")
        cb.set_ticks(lvls[1:max(1, levels÷5):end])
        cb.ax.tick_params(labelsize=18)
        fig.subplots_adjust(left = 0.01, right = 0.99, top = 0.99, wspace = 0.05)
    end
end

mainfig = figure(figsize = (14, 5))

# %% #src
# # Comparison of L
# First, compute L from ERA5
L_CERES_tz = maskedtimezonalmean(L, OCEAN_MASK, MAXDEG)

# Compute L from ERA5:
Lall_ERA5 = -ncread(ERA5_2D, "ttr")/86400 # transfrom to W/m²; divide by seconds in day
Lclr_ERA5 = -ncread(ERA5_2D, "ttrc")/86400
L_ERA5 =   Lclr_ERA5 .- Lall_ERA5
CERES_FIELDS = (Lall, Lclr, L)
ERA5_FIELDS = (Lall_ERA5, Lclr_ERA5, L_ERA5)

# Alright, plot the comparison of all fields now
compare_field_maps(CERES_FIELDS, ERA5_FIELDS)

# Next, calculate our model fit, in two ways:
# 1. On full data (as always)
# 2. On zonally+temporally averaged data (just to make case that this does not matter much)

predictors = (:Ω_std, :Ω_mean, :q700)
model_expression = "p[1]*x1 + p[2]*x2 +p[3]*x3"
eval_model_equations(model_expression, predictors)

# (1):
Ps = map(p -> getindex(field_dictionary, p), predictors)
Ymasked = ocean_masked(field_dictionary[:L], OCEAN_MASK, MAXDEG)
Psmasked = [ocean_masked(P, OCEAN_MASK, MAXDEG) for P in Ps]
Mmasked, err, p, = perform_fit(model, Ymasked, Psmasked, ones(3); n = 2)
M1 = maskedtimezonalmean(model(p, Ps...), OCEAN_MASK, MAXDEG)

# (2):
Ytz = maskedtimezonalmean(L, OCEAN_MASK, MAXDEG)
Pstz = map(P -> maskedtimezonalmean(P, OCEAN_MASK, MAXDEG), Ps)
Mtz, = perform_fit(model, Ytz, Pstz, ones(3); n = 2)

# And finally, get ERA5
L_ERA5_tz = maskedtimezonalmean(L_ERA5, OCEAN_MASK, MAXDEG)

# and plot everything:
axL = mainfig.add_subplot(121)
plot_zonal_comparison!(Ytz, L_ERA5_tz, M1, Mtz; ax = axL)


# # Comparison of C
# %% src
# We can't compute the energetically consistent effective cloud albedo of Datseris &
# Stevens 2021directly from ERA5, it requires cloud optical depth which ERA5 does not have.
# Instead, we will use cloud contribution to the atmospheric contribution to albedo,
# i.e., equation (3) of Datseris & Stevens 2021.
# This means that we can't fit temporal variability.


# In any case, we first load all fields necessary:

CERES_FIELDS = begin
    I = ncread(EBAF, "solar_mon")
    R = ncread(EBAF, "toa_sw_all_mon")
    K = ncread(EBAF, "toa_sw_clr_t_mon")
    F_s_⬆ = ncread(EBAF, "sfc_sw_up_all_mon")
    F_s_⬇ = ncread(EBAF, "sfc_sw_down_all_mon")
    F_s_⬆_K = ncread(EBAF, "sfc_sw_up_clr_t_mon")
    F_s_⬇_K = ncread(EBAF, "sfc_sw_down_clr_t_mon");
    I, R, K, F_s_⬆, F_s_⬇, F_s_⬆_K, F_s_⬇_K
end

ERA5_FIELDS = begin
    I = ncread(ERA5_2D, "tisr")/86400 # convert to W/m²
    R = I .- ncread(ERA5_2D, "tsr")/86400 # convert to W/m²
    K = I .- ncread(ERA5_2D, "tsrc")/86400 # convert to W/m²
    F_s_⬇ = ncread(ERA5_2D, "ssrd")/86400 # convert to W/m²
    F_s_⬆ = F_s_⬇ .- ncread(ERA5_2D, "ssr")/86400
    F_s_⬆_K = ncread(EBAF, "sfc_sw_up_clr_t_mon") # already in W/m²
    F_s_⬇_K = ncread(ERA5_2D, "ssrdc")/86400
    I, R, K, F_s_⬆, F_s_⬇, F_s_⬆_K, F_s_⬇_K
end

# Then calculate the albedos
α_CERES = 100donohoe_battisti_cloud_albedo(CERES_FIELDS...)
α_ERA5 = 100donohoe_battisti_cloud_albedo(ERA5_FIELDS...)

# And then, plot spatial maps of all fields just to be sure everything makes sense
CERES_FIELDS = (CERES_FIELDS..., α_CERES)
ERA5_FIELDS = (ERA5_FIELDS..., α_ERA5)
compare_field_maps(CERES_FIELDS, ERA5_FIELDS)

# ## Fit and comparison
# We can now do the full fit. Similarly with the L case, we will fit (1) over full space
# and (2) only zonal mean as well.
predictors = (:Ω_nf, :ECTEI)
model_expression = "p[1]*x1 + p[2]*x2*(1-x1)"
eval_model_equations(model_expression, predictors)

# (1):
map_to_ocean_masked = (X) -> ocean_masked(X, timemean(O) .> ocean_mask_perc, MAXDEG)
Ps = map(p -> getindex(field_dictionary, p), predictors)
Ymasked = map_to_ocean_masked(α_CERES)
Psmasked = map_to_ocean_masked.(timemean.(Ps))

Mmasked, err, p, = perform_fit(model, Ymasked, Psmasked, ones(3); n = 3)
M1 = maskedzonalmean(timemean(model(p, Ps...)), OCEAN_MASK, MAXDEG)


# (2):
Ytz = maskedzonalmean(α_CERES, OCEAN_MASK, MAXDEG)
w = cosd.(gnv(dims(Ytz, Lat))) # weights
Pstz = maskedzonalmean.(timemean.(Ps), Ref(OCEAN_MASK), Ref(MAXDEG))
Mtz, = perform_fit(model, Ytz, Pstz, ones(3); n = 50, w)

# And finally, get ERA5
α_ERA5_tz = maskedzonalmean(α_ERA5, OCEAN_MASK, MAXDEG)

# and plot everything:
axC = mainfig.add_subplot(122)
plot_zonal_comparison!(Ytz, α_ERA5_tz, M1, Mtz; ax = axC)

# Final adjustments before saving
# %% #src
mainfig.tight_layout()
axL.set_title("longwave cloud radiative effect \$L\$", fontsize = 24)
axL.legend(fontsize = 18, loc = "upper left")
axC.set_title("cloud contribution to albedo \$\\alpha^\\mathrm{CLD}\$", fontsize = 24)
axC.legend(fontsize = 18)
add_identifiers!(mainfig)
mainfig.subplots_adjust(bottom = 0.12)
wsave(papersdir("plots", "versus_era5"), mainfig)