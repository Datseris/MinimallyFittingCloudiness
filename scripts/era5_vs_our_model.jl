# This script does something rather simple:
# compares the NRMSE of our model versus that of ERA5
# For our model we have the error already in `general_model_fit.jl`

using DrWatson
@quickactivate "MinimallyFittingCloudiness"
include(scriptsdir("fields_definitions.jl"));
include(srcdir("fitting", "masking.jl"))
include(srcdir("fitting", "general.jl"))
MAXDEG = 80 # fit only within ± MAXDEG
ocean_mask_perc = 50 # points with % ≥ than this are considered "ocean"
OCEAN_MASK = O .≥ ocean_mask_perc

function plot_zonal_comparison!(CERES_tz, ERA5_tz, MODEL_tz, axZ = gca())
    latticks = [-70, -30, -0, 30, 70]
    Φlats = sind.(gnv(dims(CERES_tz, Lat)))
    axZ = (figure(); gca())
    axZ.plot(Φlats, gnv(CERES_tz); label = "CERES", lw = 3)
    axZ.plot(Φlats, gnv(ERA5_tz);  label = "ERA5", ls = "--", alpha = 0.7, lw = 2.5)
    axZ.plot(Φlats, gnv(MODEL_tz); label = "MODEL", ls = "-.", alpha = 0.7, lw = 2.5)
    axZ.set_xticks(sind.(latticks))
    axZ.set_xticklabels(string.(latticks) .* "ᵒ")
    axZ.legend(fontsize = 18)
end

# %% #src

# # Comparison of L
L_CERES_tz = maskedtimezonalmean(L, OCEAN_MASK, MAXDEG)
# Compute L from ERA5:
Lall_ERA5 = ncread(ERA5_2D, "ttr")/86400 # transfrom to W/m²; divide by seconds in day
Lclr_ERA5 = ncread(ERA5_2D, "ttrc")/86400
L_ERA5 =  Lall_ERA5 .- Lclr_ERA5
L_ERA5_tz = maskedtimezonalmean(L_ERA5, OCEAN_MASK, MAXDEG)

# %%
# Compute our model:
predictors = (:Ω_std, :Ω_mean, :q700)
model_expression = "p[1]*x1 + p[2]*x2 +p[3]*x3"
eval_model_equations(model_expression, predictors)
Ps = map(p -> getindex(field_dictionary, p), predictors)
Pstz = map(P -> maskedtimezonalmean(P, OCEAN_MASK, MAXDEG), Ps)
L_MODEL_tz, errtz, ptz, = perform_fit(model, L_CERES_tz, Pstz, rand(3); n = 100)
@show nrmse(L_CERES_tz, L_ERA5_tz)
@show nrmse(L_CERES_tz, L_MODEL_tz)

plot_zonal_comparison!(L_CERES_tz, L_ERA5_tz, L_MODEL_tz)

# # Comparison of C
# %% src
# We can't compute cloud albedo directly from ERA5, because our definition
# requires us to use cloud optical depth for computing full spatiotemporal
# cloud albedo field. So, instead, we will use the Donohoe-Battisti version.
# (We will re-fit our data to the zonal+temporal mean values of course)
CERES_FIELDS = let
    I = ncread(EBAF, "solar_mon")
    R = ncread(EBAF, "toa_sw_all_mon")
    K = ncread(EBAF, "toa_sw_clr_t_mon")
    F_s_⬆ = ncread(EBAF, "sfc_sw_up_all_mon")
    F_s_⬇ = ncread(EBAF, "sfc_sw_down_all_mon")
    F_s_⬆_K = ncread(EBAF, "sfc_sw_up_clr_t_mon")
    F_s_⬇_K = ncread(EBAF, "sfc_sw_down_clr_t_mon");
    return I, R, K, F_s_⬆, F_s_⬇, F_s_⬆_K, F_s_⬇_K
end

ERA5_FIELDS = let
    I = ?
    R = ?
    K = ?
    F_s_⬇ = ncread(ERA5_2D, "ssrd")
    F_s_⬆ = ncread(ERA5_2D, "ssr") .- F_s_⬇
    F_s_⬆_K = ncread(EBAF, "sfc_sw_up_clr_t_mon")
    F_s_⬇_K = ncread(ERA5_2D, "ssrdc");
    return I, R, K, F_s_⬆, F_s_⬇, F_s_⬆_K, F_s_⬇_K
end

α_CERES = donohoe_battisti_cloud_albedo(CERES_FIELDS...)
α_CERES = zonalmean(α_CERES, timemean(OCEAN_MASK))
α_ERA5 = donohoe_battisti_cloud_albedo(ERA5_FIELDS...)
α_ERA5 = zonalmean(α_ERA5, timemean(OCEAN_MASK))

predictors = (:Ω_nf, :ECTEI)
model_expression = "p[1]*x1 + p[2]*x2*(1-x1)"
eval_model_equations(model_expression, predictors)
Ps = map(p -> getindex(field_dictionary, p), predictors)
Pstz = map(P -> maskedtimezonalmean(P, OCEAN_MASK, MAXDEG), Ps)
α_MODEL, errtz, ptz, = perform_fit(model, α_CERES, Pstz, rand(3); n = 100)
@show nrmse(α_CERES, α_ERA5)
@show nrmse(α_CERES, α_MODEL)

plot_zonal_comparison!(α_CERES, α_ERA5, α_MODEL)



# %% src
