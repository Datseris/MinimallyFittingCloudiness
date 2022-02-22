using DrWatson
@quickactivate "MinimallyFittingCloudiness"

using ClimateBase
include(srcdir("cloudiness.jl"))
include(srcdir("plotting", "recipes.jl"))
include(srcdir("plotting", "advanced.jl"))
include(srcdir("statistics.jl"))
include(srcdir("thermodynamics.jl"))
Time = ClimateBase.Ti
GEA = 250 # Gaussian Equal Area grid spacing. 250 or 125.

EBAF = datadir("exp_pro", "CERES_EBAF_gea$(GEA).nc")
SYN1deg = datadir("exp_pro", "SYN1deg_gea$(GEA).nc")
ERA5_Ω = datadir("exp_pro", "ERA5_W_statistics_gea$(GEA).nc")
ERA5_3D = datadir("exp_pro", "ERA5_3D_gea$(GEA).nc")
ERA5_2D = datadir("exp_pro", "ERA5_2D_gea$(GEA).nc")

# CERES data provide Cloud albedo and solar zenith angle
C = effective_cloud_albedo(EBAF)
C = ClimArray(100C; name = "C") # Make albedo percentage units for more numeric similarity

# and also cloud radiative effects
Rall = ncread(EBAF, "toa_sw_all_mon")
Rclr = ncread(EBAF, "toa_sw_clr_t_mon")
Lall = ncread(EBAF, "toa_lw_all_mon")
Lclr = ncread(EBAF, "toa_lw_clr_t_mon")
F = ncread(EBAF, "cldarea_total_daynight_mon"; name = "F")
CREsw = ClimArray(Rall - Rclr; name = "CREsw")
CRElw = ClimArray(Lclr - Lall; name = "CRElw")
CRE = ClimArray(CRElw - CREsw; name = "CRE")
L = ClimArray(CRElw; name = "L")

# C[Coord(5668)] .= C[Coord(5669)] # Don't know why, don't ask why :D
I = ncread(EBAF, "solar_mon")
Z = ClimArray(acos.(I ./ 1361); name = "zenith")
O = ncread(SYN1deg, "aux_ocean_mon"; name = "ocean_fraction")
ICE = ncread(SYN1deg, "aux_snow_mon"; name = "ice_fraction")
LAND = ClimArray(100 .- O .- ICE; name = "land_fraction")

# Predictors based on wind speed
Ω_mean = ClimArray(-ncread(ERA5_Ω, "W_mean"); name = "Ω_mean") # notice multiplication with -1
Ω_std = ncread(ERA5_Ω, "W_std"; name = "Ω_std")
Ω_nf = ncread(ERA5_Ω, "W_nf"; name = "Ω_nf")
WS10 = ncread(ERA5_2D, "si10"; name = "wind10")


# Remaining predictors are thermodynamics-based
# RH = ncread(ERA5_3D, "r"; name = "relative_humidity")
T = ncread(ERA5_3D, "t"; name = "temperature")
T700 = ClimArray(T[Pre(At(700))]; name = "T700hPa")
Tsfc =  ncread(ERA5_2D, "t2m"; name= "Tsfc")
Psfc =  ncread(ERA5_2D, "sp")
Psfc = ClimArray(Psfc ./ 100; name = "p_sfc") # surface pressure is in Pascal, but we need in hPa.
LTS = ClimArray(lower_tropospheric_stability(Tsfc, T700, Psfc); name = "LTS")
V = ncread(ERA5_2D, "tcwv")

EIS = estimated_inversion_strength(Tsfc, T700; Psfc)
EIS = ClimArray(EIS; name = "EIS")

q = ncread(ERA5_3D, "q"; name= "q")
q = 1000q # multiply with 1000 to make it in units of g/kg instead of kg/kg
q.attrib["units"] = "g/kg"
q700 = ClimArray(q[Pre(At(700))]; name = "q_700hPa")
qsfc = ClimArray(q[Pre(At(1000))]; name = "q_sfc") # <- Can be improved!!!
ECTEI = estimated_cloud_top_entrainment_index(Tsfc, T700, qsfc/1000, q700/1000; Psfc)
ECTEI = ClimArray(ECTEI; name = "ECTEI")

# And we also load low-cloud cover, as this correlates with EIS/ECTEI
LCC = ncread(ERA5_2D, "lcc"; name = "cloud_cover_low")
TCC = ncread(ERA5_2D, "tcc"; name = "cloud_cover_total")
HCC = ClimArray(TCC .- LCC; name = "cloud_cover_high")


# and lastly create teh uniform predictor "U" which is just a predictor that
# has the values 1 everywhere. This is only useful for Linearized Models fit where
# it fits the intercept
U = ClimArray(ones(CRE); name = "intercept")

# Bring every field to the same time span
field_dictionary = @dict(
    CREsw, L, CRE, C, I, Z, O, Ω_mean, Ω_std, Ω_nf, WS10,
    q700, qsfc, T700, Tsfc, V, EIS, ECTEI, U, L, F, LTS
)

filed_dictionary = sametimespan(field_dictionary)
