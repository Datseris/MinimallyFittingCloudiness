using DrWatson
@quickactivate "MinimalyFittingCloudiness"

using ClimateBase
include(srcdir("cloudiness.jl"))
include(srcdir("plotting", "recipes.jl"))
include(srcdir("plotting", "advanced.jl"))
include(srcdir("statistics.jl"))
include(srcdir("thermodynamics.jl"))
Time = ClimateBase.Ti
GEA = 250 # Gaussian Equal Area grid spacing. 250 or 125.

EBAF_TOA = datadir("CERES", "EBAF_TOA_gea$(GEA).nc")
EBAF_SFC = datadir("CERES", "EBAF_SFC_gea$(GEA).nc")
SYN1deg = datadir("CERES", "SYN1deg_gea$(GEA).nc")
ERA5_Ω = datadir("ERA5", "ERA5_Wmoments_gea$(GEA).nc")
ERA5_3D = datadir("ERA5", "ERA5_3D_gea$(GEA).nc")
ERA5_2D = datadir("ERA5", "ERA5_2D_gea$(GEA).nc")

# CERES data provide Cloud albedo and solar zenith angle
C = effective_cloud_albedo(;EBAF_TOA, EBAF_SFC)
C = ClimArray(100C; name = "C") # Make albedo percentage units for more numeric similarity

# and also cloud radiative effects
Rall = ncread(EBAF_TOA, "toa_sw_all_mon")
Rclr = ncread(EBAF_TOA, "toa_sw_clr_c_mon")
Lall = ncread(EBAF_TOA, "toa_lw_all_mon")
Lclr = ncread(EBAF_TOA, "toa_lw_clr_c_mon")
F = ncread(EBAF_TOA, "cldarea_total_daynight_mon"; name = "F")

CREsw = ClimArray(Rall - Rclr; name = "CREsw")
CRElw = ClimArray(Lclr - Lall; name = "CRElw")
CRE = ClimArray(CREsw + CRElw; name = "CRE")

# C[Coord(5668)] .= C[Coord(5669)] # Don't know why, don't ask why :D
I = ncread(EBAF_TOA, "solar_mon")
Z = ClimArray(acos.(I ./ 1361); name = "zenith")
O = ncread(SYN1deg, "aux_ocean_mon"; name = "ocean_fraction")
ICE = ncread(SYN1deg, "aux_snow_mon"; name = "ice_fraction")
LAND = ClimArray(100 .- O .- ICE; name = "land_fraction")

# Predictors based on wind speed
Ω_mean = ClimArray(-ncread(ERA5_Ω, "Ω_mean"); name = "Ω_mean") # notice multiplication with -1
Ω_std = ncread(ERA5_Ω, "Ω_std")
Ω_nf = ncread(ERA5_Ω, "Ω_negative_fraction"; name = "Ω_neg_frac")
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

# Now the problem is, Z is in **complete antiphase with** C...
# However, our basic assumption is that C in extratropics is larger because
# of larger zenith angle.
# But if we use a timevarying zenith angle we have a problem. Because on one hand
# We expect C ~ Z correlation, but in the _time domain_ we have C ~ -Z anticorrelation.

# I'll do a dirty hack at the moment: I will replace temporal dependence of zenith
# angle its time-averaged version. This way we keep the latitudinal dependence
# but we remove the zenith dependence.

# ztimemeaned = timemean(ζ)
# for i in 1:length(dims(ζ, Time))
#     ζ[Time(i)] .= ztimemeaned
# end

# Aaaanyways, bring every field to the same time span

CREsw, CRElw, CRE, C, I, Z, O, Ω_mean, Ω_std, Ω_nf, WS10, q700, qsfc, T700, Tsfc, V, EIS, ECTEI, F, LTS = 
sametimespan(CREsw, CRElw, CRE, C, I, Z, O, Ω_mean, Ω_std, Ω_nf, WS10, q700, qsfc, T700, Tsfc, V, EIS, ECTEI, F, LTS)

# and lastly create teh uniform predictor "U" which is just a predictor that
# has the values 1 everywhere. This is only useful for Linearized Models fit where
# it gives offset

U = ClimArray(ones(CRE); name = "intercept")
L = CRElw

field_dictionary = @dict(CREsw, CRElw, CRE, C, I, Z, O, Ω_mean, Ω_std, Ω_nf, WS10, q700, qsfc, T700, Tsfc, V, EIS, ECTEI, U, L, F, LTS)