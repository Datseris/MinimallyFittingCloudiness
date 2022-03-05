# In this file we provide functions that create definitions of
# cloudiness and return them as spatiotemporal data

"""
    total_toa_albedo(a, s, t) = a + s*t^2/(1-a*s)
Combine given atmosphere albedo `a`, surface albedo `s` and atmosphere transmittance `t`
into a total top-of-the-atmosphere albedo `α` according to the model of Donohoe & Battisti (2011).
"""
total_toa_albedo(a, s, t) = a + s*t^2/(1-a*s)


"""
According to paper by Datseris & Stevens.
Also return the (time averaged)
"""
function effective_cloud_albedo(EBAF)

    R = ncread(EBAF, "toa_sw_all_mon")
    K = ncread(EBAF, "toa_sw_clr_t_mon")
    I = ncread(EBAF, "solar_mon")
    F = ncread(EBAF, "cldarea_total_daynight_mon") ./ 100
    T = ncread(EBAF, "cldtau_total_day_mon")
    Tcorrected = sinusoidal_continuation(T, [1, 2.0]; Tmin = 0)
    # C1 = _cloud_effective_albedo(F, Tcorrected, Float32(0.9))
    # C1_timeavg = timemean(C1)

    # use `\:arrow_down:` and tab for arrows
    F_s_⬆ = ncread(EBAF, "sfc_sw_up_all_mon")
    F_s_⬇ = ncread(EBAF, "sfc_sw_down_all_mon")
    F_s_⬆_K = ncread(EBAF, "sfc_sw_up_clr_t_mon")
    F_s_⬇_K = ncread(EBAF, "sfc_sw_down_clr_t_mon");

    argsall = timemean.((I, R, F_s_⬆, F_s_⬇))
    argsclr = timemean.((I, K, F_s_⬆_K, F_s_⬇_K))
    α_ATM, α_SFC = surface_atmosphere_contributions(argsall...) # all sky
    α_ATM_K, α_K_SFC = surface_atmosphere_contributions(argsclr...) # clear sky

    # This is the normalization done to have energetically consistent albedo value
    # according to Donohoe&Battisti model
    C2_timeavg = α_ATM .- α_ATM_K
    a = C2_timeavg; τ = timemean(Tcorrected); f = timemean(F)
    _x = -2a ./ (τ .* a - f .* τ)
    g = clamp.(1 .- _x ./ Float32(√3), 0, 1)

    C = _cloud_effective_albedo(F, Tcorrected, g)
end

function donohoe_battisti_cloud_albedo(I, R, K, F_s_⬆, F_s_⬇, F_s_⬆_K, F_s_⬇_K)
    argsall = timemean.((I, R, F_s_⬆, F_s_⬇))
    argsclr = timemean.((I, K, F_s_⬆_K, F_s_⬇_K))
    α_ATM, α_SFC = surface_atmosphere_contributions(argsall...) # all sky
    α_ATM_K, α_K_SFC = surface_atmosphere_contributions(argsclr...) # clear sky
    return α_ATM .- α_ATM_K
end


"""
    surface_atmosphere_contributions(I, F_toa_⬆, F_s_⬆, F_s_⬇) → α_ATM, α_SFC
Calculate the atmospheric and surface _contributions_ of the planetary albedo, so that
the TOA albedo is `α = α_ATM + α_SFC`, using the
simple 1-layer radiative transfer model by Donohoe & Battisti (2011) or G. Stephens (2015).
Stephens' formulas are incorrect and I have corrected them!
"""
function surface_atmosphere_contributions(I, F_toa_⬆, F_s_⬆, F_s_⬇)
    R = F_toa_⬆ ./ I   # planetary albedo (system reflectance)
    T = F_s_⬇ ./ I     # system transmisttance
    α = F_s_⬆ ./ F_s_⬇ # surface albedo

    # Formulas by Graeme, which are wrong!
    # τ = @. (1 - α*R) / (1 - (α^2) * (T^2)) # atmosphere transmittance
    # r = @. R - τ*α*T # atmospheric contribution to albedo

    # My calculations:
    r = @. (R - α*T^2) / (1 - (α*T)^2) # atmospheric contribution to albedo
    t = @. T*(1 - r*α)                 # atmospheric transmittance
    s = @. (α*t^2) / (1 - r*α)         # surface contribution to albedo
    return r, s
end


# This actually implements the formula from Lacis
function _cloud_effective_albedo(f, τ, g = 0.9; frequencies = [1.0, 2.0])
    T = τ
    p = @. Float32(sqrt(3))*(1 - g)
    R = p .* T ./ (2 .+ p .* T)
    cloud_frac = any(x -> x > 1, f) ? f ./ 100 : f
    eca = cloud_frac .* R
    gstr = g isa Real ? string(g) : "array"
    return ClimArray(Float32.(eca.data), eca.dims, "effective cloud albedo, g = $(gstr)")
end

function cloud_albedo_lacis(f, τ, g)
    p = eltype(g)(sqrt(3))*(1 .- g)
    r = p .* τ ./ (2 .+ p .* τ)
    return f .* r
end
