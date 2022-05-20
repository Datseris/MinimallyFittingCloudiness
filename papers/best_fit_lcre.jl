# Due to the hack done in meta-programmically creating the model from
# strings, one needs to run this script twice for it to work.
# (I don't know how to fix this bug yet)
using DrWatson
@quickactivate "MinimallyFittingCloudiness"
include(scriptsdir("fields_definitions.jl"));
include(srcdir("fitting", "masking.jl"))
include(srcdir("fitting", "general.jl"))

predictors = (:Ω_std, :Ω_mean, :Tsfc)
predicted = :L
Φname = "\$L\$"
expression = "p[1]*x1 + p[2]*x2 + p[3]*x3"
params = [42.67918994446545, 208.88564817976484, 0.06557681258808708]
MAXDEG = 70 # fit only within ± MAXDEG
ocean_mask_perc = 50 # points with % ≥ than this are considered "ocean"
limits = (0, 70)
seasonal_limits = (-8, 25)
seasonal_offset = (seasonal_limits[2] - seasonal_limits[1])/2

function generate_contributions(p, Ps)
    c1 = p[1]*Ps[1]
    c2 = p[2]*Ps[2]
    c3 = p[3]*Ps[3]
    return c1, c2, c3
end
contribution_titles = ("\$p_1 \\omega_\\mathrm{std}\$", "\$p_2 \\omega_{500}\$",  "\$p_3 \\mathrm{SST}\$",)
contrib_cmaps = (:inferno, :BrBG, :inferno, )

# The rest of this will be made on script
include(papersdir("best_fit_plot.jl"))

# %% Total longwave CRE
total_lwCRE = spacemean(timemean(Y), timemean(OCEAN_MASK))
total_lwCRE_model = spacemean(timemean(M), timemean(OCEAN_MASK))

total_lwCRE, total_lwCRE_model
