using DrWatson
@quickactivate "AlbedoBounds"
include(scriptsdir("predictors", "cloudiness_predictors_definition.jl"));

# %%
ocean_mask_frac = 50
OCEAN_MASK = O .> ocean_mask_frac # do we want ice or not in our masking?
MAXDEG = 90
latzones = (-90, -30, 0, 30, 90)

function maxdegsel(F, MAXDEG = MAXDEG)
    MAXDEGSEL = Lat((-MAXDEG)..(MAXDEG))
    if hasdim(F, Lat)
        return F[MAXDEGSEL]
    elseif hasdim(F, Coord)
        return F[Coord(MAXDEGSEL)]
    else
        return F
    end
end

function mean_over_latzones(_F, _OCEAN_MASK = OCEAN_MASK, latzones = latzones)
    F = maxdegsel(_F)
    OCEAN_MASK = maxdegsel(_OCEAN_MASK)
    means = zeros(length(latzones)-1)
    for j in 1:length(means)
        l1, l2 = latzones[j], latzones[j+1]
        Fsel = F[Coord(Lat(l1..l2))]
        Wsel = OCEAN_MASK[Coord(Lat(l1..l2))]
        means[j] = mean(vec(Fsel), weights(vec(Wsel)))
    end
    return means
end

fields_to_use = (C, L, F)
D = length(fields_to_use)
fig, axs = subplots(D, 2)

t = dims(C, Time).val
t = 1:length(t)

names = ["\$C\$", "\$L\$"]
for (i, F) in enumerate(fields_to_use)
    # Timeseries plot:
    means = zeros(length(latzones)-1)
    for j in 1:length(latzones)-1
        l1, l2 = latzones[j], latzones[j+1]
        Fsel = F[Coord(Lat(l1..l2))]
        Wsel = OCEAN_MASK[Coord(Lat(l1..l2))]
        out = spacemean(Fsel, Wsel)
        # axs[i, 2].plot(t[100:200], out[100:200]; lw = 1.0, label = "$j")
        dates, vals = seasonality(out)
        m = mean.(vals)
        m = m .- mean(m)
        push!(m, m[1])
        v = std.(vals)
        push!(v, v[1])
        axs[i, 2].plot(1:13, m; color = "C$(j-1)", lw = 2, label = "($(l1),$(l2))")
        axs[i, 2].fill_between(1:13, m.-v, m.+v; color = "C$(j-1)", alpha = 0.5)
        # axs[i, 2].plot(mean.(vals); color = "C$j", lw = 2)   
    end
    # Fn, Fs = hemispheric_means(maxdegsel(F), maxdegsel(OCEAN_MASK))
    # axs[i,2].plot(t[100:200], Fn[100:200], label = "\$\\mathcal{N}\$")
    # axs[i,2].plot(t[100:200], Fs[100:200]; linestyle = :dashed, label = "\$\\mathcal{S}\$")
    # Zonally average plot:
    Fz = maxdegsel(timemean(zonalmean(F, OCEAN_MASK)))
    lats = sind.(dims(Fz, Lat).val)
    axs[i,1].plot(lats, Fz; color = "C4", linewidth = 3)
    means = mean_over_latzones(F)
    for j in 1:length(latzones)-1
        l1, l2 = latzones[j], latzones[j+1]
        x = (sind(l1) + sind(l2))/2
        s = string(rdspl(means[j], 3))
        axs[i, 1].text(x, 0.5, s; ha = :center, transform = axs[i, 1].get_xaxis_transform())
        # axs[i, 1].axvline(sind(l1); color = "k")
        axs[i, 1].axvline(sind(l2); color = :k, lw = 1)
        axs[i, 1].axvspan(sind(l1), sind(l2); color = "C$(j-1)", alpha = 0.2)
    end
    axs[i,1].set_xlim(-1,1)
    axs[i,2].set_xlim(1,13)
end
for i in 1:D
    axs[i,1].set_xticks(sind.(-70:20:70))
    axs[i,2].set_xticks(1:2:13)
    axs[i,1].set_xticklabels(-70:20:70)
end

setp(axs[1,1].get_xticklabels(), visible=false)
# axs[1,2].legend(loc = "upper left", framealpha = 0.95)
axs[1,1].set_ylabel("\$C\$ (%)")
axs[2,1].set_ylabel("\$L\$ (W/m²)")
axs[3,1].set_ylabel("\$F\$ (%)")
axs[D,1].set_xlabel("lat (ᵒ)")
axs[D,2].set_xlabel("time (month)")
for i in 1:D-1
    setp(axs[i,1].get_xticklabels(), visible=false)
    setp(axs[i,2].get_xticklabels(), visible=false)
end
add_identifiers!(fig)
fig.tight_layout(pad=0.3)

wsave(papersdir("plots", basename(@__FILE__)[1:end-3]), fig)
