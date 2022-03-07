include("recipes.jl")
include("../statistics.jl")

function plot_zonal_averages(Xs, latzones = (-90, -30, 0, 30, 90))
    labels = [String(X.name) for X in Xs]
    # Zonal averages
    fig, axs = subplots(length(Xs), 2; figsize = (20, 18))
    for (i, X) in enumerate(Xs)
        Xz = to_zonalmean(X)

        Xzsouth = reverse(Xz[Lat((-90)..(0))])
        Xznorth = Xz[Lat(Between(0, 90))]
        dl = length(Xznorth) - length(Xzsouth)
        Xznorth = Xznorth[1:end-dl]
        Diff_Xz = Xznorth .- Xzsouth

        latd = sind.(dims(Xz, Lat).val)
        diff_latd = sind.(dims(Diff_Xz, Lat).val)
        axs[i, 1].plot(latd, Array(Xz); color = "C$(i-1)")
        axs[i, 1].set_ylabel(labels[i]; color = "C$(i-1)")
        mean_val = latmean(Xz)
        axs[i, 1].text(0.5, 0.8, string(rdspl(mean_val)), ha = :center, va = :center, transform = axs[i,1].transAxes)
        axs[i, 2].plot(diff_latd, Array(Diff_Xz); color = "C$(i-1)")
        mean_diff = latmean(Diff_Xz)
        axs[i, 2].text(0.5, 0.8, string(rdspl(mean_diff/mean_val)), ha = :center, va = :center, transform = axs[i,2].transAxes)

        for j in 1:length(latzones)-1
            l1, l2 = latzones[j], latzones[j+1]
            mv = latmean(Xz[Lat(Between(l1, l2))])
            x = (sind(l1) + sind(l2))/2
            s = string(rdspl(mv, 3))
            axs[i, 1].text(x, 0.1, s; ha = :center, transform = axs[i, 1].get_xaxis_transform())
            axs[i, 1].axvspan(sind(l1), sind(l2); alpha = 0.05, color = "C$(j+1)")
            if i == length(Xs)
                axs[i,1].set_xticks(sind.(-70:20:70))
                axs[i,1].set_xticklabels(-70:20:70)
                axs[i,2].set_xticks(sind.(0:10:80))
                axs[i,2].set_xticklabels(0:10:80)
            else
                axs[i,1].set_xticks(sind.(-70:20:70))
                axs[i,2].set_xticks(sind.(0:10:80))
                axs[i,1].set_xticklabels([])
                axs[i,2].set_xticklabels([])
            end
        end
    end
    axs[1,1].set_title("Zonal average")
    axs[1,2].set_title("Hemispheric difference")
    fig.tight_layout(pad = 0.5)
    fig.subplots_adjust(hspace = 0.1)
end

# %% Linear correlation map:
function plot_correlation_maps(X, Y)
    M = copy(timemean(Y))
    for i in spatialidxs(Y)
        M[i] = Statistics.cor(X[i...], Y[i...])
    end
    aaa = spacemean(dropnan(M))
    fig, ax, cb = earthscatter(M; s = 30, vmin = -1, vmax = 1, cmap = "PRGn")
    ax.set_title("corr. ($(X.name), $(Y.name)), avg. corr = $(rdspl(aaa))")
    return aaa
end

function plot_spacetime_comparison(X, M)
    fig, axs = subplots(3, 1, figsize = (12,10), sharex = true)
    xcoord = dims(X, Lat).val
    ycoord = 1:size(X, Time)
    vmin, vmax = extrema(X)
    if isnan(vmin) || isnan(vmax)
        z = vec(X)
        j = findall(!isnan, z)
        vmin, vmax = extrema(z[j])
    end
    im = axs[1].pcolormesh(ycoord, xcoord, X.data; vmin, vmax)
    fig.colorbar(im, ax=axs[1])
    im = axs[2].pcolormesh(ycoord, xcoord, M.data; vmin, vmax)
    fig.colorbar(im, ax=axs[2])
    D = X.data .- M.data
    v = maximum(abs.(D))
    isnan(v) && (v = vmax/3)
    im = axs[3].pcolormesh(ycoord, xcoord, D; cmap = "PRGn", vmax = v, vmin = -v)
    fig.colorbar(im, ax=axs[3], extend="both")
    axs[1].set_ylabel(string(C.name))
    axs[2].set_ylabel("model")
    axs[3].set_ylabel("diff")
    for ax in axs; ax.set_aspect("auto"); end
    axs[1].set_title("NRMSE = $(round(nrmse(X, M);sigdigits=4))")
    return fig, axs
end

function plot_spacedependent_comparison(X, Y, OCEAN_MASK, pfit)
    titlestr = join(rdspl.(pfit, 2), ", ")
    C = copy(X)
    vmin, vmax = extrema(C)
    C[.!(OCEAN_MASK)] .= NaN
    M = copy(X)
    M[.!(OCEAN_MASK)] .= NaN
    M[OCEAN_MASK] .= Y
    vmin, vmax = extrema(timemean(X))
    fig, ax = earthscatter(timemean(C); vmin, vmax)
    ax.set_title("real data")
    fig, ax = earthscatter(timemean(M); vmin, vmax)
    ax.set_title("model fit")
    err = begin
        x, y = timemean(C), timemean(M)
        is = findall(!isnan, x)
        nrmse(x[is], y[is])
    end
    println("Error of spatial map: ", err)
    return err
end

using PyCall, Statistics
seaborn = pyimport("seaborn")

function plot_scatter_pdf(ω, c; color = "C0", xlabel = "x", ylabel="y")
    # axs[1, i].hist(Array(c); bins = 20)
    # axs[2, i].hist(Array(ω); bins = 20)
    ret = seaborn.jointplot(; y = Array(c), x = Array(ω),
        height=10, ratio=3, kind = "hist",
        # color, alpha = 0.1,
        cmap = "viridis"
    )
    scax = ret.ax_joint
    scax.set_xlabel(xlabel)
    scax.set_ylabel(ylabel)

    slope, _, fit = linreg(ω,c)
    scax.plot(ω, fit; color="red")
    scax.text(0.95, 0.95, "s = $(rdspl(slope))", transform = scax.transAxes, ha = :right, color = "red")

    ωax = ret.ax_marg_x
    ωax.axvline(mean(ω); ls = "dashed", color="red")
    cax = ret.ax_marg_y
    cax.axhline(mean(c); ls = "dashed", color="red")
    ret.fig.tight_layout()
    return ret.fig, (scax, ωax, cax)
end

function plot_spatial_map_with_mask_lat(C, OCEAN_MASK, MAXDEG = 91;
    vmin = nothing, vmax = nothing, cmap = :viridis,
    levels = 41)
    if MAXDEG < 90
        C = C[Coord(Lat(Between(-MAXDEG, MAXDEG)))]
        OCEAN_MASK = OCEAN_MASK[Coord(Lat(Between(-MAXDEG, MAXDEG)))]
    end
    tmXX = timemean(C, OCEAN_MASK)
    if isnothing(vmin)
        vmin, vmax = extrema(tmXX)
    end
    if isnan(vmin) || isnan(vmax)
        z = vec(tmXX)
        j = findall(!isnan, z)
        vmin, vmax = extrema(z[j])
    end
    fig = figure()
    axC = subplot(111, projection=DEFPROJ)
    earthscatter!(axC, tmXX; vmin, vmax, levels, cmap)
    return fig, axC
end

function plot_total_model_comparison(C, Xs, model, p, OCEAN_MASK; MAXDEG = 91)
    name = string(C.name)

    println("Creating full model field...")
    M = ClimArray(model(p, Xs...); name = "Model")

    if MAXDEG < 90
        C = C[Coord(Lat(Between(-MAXDEG, MAXDEG)))]
        M = M[Coord(Lat(Between(-MAXDEG, MAXDEG)))]
        OCEAN_MASK = OCEAN_MASK[Coord(Lat(Between(-MAXDEG, MAXDEG)))]
        Xs = map(X -> X[Coord(Lat(Between(-MAXDEG, MAXDEG)))], Xs)
    end
    full_nrmse = nrmse(C[OCEAN_MASK], M[OCEAN_MASK])
    println("Full NRMSE: ", full_nrmse)


    # Compare timemeans (spatial maps)
    fig = figure()
    axC = subplot(121, projection=DEFPROJ)
    axM = subplot(122, projection=DEFPROJ)
    Cmap = timemean(C, OCEAN_MASK)
    Mmap = timemean(M, OCEAN_MASK)
    _, cb = earthscatter!(axC, Cmap; levels = 21)
    earthscatter!(axM, Mmap; vmin = cb.vmin, vmax = cb.vmax, levels = 21)
    axC.set_title(name)
    axM.set_title("M")
    fig.suptitle("NRMSE=$(round(full_nrmse;sigdigits=3)),\nMAXDEG=$(MAXDEG), OCEAN_FRAC=$(ocean_mask_perc)")
    timemean_error = nrmse(Cmap, Mmap)
    println("Time-averaged NRMSE: ", timemean_error)

    # plot difference map
    fig = figure()
    axD = subplot(111, projection=DEFPROJ)
    dmap = Cmap .- Mmap
    vdif = cb.vmax/3
    earthscatter!(axD, dmap; vmin=-vdif, vmax=vdif, levels = 21, cmap = "PuOr")
    axD.set_title("$name - M")

    # Time variability stuff:

    # spatial correlation plot
    CC = copy(C); MM = copy(M)
    CC[.!OCEAN_MASK] .= NaN
    MM[.!OCEAN_MASK] .= NaN
    mean_correlation = plot_correlation_maps(CC, MM)
    println("Average pearson coefficien: ", mean_correlation)

    # latitude zones timeseries (new)
    timeseries_errors = Float64[]
    latzones = (-90, -30, 0, 30, 90)
    fig, axs = subplots(4,1; sharex = true)
    for j in 1:length(latzones)-1
        l1, l2 = latzones[j], latzones[j+1]
        Csel = C[Coord(Lat(l1..l2))]
        Wsel = OCEAN_MASK[Coord(Lat(l1..l2))]
        Msel = M[Coord(Lat(l1..l2))]
        Csel = spacemean(Csel, Wsel)
        Msel = spacemean(Msel, Wsel)

        demeaned = []

        # axs[i, 2].plot(t[100:200], out[100:200]; lw = 1.0, label = "$j")
        for (n, out) in enumerate((Csel, Msel))
            dates, vals = seasonality(out)
            m = mean.(vals)
            m = m .- mean(m)
            push!(demeaned, m)
            push!(m, m[1])
            v = std.(vals)
            push!(v, v[1])
            axs[j].plot(1:13, m; color = "C$(n-1)", lw = 2, label = n == 1 ? "C" : "M",
            ls = n == 1 ? "-" : "--")
            axs[j].fill_between(1:13, m.-v, m.+v; color = "C$(n-1)", alpha = 0.5)
        end
        # push!(timeseries_errors, nrmse(Csel, Msel))
        push!(timeseries_errors, nrmse(demeaned[1], demeaned[2]))
        axs[j].set_ylabel("($(l1),$(l2))")
    end

    println("Timeseries error: ", median(timeseries_errors))


    # # hemispheric timeseries (old)
    # Xnh, Xsh = hemispheric_means(C, OCEAN_MASK)
    # Mnh, Msh = hemispheric_means(M, OCEAN_MASK)
    # t = dims(Xnh, Time).val
    # fig, axs = subplots(2,1; sharex = true)
    # axs[1].plot(t, Xnh; label = "\$\\mathcal{T}($name)\\approx $(rdspl(timemean(Xnh)))\$")
    # axs[1].plot(t, Mnh; label = "\$\\mathcal{T}(M)\\approx $(rdspl(timemean(Mnh)))\$", ls = "--")
    # axs[1].set_ylabel("NH")
    # axs[2].plot(t, Xsh; label = "\$\\mathcal{T}($name)\\approx $(rdspl(timemean(Xsh)))\$")
    # axs[2].plot(t, Msh; label = "\$\\mathcal{T}(M)\\approx $(rdspl(timemean(Msh)))\$", ls = "--")
    # axs[2].set_ylabel("SH")
    # axs[1].set_xlim(Date(2004), Date(2010))
    # for ax in axs; ax.legend(); end
    # nh_error = nrmse(Xnh, Mnh)
    # sh_error = nrmse(Xsh, Msh)
    # println("NH timeseries NRMSE: ", nh_error)
    # println("SH timeseries NRMSE: ", sh_error)

    # Zonal mean spacetime plots
    tozonalmean(X) = zonalmean(X, OCEAN_MASK)
    zC, zM = tozonalmean(C), tozonalmean(M)
    println("Zonally-avg. NRMSE: ", nrmse(zC, zM))
    plot_spacetime_comparison(zC, zM)
    timezonal_full_error = nrmse(timemean(zC), timemean(zM))
    println("Zonal+time-avg NRMSE: ", timezonal_full_error)

    oceany_zonal(X) = to_non_nan_lats(timemean(zonalmean(X, OCEAN_MASK)))
    zP = [oceany_zonal(P) for P in Xs]
    zX = oceany_zonal(C)
    fig, axs = plot_step_by_step_comparison(zX, zP, model, p)
    axs[1].set_title("fit using full data\nMAXDEG=$(MAXDEG), OCEAN_FRAC=$(ocean_mask_perc)")
    fig.tight_layout(pad=0.3)
    return timemean_error, mean_correlation, timeseries_errors, timezonal_full_error
end

function plot_step_by_step_comparison(X, Ps, model, p; latzones = (-90, -30, 0, 30, 90))
    M = model(p, Ps...)
    R = X .- M # residuals of first fit
    fig, axs = subplots(2, 1, figsize = (12,8), sharex = true)

    # First plot full model comparison with zonal means
    ax = axs[1]
    err = nrmse(X.data, M.data)
    ax.plot(sind.(dims(X, Lat).val), X.data; label = X.name)
    ax.plot(sind.(dims(M, Lat).val), M.data; label = "M, nrmse = $(round(err;sigdigits=3))", ls = "-.")

    for i in 1:length(latzones)-1
        l1, l2 = latzones[i], latzones[i+1]
        mv = latmean(M[Lat(Between(l1, l2))])
        x = (sind(l1) + sind(l2))/2
        s = string(rdspl(mv, 3))
        ax.text(x, 0.05, s, ha = :center, transform = ax.get_xaxis_transform())
        mv = latmean(X[Lat(Between(l1, l2))])
        x = (sind(l1) + sind(l2))/2
        s = string(rdspl(mv, 3))
        ax.text(x, 0.15, s, ha = :center, transform = ax.get_xaxis_transform(); color = "C0")
        ax.axvspan(sind(l1), sind(l2); alpha = 0.05, color = "C$(i+1)")
    end
    ax.legend()

    # Now plot residuals of first fit and their fit
    err2 = nrmse(R.data, M.data)
    axs[2].plot(sind.(dims(R, Lat).val), R; label = "residuals", color = "C2")
    axs[2].legend()
    v = latmean(R)
    axs[2].text(0.5, 0.15, string(rdspl(v)), ha = :center, transform = axs[2].transAxes; color = "C2")
    latitudinal_axis!(axs[2])
    axs[2].set_xlabel("sin(φ)")
    fig.tight_layout(pad=0.3)
    return fig, axs
end


function full_plot_field(C, tit = string(C.name);
        usescatter = spacestructure(C) == CoordinateSpace(), cmap = "YlGnBu_r",
        save = false
    )
    ff(x) = rdspl(x, 3)
    fig = figure(figsize = (13, 16))
    ax1 = fig.add_subplot(2,1,1; projection = DEFPROJ)

    Ct = timemean(C)
    vmin, vmax = Statistics.quantile(vec(Ct.data), [0.025, 0.975])
    if usescatter
        earthscatter!(ax1, Ct; vmin, vmax, levels = 21, cmap)
    else
        earthsurface!(ax1, Ct; vmin, vmax, levels = 21, cmap)
    end
    ax1.set_title(tit)

    ax2 = fig.add_subplot(4,1,4) # timeseries
    Cn, Cs = hemispheric_means(C)
    t = dims(C, Time).val

    # Cn, Cs = hemispheric_means(C)
    # ax.plot(t, Cn, label = "\$\\mathcal{N}(C)\$")
    # ax.plot(t, Cs; ls = "dashed", label = "\$\\mathcal{S}(C)\$")

    ax2.plot(t, Cn; label = "\$\\mathcal{N}\$")
    ax2.plot(t, Cs; ls = "dashed", label = "\$\\mathcal{S}\$")
    Cg = spacemean(C)
    ax2.plot(t, Cg; ls = "-.", label = "\$\\mathcal{G}\$")
    ax2.legend(ncol = 3)
    # ax2.set_xlim(t[50], t[150])
    ax2.set_ylabel("timeseries")

    ax3 = fig.add_subplot(4,1,3) # zonal plot
    Cl = zonalmean(timemean(C))
    nice_lat_plot!(ax3, Cl)
    ax3.set_ylabel("zonal mean")

    # add stats statements
    g = spacemean(Ct) # total mean
    nh, sh = hemispheric_means(Ct)
    gstd, nstd, sstd = timeagg.(std, (Cg, Cn, Cs))


    ax3.set_title(
    "Stats: \$\\mathcal{GT} = $(ff(g)), \\mathcal{NT} = $(ff(nh)), \\mathcal{ST} = $(ff(sh)), \$\n\$"*
    "\\mathcal{G}_\\mathrm{std} = $(ff(gstd)), \\mathcal{N}_\\mathrm{std} = $(ff(nstd)), \\mathcal{S}_\\mathrm{std} = $(ff(sstd)), "*
    "\\mathcal{T}_\\mathrm{std} = $(ff(std(Ct)))\$",
    size = 20)

    save && wsave(plotsdir("info", tit*".png"), fig)
    return fig
end

function nice_lat_plot!(ax, Cl, latzones = (-90, -30, 0, 30, 90))
    ax.plot(sind.(dims(Cl, Lat).val), Cl.data; color = "C1")

    for i in 1:length(latzones)-1
        l1, l2 = latzones[i], latzones[i+1]
        mv = latmean(Cl[Lat(Between(l1, l2))])
        x = (sind(l1) + sind(l2))/2
        s = string(rdspl(mv, 3))
        ax.text(x, 0.15, s, ha = :center, transform = ax.get_xaxis_transform(); color = "C1")
        ax.axvspan(sind(l1), sind(l2); alpha = 0.05, color = "C$(i+1)")
    end

    latitudinal_axis!(ax)
end
