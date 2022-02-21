import LsqFit

# Function start:
function hemispheric_model_fit(predicting_model, C, Xs, p0; latzones = (-90, -30, 0, 30, 90))
    L = length(latzones)-1
    M = ClimArray(copy(C); name = "Model fit", attrib = [copy(p0) for j in 1:L÷2])
    M .= 0
    Xs_split = [split_to_latzones(p, latzones) for p in Xs]
    Csplit = split_to_latzones(C, latzones)
    ## this is information necessary if the fields have spatial distribution,
    ## and not have been zonally averaged
    if hasdim(Csplit[1], Coord)
        points_per_latzone = [size(Csplit[i], Coord) for i in 1:L-1]
        pushfirst!(points_per_latzone, 0)
        additional_points = cumsum(points_per_latzone)
    else
        additional_points = Int[]
    end

    # Main fit loop over the latitude zones
    for i in 1:L÷2
        ## fitting field north and south hemisphere datapoints
        cs, cn = Csplit[i], Csplit[L-i+1]

        ## predictors north and south hemisphere datapoints
        ps = [p[i] for p in Xs_split]
        pn = [p[L-i+1] for p in Xs_split]

        ## training data: combine both hemispheres in the appropriate latitude zone
        ydata = vcat(vec(Array(cs)), vec(Array(cn)))
        xdata = hcat([vcat(vec(fs), vec(fn)) for (fs, fn) in zip(ps, pn)]...)

        ## Perform the model fit
        modelfit = LsqFit.curve_fit(predicting_model, xdata, ydata, p0)
        pfit = modelfit.param
        M.attrib[i] .= pfit

        ## Apply the model fit, now at each individual latzone uniquely
        apply_model_fit_per_latzone!(M, predicting_model, i, cs, cn, ps, pn, pfit, additional_points)
    end
    return M, nrmse(C.data, M.data)
end

# Helper functions:
function apply_model_fit_per_latzone!(M, predicting_model, i, cs, cn, ps, pn, pfit, additional_points)
    if hasdim(M, Coord) # version without zonal average
        coords = dims(M, Coord).val
        for j in 1:size(cs, Coord)
            rj = j + additional_points[i]
            coord_c = dims(cs, Coord)[j]
            coord_m = dims(M, Coord)[rj]
            @assert coord_c == coord_m # this must give true always
            ps_vals = [p[Coord(j)] for p in ps]
            if length(dims(M)) > 1
                M[Coord(rj)] .= predicting_model(hcat(ps_vals...), pfit)
            else
                M[rj] = predicting_model(hcat(ps_vals...), pfit)
            end
        end
        for j in 1:size(cs, Coord)
            rj = j + additional_points[L-i+1]
            coord_c = dims(cn, Coord)[j]
            coord_m = dims(M, Coord)[rj]
            @assert coord_c == coord_m # this must give true always
            pn_vals = [p[Coord(j)] for p in pn]
            M[rj] = multimodel(pfit, pn_vals...)
            if length(dims(M)) > 1
                M[Coord(rj)] .= predicting_model(hcat(pn_vals...), pfit)
            else
                M[rj] = predicting_model(hcat(pn_vals...), pfit)
            end
        end

        # TODO: I think we have a problem here: the latitudes are always
        # ordered from decreasing to increasing,
    elseif hasdim(M, Lat) # latitudinally averaged version
        for j in 1:size(cs, Lat) # south hemisphere
            l = dims(cs, Lat)[j]
            ps_vals = [p[Lat(At(l))] for p in ps]
            if length(dims(M)) > 1
                M[Lat(At(l))] .= predicting_model(hcat(ps_vals...), pfit)
            else
                M[Lat(At(l))] = predicting_model(ps_vals, pfit)
            end
        end
        for j in 1:size(cs, Lat)
            l = dims(cn, Lat)[j]
            pn_vals = [p[Lat(At(l))] for p in pn]
            if length(dims(M)) > 1
                M[Lat(At(l))] .= predicting_model(hcat(pn_vals...), pfit)
            else
                M[Lat(At(l))] = predicting_model(pn_vals, pfit)
            end
        end
    end
end

function split_to_latzones(A, latzones)
    if hasdim(A, Coord)
        map(i -> A[Coord(Between(latzones[i], latzones[i+1]))], 1:length(latzones)-1)
    elseif hasdim(A, Lat)
        map(i -> A[Lat(Between(latzones[i], latzones[i+1]))], 1:length(latzones)-1)
    else
        error("No latitudinal information present in `A`")
    end
end
