const _bss_assets = @path joinpath(ASSETS, "best_sum_spacings")

get_cdf_bss_n_k = let data=load(joinpath(_bss_assets, "data_min_sum_ordered_spacings.jld2"))
    (n::Integer, k::Integer) -> begin

        _present = false
        if n in data["n"]
            _idx_n = findfirst(x -> x==n, data["n"])
            if k in data["k"][_idx_n]
                _present = true
            end
        end

        if _present
            _idx_n = findfirst(x -> x==n, data["n"])
            _idx_k = findfirst(x -> x==k, data["k"][_idx_n])
            cdf_model = Interpolations.interpolate(data["t_val_single"][_idx_n][_idx_k], data["p_val_single"], FritschButlandMonotonicInterpolation())
            return cdf_model
        else

            idx_n = findall(x -> x >= 10, data["n"])

            aux_t_vals = zeros(length(idx_n))
            t_vals = zeros(length(data["p_val_single"]) - 2)

            aux_k = data["n"][idx_n] .* ((k-1) / (n-1)) .+ ((n-k) / (n-1))

            for j in 1:(length(data["p_val_single"]) - 2)
                aux_t_vals .= 0
                for (i, _idx_n) in enumerate(idx_n)
                    _vals = [_vec[j+1] for _vec in data["t_val_single"][_idx_n]]
                    _model_1 = Interpolations.interpolate(data["k"][_idx_n], _vals, FritschButlandMonotonicInterpolation())
                    aux_t_vals[i] = _model_1(aux_k[i])
                end
                _model_2 = Interpolations.interpolate(data["n"][idx_n], aux_t_vals, FritschButlandMonotonicInterpolation())
                t_vals[j] = _model_2(n)
            end

            max_val = 1 / floor((n + 1) / k)
            
            idx_res = findall(y -> y <= max_val, t_vals)
            knots_res = vcat([0], t_vals[idx_res], [max_val])
            sort!(knots_res)
            p_vals_res = vcat([0], data["p_val_single"][idx_res .+ 1], [1])
            cdf_model = Interpolations.interpolate(knots_res, p_vals_res, FritschButlandMonotonicInterpolation())

            return cdf_model
        end
    end
end



function get_cdf_bss_n(n::Integer)

    cdfs = [get_cdf_bss_n_k(n, k) for k in 1:n]
    return cdfs
end






get_cdf_bss_min_n = let data=load(joinpath(_bss_assets, "data_min_sum_ordered_spacings.jld2"))
    (n::Integer) -> begin

        t_vals = zeros(length(data["p_val_best"]) - 2)

        for j in 1:(length(data["p_val_best"]) - 2)
            _vals = [_vec[j+1] for _vec in data["t_val_best"]]
            _model = Interpolations.interpolate(data["n"], _vals, FritschButlandMonotonicInterpolation())
            t_vals[j] = _model(n)
        end


        knots_res = vcat([0], t_vals, [1])
        cdf_model = Interpolations.interpolate(knots_res, data["p_val_best"], FritschButlandMonotonicInterpolation())

        return cdf_model
    end
end


function get_bss_min_TS(n::Integer)

    cdfs = get_cdf_bss_n(n)
    cdf_min = get_cdf_bss_min_n(n)

    TS = let cdfs=cdfs, cdf_min=cdf_min
        (events::AbstractVector{<:Real}) -> begin

            n = length(events)

            _min_spacing = 0.0
            _min_p_val = 2.0
            
            for k in eachindex(cdfs)
                _min_spacing = events[k]
                for i in (k+1):n
                    _min_spacing = min(_min_spacing, events[i] - events[i-k])
                end
                _min_spacing = min(_min_spacing, 1 - events[n+1-k])
                _min_p_val = min(_min_p_val, cdfs[k](_min_spacing))
            end

            res = cdf_min(_min_p_val)
            
            return res
        end
    end

    return TS
end

function get_bss_p_value(events::AbstractVector{<:Real})

    cdfs = get_cdf_bss_n(n)
    cdf_min = get_cdf_bss_min_n(n)

    _min_spacing = 0.0
    _min_p_val = 2.0
    
    for k in eachindex(cdfs)
        _min_spacing = events[k]
        for i in (k+1):n
            _min_spacing = min(_min_spacing, events[i] - events[i-k])
        end
        _min_spacing = min(_min_spacing, 1 - events[n+1-k])
        _min_p_val = min(_min_p_val, cdfs[k](_min_spacing))
    end

    res = cdf_min(_min_p_val)
    
    return res
end


