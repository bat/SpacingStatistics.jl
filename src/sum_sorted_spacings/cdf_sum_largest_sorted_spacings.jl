function cdf_sum_largest_sorted_spacings_poisson(μ::Real, k::Integer, x::Real; N_max::Integer = 750, ϵ::Real=eps(Float64))

    p = 0.0
    dist = Poisson(μ)
    m = ceil(Int, μ)

    for n in max(k, m):N_max
        tmp_p = pdf(dist, n) * cdf_sum_max_var_precision_float(n, k, x; edge=true, poisson=false, BF_bits_max=32768, xtol=eps(Float64))

        if tmp_p < p*ϵ && n>m
            p += tmp_p
            break
        end
        p += tmp_p
    end

    for n in (max(k, m)-1):-1:k
        tmp_p = pdf(dist, n) * cdf_sum_max_var_precision_float(n, k, x; edge=true, poisson=false, BF_bits_max=32768, xtol=eps(Float64))

        if tmp_p < p*ϵ
            p += tmp_p
            break
        end
        p += tmp_p
    end

    return clamp(p, 0, 1)
end


function max_cdf_sum_largest_sorted_spacings_poisson(μ::Real, events::AbstractVector{<:Real}; N_max::Integer = 750)

    N_observed = length(events)
    @argcheck N_observed <= N_max

    gaps = reverse(sort(diff(vcat([0], events, [1]))))

    p = 0.0

    for k in Base.OneTo(N_observed)
        tmp_p = cdf_sum_largest_sorted_spacings_poisson(μ, k, sum(gaps[1:k]); N_max=N_max)
        p = max(tmp_p, p)
    end

    return p
end


function get_model_SSLS_poisson(p::Real; n_dims::Integer = 1)

    data_mean = CSV.File(joinpath(_sum_spacings_assets, "SLSS__$(n_dims)D__sampling_table__mean__p_to_x.csv"))
    N = parse.(Float64, string.(keys(data_mean[1])[2:end]))
    p_vals = data_mean.p_vals

    @argcheck p <= p_vals[end]

    idx = searchsortedfirst(p_vals, p)
    n0 = -log(1 - p)
    idxs = findall(x -> x > n0, N)

    x = nothing

    if p in p_vals
        x = values(data_mean[idx])[idxs]
    else
        x_up = values(data_mean[idx])[idxs]
        x_down = values(data_mean[idx-1])[idxs]
        x = x_down .+ (x_up .- x_down) .* (p - p_vals[idx-1]) ./ (p_vals[idx] - p_vals[idx-1])
    end

    model = LinearInterpolation(vcat([n0], N[idxs]), vcat([p], x))

    data_mean = []
    GC.gc()

    return model
end




