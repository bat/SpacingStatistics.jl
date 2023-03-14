__f = open(joinpath(_pcs_assets, "product_complementary_spacings__tot_param_lengths.bin"), "r")
__length_n = ntoh(read(__f, Int))
__length_A = ntoh(read(__f, Int))
__size_knots_1 = ntoh(read(__f, Int))
__size_knots_2 = ntoh(read(__f, Int))
__size_knots_3 = ntoh(read(__f, Int))
close(__f)


function pcs_ts(events::AbstractVector{<:Real})
     
    n = length(events)
    
    res = -log(1 - events[1])
    for i in Base.OneTo(n-1)
        res += -log(1 - events[i+1] + events[i])
    end
    res += -log(events[end])
    
    return res
end


min_pcs_ts(n::Integer, _n_dims::Integer) = _n_dims * (n+1) * (log(n+1) - log(n))


cdf_pcs_inf = let L=__length_A, A=ntoh.(Mmap.mmap(open(joinpath(_pcs_assets, "product_complementary_spacings__A.bin"), "r"), Vector{Float64}, (__length_A))), knots_arr=Mmap.mmap(open(joinpath(_pcs_assets, "product_complementary_spacings__knots.bin"), "r"), Array{Float64, 3}, (__size_knots_1, __size_knots_2, __size_knots_3)), m_arr=Mmap.mmap(open(joinpath(_pcs_assets, "product_complementary_spacings__m.bin"), "r"), Array{Float64, 3}, (__size_knots_1, __size_knots_2, __size_knots_3)), c_arr=Mmap.mmap(open(joinpath(_pcs_assets, "product_complementary_spacings__c.bin"), "r"), Array{Float64, 3}, (__size_knots_1 - 1, __size_knots_2, __size_knots_3)), d_arr=Mmap.mmap(open(joinpath(_pcs_assets, "product_complementary_spacings__d.bin"), "r"), Array{Float64, 3}, (__size_knots_1 - 1, __size_knots_2, __size_knots_3))
    (_x::Real, n::Integer, n_dims::Integer) -> begin

        x = 1 / (_x - min_pcs_ts(n, n_dims) + 1)
        
        knots = mappedarray(ntoh, view(knots_arr, 1:L, n, n_dims))
        m = mappedarray(ntoh, view(m_arr, 1:L, n, n_dims))
        c = mappedarray(ntoh, view(c_arr, 1:L-1, n, n_dims))
        d = mappedarray(ntoh, view(d_arr, 1:L-1, n, n_dims))

        res = 0.0

        if x <= knots[begin]
            res = 0.0
        elseif x >= knots[end]
            res = 1.0
        else
            j = searchsortedfirst(knots, x)
            if j > 1
                j -= 1
            end
            xdiff = x - knots[j]
            res = A[j] + m[j]*xdiff + c[j]*xdiff*xdiff + d[j]*xdiff*xdiff*xdiff
        end
        
        return 1 - res
    end
end


function cdf_pcs_inf_poisson(x::Real, μ::Real, n_dims::Integer; ϵ::Real=eps(Float64))
    
    n_min = findlast(a -> min_pcs_ts(a, n_dims) >= x, 1:1000)
    n_min = isnothing(n_min) ? 0 : n_min
    m = ceil(Int, μ)
    p_val = 0
    
    for n in max(n_min+1, m):1000
        tmp_p = pdf(Poisson(μ), n) * cdf_pcs_inf(x, n, n_dims)
        if tmp_p < p_val*ϵ && n>m
            p_val += tmp_p
            break
        end
        p_val += tmp_p
    end
    for n in (max(n_min+1, m)-1):-1:(n_min+1)
        tmp_p = pdf(Poisson(μ), n) * cdf_pcs_inf(x, n, n_dims)
        if tmp_p < p_val*ϵ
            p_val += tmp_p
            break
        end
        p_val += tmp_p
    end
    return clamp(p_val, 0, 1)
end


function cdf_pcs_inf_poisson_max(x::Real, μ::Real, n_dims::Integer; ϵ::Real=eps(Float64))

    n_min = findlast(a -> min_pcs_ts(a, 1) >= x, 1:1000)
    n_min = isnothing(n_min) ? 0 : n_min
    m = ceil(Int, μ)
    p_val = 0
    
    for n in max(n_min+1, m):1000
        tmp_p = pdf(Poisson(μ), n) * cdf(Beta(n_dims, 1), cdf_pcs_inf(x, n, 1))
        if tmp_p < p_val*ϵ && n>m
            p_val += tmp_p
            break
        end
        p_val += tmp_p
    end
    for n in (max(n_min+1, m)-1):-1:(n_min+1)
        tmp_p = pdf(Poisson(μ), n) * cdf(Beta(n_dims, 1), cdf_pcs_inf(x, n, 1))
        if tmp_p < p_val*ϵ
            p_val += tmp_p
            break
        end
        p_val += tmp_p
    end
    return clamp(p_val, 0, 1)
end
