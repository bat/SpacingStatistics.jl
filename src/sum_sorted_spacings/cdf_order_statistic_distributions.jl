function heaviside(x::Real)
    return x>=0 ? x : 0
end

heaviside(x::Interval) = interval(heaviside(x.lo), heaviside(x.hi))

gamma_int(x::Interval) = interval(gamma(x.lo), gamma(x.hi))

function pdf_poisson(λ::Interval, x::Interval)
    return pow(λ, x) * exp(-λ) / gamma_int(x + 1)
end


### SUM OF MINIMUM SORTED GAPS ###


function _cdf_sum_min_with_edge(N::T, k::T, s::T) where {T<:Interval{BigFloat}}

    B0 = interval(BigFloat(0))
    B1 = interval(BigFloat(1))
    B2 = interval(BigFloat(2))
    B3 = interval(BigFloat(3))
    B_1 = interval(BigFloat(-1))

    @argcheck B1 <= k <= N

    if s <= B0
        res = B0
    elseif s < k / (N + B1)
        A = (gamma_int(N + B2) / gamma_int(N + B2 - k)) / pow(N + B1 - k, k - B1)
        res = B0
        M = floor(Int, min(k.lo, ((k + B1 - s * (N + B2)) / (B1 - s)).lo))
        for i in interval.(BigFloat.(Base.OneTo(M)))
            res += pow(B_1, i - B1) * pow(k + B1 - i, k - B1) / ((N + B2 - i) * gamma_int(k + B1 - i) * gamma_int(i)) *
                   (B1 - pow(B1 - s * (N + B2 - i) / (k + B1 - i), N))
        end
        for i in interval.(BigFloat.((M + 1):1:Int(k.lo)))
            res += pow(B_1, i - B1) * pow(k + B1 - i, k - B1) / ((N + B2 - i) * gamma_int(k + B1 - i) * gamma_int(i))
        end
        res *= A
    else
        res = B1
    end

    return res
end


function _cdf_sum_min_no_edge(N::T, k::T, s::T) where {T<:Interval{BigFloat}}

    B0 = interval(BigFloat(0))
    B1 = interval(BigFloat(1))
    B2 = interval(BigFloat(2))
    B3 = interval(BigFloat(3))
    B_1 = interval(BigFloat(-1))

    @argcheck B1 <= k <= N - B2

    if s <= B0
        res = B0
    elseif s < k / (N - B1)
        A = (gamma_int(N) / gamma_int(N - k)) / pow(N - B1 - k, k - B1)
        res = B0
        M = floor(Int, min(k.lo, ((k + B1 - s * N) / (B1 - s)).lo))
        for i in interval.(BigFloat.(Base.OneTo(M)))
        res += pow(B_1, i - B1) * pow(k + B1 - i, k - B1) / ((N - i) * gamma_int(k + B1 - i) * gamma_int(i)) *
               (B1 - pow(B1 - s * (N - i) / (k + B1 - i), N))
        end
        for i in interval.(BigFloat.((M + 1):1:Int(k.lo)))
            res += pow(B_1, i - B1) * pow(k + B1 - i, k - B1) / ((N - i) * gamma_int(k + B1 - i) * gamma_int(i))
        end
        res *= A
    else
        res = B1
    end

    return res
end


function _cdf_sum_min(N::Interval{BigFloat}, k::Interval{BigFloat}, s::Interval{BigFloat}; edge::Bool = false)
    return edge ? _cdf_sum_min_with_edge(N, k, s) : _cdf_sum_min_no_edge(N, k, s)
end


function _cdf_sum_min_poisson(λ::Interval{BigFloat}, k::Interval{BigFloat}, s::Interval{BigFloat}; edge::Bool = false, ϵ::BigFloat = BigFloat(1e-5))

    B0 = interval(BigFloat(0))
    B1 = interval(BigFloat(1))
    B2 = interval(BigFloat(2))
    B3 = interval(BigFloat(3))
    B_1 = interval(BigFloat(-1))

    N_max = cquantile(Poisson(Float64(λ.hi)), Float64(ϵ))

    offset = edge ? B0 : B2
    N = interval.((k + offset).lo:B1.lo:BigFloat(N_max))
    norm = B1 / (B1 - sum(pdf_poisson.(λ, interval.(B0.lo:B1.lo:max(B0, k + offset - B1).lo))))
    ϵ_cor = (ϵ / norm).lo

    res_old = B0
    for n in N
        res_old += pdf_poisson(λ, n) * _cdf_sum_min(n, k, s; edge = edge)
        if !(res_old * norm ⊆ interval(0.0, 1.0))
            break
        end
    end

    res = res_old * norm
    return res
end


function cdf_sum_min(N::Interval{BigFloat}, k::Interval{BigFloat}, s::Interval{BigFloat}; edge::Bool = false, poisson::Bool = false, ϵ::BigFloat = BigFloat(1e-5))
    return poisson ? _cdf_sum_min_poisson(N, k, s; edge = edge, ϵ = ϵ) : _cdf_sum_min(N, k, s; edge = edge)
end


function cdf_sum_min(N::Real, k::Real, s::Real; edge::Bool = false, poisson::Bool = false, ϵ::Real = 1e-5)
    return cdf_sum_min(interval(BigFloat(N)), interval(BigFloat(k)), interval(BigFloat(s)); edge=edge, poisson=poisson, ϵ=BigFloat(ϵ))
end


function cdf_sum_min_set_precision(N::Real, k::Real, s::Real; edge::Bool = false, poisson::Bool = false, ϵ::Real = 1e-5, BF_bits::Integer = 64)
    setprecision(BigFloat, BF_bits)
    return cdf_sum_min(N, k, s; edge=edge, poisson=poisson, ϵ=ϵ)
end


function cdf_sum_min_var_precision(N::Real, k::Real, s::Real; edge::Bool = false, poisson::Bool = false, ϵ::Real = 1e-5, BF_bits_min::Integer = 64, BF_bits_max::Integer = 8192, xtol::Real = 1e-20)
    prc = BF_bits_min
    res = cdf_sum_min_set_precision(N, k, s; edge = edge, poisson = poisson, ϵ = ϵ, BF_bits = prc)
    while (!(res ⊆ interval(0.0, 1.0)) || diam(res)>xtol) && prc <= BF_bits_max
        prc *= 2
        res = cdf_sum_min_set_precision(N, k, s; edge = edge, poisson = poisson, ϵ = ϵ, BF_bits = prc)
    end
    if (!(res ⊆ interval(0.0, 1.0)) || diam(res)>xtol)
        res = interval(NaN)
    end
    return res
end


function cdf_sum_min_var_precision_float(N::Real, k::Real, s::Real; edge::Bool = false, poisson::Bool = false, ϵ::Real = 1e-5, BF_bits_min::Integer = 64, BF_bits_max::Integer = 8192, xtol::Real = 1e-20)
    return Float64(mid(cdf_sum_min_var_precision(N, k, s; edge=edge, poisson=poisson, ϵ=ϵ, BF_bits_min=BF_bits_min, BF_bits_max=BF_bits_max, xtol=xtol)))
end



### SUM OF MAXIMUM SORTED GAPS ###



function _cdf_sum_max_with_edge(N::T, k::T, s::T) where {T<:Interval{BigFloat}}

    B0 = interval(BigFloat(0))
    B1 = interval(BigFloat(1))
    B2 = interval(BigFloat(2))
    B3 = interval(BigFloat(3))
    B_1 = interval(BigFloat(-1))

    @argcheck B1 <= k <= N

    if s <= k / (N + B1)
        res = B0
    elseif s < B1
        A = (gamma_int(N + B2) / gamma_int(k + B1)) / pow(k, N - k)
        res = B0
        M = floor(Int, min((N + B1 - k).lo, ((N + B2) - k / s).lo))
        for i in interval.(BigFloat.(Base.OneTo(M)))
        res += pow(B_1, i - B1) * pow(N + B2 - k - i, N - k) / ((N + B2 - i) * gamma_int(N + B2 - k - i) * gamma_int(i)) *
               pow(B1 - (B1 - s) * (N + B2 - i) / (N + B2 - k - i), N)
        end
        res *= A
    else
        res = B1
    end

    return res
end


function _cdf_sum_max_no_edge(N::T, k::T, s::T) where {T<:Interval{BigFloat}}

    B0 = interval(BigFloat(0))
    B1 = interval(BigFloat(1))
    B2 = interval(BigFloat(2))
    B3 = interval(BigFloat(3))
    B_1 = interval(BigFloat(-1))

    @argcheck B1 <= k <= N - B2

    if s <= B0
        res = B0
    elseif s < B1
        A = N * (gamma_int(N) / gamma_int(k + B1)) / pow(k, N - k)
        res = B0
        M = floor(Int, min((N - B1 - k).lo, max(0, (N - k / s).lo)))
        for i in interval.(BigFloat.(Base.OneTo(M)))
            res += pow(B_1, i - B1) * pow(N - k - i, N - B3 - k) / gamma_int(N - k - i) / gamma_int(i) *
                (
                    (pow(s, N - B1) * (k - s * (N * (k + B1) - i - B2 * k) / N) * (N - k - i)) +
                    (pow(s * (N - i) - k, N) / (pow(N - i - k, N - B3) * N * (N - i)))
                )
        end
        for i in interval.(BigFloat.((M + 1):1:Int((N - B1 - k).lo)))
            res += pow(B_1, i - B1) * pow(N - k - i, N - B3 - k) / gamma_int(N - k - i) / gamma_int(i) *
                   (pow(s, N - B1) * (k - s * (N * (k + B1) - i - B2 * k) / N) * (N - k - i))
        end
        res *= A
    else
        res = B1
    end

    return res
end


function _cdf_sum_max(N::Interval{BigFloat}, k::Interval{BigFloat}, s::Interval{BigFloat}; edge::Bool = false)
    return edge ? _cdf_sum_max_with_edge(N, k, s) : _cdf_sum_max_no_edge(N, k, s)
end


function _cdf_sum_max_poisson(λ::Interval{BigFloat}, k::Interval{BigFloat}, s::Interval{BigFloat}; edge::Bool = false, ϵ::BigFloat = BigFloat(1e-5))

    B0 = interval(BigFloat(0))
    B1 = interval(BigFloat(1))
    B2 = interval(BigFloat(2))
    B3 = interval(BigFloat(3))
    B_1 = interval(BigFloat(-1))

    N_max = cquantile(Poisson(Float64(λ.hi)), Float64(ϵ))

    offset = edge ? B0 : B2
    N = interval.((k + offset).lo:B1.lo:BigFloat(N_max))
    norm = B1 / (B1 - sum(pdf_poisson.(λ, interval.(B0.lo:B1.lo:max(B0, k + offset - B1).lo))))
    ϵ_cor = (ϵ / norm).lo

    res = B0
    for n in N
        res += pdf_poisson(λ, n) * _cdf_sum_max(n, k, s; edge = edge) * norm
        if !(res ⊆ interval(0.0, 1.0))
            break
        end
    end
    return res
end


function cdf_sum_max(N::Interval{BigFloat}, k::Interval{BigFloat}, s::Interval{BigFloat}; edge::Bool = false, poisson::Bool = false, ϵ::BigFloat = BigFloat(1e-5))
    return poisson ? _cdf_sum_max_poisson(N, k, s; edge = edge, ϵ = ϵ) : _cdf_sum_max(N, k, s; edge = edge)
end


function cdf_sum_max(N::Real, k::Real, s::Real; edge::Bool = false, poisson::Bool = false, ϵ::Real = 1e-5)
    return cdf_sum_max(interval(BigFloat(N)), interval(BigFloat(k)), interval(BigFloat(s)); edge=edge, poisson=poisson, ϵ=BigFloat(ϵ))
end


function cdf_sum_max_set_precision(N::Real, k::Real, s::Real; edge::Bool = false, poisson::Bool = false, ϵ::Real = 1e-5, BF_bits::Integer = 64)
    setprecision(BigFloat, BF_bits)
    return cdf_sum_max(N, k, s; edge=edge, poisson=poisson, ϵ=ϵ)
end


function cdf_sum_max_var_precision(N::Real, k::Real, s::Real; edge::Bool = false, poisson::Bool = false, ϵ::Real = 1e-5, BF_bits_min::Integer = 64, BF_bits_max::Integer = 8192, xtol::Real = 1e-20)
    prc = BF_bits_min
    res = cdf_sum_max_set_precision(N, k, s; edge = edge, poisson = poisson, ϵ = ϵ, BF_bits = prc)
    while (!(res ⊆ interval(0.0, 1.0)) || (mid(res) != zero(typeof(mid(res))) ? (diam(res)/mid(res))>xtol : diam(res)>xtol)) && prc <= BF_bits_max
        prc *= 2
        res = cdf_sum_max_set_precision(N, k, s; edge = edge, poisson = poisson, ϵ = ϵ, BF_bits = prc)
    end
    if (!(res ⊆ interval(0.0, 1.0)) || (mid(res) != zero(typeof(mid(res))) ? (diam(res)/mid(res))>xtol : diam(res)>xtol))
        res = interval(NaN)
    end
    return res
end


function cdf_sum_max_var_precision_float(N::Real, k::Real, s::Real; edge::Bool = false, poisson::Bool = false, ϵ::Real = 1e-5, BF_bits_min::Integer = 64, BF_bits_max::Integer = 8192, xtol::Real = eps(Float64))
    return Float64(mid(cdf_sum_max_var_precision(N, k, s; edge=edge, poisson=poisson, ϵ=ϵ, BF_bits_min=BF_bits_min, BF_bits_max=BF_bits_max, xtol=xtol)))
end



### SUM OF MINIMUM AND MIDDLE SORTED GAPS ###


function _cdf_sum_middle_with_edge(N::T, r::T, s::T, x::T) where {T<:Interval{BigFloat}}

    B0 = interval(BigFloat(0))
    B1 = interval(BigFloat(1))
    B2 = interval(BigFloat(2))
    B3 = interval(BigFloat(3))
    B_1 = interval(BigFloat(-1))

    @argcheck B1 <= s <= N
    @argcheck B1 <= r < s

    if x <= B0
        res = B0
    elseif x < (s - r + B1) / (N + B2 - r)
        A = gamma_int(N + B2) / gamma_int(N + B2 - s)/ gamma_int(s - r + B1)/ gamma_int(r) / pow(N + B1 - s, s - r - B1)
        res = B0
        for i in interval.(BigFloat.(0:Int(r.lo - 1)))
            for j in interval.(BigFloat.(0:Int(s.lo - r.lo - 1)))
                F = (j + B1) * (N + B1 - s) - i * (s - r - j)
                H_i = B1 - x * (N + B2 - r + i) / (s - r + B1)
                H_j = B1 - x * (N + B1 - r - j) / (s - r - j)

                if F == B0
                    H = (B1 - pow(heaviside(H_i), N - B1) * ((N + B1 - r - j) / (s - r - j) * (N - B1) * x + B1)) / pow((N + B1 - r - j) / (s - r - j), B2)
                    res += pow(B_1, i+j) * gamma_int(r) / gamma_int(i + B1) / gamma_int(r - i) * gamma_int(s - r + B1) / gamma_int(j + B1) / gamma_int(s - r - j + B1) * pow(s - r - j, s - r - B1) / (s - r + B1) * H
                else
                    H = (s - r - j) / (N - r - j + B1) * pow(heaviside(H_j), N) - (s - r + B1) / (N - r + i + B2) * pow(heaviside(H_i), N) + F / (N - r - j + B1) / (N - r + i + B2)
                    res += pow(B_1, i+j) * gamma_int(r) / gamma_int(i + B1) / gamma_int(r - i) * gamma_int(s - r + B1) / gamma_int(j + B1) / gamma_int(s - r - j + B1) * pow(s - r - j, s - r) * H / F
                end
            end
        end
        res *= A
    else
        res = B1
    end

    return res
end


function _cdf_sum_middle_no_edge(N::T, r::T, s::T, x::T) where {T<:Interval{BigFloat}}

    B0 = interval(BigFloat(0))
    B1 = interval(BigFloat(1))
    B2 = interval(BigFloat(2))
    B3 = interval(BigFloat(3))
    B_1 = interval(BigFloat(-1))

    @argcheck B1 <= s <= N - B2
    @argcheck B1 <= r < s

    if x <= B0
        res = B0
    elseif x < (s - r + B1) / (N - r)
        A = gamma_int(N) / gamma_int(N - s)/ gamma_int(s - r + B1)/ gamma_int(r) / pow(N - B1 - s, s - r - B1)
        res = B0
        for i in interval.(BigFloat.(0:Int(r.lo - 1)))
            for j in interval.(BigFloat.(0:Int(s.lo - r.lo - 1)))
                F = (j + B1) * (N - B1 - s) - i * (s - r - j)
                H_i = B1 - x * (N - r + i) / (s - r + B1)
                H_j = B1 - x * (N - B1 - r - j) / (s - r - j)
                
                if F == B0
                    H = (B1 - pow(heaviside(H_i), N - B1) * ((N - B1 - r - j) / (s - r - j) * (N - B1) * x + B1)) / pow((N - B1 - r - j) / (s - r - j), B2)
                    res += pow(B_1, i + j) * gamma_int(r) / gamma_int(i + B1) / gamma_int(r - i) * gamma_int(s - r + B1) / gamma_int(j + B1) / gamma_int(s - r - j + B1) * pow(s - r - j, s - r - B1) / (s - r + B1) * H
                else
                    H = (s - r - j) / (N - r - j - B1) * pow(heaviside(H_j), N) - (s - r + B1) / (N - r + i) * pow(heaviside(H_i), N) + F / (N - r - j - B1) / (N - r + i)
                    res += pow(B_1, i + j) * gamma_int(r) / gamma_int(i + B1) / gamma_int(r - i) * gamma_int(s - r + B1) / gamma_int(j + B1) / gamma_int(s - r - j + B1) * pow(s - r - j, s - r) * H / F
                end
            end
        end
        res *= A
    else
        res = B1
    end

    return res
end


function _cdf_sum_middle(N::Interval{BigFloat}, r::Interval{BigFloat}, s::Interval{BigFloat}, x::Interval{BigFloat}; edge::Bool = false)
    return edge ? _cdf_sum_middle_with_edge(N, r, s, x) : _cdf_sum_middle_no_edge(N, r, s, x)
end


function _cdf_sum_middle_poisson(λ::Interval{BigFloat}, r::Interval{BigFloat}, s::Interval{BigFloat}, x::Interval{BigFloat}; edge::Bool = false, ϵ::BigFloat = BigFloat(1e-5))

    B0 = interval(BigFloat(0))
    B1 = interval(BigFloat(1))
    B2 = interval(BigFloat(2))
    B3 = interval(BigFloat(3))
    B_1 = interval(BigFloat(-1))

    N_max = cquantile(Poisson(Float64(λ.hi)), Float64(ϵ))

    offset = edge ? B0 : B2
    N = interval.((s + offset).lo:B1.lo:BigFloat(N_max))
    norm = B1 / (B1 - sum(pdf_poisson.(λ, interval.(B0.lo:B1.lo:max(B0, s + offset - B1).lo))))
    ϵ_cor = (ϵ / norm).lo

    res = B0
    for n in N
        res += pdf_poisson(λ, n) * _cdf_sum_middle(n, r, s, x; edge = edge) * norm
        if !(res ⊆ interval(0.0, 1.0))
            break
        end
    end
    return res
end


function cdf_sum_middle(N::Interval{BigFloat}, r::Interval{BigFloat}, s::Interval{BigFloat}, x::Interval{BigFloat}; edge::Bool = false, poisson::Bool = false, ϵ::BigFloat = BigFloat(1e-5))
    return poisson ? _cdf_sum_middle_poisson(N, r, s, x; edge = edge, ϵ = ϵ) : _cdf_sum_middle(N, r, s, x; edge = edge)
end


function cdf_sum_middle(N::Real, r::Real, s::Real, x::Real; edge::Bool = false, poisson::Bool = false, ϵ::Real = 1e-5)
    return cdf_sum_middle(interval(BigFloat(N)), interval(BigFloat(r)), interval(BigFloat(s)), interval(BigFloat(x)); edge=edge, poisson=poisson, ϵ=BigFloat(ϵ))
end


function cdf_sum_middle_set_precision(N::Real, r::Real, s::Real, x::Real; edge::Bool = false, poisson::Bool = false, ϵ::Real = 1e-5, BF_bits::Integer = 64)
    setprecision(BigFloat, BF_bits)
    return cdf_sum_middle(N, r, s, x; edge=edge, poisson=poisson, ϵ=ϵ)
end


function cdf_sum_middle_var_precision(N::Real, r::Real, s::Real, x::Real; edge::Bool = false, poisson::Bool = false, ϵ::Real = 1e-5, BF_bits_min::Integer = 64, BF_bits_max::Integer = 8192, xtol::Real = 1e-20)
    prc = BF_bits_min
    res = cdf_sum_middle_set_precision(N, r, s, x; edge = edge, poisson = poisson, ϵ = ϵ, BF_bits = prc)
    while (!(res ⊆ interval(0.0, 1.0)) || diam(res)>xtol) && prc <= BF_bits_max
        prc *= 2
        res = cdf_sum_middle_set_precision(N, r, s, x; edge = edge, poisson = poisson, ϵ = ϵ, BF_bits = prc)
    end
    if (!(res ⊆ interval(0.0, 1.0)) || diam(res)>xtol)
        res = interval(NaN)
    end
    return res
end


function cdf_sum_middle_var_precision_float(N::Real, r::Real, s::Real, x::Real; edge::Bool = false, poisson::Bool = false, ϵ::Real = 1e-5, BF_bits_min::Integer = 64, BF_bits_max::Integer = 8192, xtol::Real = 1e-20)
    return Float64(mid(cdf_sum_middle_var_precision(N, r, s, x; edge=edge, poisson=poisson, ϵ=ϵ, BF_bits_min=BF_bits_min, BF_bits_max=BF_bits_max, xtol=xtol)))
end



### MINIMUM SORTED GAPS ###


function _cdf_minimum_with_edge(N::T, k::T, x::T) where {T<:Interval{BigFloat}}

    B0 = interval(BigFloat(0))
    B1 = interval(BigFloat(1))
    B2 = interval(BigFloat(2))
    B3 = interval(BigFloat(3))
    B_1 = interval(BigFloat(-1))

    @argcheck B1 <= k <= N + B1

    if x <= B0
        res = B0
    elseif x < B1 / (N + B2 - k)
        # M = max(0, Int(floor(((x * (N + B2) - B1) / x).lo)))
        M = floor(Int, max(0, (N + B2 - B1 / x).lo))
        A = (N + B1) * gamma_int(N + B1) / gamma_int(k) / gamma_int(N - k + B2)
        res = B0
        # for i in interval.(BigFloat.(Base.OneTo(Int(k.lo))))
        #     res += pow(B_1, k - i) * gamma_int(k) / gamma_int(i) / gamma_int(k - i + B1) / (N + B2 - i) * (B1 - pow(heaviside(B1 - x * (N + B2 - i)), N))
        # end

        for i in interval.(BigFloat.(Base.OneTo(M)))
            res += pow(B_1, k - i) * gamma_int(k) / gamma_int(i) / gamma_int(k - i + B1) / (N + B2 - i)
        end

        for i in interval.(BigFloat.((M + 1):1:Int(k.lo)))
            res += pow(B_1, k - i) * gamma_int(k) / gamma_int(i) / gamma_int(k - i + B1) / (N + B2 - i) * (B1 - pow(B1 - x * (N + B2 - i), N))
        end
        res *= A
    else
        res = B1
    end
    
    return res
end


function _cdf_minimum_no_edge(N::T, k::T, x::T) where {T<:Interval{BigFloat}}

    B0 = interval(BigFloat(0))
    B1 = interval(BigFloat(1))
    B2 = interval(BigFloat(2))
    B3 = interval(BigFloat(3))
    B_1 = interval(BigFloat(-1))

    @argcheck B1 <= k <= N - B1
    
    if x <= B0
        res = B0
    elseif x < B1 / (N - k)
        M = floor(Int, max(0, (N - B1 / x).lo))
        A = (N - B1) * gamma_int(N - B1) / gamma_int(k) / gamma_int(N - k)
        res = B0
        # for i in interval.(BigFloat.(Base.OneTo(Int(k.lo))))
        #     res += pow(B_1, k - i) * gamma_int(k) / gamma_int(i) / gamma_int(k - i + B1) / (N - i) * (B1 - pow(heaviside(B1 - x * (N - i)), N))
        # end

        for i in interval.(BigFloat.(Base.OneTo(M)))
            res += pow(B_1, k - i) * gamma_int(k) / gamma_int(i) / gamma_int(k - i + B1) / (N - i)
        end

        for i in interval.(BigFloat.((M + 1):1:Int(k.lo)))
            res += pow(B_1, k - i) * gamma_int(k) / gamma_int(i) / gamma_int(k - i + B1) / (N - i) * (B1 - pow(B1 - x * (N - i), N))
        end
        res *= A
    else
        res = B1
    end
    
    return res
end


function _cdf_minimum(N::Interval{BigFloat}, k::Interval{BigFloat}, x::Interval{BigFloat}; edge::Bool = false)
    return edge ? _cdf_minimum_with_edge(N, k, x) : _cdf_minimum_no_edge(N, k, x)
end


function _cdf_minimum_poisson(λ::Interval{BigFloat}, k::Interval{BigFloat}, x::Interval{BigFloat}; edge::Bool = false, ϵ::BigFloat = BigFloat(1e-5))

    B0 = interval(BigFloat(0))
    B1 = interval(BigFloat(1))
    B2 = interval(BigFloat(2))
    B3 = interval(BigFloat(3))
    B_1 = interval(BigFloat(-1))

    N_max = cquantile(Poisson(Float64(λ.hi)), Float64(ϵ))

    offset = edge ? B_1 : B1
    N = interval.((k + offset).lo:B1.lo:BigFloat(N_max))
    norm = B1 / (B1 - sum(pdf_poisson.(λ, interval.(B0.lo:B1.lo:max(B0, k + offset - B1).lo))))
    ϵ_cor = (ϵ / norm).lo

    res = B0
    for n in N
        res += pdf_poisson(λ, n) * _cdf_minimum(n, k, x; edge = edge) * norm
        if !(res ⊆ interval(0.0, 1.0))
            break
        end
    end
    return res
end


function cdf_minimum(N::Interval{BigFloat}, k::Interval{BigFloat}, x::Interval{BigFloat}; edge::Bool = false, poisson::Bool = false, ϵ::BigFloat = BigFloat(1e-5))
    return poisson ? _cdf_minimum_poisson(N, k, x; edge = edge, ϵ = ϵ) : _cdf_minimum(N, k, x; edge = edge)
end


function cdf_minimum(N::Real, k::Real, x::Real; edge::Bool = false, poisson::Bool = false, ϵ::Real = 1e-5)
    return cdf_minimum(interval(BigFloat(N)), interval(BigFloat(k)), interval(BigFloat(x)); edge=edge, poisson=poisson, ϵ=BigFloat(ϵ))
end


function cdf_minimum_set_precision(N::Real, k::Real, x::Real; edge::Bool = false, poisson::Bool = false, ϵ::Real = 1e-5, BF_bits::Integer = 64)
    setprecision(BigFloat, BF_bits)
    return cdf_minimum(N, k, x; edge=edge, poisson=poisson, ϵ=ϵ)
end


function cdf_minimum_var_precision(N::Real, k::Real, x::Real; edge::Bool = false, poisson::Bool = false, ϵ::Real = 1e-5, BF_bits_min::Integer = 64, BF_bits_max::Integer = 8192, xtol::Real = 1e-20)
    prc = BF_bits_min
    res = cdf_minimum_set_precision(N, k, x; edge = edge, poisson = poisson, ϵ = ϵ, BF_bits = prc)
    while (!(res ⊆ interval(0.0, 1.0)) || diam(res)>xtol) && prc <= BF_bits_max
        prc *= 2
        res = cdf_minimum_set_precision(N, k, x; edge = edge, poisson = poisson, ϵ = ϵ, BF_bits = prc)
    end
    if (!(res ⊆ interval(0.0, 1.0)) || diam(res)>xtol)
        res = interval(NaN)
    end
    return res
end


function cdf_minimum_var_precision_float(N::Real, k::Real, x::Real; edge::Bool = false, poisson::Bool = false, ϵ::Real = 1e-5, BF_bits_min::Integer = 64, BF_bits_max::Integer = 8192, xtol::Real = 1e-20)
    return Float64(mid(cdf_minimum_var_precision(N, k, x; edge=edge, poisson=poisson, ϵ=ϵ, BF_bits_min=BF_bits_min, BF_bits_max=BF_bits_max, xtol=xtol)))
end


### SORTED ORDER STATISTICS ###


function cdf_sorted_gaps(N::Real, r::Real, s::Real, x::Real; edge::Bool = false, poisson::Bool = false, ϵ::Real = 1e-5, BF_bits_min::Integer = 64, BF_bits_max::Integer = 8192, xtol::Real = 1e-20)
    
    res = 0
    offset = edge ? 1 : -1

    if r == s && s !== N + offset
        res = cdf_minimum_var_precision_float(N, s, x; edge=edge, poisson=poisson, ϵ=ϵ, BF_bits_min=BF_bits_min, BF_bits_max=BF_bits_max, xtol=xtol)
    elseif !(edge) && r==1 && s==(N-1)
        res = beta_inc(s-r+1, N-s+r, x, 1-x)[1]
    elseif r == 1
        res = cdf_sum_min_var_precision_float(N, s, x; edge=edge, poisson=poisson, ϵ=ϵ, BF_bits_min=BF_bits_min, BF_bits_max=BF_bits_max, xtol=xtol)
    elseif s == N + offset
        res = cdf_sum_max_var_precision_float(N, s-r+1, x; edge=edge, poisson=poisson, ϵ=ϵ, BF_bits_min=BF_bits_min, BF_bits_max=BF_bits_max, xtol=xtol)
    else
        res = cdf_sum_middle_var_precision_float(N, r, s, x; edge=edge, poisson=poisson, ϵ=ϵ, BF_bits_min=BF_bits_min, BF_bits_max=BF_bits_max, xtol=xtol)
    end

    return res
end

function limits_sorted_gaps(N::Real, r::Real, s::Real, edge::Bool = false)
    
    lower = 0
    upper = 1
    offset = edge ? 1 : -1
    
    if s == N + offset
        lower = edge ? (s - r + 1) / (N + 1) : 0
        upper = 1
    else
        lower = 0
        upper = (s - r + 1) / (N + offset - r + 1)
    end
    return (lower, upper)
end





### ORDER STATISTICS ###


function _cdf_ordered_gaps_with_edge(N::Integer, r::Integer, s::Integer, x::Real)

    @argcheck 1 <= s <= N+1
    @argcheck 1 <= r <= s

    res = beta_inc(s-r+1, N-s+r, x, 1-x)[1]
    return res
end

function _cdf_ordered_gaps_no_edge(N::Integer, r::Integer, s::Integer, x::Real)

    @argcheck 1 <= s <= N-1
    @argcheck 1 <= r <= s

    res = beta_inc(s-r+1, N-s+r, x, 1-x)[1]
    return res
end


function cdf_ordered_gaps(N::Integer, r::Integer, s::Integer, x::Real; edge::Bool = false)
    return edge ? _cdf_ordered_gaps_with_edge(N, r, s, x) : _cdf_ordered_gaps_no_edge(N, r, s, x)
end