# This file is a part of SpacingStatistics.jl, licensed under the MIT License (MIT).


using FileIO
using Mmap
using Interpolations
using SciPy
using Roots
using ArgCheck
using Distributions

const _rps_assets = @path joinpath(ASSETS, "rps")


function minimum_rps_ts(N::Integer)
    return sum(j * log(j) for j in Base.OneTo(N+1))
end

function prod_spacings(events::AbstractVector{<:Real})
    TS::float(eltype(events)) = 0
    @inbounds @simd for i in eachindex(events)[begin:end-1]
        TS -= log(events[i+1] - events[i])
    end
    return TS
end

function _reduce!(events::AbstractVector{<:Real}, n::Integer)
    idxs = eachindex(events)
    @inbounds for i in idxs[begin:begin+n-2]
        events[i] = (events[i] + events[i+1]) / 2
    end
    a = events[first(idxs)]
    b = events[first(idxs)+n-2]
    inv_c = inv(b - a)
    inv_d = - a * inv_c
    @inbounds for i in idxs[begin:begin+n-2]
        events[i] = muladd(events[i], inv_c, inv_d)
    end
    nothing
end

function prod_spacings_layers(eventos::AbstractVector{<:Real})
    
    events = vcat([0], sort(eventos), [1])
    n = length(events)
    TS::float(eltype(events)) = 0
    idxs = eachindex(events)
    
    for i in idxs[begin:end-1]
        TS += prod_spacings(view(events, idxs[begin:begin+n-1]))
        _reduce!(events, n)
        n -= 1
    end

    return TS
end



function _get_rps_models()

    degree = 3

    s_knots = open(joinpath(_rps_assets, "RPS_table_fit_Bspline_knots_hton.bin"), "r")
    n_rows_knots = ntoh(read(s_knots, Int64))
    n_cols_knots = ntoh(read(s_knots, Int64))
    data_knots = ntoh.(Mmap.mmap(s_knots, Matrix{Float64}, (n_rows_knots, n_cols_knots)))

    s_coefficients = open(joinpath(_rps_assets, "RPS_table_fit_Bspline_coefficients_hton.bin"), "r")
    n_rows_coefficients = ntoh(read(s_coefficients, Int64))
    n_cols_coefficients = ntoh(read(s_coefficients, Int64))
    data_coefficients = ntoh.(Mmap.mmap(s_coefficients, Matrix{Float64}, (n_rows_coefficients, n_cols_coefficients)))

    models = [SciPy.interpolate.BSpline(filter(a -> a >= 0.0, data_knots[i, 2:end]), filter(a -> a >= 0.0, data_coefficients[i, 2:end]), degree) for i in Base.OneTo(n_rows_knots)]

    return models
end

_get_rps_x_vals = let m_N = _get_rps_models()
    (N::Integer) -> begin
        return [m_N[i](log10(N))[] for i in eachindex(m_N)]
    end
end


_rps_s_p_values = open(joinpath(_rps_assets, "RPS_table_fit_p_values_hton.bin"), "r")
_rps_n_rows_p_values = ntoh(read(_rps_s_p_values, Int64))
_rps_p = ntoh.(Mmap.mmap(_rps_s_p_values, Vector{Float64}, (_rps_n_rows_p_values));)

specific_rps_cdf = let p = _rps_p
    (N::Integer) -> begin
        x = _get_rps_x_vals(N)
        specific_cdf = Interpolations.interpolate(vcat([0], x, [1]), vcat([0], p, [1]), FritschButlandMonotonicInterpolation())
        return specific_cdf
    end
end

function specific_rps_quantile(N::Integer)
    specific_cdf = specific_rps_cdf(N)
    specific_quantile = let model = specific_cdf
        (q::Real) -> begin
            @argcheck 0<=q<=1
            f(x) = model(x) - q
            x = find_zero(f, (0,1))
            return x
        end
    end
    return specific_quantile
end

function rps_cdf(
    x::Real, 
    N::Integer
) 
    return specific_rps_cdf(N)(x)
end


_rps_s_N_values = open(joinpath(_rps_assets, "RPS_table_fit_N_values_hton.bin"), "r")
_rps_n_rows_N_values = ntoh(read(_rps_s_N_values, Int64))
_rps_N_max = Int(maximum(ntoh.(Mmap.mmap(_rps_s_N_values, Vector{Int64}, (_rps_n_rows_N_values)))))


rps_cdf_fast = let specific_cdf = [specific_rps_cdf(i) for i in Base.OneTo(_rps_N_max)]
    (x::Real, N::Integer) -> begin
        return specific_cdf[N](x)
    end
end

function rps_quantile(
    q::Real, 
    N::Integer
)
    return specific_rps_quantile(N)(q)
end

rps_quantile_fast = let specific_quantile = [specific_rps_quantile(i) for i in Base.OneTo(_rps_N_max)]
    (x::Real, N::Integer) -> begin
        return specific_quantile[N](x)
    end
end

function rps_ts(events::AbstractVector{<:Real})
    return minimum_rps_ts(length(events)) / prod_spacings_layers(events)
end

function rps_p_value(events::AbstractVector{<:Real})
    return rps_cdf(rps_ts(events), length(events))
end

function rps_p_value_fast(events::AbstractVector{<:Real})
    return rps_cdf_fast(rps_ts(events), length(events))
end


"""
    SpacingStatistics.rps(events::AbstractVector{<:Real}, dist::UnivariateDistribution)

Returns the RPS test statistic and associated p-value for a list of events assuming Distribution dist as the null-hypotheis.

```jldoctest
using SpacingStatistics
using Distributions

events = [0.1, 0.4, 0.56, 0.3, 0.12]
SpacingStatistics.rps(events, Unifrom())

# output

(0.901005261912951, 0.5912765565468424)
```
"""
function rps(events::AbstractVector{<:Real}, dist::UnivariateDistribution)
    unit_events = sort(cdf.(dist, events))
    rps_ts_val = rps_ts(unit_events)
    rps_p_val = rps_cdf_fast(rps_ts_val, length(events))
    return (rps_ts_val, rps_p_val)
end






