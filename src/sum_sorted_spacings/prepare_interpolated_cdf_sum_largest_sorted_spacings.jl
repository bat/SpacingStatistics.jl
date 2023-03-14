using FileIO
using Mmap
using Interpolations
using Distributions
using MappedArrays
using ArgCheck
using CSV

###########

function get_index__sum_max__k_ordering(
    n::Integer,
    k::Integer;
    N_max::Integer = 750,
)

    res = Int((N_max+1)*(k-1) - (k+1)*k/2 + n+1)
    return res
end


function get_m_c_d(knots, A, it)
    
    TWeights = Interpolations.tweight(A)
    TCoeffs1 = typeof(oneunit(eltype(A)) / oneunit(eltype(knots)))
    TCoeffs2 = typeof(oneunit(eltype(A)) / oneunit(eltype(knots))^2)
    TCoeffs3 = typeof(oneunit(eltype(A)) / oneunit(eltype(knots))^3)
    
    n = length(knots)
    
    m, Δ = Interpolations.calcTangents(TCoeffs1, knots, A, it)
    
    c = Vector{TCoeffs2}(undef, n-1)
    
    d = Vector{TCoeffs3}(undef, n-1)
    
    for k ∈ eachindex(c)
        if it == LinearMonotonicInterpolation
            c[k] = zero(TCoeffs2)
            d[k] = zero(TCoeffs3)
        else
            xdiff = knots[k+1] - knots[k]
            c[k] = (3*Δ[k] - 2*m[k] - m[k+1]) / xdiff
            d[k] = (m[k] + m[k+1] - 2*Δ[k]) / (xdiff * xdiff)
        end
    end
    
    return m, c, d
end 

__f = open("/Users/loliansh/Documents/transient/sum_sorted_gaps_interpolation_parameters_v2/stored_cdf_coeff_sum_sorted_max_edge_1_1_750__partial_param_lengths.bin", "r")
__length_L = ntoh(read(__f, Int))
__length_knots = ntoh(read(__f, Int))
__length_A = ntoh(read(__f, Int))
close(__f)



partial_param_lengths = [
    __length_L,
    __length_knots,
    __length_A,
]
s = open("/Users/loliansh/Documents/transient/sum_sorted_gaps_interpolation_parameters_v2/stored_cdf_coeff_sum_sorted_max_edge_1_1_750__partial_param_lengths.bin", "w+")
write(s, hton.(partial_param_lengths))
close(s)




L_vec = Mmap.mmap(open("/Users/loliansh/Documents/transient/sum_sorted_gaps_interpolation_parameters_v2/stored_cdf_coeff_sum_sorted_max_edge_1_1_750__L.bin", "r"), Vector{Int64}, (__length_L))
knots_vec = Mmap.mmap(open("/Users/loliansh/Documents/transient/sum_sorted_gaps_interpolation_parameters_v2/stored_cdf_coeff_sum_sorted_max_edge_1_1_750__knots.bin", "r"), Vector{Float64}, (__length_knots))
A_vec = Mmap.mmap(open("/Users/loliansh/Documents/transient/sum_sorted_gaps_interpolation_parameters_v2/stored_cdf_coeff_sum_sorted_max_edge_1_1_750__A.bin", "r"), Vector{Float64}, (__length_A))

if !isfile("/Users/loliansh/Documents/transient/sum_sorted_gaps_interpolation_parameters_v2/stored_cdf_coeff_sum_sorted_max_edge_1_1_750__m.bin") || !isfile("/Users/loliansh/Documents/transient/sum_sorted_gaps_interpolation_parameters_v2/stored_cdf_coeff_sum_sorted_max_edge_1_1_750__c.bin") || !isfile("/Users/loliansh/Documents/transient/sum_sorted_gaps_interpolation_parameters_v2/stored_cdf_coeff_sum_sorted_max_edge_1_1_750__d.bin")

    get_m_c_d_vec = let L_vec = Mmap.mmap(open("/Users/loliansh/Documents/transient/sum_sorted_gaps_interpolation_parameters_v2/stored_cdf_coeff_sum_sorted_max_edge_1_1_750__L.bin", "r"), Vector{Int64}, (__length_L)), knots_vec = Mmap.mmap(open("//Users/loliansh/Documents/transient/sum_sorted_gaps_interpolation_parameters_v2/stored_cdf_coeff_sum_sorted_max_edge_1_1_750__knots.bin", "r"), Vector{Float64}, (__length_knots)), A_vec = Mmap.mmap(open("/Users/loliansh/Documents/transient/sum_sorted_gaps_interpolation_parameters_v2/stored_cdf_coeff_sum_sorted_max_edge_1_1_750__A.bin", "r"), Vector{Float64}, (__length_A))

        () -> begin

            m_vec = zeros(length(A_vec))
            c_vec = zeros(length(A_vec) - length(L_vec))
            d_vec = zeros(length(A_vec) - length(L_vec))

            for idx in eachindex(L_vec)

                L_start_long = idx==1 ? 1 : (ntoh(L_vec[idx - 1]) + 1)
                L_start_short = idx==1 ? 1 : (L_start_long - idx + 1)
                L_stop_long = ntoh(L_vec[idx])
                L_stop_short = L_stop_long - idx
                
                knots = mappedarray(ntoh, view(knots_vec, L_start_long:L_stop_long))
                A = mappedarray(ntoh, view(A_vec, L_start_long:L_stop_long))

                _m, _c, _d = get_m_c_d(knots, A, FritschCarlsonMonotonicInterpolation())


                m_vec[L_start_long:L_stop_long] .= _m
                c_vec[L_start_short:L_stop_short] .= _c
                d_vec[L_start_short:L_stop_short] .= _d
            end

            s = open("/Users/loliansh/Documents/transient/sum_sorted_gaps_interpolation_parameters_v2/stored_cdf_coeff_sum_sorted_max_edge_1_1_750__m_lol.bin", "w+")
            write(s, hton.(m_vec))
            close(s)

            s = open("/Users/loliansh/Documents/transient/sum_sorted_gaps_interpolation_parameters_v2/stored_cdf_coeff_sum_sorted_max_edge_1_1_750__c_lol.bin", "w+")
            write(s, hton.(c_vec))
            close(s)

            s = open("/Users/loliansh/Documents/transient/sum_sorted_gaps_interpolation_parameters_v2/stored_cdf_coeff_sum_sorted_max_edge_1_1_750__d_lol.bin", "w+")
            write(s, hton.(d_vec))
            close(s)

            tot_param_lengths = [
                length(L_vec),
                length(knots_vec),
                length(A_vec),
                length(m_vec),
                length(c_vec),
                length(d_vec),
            ]
            s = open("/Users/loliansh/Documents/transient/sum_sorted_gaps_interpolation_parameters_v2/stored_cdf_coeff_sum_sorted_max_edge_1_1_750__tot_param_lengths.bin", "w+")
            write(s, hton.(tot_param_lengths))
            close(s)

            nothing

        end
    end

    @time get_m_c_d_vec()

end



