function get_m_c_d(knots, A; it=FritschCarlsonMonotonicInterpolation())

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


if !isfile(joinpath(_pcs_assets, "product_complementary_spacings__m.bin")) || !isfile(joinpath(_pcs_assets, "product_complementary_spacings__c.bin")) || !isfile(joinpath(_pcs_assets, "product_complementary_spacings__d.bin"))

    __f = open(joinpath(_pcs_assets, "product_complementary_spacings__tot_param_lengths.bin"), "r")
    __length_n = ntoh(read(__f, Int))
    __length_A = ntoh(read(__f, Int))
    __size_knots_1 = ntoh(read(__f, Int))
    __size_knots_2 = ntoh(read(__f, Int))
    __size_knots_3 = ntoh(read(__f, Int))
    close(__f)

    A = mappedarray(ntoh, view(Mmap.mmap(open(joinpath(_pcs_assets, "product_complementary_spacings__A.bin"), "r"), Vector{Float64}, (__length_A)), 1:__length_A))
    knots_arr = Mmap.mmap(open(joinpath(_pcs_assets, "product_complementary_spacings__knots.bin"), "r"), Array{Float64, 3}, (__size_knots_1, __size_knots_2, __size_knots_3))

    m = zeros(__size_knots_1, __size_knots_2, __size_knots_3)
    c = zeros(__size_knots_1 - 1, __size_knots_2, __size_knots_3)
    d = zeros(__size_knots_1 - 1, __size_knots_2, __size_knots_3)

    for i in 1:__size_knots_3
        for j in 1:__size_knots_2
            knots = ntoh.(knots_arr[:, j, i])
            _m, _c, _d = get_m_c_d(knots, A)
            m[:, j, i] .= _m
            c[:, j, i] .= _c
            d[:, j, i] .= _d
        end
    end

    s = open(joinpath(_pcs_assets, "product_complementary_spacings__m.bin"), "w+")
    write(s, hton.(m))
    close(s)

    s = open(joinpath(_pcs_assets, "product_complementary_spacings__c.bin"), "w+")
    write(s, hton.(m))
    close(s)

    s = open(joinpath(_pcs_assets, "product_complementary_spacings__d.bin"), "w+")
    write(s, hton.(m))
    close(s)

end



