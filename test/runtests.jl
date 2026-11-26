using Test
using StandingColumnWell

@testset "StandingColumnWell.jl" begin
    # SCW of the Beaudry validation case (script/script_validation_GBeaudry.jl)
    t = collect(3600.0:3600.0:(3600.0 * 24 * 365))
    V, β, T₀ = 0.001633, 0.1, 10.0
    ΔHₚ, rₚ = 0.22, 100.0
    rb, ri, ro, H = 0.0855, 0.0247, 0.0302, 213.0
    kp, ϵb, ϵp = 0.4, 0.0, 0.0
    k = [73.0 2.54 2.17e6; 71.0 2.74 2.27e6; H - 144.0 2.94 2.37e6]
    K = [7.0 1e-5; H - 7.0 1e-7]

    @testset "borehole_resistance_scw" begin
        Rfb, Rv, f = borehole_resistance_scw(V, H, rb, ro, ri, kp, T₀, ϵb, ϵp)
        @test Rfb > 0 && Rv > 0 && f > 0
        # Hydraulic resistance Rv = H / (ρf cf V), Jacques et al. Eq. 18
        @test Rv ≈ H / (4.19e6 * V) rtol = 0.01
    end

    @testset "convergence_flow" begin
        # Homogeneous aquifer: balance point at mid-depth, effective properties unchanged
        H̃₁, H̃₂, K̃₁, K̃₂ = convergence_flow([H 1e-6])
        @test H̃₁ ≈ H / 2 && H̃₂ ≈ H / 2
        @test K̃₁ ≈ 1e-6 && K̃₂ ≈ 1e-6
        # Layered aquifer: the lengths still sum to H
        H̃₁, H̃₂, _, _ = convergence_flow([20.0 1e-5; 80.0 1e-7])
        @test H̃₁ > 0 && H̃₂ > 0
    end

    @testset "βils and βils_outlet" begin
        Rfb, Rv, f = borehole_resistance_scw(V, H, rb, ro, ri, kp, T₀, ϵb, ϵp)
        gb = βils(t, V, β, H, rb, ro, k, K, f, ΔHₚ, rₚ)
        @test length(gb) == length(t)
        @test all(gb .>= 0)

        # The three outlet methods give the same g-function, never negative
        g1 = βils_outlet(gb, Rfb, Rv, β)
        g2 = βils_outlet(t, V, β, H, rb, ro, k, K, Rfb, Rv, f, ΔHₚ, rₚ)
        g3 = βils_outlet(t, V, β, H, rb, ro, ri, k, kp, K, T₀, ϵb, ϵp, ΔHₚ, rₚ)
        @test g1 ≈ g2 ≈ g3
        @test all(g1 .>= 0)

        # Input checks on the layer matrices
        @test_throws ErrorException βils(t, V, β, H, rb, ro, [H 2.5], K, f, ΔHₚ, rₚ)
        @test_throws ErrorException βils(t, V, β, H, rb, ro, k, [H 1e-6 0.0], f, ΔHₚ, rₚ)
    end
end
