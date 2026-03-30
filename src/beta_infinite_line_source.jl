using SpecialFunctions: erf

"""
This file also uses the infinite line source (ILS) model and the convolution product implemented
in the GroundHeatExchanger.jl
"""

"""
    βils_outlet(t, k, kp, K, rb, ro, ri, H, Hp, V, β, T₀, dT, ϵb, ϵp)

Compute the borehole outlet transfer function from the β-ILS model and the thermal resistance
of the SCW. The function computes the borehole wall (g-function) using the β-ILS model, then
adjusts it to obtain the outlet transfer function using the development from Jacques et al. (2025).
The obtained transfer function is for an impulse of 1 W/m.
# Arguments
    - `t`: Time vector (nₜ x 1) [s]
    - `k`: A (ix3) matrix of "i" geological layers and their thermal properties:
        - Column 1: Thickness [m]
        - Column 2: Thermal conductivity [W/mK]
        - Column 3: Volumetric specific heat [J/m³K]
        - Example: [20 k₁ C₁; 80 K₂ C₂] for 2 layers of 20m and 80m thicknesses
        - Note: For a single layer ground, the first value is equal to `H`.
    - `kp`: Pipe thermal conductivity [W/mK]
    - `K`: A (jx2) matrix of "j" hydrogeological layers and their hydraulic properties:
        - Column 1: Thickness [m]
        - Column 2: Hydraulic conductivitiy [m/s]
        - Example: [20 K₁; 80 K₂] for 2 layers of 20m and 80m thicknesses
    - `rb`: Borehole radius [m]
    - `ro`: Outer pipe radius [m]
    - `ri`: Inner pipe radius [m]
    - `H`: Borehole depth [m]
    - `Hp`: Pipe length [m]
    - `V`: Circulating fluid flow rate [m³/s]
    - `β`: Fluid bleed rate ratio (between 0.01 and 1) [-]
    - `T₀`: Initial ground temperature [°C] (used for water properties)
    - `dT`: Temperature difference at the heat exchanger (for fluid properties) [K]
        - Positive value: Ground heating -> the fluid is being cooled (Tb<Tf)
        - Negative value: Ground cooling -> the fluid is being heating (Tb>Tf)
    - `ϵb`: Borehole roughness height [m]
    - `ϵp`: Pipe roughness height [m]
# Output
    - `g̃`: A transfer function for the outlet fluid temperature of SCW [°Cm/W]
# Reference
    - Nguyen, A., Jacques, L., & Pasquier, P. (2025). An easy-to-use analytical model for 
        standing column wells operating with bleed. Applied Thermal Engineering.
        https://doi.org/10.1016/j.applthermaleng.2024.124543
    - Jacques, L., Pasquier, P., Nguyen, A., & Beaudry, G. (2025). Improvement and experimental 
        validation of an analytical model for standing column wells operated with bleed. Applied 
        Thermal Engineering, 279, 127620.
"""
function βils_outlet(t, k::AbstractMatrix{<:Real}, kp::Real, K::AbstractMatrix{<:Real}, rb::Real,
    ro::Real, ri::Real, H::Real, Hp::Real, V::Real, β::Real, T₀::Real, dT::Real, ϵb::Real, ϵp::Real)
    # 1. Compute the thermal resistance of the SCW
    Rfb, Rv, f = borehole_resistance_scw(V, kp, rb, ro, ri, H, Hp, T₀, dT, ϵb, ϵp)

    # 2. Compute the borehole wall g-function using the β-ILS model
    gb = βils(t, k, K, rb, ro, H, V, β, f) # Jacques et al. Eq. 8

    # 3. Adjust the borehole wall g-function to obtain the outlet g-function for 1 W/m (g̃ᵢ ≥ 0)
    return @. max(gb - (1 - β) * (Rv / (2 - β) - Rfb), 0.0) # Jacques et al. Eq. 27
end

"""
    βils(t, k, K, rb, ro, H, V, β, f)

Compute the analytical model developped by Nguyen et al. (2025) for standing column
wells (SCW). The output is a g-function that requires a heat load per unit of borehole 
length [W/m] to provide the borehole wall temperature.
# Arguments
    - `t`: Time vector (nₜ x 1) [s]
    - `k`: A (ix3) matrix of "i" geological layers and their thermal properties:
        - Column 1: Thickness [m]
        - Column 2: Thermal conductivity [W/mK]
        - Column 3: Volumetric specific heat [J/m³K]
        - Example: [20 k₁ C₁; 80 K₂ C₂] for 2 layers of 20m and 80m thicknesses
        - Note: For a single layer ground, the first value is equal to `H`.
    - `K`: A (jx2) matrix of "j" hydrogeological layers and their hydraulic properties:
        - Column 1: Thickness [m]
        - Column 2: Hydraulic conductivitiy [m/s]
        - Example: [20 K₁; 80 K₂] for 2 layers of 20m and 80m thicknesses
    - `V`: Circulating fluid flow rate ([m³/s])
    - `β`: Fluid bleed rate ratio (between 0.01 and 1) [-]
    - `H`: Borehole depth [m]
    - `rb`: Borehole radius [m]
    - `f`: Borehole (not pipe) friction factor from the thermal resistance (Rb_SCW()) [-]
# Output
    - `g`: A g-function corresponding to the borehole wall temperature of the SCW [°Cm/W]
# Reference
    - Nguyen, A., Jacques, L., & Pasquier, P. (2025). An easy-to-use analytical model for 
        standing column wells operating with bleed. Applied Thermal Engineering.
        https://doi.org/10.1016/j.applthermaleng.2024.124543
    - Jacques, L., Pasquier, P., Nguyen, A., & Beaudry, G. (2025). Improvement and experimental 
        validation of an analytical model for standing column wells operated with bleed. Applied 
        Thermal Engineering, 279, 127620.
"""
function βils(t, k::AbstractMatrix{<:Real}, K::AbstractMatrix{<:Real}, rb::Real, ro::Real, H::Real,
    V::Real, β::Real, f::Real)
    # 0. Check inputs
    if size(k, 2) != 3
        error("k must be a (ix3) matrix of geological layers and their thermal properties.")
    end
    if size(K, 2) != 2
        error("K must be a (jx2) matrix of hydrogeological layers and their hydraulic properties.")
    end
    # Replace β with a small value if it's 0
    β = β == 0 ? 1e-9 : β

    # 1. Effective thermal properties (arithmetic weighted average)
    ks = sum(k[:, 2] .* k[:, 1]) / sum(k[:, 1])
    Cs = sum(k[:, 3] .* k[:, 1]) / sum(k[:, 1])

    # 2. Induced convergent flow consideration
    H̃₁, _, K̃₁, _ = convergence_flow(K)      # Compute the effective hydraulic properties
    ΔHa = 0.22                              # Assumption: Drawdown caused by recirculation
    rₘ = 100                                # Assumption: Drawdown horizontal influence radius
    # Jacques et al. Eq. 35
    ΔH̃ = ΔHa - (f * H * V^2 / (4 * π^2 * 9.81 * (rb - ro) * (rb^2 - ro^2)^2)) # Succion drawdown
    ΔH̃ < 0 ? ΔH̃ = 0 : ΔH̃ = ΔH̃               # Avoid negative value of ΔH̃
    # Jacques et al. Eq. 34
    Ṽ = 2 * π * K̃₁* ΔH̃ * H̃₁ / (log(rₘ / rb)) # Rate of convergence flow caused by depressurization

    # 3. Layered heterogeneity - Effective length (if present)
    Vb = V * β                              # Bleed flow rate [m³/s]
    Kₕ = sum(K[:, 2] .* K[:, 1]) / H        # Arithmetic weighted average
    Kᵥ = H / sum(K[:, 1] ./ K[:, 2])        # Harmonic weighted average
    # Jacques et al. Eq. 31
    K̃ = sqrt(Kₕ * Kᵥ)                       # Effective hydraulic conductivity
    He = H * K̃ / Kₕ                         # Effective length only considering heterogeneity

    # 3. Compute the Peclet number and dimensionless time
    ta = @. (Ṽ / H̃₁ + Vb / He) * t / (π * rb^2) # Dimensionless time, Jacques et al. Eq. 33
    Pe = (Ṽ / H̃₁ + Vb / He) / (2 * π * (ks / Cs)) # Effective Peclet number, Jacques et al. Eq. 32

    # 4. Range of validity of the inputs based on scenarios in Nguyen et al. (2025)
    # values_validity(t, ks, Cs, rb, H, V, β, Vb, Pe)

    # 5. Compute the scaling function h (see Table 1 in Jacques et al.)
    coefs₁ = [-0.4478, 0.1288, 0.6458, 0.1483]
    coefs₂ = [0.7110, -0.7595, 0.7838, 0.5114]
    coefs₃ = [1.0011, 0.7339, 0.2383, -0.0100]
    ϵ = @. coefs₁ * Pe^coefs₂ + coefs₃      # Jacques et al. Eq. 7
    # Jacques et al. Eq. 4
    h = @. ϵ[1] * (1 - ϵ[2] * β) * (1 + erf((sqrt(Pe) / ϵ[3]) * ((1 - ta^ϵ[4]) / (ta^ϵ[4]))))
    
    # 6. Initial transfer function g_0 based on the SLI
    g₀ = ils(t, rb, ks, Cs)                 # Jacques et al. Eq. 3

    # 7. Combine initial and scaling function
    gb = convolution(g₀, h)      # Jacques et al. Eq. 8
    return gb
end

"""
    convergence_flow(K)

Evaluate the balance point "m" and the hydraulic conductivities over and under it in a SCW caused by
recirculation flow. This is used to compute the effective length in the β-ILS model that considers
induced convergent flow in a SCW, as defined by Jacques et al. (2025).
# Arguments
    - `K`: A (jx2) matrix of "j" hydrogeological layers and their hydraulic properties:
        - Column 1: Thickness [m]
        - Column 2: Hydraulic conductivitiy [m/s]
# Output
    - `H̃₁`: Depressurization length in the borehole [m]
    - `H̃₂`: Surpressurization length in the borehole [m]
    - `K̃₁`: Hydraulic conductivity above the balance point [m/s]
    - `K̃₂`: Hydraulic conductivity above the balance point [m/s]
# Reference
    - Jacques, L., Pasquier, P., Nguyen, A., & Beaudry, G. (2025). Improvement and experimental 
        validation of an analytical model for standing column wells operated with bleed. Applied 
        Thermal Engineering, 279, 127620.
"""
function convergence_flow(K::AbstractMatrix{<:Real})
    # 0. Constants and initialization
    nlayer = size(K, 1)                 # Number of hydrogeologic layers
    ΔH = sum(@view K[:, 1])             # Total zone thickness where the property are calculated
    mᵢ = ΔH / 2                         # Initial guess of balance point m (dimensionless)
    i = 1                               # Iterator
    iₘ = 100                            # Max number of iteration to convergence
    Δm = 1.0                            # Middle point convergence
    tol = 1e-2                          # Tolerance criteria before convergence in meters
    
    # Preallocation
    T = float(promote_type(eltype(K)))
    m, Kh1, Kh2, Kv1, Kv2 = zero(T), zero(T), zero(T), zero(T), zero(T), zero(T), zero(T)
    K̃₁, K̃₂ = zero(T), zero(T)

    # 1. Iteration loop
    while i ≤ iₘ && Δm ≥ tol
        # Depth of the middle point
        zm = ΔH / mᵢ # Weighted average depth of middle point

        # Initialize variables for accumulation
        H_up, H_low = zero(T), zero(T)
        ΣKh1, ΣKh2, ΣKv1, ΣKv2 = zero(T), zero(T), zero(T), zero(T)
        zbot, ztop = zero(T), zero(T)
        
        # Loop for each layer
        @inbounds for j in 1:nlayer
            thickness = K[j, 1]
            Kval      = K[j, 2]
            zbot = ztop + thickness
            # Upper zone
            if ztop < zm
                dz = min(zbot, zm) - ztop
                if dz > 0
                    H_up  += dz
                    ΣKh1 += dz * Kval
                    ΣKv1 += dz / Kval
                end
            end
            # Lower zone
            if zbot > zm
                dz = zbot - max(ztop, zm)
                if dz > 0
                    H_low  += dz
                    ΣKh2 += dz * Kval
                    ΣKv2 += dz / Kval
                end
            end
            ztop = zbot
        end

        # Effective conductivities
        Kh1 = ΣKh1 / H_up
        Kv1 = H_up / ΣKv1
        K̃₁ = sqrt(Kh1 * Kv1)            # Jacques et al. Eq. 31

        Kh2 = ΣKh2 / H_low
        Kv2 = H_low / ΣKv2
        K̃₂ = sqrt(Kh2 * Kv2)            # Jacques et al. Eq. 31

        # Update balance point
        zmp = Kh2 * ΔH / (Kh1 + Kh2)
        m   = ΔH / zmp

        # Convergence criteria
        Δm = abs(m - mᵢ)
        mᵢ = m
        i += 1
    end
    # Obtain effective length in depressurization and surpressurization zones
    H₁ = H / m                          # Depressurization lenth
    H₂ = H - H₁                         # Surpressurization length
    H̃₁ = H₁ * K̃₁ / Kh1                  # Jacques et al. Eq. 31
    H̃₂ = H₂ * K̃₂ / Kh2                  # Jacques et al. Eq. 31
    return H̃₁, H̃₂, K̃₁, K̃₂#, Kh1, Kv1, Kh2, Kv2
end

"""
Validate the values of the parameters used to generate the transfer function. See Table 1 in Nguyen 
et al. (2025).
"""
function values_validity(t, ks, Cs, rb, H, V, β, B, Pe)
    if any(t -> t < 0 || t > 50 * 365 * 24 * 3600, t)
        @warn "t is out of bounds"
    end
    if ks < 1 || ks > 5
        @warn "ks is out of bounds"
    end
    if Cs < 2.1e6 || Cs > 2.6e6
        @warn "Cs is out of bounds"
    end
    if rb < 0.076 || rb > 0.102
        @warn "rb is out of bounds"
    end
    if H < 100 || H > 500
        @warn "L is out of bounds"
    end
    if any(V -> V < 54 / 60000 || V > 600 / 60000, V)
        @warn "V is out of bounds"
    end
    if any(β -> β < 0.0096 || β > 1, β)
        @warn "β is out of bounds"
    end
    if any(B -> B < 0.72 / 60000 || B > 449 / 60000, B)
        @warn "B is out of bounds"
    end
    if Pe < 0.05 || Pe > 1
        @warn "Pe is out of bounds"
    end
end