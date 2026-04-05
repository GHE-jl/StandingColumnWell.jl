using SpecialFunctions: erf

"""
This file also uses the infinite line source (ILS) model and the convolution product implemented
in the GroundHeatExchanger.jl
"""

"""
    βils_outlet(gb, Rfb, Rv, β)
    βils_outlet(t, V, β, H, rb, ro, k, K, Rfb, Rv, f, ΔHₚ, rₚ)
    βils_outlet(t, V, β, H, rb, ro, ri, k, kp, K, T₀, ϵb, ϵp, ΔHₚ, rₚ)

Compute the borehole outlet transfer function from the β-ILS model and the thermal resistance
of the SCW. The function computes the borehole wall (g-function) using the β-ILS model, then
adjusts it to obtain the outlet transfer function using the development from Jacques et al. (2025).
The obtained transfer function is for an impulse of 1 W/m.
# Arguments
    - `t`: Time vector (nₜ x 1) [s]
    - `V`: Circulating fluid flow rate [m³/s]
    - `β`: Fluid bleed rate ratio (between 0.01 and 1) [-]
    - `H`: Borehole depth [m]
    - `rb`: Borehole radius [m]
    - `ro`: Outer pipe radius [m]
    - `ri`: Inner pipe radius [m]
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
    - `T₀`: Initial ground temperature [°C] (used for water properties)
    - `ϵb`: Borehole roughness height [m]
    - `ϵp`: Pipe roughness height [m]
    - `ΔHₚ`: Drawdown caused by recirculation flow in the SCW [m]
    - `rₚ`: Horizontal influence radius of the recirculation flow in the SCW [m]
    - `gb`: A g-function corresponding to the borehole wall temperature of the SCW [°Cm/W]
    - `Rfb`: Subsurface effective borehole thermal resistance in a SCW [K·m/W]
    - `Rv`: Hydraulic resistance due to the fluid flow [K·m/W]
    - `f`: Friction factor of the borehole [-]
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
function βils_outlet(gb::AbstractArray{<:Real}, Rfb::Real, Rv::Real, β::Real)
    # Adjust the borehole wall g-function to obtain the outlet g-function for 1 W/m (g̃ᵢ ≥ 0)
    return @. max(gb - (1 - β) * (Rv / (2 - β) - Rfb), 0.0) # Jacques et al. Eq. 27
end
function βils_outlet(t, V::Real, β::Real, H::Real, rb::Real, ro::Real, k::AbstractMatrix{<:Real},
    K::AbstractMatrix{<:Real}, Rfb::Real, Rv::Real, f::Real, ΔHₚ::Real, rₚ::Real)
    # Compute the borehole wall g-function using the β-ILS model
    gb = βils(t, V, β, H, rb, ro, k, K, f, ΔHₚ, rₚ) # Jacques et al. Eq. 8

    return βils_outlet(gb, Rfb, Rv, β)
end
function βils_outlet(t, V::Real, β::Real, H::Real, rb::Real, ro::Real, ri::Real, 
    k::AbstractMatrix{<:Real}, kp::Real, K::AbstractMatrix{<:Real}, T₀::Real,
    ϵb::Real, ϵp::Real, ΔHₚ::Real, rₚ::Real)
    # Compute the thermal resistance of the SCW
    Rfb, Rv, f = borehole_resistance_scw(V, H, rb, ro, ri, kp, T₀, ϵb, ϵp) # Jacques et al. Eq. 15

    # Compute the borehole wall g-function using the β-ILS model
    gb = βils(t, V, β, H, rb, ro, k, K, f, ΔHₚ, rₚ) # Jacques et al. Eq. 8

    return βils_outlet(gb, Rfb, Rv, β)
end

"""
    βils(t, V, β, H, rb, ro, k, K, f, ΔHₚ, rₚ)
    βils(t, V, β, H, rb, ro, k, K, Re, ϵb, ΔHₚ, rₚ)
    βils(t, V, β, H, rb, ro, k, K, ϵb, ρ, μ, ΔHₚ, rₚ)

Compute the analytical model developped by Nguyen et al. (2025) for standing column
wells (SCW). The output is a g-function that requires a heat load per unit of borehole 
length [W/m] to provide the borehole wall temperature.
# Arguments
    - `t`: Time vector (nₜ x 1) [s]
    - `V`: Circulating fluid flow rate ([m³/s])
    - `β`: Fluid bleed rate ratio (between 0.01 and 1) [-]
    - `H`: Borehole depth [m]
    - `rb`: Borehole radius [m]
    - `ro`: Outer pipe radius [m]
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
    - `f`: Borehole (not pipe) friction factor from the thermal resistance (Rb_SCW()) [-]
    - `ΔHₚ`: Drawdown caused by starting the recirculation flow in the SCW [m]
    - `rₚ`: Horizontal influence radius of the recirculation flow in the SCW [m]
    - `Re`: Reynolds number of the flow in the borehole (rb - ro) [-]
    - `ϵb`: Borehole roughness height [m]
    - `ρ`: Fluid density [kg/m³]
    - `μ`: Fluid dynamic viscosity [Pa·s]
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
function βils(t, V::Real, β::Real, H::Real, rb::Real, ro::Real, k::AbstractMatrix{<:Real},
    K::AbstractMatrix{<:Real}, f::Real, ΔHₚ::Real, rₚ::Real)
    # Check inputs
    if size(k, 2) != 3
        error("k must be a (ix3) matrix of geological layers and their thermal properties.")
    end
    if size(K, 2) != 2
        error("K must be a (jx2) matrix of hydrogeological layers and their hydraulic properties.")
    end
    # Replace β with a small value if it's 0
    β = β == 0 ? 1e-9 : β

    # Effective thermal properties (arithmetic weighted average)
    ks = sum(k[:, 2] .* k[:, 1]) / sum(k[:, 1])
    Cs = sum(k[:, 3] .* k[:, 1]) / sum(k[:, 1])

    # Induced convergent flow consideration
    H̃₁, _, K̃₁, _ = convergence_flow(K)      # Compute the effective hydraulic properties
    # Jacques et al. Eq. 35
    ΔH̃ = ΔHₚ - (f * H * V^2 / (4 * π^2 * 9.81 * (rb - ro) * (rb^2 - ro^2)^2)) # Succion drawdown
    ΔH̃ < 0 ? ΔH̃ = 0 : ΔH̃ = ΔH̃               # Avoid negative value of ΔH̃
    # Jacques et al. Eq. 34
    Ṽ = 2 * π * K̃₁* ΔH̃ * H̃₁ / (log(rₚ / rb)) # Rate of convergence flow caused by depressurization

    # Layered heterogeneity - Effective length (if present)
    Vb = V * β                              # Bleed flow rate [m³/s]
    Kₕ = sum(K[:, 2] .* K[:, 1]) / H        # Arithmetic weighted average
    Kᵥ = H / sum(K[:, 1] ./ K[:, 2])        # Harmonic weighted average
    # Jacques et al. Eq. 31
    K̃ = sqrt(Kₕ * Kᵥ)                       # Effective hydraulic conductivity
    He = H * K̃ / Kₕ                         # Effective length only considering heterogeneity

    # Compute the Peclet number and dimensionless time
    ta = @. (Ṽ / H̃₁ + Vb / He) * t / (π * rb^2) # Dimensionless time, Jacques et al. Eq. 33
    Pe = (Ṽ / H̃₁ + Vb / He) / (2 * π * (ks / Cs)) # Effective Peclet number, Jacques et al. Eq. 32

    # Compute the scaling function h (see Table 1 in Jacques et al.)
    coefs₁ = [-0.4478, 0.1288, 0.6458, 0.1483]
    coefs₂ = [0.7110, -0.7595, 0.7838, 0.5114]
    coefs₃ = [1.0011, 0.7339, 0.2383, -0.0100]
    ϵ = @. coefs₁ * Pe^coefs₂ + coefs₃      # Jacques et al. Eq. 7
    # Jacques et al. Eq. 4
    h = @. ϵ[1] * (1 - ϵ[2] * β) * (1 + erf((sqrt(Pe) / ϵ[3]) * ((1 - ta^ϵ[4]) / (ta^ϵ[4]))))
    
    # Initial transfer function g_0 based on the SLI
    g₀ = ils(t, rb, ks, Cs)                 # Jacques et al. Eq. 3

    # Combine initial and scaling function
    gb = convolution(g₀, h)                 # Jacques et al. Eq. 8
    return gb
end

function βils(t, V::Real, β::Real, H::Real, rb::Real, ro::Real, k::AbstractMatrix{<:Real},
    K::AbstractMatrix{<:Real}, Re::Real, ϵb::Real, ΔHₚ::Real, rₚ::Real)
    # Friction factor using the Colebrook-White equation
    f = friction_factor_Colebrook_White(Re, rb - ro, ϵb)

    return βils(t, V, β, H, rb, ro, k, K, f, ΔHₚ, rₚ)
end
function βils(t, V::Real, β::Real, H::Real, rb::Real, ro::Real, k::AbstractMatrix{<:Real},
    K::AbstractMatrix{<:Real}, ϵb::Real, ρ::Real, μ::Real, ΔHₚ::Real, rₚ::Real)
    # Reynold calculation
    Re = Reynold(V / ((π * rb^2) - (π * ro^2)), rb - ro, ρ, μ)
    # Friction factor using the Colebrook-White equation
    f = friction_factor_Colebrook_White(Re, rb - ro, ϵb)
    return βils(t, V, β, H, rb, ro, k, K, f, ΔHₚ, rₚ)
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
    - `K̃₂`: Hydraulic conductivity below the balance point [m/s]
# Reference
    - Jacques, L., Pasquier, P., Nguyen, A., & Beaudry, G. (2025). Improvement and experimental 
        validation of an analytical model for standing column wells operated with bleed. Applied 
        Thermal Engineering, 279, 127620.
"""
function convergence_flow(K::AbstractMatrix{<:Real})
    # Constants and initialization
    nlayer = size(K, 1)                 # Number of hydrogeologic layers
    ΔH = sum(@view K[:, 1])             # Total zone thickness where the property are calculated
    i = 1                               # Iterator
    iₘ = 100                            # Max number of iteration to convergence
    Δm = 1.0                            # Middle point convergence
    tol = 1e-2                          # Tolerance criteria before convergence in meters
    
    # Preallocation (use Float64 for intermediate calculations)
    T = Float64
    K̃₁, K̃₂, Kh1, Kh2 = zero(T), zero(T), zero(T), zero(T)


    # Iteration loop: iterate on the middle depth `zm`
    zm = ΔH / 2.0
    for i in 1:iₘ

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

        # Update balance point (depth) and check convergence
        zmp = Kh2 * ΔH / (Kh1 + Kh2)
        Δm = abs(zmp - zm)
        zm = zmp
        if Δm < tol
            break
        end
    end

    # Obtain effective length in depressurization and surpressurization zones
    H₁ = zm                           # Depressurization length (upper zone)
    H₂ = ΔH - H₁                      # Surpressurization length (lower zone)
    H̃₁ = H₁ * K̃₁ / Kh1               # Jacques et al. Eq. 31
    H̃₂ = H₂ * K̃₂ / Kh2               # Jacques et al. Eq. 31
    return H̃₁, H̃₂, K̃₁, K̃₂
end