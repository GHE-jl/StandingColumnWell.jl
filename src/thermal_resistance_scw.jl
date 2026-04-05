"""
    borehole_resistance_scw(V, H, rb, ro, ri, kp, T₀, ϵb, ϵp)

Computation of theoretical borehole thermal resistance in a SCW based on Jacques et al. 2025. This 
formulation assumes uniform flux along the borehole wall.
# Arguments
    - `V`: Circulating fluid flow rate [m³/s]
    - `H`: Borehole depth [m]
    - `rb`: Borehole radius [m]
    - `ro`: Outer pipe radius [m]
    - `ri`: Inner pipe radius [m]
    - `kp`: Pipe thermal conductivity [W/mK]
    - `T₀`: Fluid temperature [°C]
    - `ϵb`: (Optional) Borehole roughness height [m]
    - `ϵp`: (Optional) Pipe roughness height [m]
# Output
    - `Rfb`: Subsurface effective borehole thermal resistance in a SCW [K·m/W]
    - `Rv`: Hydraulic resistance due to the fluid flow [K·m/W]
    - `f`: Friction factor of the borehole [-]
# Reference
    - Jacques, L., Pasquier, P., Nguyen, A., & Beaudry, G. (2025). Improvement and experimental 
        validation of an analytical model for standing column wells operated with bleed. Applied 
        Thermal Engineering, 279, 127620.
    - Todorov, O., Alanne, K., Virtanen, M., & Kosonen, R. (2021). Different Approaches for 
        Evaluation and Modeling of the Effective Thermal Resistance of Groundwater-Filled Boreholes.
    - Lamarche, L. (2023). Fundamentals of Geothermal Heat Pump Systems: Design and Application. 
        Springer Nature Switzerland.
"""
function borehole_resistance_scw(V::Real, H::Real, rb::Real, ro::Real, ri::Real, kp::Real, T₀::Real,
    ϵb::Real, ϵp::Real)
    # Initialize parameters
    kf = water_k(T₀)                    # Water thermal conductivity at working fluid temperature
    ρf = water_ρ(T₀)                    # Water density at working fluid temperature
    cpf = water_cp(T₀)                  # Water specific heat at working fluid temperature
    μf = water_μ(T₀)                    # Water viscosity at working fluid temperature
    ϵ = (ϵb * rb + ϵp * ro) / (rb + ro) # Equivalent roughness [m] Jacques et al. Eq. 47

    # Preallocation of properties for pipe and borehole
    rₑ, V̇, Re, f, Nu, h, Rf = zeros(2), zeros(2), zeros(2), zeros(2), zeros(3),
        zeros(2), zeros(2)
    
    # Equivalent radius
    rₑ[1] = ri                          # Equivalent radius in the pipe
    rₑ[2] = (rb - ro)                   # Equivalent radius in the borehole

    # Reynolds numbers
    V̇[1] = V / (π * ri^2)               # Fluid velocity in the pipe [m/s]
    V̇[2] = V / ((π * rb^2) - (π * ro^2)) # Fluid velocity in the borehole [m/s]
    Re[1] = Reynold(V̇[1], rₑ[1], ρf, μf) # Reynolds number in the pipe
    Re[2] = Reynold(V̇[2], rₑ[2], ρf, μf) # Reynolds number in the borehole

    # Prandtl numbers (same as in the pipe since same fluid properties)
    Pr = Prandtl(kf, cpf, μf)

    # Friction factors using the Colebrook-White equation
    # f[1] = friction_factor_Colebrook_White(Re[1], rₑ[1], ϵp) # Friction factor in the pipe
    f[2] = friction_factor_Colebrook_White(Re[2], rₑ[2], ϵ) # Friction factor in the borehole
    
    # Nusselt number
    Nu[1] = Nusselt(Re[1], Pr, rₑ[1], ϵp) # Nusselt number in the pipe
    Nu[2] = Nusselt(Re[2], Pr, rₑ[2], ϵ) # Nusselt number in the borehole
    # Nu[2] = Nusselt_annulus(Re[2], Pr, rₑ[2], ϵb, ϵp) # Nusselt number in the borehole

    # Film coefficient (h)
    # h[1] = Nu[1] * kf / (2 * rₑ[1])
    # h[2] = Nu[2] * kf / (2 * rₑ[2])

    # Fluid resistance (Rf)
    # Rf[1] = 1 / (2 * π * rₑ[1] * h[1])
    # Rf[2] = 1 / (2 * π * rₑ[2] * h[2])

    # Convective thermal resistances for the fluid in pipe and borehole
    RCi = 1 / (Nu[1] * kf * π)          # In pipe at ri, Jacques et al. Eq. A. 41
    RCo = 1 / (Nu[2] * kf * π) * rₑ[2] / ro # In borehole at ro, Jacques et al. Eq. 42
    RCb = 1 / (Nu[2] * kf * π) * rₑ[2] / rb # In borehole at rb, Jacques et al. Eq. 43

    # Pipe thermal resistance
    Rp = resistance_pipe(ro, ri, kp)

    # Total borehole thermal resistance in SCW
    Rv = H / (cpf * ρf * V)             # Fluid flow hydraulic resistance, Jacques et al. Eq. 18
    R12 = RCi + RCo + Rp                # Internal thermal resistance, Jacques et al. Eq. 17
    Rsb = RCb * (1 + (Rv^2 / (3 * RCb * R12))) # Subsurface borehole resist., Jacques et al. Eq. 15
    Rfb = Rsb - RCb                     # Fluid thermal resistance of the borehole in a SCW
    return Rfb, Rv, f[2]                # Only friction factor in borehole is relevant in β-ILS
end