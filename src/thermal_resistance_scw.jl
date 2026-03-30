"""
    borehole_resistance_scw(V, kp, rb, ro, ri, H, Hp, Tf, dT=2.0, ϵb=5e-4, ϵp=5e-4)

Computation of theoretical borehole thermal resistance in a SCW based on Jacques et al. 2025. This 
formulation assumes uniform flux along the borehole wall.
# Arguments
    - `V`: Circulating fluid flow rate [m³/s]
    - `kp`: Pipe thermal conductivity [W/mK]
    - `rb`: Borehole radius [m]
    - `ro`: Outer pipe radius [m]
    - `ri`: Inner pipe radius [m]
    - `H`: Borehole depth [m]
    - `Hp`: Pipe length [m]
    - `Tf`: Fluid temperature [°C]
    - `dT`: (Optional) Temperature difference at the heat exchanger (for fluid properties) [K]
        - Positive value: Ground heating -> the fluid is being cooled (Tb<Tf)
        - Negative value: Ground cooling -> the fluid is being heating (Tb>Tf)
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
    - Seol, H., Jeong, S., Cho, C., & You, K. (2008). Shear load transfer for rock-
        socketed drilled shafts based on borehole roughness and geological strength index (GSI). 
        International Journal of Rock Mechanics and Mining Sciences, 45(6), 848–861. 
        https://doi.org/10.1016/j.ijrmms.2007.09.008
    - Todorov, O., Alanne, K., Virtanen, M., & Kosonen, R. (2021). Different Approaches for 
        Evaluation and Modeling of the Effective Thermal Resistance of Groundwater-Filled Boreholes.
    - Lamarche, L. (2023). Fundamentals of Geothermal Heat Pump Systems: Design and Application. 
        Springer Nature Switzerland.
"""
function borehole_resistance_scw(V::Real, kp::Real, rb::Real, ro::Real, ri::Real, H::Real, Hp::Real,
    Tf::Real, dT::Real = 2.0, ϵb::Real=5e-4, ϵp::Real=5e-4)
    # Initialize parameters
    kf = water_k(Tf)                    # Water thermal conductivity at working fluid temperature
    ρf = water_ρ(Tf)                    # Water density at working fluid temperature
    cpf = water_cp(Tf)                  # Water specific heat at working fluid temperature
    μf = water_μ(Tf)                    # Water viscosity at working fluid temperature
    μf_dT = water_μ(Tf + dT)            # Water viscosity at working fluid temperature + dT
    np = 1                              # number of pipes (1 in a SCW)
    L = [Hp, H]                         # Lengths for (1) pipe and (2) borehole for [m]

    # Preallocation of properties for pipe and borehole
    rₑ, V̇, Re, Pr, f, Nu, h, Rf = zeros(2), zeros(2), zeros(2), zeros(2), zeros(2), zeros(3),
        zeros(2), zeros(2)
    
    # Equivalent radius
    rₑ[1] = ri                          # Equivalent radius in the pipe
    rₑ[2] = (rb - ro)                   # Equivalent radius in the borehole

    # Reynolds numbers
    V̇[1] = V / (π * ri^2)                 # Fluid velocity in the pipe [m/s]
    V̇[2] = V / ((π * rb^2) - (π * ro^2))  # Fluid velocity in the borehole [m/s]
    Re[1] = Reynold(V̇[1], rₑ[1], ρf, μf)  # Reynolds number in the pipe
    Re[2] = Reynold(V̇[2], rₑ[2], ρf, μf)  # Reynolds number in the borehole

    # Prandtl numbers
    Pr[1] = Prandtl(kf, cpf, μf)
    Pr[2] = Prandtl(kf, cpf, μf)

    # Friction factors using the Colebrook-White equation
    ϵ = [ϵp, (ϵb * rb + ϵp * ro) / (rb + ro)] # Equivalent roughness [m] Jacques et al. Eq. 47
    f[1] = friction_factor_Colebrook_White(Re[1], rₑ[1], ϵ[1]) # Friction factor in the pipe
    f[2] = friction_factor_Colebrook_White(Re[2], rₑ[2], ϵ[2]) # Friction factor in the borehole

    # Nusselt number
    Nu[1] = Nusselt(Re[1], Pr[1], rₑ[1], ϵp) # Nusselt number in the pipe
    Nu[2] = Nusselt(Re[2], Pr[2], rₑ[2], ϵ[2]) # Nusselt number in the borehole
    # Nu[2] = Nusselt_annulus(Re[2], Pr[2], rₑ[2], ϵb, ϵp) # Nusselt number in the borehole

    # Film coefficient (h)
    h[1] = Nu[1] * kf / (2 * rₑ[1])
    h[2] = Nu[2] * kf / (2 * rₑ[2])

    # Fluid resistance (Rf)
    Rf[1] = 1 / (2 * np * π * rₑ[1] * h[1])
    Rf[2] = 1 / (2 * np * π * rₑ[2] * h[2])

    # Nusselt number and fluid resistance, loop for regions i = 1 (pipe), i = 2 (borehole)
    #=
    # Preallocation
    f1, Nu1, h1, Rf1 = zeros(2), zeros(2), zeros(2), zeros(2)

    for i in 1:2
        if Re[i] < 2300                 # Laminar flow in pipe or borehole
            term1 = (Re[i] * Pr[i] * 2 * rₑ[i] / L[i])^0.333 * (μf / μf_dT)^0.14
            if i == 1 && term1 >= 2
                Nu1[i] = 4.36    # Lamarche Eq. 2.42
            elseif i == 1
                Nu1[i] = 1.86 * term1
            elseif i == 2
                f = 64 / Re[i]
                Nu1[i] = 3.66
            else
                error("Unexpected case in Nusselt number calculation. Laminar flow.")
            end
        elseif Re[i] > 10000 || i == 2  # Turbulent flow (Re>2300 in borehole, Todorov et al. 2021)
            if i == 1                   # Turbulent flow in pipe
                # Dittus-Boelter equation to solve for Nu
                # dT > 0 ? n = 0.3 : n = 0.4 # Define fluid heating or cooling
                # Nu1[i] = 0.023 * Re[i]^0.8 * Pr[i]^n # Jacques et al. Eq. A.44
                # Colebrook-White equation to solve for f in borehole
                f1[i] = friction_factor_Colebrook_White(Re[i], rₑ[i], ϵp)
                # Nusselt number in borehole (see Todorov et al. 2021)
                Nu1[i] = ((f1[i] / 8) * (Re[i] - 1000.0) * Pr[i]) / (1 + 12.7 * (f1[i] / 8)^0.5 * 
                    (Pr[i]^(2 / 3) - 1)) # Jacques et al. Eq. A.46
            elseif i == 2               # Turbulent flow in borehole
                # Colebrook-White equation to solve for f in borehole
                ϵ = (ϵb * rb + ϵp * ro) / (rb + ro) # Equivalent roughness [m] Jacques et al. Eq. 47
                f1[i] = friction_factor_Colebrook_White(Re[i], rₑ[i], ϵ)
                # Nusselt number in borehole (see Todorov et al. 2021)
                Nu1[i] = ((f1[i] / 8) * (Re[i] - 1000.0) * Pr[i]) / (1 + 12.7 * (f1[i] / 8)^0.5 * 
                    (Pr[i]^(2 / 3) - 1)) # Jacques et al. Eq. A.46
            else
                error("Unexpected case in Nusselt number calculation. Turbulent flow.")
            end
        else    # Transitional flow in pipe only (i = 1 and 2300 < Re < 10000)
            f1[i] = friction_factor_Colebrook_White(ϵp, rₑ[i], Re[i])
            Nu1[i] = ((f1[i] / 2) * (Re[i] - 1000.0) * Pr[i]) / (1 + 12.7 * (f1[i] / 2)^0.5 * 
                (Pr[i]^(2 / 3) - 1)) # Jacques et al. Eq. A.46
        end
        # Film coefficient (h) and fluid resistance (Rf)
        h1[i] = Nu1[i] * kf / (2 * rₑ[i])
        Rf1[i] = 1 / (2 * np * π * rₑ[i] * h1[i]) # Fluid thermal resistance
    end=#

    # Convective thermal resistances for the fluid in pipe and borehole
    RCi = 1 / (Nu[1] * kf * π)          # In pipe at ri, Jacques et al. Eq. A. 41
    RCo = 1 / (Nu[2] * kf * π) * rₑ[2] / ro # In borehole at ro, Jacques et al. Eq. 42
    RCb = 1 / (Nu[2] * kf * π) * rₑ[2] / rb # In borehole at rb, Jacques et al. Eq. 43

    # Pipe thermal resistance
    # Rp = log(ro / ri) / (2 * π * np * kp * Hp) # Jacques et al. Eq. A.40
    Rp = resistance_pipe(ro, ri, kp)

    # Total borehole thermal resistance in SCW
    Rv = H / (cpf * ρf * V)             # Fluid flow hydraulic resistance, Jacques et al. Eq. 18
    R12 = RCi + RCo + Rp                # Internal thermal resistance, Jacques et al. Eq. 17
    Rsb = RCb * (1 + (Rv^2 / (3 * RCb * R12))) # Subsurface borehole resist., Jacques et al. Eq. 15
    Rfb = Rsb - RCb                     # Fluid thermal resistance of the borehole in a SCW
    # Rfb2 = R12 - RCb                    # As implemented in the codes of Jacques et al. (2025)
    return Rfb, Rv, f[2]                # Only friction factor in borehole is relevant in β-ILS
end