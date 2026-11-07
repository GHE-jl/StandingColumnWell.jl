```@meta
CurrentModule = StandingColumnWell
```

# StandingColumnWell.jl

*Thermal response of standing column wells (SCW), part of the GHE.jl ecosystem.*

The package implements the β-ILS analytical model of Nguyen et al. (2025), improved by Jacques et
al. (2025), for a standing column well operated with bleed. It depends on
[`GroundHeatExchanger.jl`](https://github.com/GHE-jl/GroundHeatExchanger.jl) for the infinite line
source, the temporal superposition and the fluid properties/correlations.

## What it does

- **SCW thermal resistance** ([`borehole_resistance_scw`](@ref)) — subsurface effective borehole
  thermal resistance `Rfb`, hydraulic resistance `Rv` and borehole friction factor `f` of a SCW
  (Jacques et al. 2025, Eqs. 15-18).
- **Induced convergent flow** ([`convergence_flow`](@ref)) — balance point and effective hydraulic
  properties of a layered aquifer under the recirculation flow (Jacques et al. 2025, Eq. 31).
- **β-ILS g-function** ([`βils`](@ref)) — borehole wall g-function of the SCW, from the ILS and a
  scaling function of the Péclet number and the bleed ratio `β` (Jacques et al. 2025, Eq. 8).
- **Outlet transfer function** ([`βils_outlet`](@ref)) — outlet fluid temperature response of the
  SCW for 1 W/m (Jacques et al. 2025, Eq. 27), to be used with `convolution` or `convolution_ns`
  from `GroundHeatExchanger.jl`.

## References

- Nguyen, A., Jacques, L., & Pasquier, P. (2025). An easy-to-use analytical model for standing
  column wells operating with bleed. *Applied Thermal Engineering*.
  <https://doi.org/10.1016/j.applthermaleng.2024.124543>
- Jacques, L., Pasquier, P., Nguyen, A., & Beaudry, G. (2025). Improvement and experimental
  validation of an analytical model for standing column wells operated with bleed. *Applied
  Thermal Engineering*, 279, 127620.
