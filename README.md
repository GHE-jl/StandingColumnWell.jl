# StandingColumnWell.jl

[![CI](https://github.com/GHE-jl/StandingColumnWell.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/GHE-jl/StandingColumnWell.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Documentation](https://github.com/GHE-jl/StandingColumnWell.jl/actions/workflows/Documenter.yml/badge.svg?branch=main)](https://GHE-jl.github.io/StandingColumnWell.jl)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPLv3-blue.svg)](LICENSE)

Thermal response of standing column wells (SCW), part of the
[GHE.jl](https://GHE-jl.github.io/GHE.jl-docs) ecosystem. It implements the β-ILS analytical model
of Nguyen et al. (2025), improved by Jacques et al. (2025), for a SCW operated with bleed, and
depends on [`GroundHeatExchanger.jl`](https://github.com/GHE-jl/GroundHeatExchanger.jl) for the
infinite line source, the temporal superposition and the fluid properties/correlations.

## Quick start

```julia
using StandingColumnWell

t = collect(3600.0:3600.0:(3600.0 * 24 * 365)) # time [s]
V, β, T₀ = 0.001633, 0.1, 10.0                 # flow rate [m³/s], bleed ratio [-], T [°C]
rb, ri, ro, H = 0.0855, 0.0247, 0.0302, 213.0   # borehole/pipe radii [m], saturated length [m]
kp, ϵb, ϵp = 0.4, 0.0, 0.0                      # pipe conductivity [W/mK], roughness [m]
ΔHₚ, rₚ = 0.22, 100.0                           # drawdown [m], influence radius [m]
k = [H 2.7 2.2e6]                               # thermal layers: thickness, ks, Cs
K = [H 1e-6]                                    # hydraulic layers: thickness, K [m/s]

# Outlet transfer function of the SCW for 1 W/m
g̃ = βils_outlet(t, V, β, H, rb, ro, ri, k, kp, K, T₀, ϵb, ϵp, ΔHₚ, rₚ)
```

## What it does

- **SCW thermal resistance** (`borehole_resistance_scw`): subsurface effective borehole thermal
  resistance `Rfb`, hydraulic resistance `Rv` and borehole friction factor `f`.
- **Induced convergent flow** (`convergence_flow`): balance point and effective hydraulic
  properties of a layered aquifer under the recirculation flow.
- **β-ILS g-function** (`βils`): borehole wall g-function of the SCW, from the ILS and a scaling
  function of the Péclet number and the bleed ratio `β`.
- **Outlet transfer function** (`βils_outlet`): outlet fluid temperature response of the SCW for
  1 W/m, to be used with `convolution` or `convolution_ns` from `GroundHeatExchanger.jl`.

See the [documentation](https://GHE-jl.github.io/StandingColumnWell.jl) for the full API reference.

## Scripts

Run from the package root with `julia --project=script/ script/<name>.jl`. First-time setup:
```
julia --project=script/ -e 'using Pkg; Pkg.instantiate()'
```

| Script | What it shows |
|---|---|
| `script_validation_GBeaudry.jl` | Borehole wall and outlet β-ILS g-functions of a SCW in a layered ground |
| `script_validation_LJacques.jl` | Non-stationary simulation of a SCW with variable pumping and bleed flow rates |

## Installation

This package is not yet registered in the Julia General registry. Install directly from GitHub:

```julia
using Pkg
Pkg.add(url = "https://github.com/GHE-jl/StandingColumnWell.jl")
```

## Dependencies

### Library

| Package | Used for |
|---|---|
| [GroundHeatExchanger.jl](https://github.com/GHE-jl/GroundHeatExchanger.jl) | Infinite line source, temporal superposition, fluid properties and convection/friction correlations |
| [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl) | Error function in the β-ILS scaling function |

### Scripts only

| Package | Used in |
|---|---|
| [CairoMakie.jl](https://github.com/MakieOrg/Makie.jl) | Plotting in all `script_*.jl` examples |

## References

- Nguyen, A., Jacques, L., & Pasquier, P. (2025). An easy-to-use analytical model for standing
  column wells operating with bleed. *Applied Thermal Engineering*.
  <https://doi.org/10.1016/j.applthermaleng.2024.124543>
- Jacques, L., Pasquier, P., Nguyen, A., & Beaudry, G. (2025). Improvement and experimental
  validation of an analytical model for standing column wells operated with bleed. *Applied
  Thermal Engineering*, 279, 127620.

## License

AGPL-3.0, see [LICENSE](LICENSE).
