"""
Script to test the implementation of the β-ILS model for a standing column well (SCW) with
heterogeneous ground properties.
"""

using Pkg; Pkg.activate(joinpath(@__DIR__, "..", "script"))

using StandingColumnWell
using GroundHeatExchanger: fluid_property, Reynolds
using CairoMakie

# Operational parameters
t = 0.0:3600:3600*24*365    # Time vector (s)
V = 0.001633                # Pumping rate (m3/s)
B = 0.0/1000/60             # Bleed rate (m3/s)
β = B / V                   # Bleed ratio [-]
T₀ = 10.0                   # Average fluid temperature for property evaluation (°C)
ΔHₚ = 0.22                  # Observed drawdown due to pumping [m]
rₚ = 100.0                  # Influence radius of the borehole [m]
_, _, ρ, μ = fluid_property(T₀, :water) # Water density [kg/m³], viscosity [Pa·s]

# Geometric parameters
rb = 0.0855                 # SCW radius (m)
z = 215.0                   # SCW depth (m)
zw = 2.0                    # Static water level depth (m)
ri = 0.0247                 # Inner radius of the injection pipe (m)
ro = 0.0302                 # Outer radius of the injection pipe (m)
H = z - zw                  # Length of the SCW in the saturated zone (m)
Hp = H - 10.0               # Pipe length (m)

# Properties for resistance calculations
kp = 0.4                    # Thermal conductivity of the pipe (W/mK)
ϵb = 0.0                    # Roughness height at the well wall (m)
ϵp = 0.0                    # Roughness height at the pipe wall (m)

# Heterogeneous ground properties
k = [73.0 2.54 2.17e6 ; 71.0 2.74 2.27e6 ; H - 144.0 2.94 2.37e6] # Thermal properties
K = [7.0 1e-5; H - 7.0 1e-7] # Hydraulic properties

# Properties
Re = Reynolds(V / ((π * rb^2) - (π * ro^2)), rb - ro, ρ, μ)
Rfb, Rv, f = borehole_resistance_scw(V, H, rb, ro, ri, kp, T₀, ϵb, ϵp)

# Methods for the β-ILS model
gb1 = βils(t, V, β, H, rb, ro, k, K, f, ΔHₚ, rₚ)
gb2 = βils(t, V, β, H, rb, ro, k, K, Re, ϵb, ΔHₚ, rₚ)
gb3 = βils(t, V, β, H, rb, ro, k, K, ϵb, ρ, μ, ΔHₚ, rₚ)

@assert gb1 ≈ gb2 ≈ gb3 "The different methods for the β-ILS model should yield the same 
    g-function values."

# Method for the β-ILS outlet g-function

g1 = βils_outlet(gb1, Rfb, Rv, β)
g2 = βils_outlet(t, V, β, H, rb, ro, k, K, Rfb, Rv, f, ΔHₚ, rₚ)
g3 = βils_outlet(t, V, β, H, rb, ro, ri, k, kp, K, T₀, ϵb, ϵp, ΔHₚ, rₚ)

@assert g1 ≈ g2 ≈ g3 "The different methods for the β-ILS outlet g-function should yield the same 
    g-function values."

fig = Figure()
ax1 = Axis(fig[1, 1], xlabel="Time (h)", ylabel="β-ILS")
lines!(ax1, t / 3600, gb1, color=:blue, linestyle=:solid, label="β-ILS (1)")
lines!(ax1, t / 3600, gb2, color=:red, linestyle=:dash, label="β-ILS (2)")
lines!(ax1, t / 3600, gb3, color=:green, linestyle=:dot, label="β-ILS (3)")
axislegend(ax1, position=:rc)

ax2 = Axis(fig[2, 1], xlabel="Time (h)", ylabel="Outlet β-ILS")
lines!(ax2, t / 3600, g1, color=:blue, linestyle=:solid, label="β-ILS outlet (1)")
lines!(ax2, t / 3600, g2, color=:red, linestyle=:dash, label="β-ILS outlet (2)")
lines!(ax2, t / 3600, g3, color=:green, linestyle=:dot, label="β-ILS outlet (3)")
axislegend(ax2, position=:rc)
display(fig)