"""
Validation using the inputs provided by Gabrielle Beaudry.
"""

using CairoMakie

# includet("../../GroundHeatExchanger.jl/src/GroundHeatExchanger.jl")
# using .GroundHeatExchanger

includet("../src/StandingColumnWell.jl")
using .StandingColumnWell

# Operational parameters
t = 0.0:3600:3600*24*365    # Time vector (s)
Vp = 0.001633               # Pumping rate (m3/s)
Vb = 0.0/1000/60            # Bleed rate (m3/s)
T₀ = 10.0                   # Average fluid temperature for property evaluation (°C)

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
dT = 2.0                    # Temperature difference at the heat exchanger (K)

# Heterogeneous ground properties
k = [73.0 2.54 2.17e6 ; 71.0 2.74 2.27e6 ; H - 144.0 2.94 2.37e6] # Thermal properties
K = [7.0 1e-5; H - 7.0 1e-7] # Hydraulic properties
@assert sum(k[:, 1]) == H "The sum of the layer in k must equal the total depth of the SCW."
@assert sum(K[:, 1]) == H "The sum of the layer in K must equal the total depth of the SCW."

# Compute the outputs of the β-ILS model
 # Thermal resistance of the SCW
Rb, Rv, f = borehole_resistance_scw(Vp, kp, rb, ro, ri, H, Hp, T₀, dT, ϵb, ϵp)
 # Effective length of the SCW considering recirculation flow
H̃₁, H̃₂, K̃₁, K̃₂ = convergence_flow(K)
 # Borehole wall g-function
gb = βils(t, k, K, rb, ro, H, Vp, Vb, f)
# Outlet g-function
gout = βils_outlet(t, k, kp, K, rb, ro, ri, H, Hp, Vp, Vb, T₀, dT, ϵb, ϵp)

fig = Figure()
ax1 = Axis(fig[1, 1], xlabel="Time (s)", ylabel="g-function")
lines!(ax1, t, gb, label="Borehole wall g-function")
lines!(ax1, t, gout, label="Outlet g-function")
axislegend(ax1, position=:rc)
display(fig)