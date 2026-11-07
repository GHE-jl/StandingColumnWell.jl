# Generation of SCW thermal resistance values for different depth
# Code to help in the CMM project

using Pkg; Pkg.activate(joinpath(@__DIR__, "..", "script"))

using StandingColumnWell

# Define variables
V = 100.0 / 1000 / 60               # Pumping rate (m3/s)
H = [200.0, 300.0, 500.0, 700.0]    # SCW depths (m)
rb = 0.0855                         # SCW radius (m)
ro = 0.0302                         # Outer radius of the injection pipe (m)
ri = 0.0247                         # Inner radius of the injection pipe (m)
kp = 0.4                            # Thermal conductivity of the pipe (W/mK)
T₀ = 10.0                           # Average fluid temperature for property evaluation (°C)
ϵb = 5e-3                           # Roughness height at the well wall (m)
ϵp = 0.0                            # Roughness height at the pipe wall (m)

# Compute the SCW thermal resistance for different depths
Rb = zeros(length(H), 1)  # Preallocate array for results
for (i, h) in enumerate(H)
    tmp = borehole_resistance_scw(V, h, rb, ro, ri, kp, T₀, ϵb, ϵp)
    Rb[i] = tmp[1]  # Extract the subsurface borehole resistance (Rfb)
    println("Depth: $h m, Borehole Resistance: $(round(Rb[i], digits=4)) K·m/W")
end