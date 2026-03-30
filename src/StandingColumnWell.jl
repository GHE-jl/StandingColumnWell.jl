module StandingColumnWell

using Revise

include("../../GroundHeatExchanger.jl/src/GroundHeatExchanger.jl")    # Load the GroundHeatExchanger
using .GroundHeatExchanger

includet("beta_infinite_line_source.jl")
includet("thermal_resistance_scw.jl")

export βils_outlet, βils
export borehole_resistance_scw

end
