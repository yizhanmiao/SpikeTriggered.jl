module SpikeTriggered

import Reexport: @reexport

include("Waveforms/Waveforms.jl") # ==> simple event detection and removal
include("Stats/root.jl") # ==> spike related statistics: sta, stc, burst, power spectrum, etc.
# include("PSTH/PSTH.jl") # ==> binned spikes

@reexport using .Stats

end # module
