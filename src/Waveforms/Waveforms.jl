module Waveforms
using GSL: GSL
using DSP: conv, Lowpass, Bandpass, Highpass, Butterworth, filtfilt, digitalfilter
using Statistics: median

include("events.jl")
include("filter.jl")

end
