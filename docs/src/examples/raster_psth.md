```@meta
EditURL = "raster_psth.jl"
```

# Raster and PSTH

Rasters and [peristimulus time histograms](https://en.wikipedia.org/wiki/Peristimulus_time_histogram)(PSTHs)
are very useful for inspecting neuronal activities to periodic stimuli or repeated trials.
Patterns of firing can be easily seen and further analyized
by aligning spike times to onsets of the stimulus.

Spike trains can be represented as `Vector{T}` in Julia.
Mathematically, it can be represented as series of [Dirac delta function](https://en.wikipedia.org/wiki/Dirac_delta_function):

```math
    \operatorname{s}(t) = \sum_i \delta(t - t'_i)
```
where $t'_i$ is the ith spike in the spike train vector.

!!! warning "Important"
    Spike trains sholud be monotonically increasing;
    and most of the analysis assume that values are sorted.

!!! tip
    in mose cases, `Float32` should be good enough for spike times;
    for a 20min recording, the epsilon for `Float32` is around 0.1 ms.
    You can use `Float64` if you have longer recordings or need higher precision.

````@example raster_psth
using SpikeTriggered # hide
````

## Getting Rasters
Rasters are created by calling [`spike_raster`](@ref) as `Vector{Vector{T}}`.
By default, each element of rasters are spike times relative the *onset* time of each trial.

!!! note
    So far, `spike_raster` assumes that trials are equal duration


!!! note
    You can also split rasters into different subgroups (e.g., directions of your drifting grating),
    by using [`spike_raster_groupby`](@ref).

````@example raster_psth
spike_train = sort(rand(1000) .* 50 .- 0.5) # make some fake spikes
marker_onset = 0:10:40 # fake trial onset time
spike_raster(spike_train, marker_onset;
    head=0.5, # include 0.5 seconds before trial onset
    duration=2.0, # each trial lasts for 2 seconds
    tail=0.5, # include 0.5 seconds after each trial
)
````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

