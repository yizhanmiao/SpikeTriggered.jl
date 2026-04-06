# Raster and PSTH

Raster plot and [peristimulus time histogram](https://en.wikipedia.org/wiki/Peristimulus_time_histogram)(PSTH) are very useful for neural activities to periodic stimuli or repeated trials. Patterns of firing can be easily seen and further analyized when spikes are aligned with the onsets of the stimulus.

In most cases, spike trains are stored as `Vector{T}` where `T <: Real`.
Mathematically, it can be represented as series of [Dirac delta function](https://en.wikipedia.org/wiki/Dirac_delta_function):

```math
    \operatorname{s}(t) = \sum_i \delta(t - t'_i)
```

where $t'_i$ is the spike times in the spike train vector.

## Getting Rasters

Rasters will be created as `Vector{Vector{T}}` using [`spike_raster`](@ref), where each subvector will be spike times relative to the **stimulus onset**. Besides the spike train, you should at least have the onset time of the stimulus (`marker::AbstractVector`), and specify the length of time before the onset (`head`) and `duration` of the stimulus and time after the stimulus (`tail`).

```julia
using SpikeTriggered
spike_train = rand(1000) * 50
marker_onset = 0:10:40
raster = spike_raster(spike_train, marker_onset; head=0.5, duration=1.0, tail=0.5)
...
```

## Getting PSTHs

Histogram style PSTHs can be generated with [`spike_histogram`](@ref), and `FHist` is used under the hook.
In order to have finer control and a nonlinear timescale
(e.g., when you want to group each trial into 0-100ms and 100ms-500ms as triggering phase and steady phase), `edges` of the histogram is taken. So, don't forget to append the last right edge.

```julia
edges = 0:0.01:50.0
psth = spike_histogram(spike_or_raster, edges)
```

!!! note
    If you need PSTHs for each single trial, you can organize the edges into 2d Matrix
    as `[edges x nTrials]`, which results into 2d raster as `[n x nTrials]`;
    an easy way to generate such matrix:
    ```julia
    projection = 0:0.001:2
    trial_onset = [0,5,10,15,20]
    edge_2d = projection .+ trial_onset' # 2d edge matrix
    psth_2d = spike_histogram(spike_train, edge_2d)
    ```

!!! note
    the PSTH from raster has an additional kwargs `norm=true` which will normalized
    the PSTH by the number of trials.

## Getting smoothed PSTHs

Sometimes, it's arbitary to find a good bin size, or you just want to smooth the histogram;
you can use a kernel function to generate a smooth PSTH using [`spike_histogram_smoothed`](@ref).
Essentially, it is a convolution over time steps provided by `projection::AbstractVector`.

```math
\operatorname{PSTH}(t) = (\operatorname{h} * \operatorname{s})(t) = \sum_i \delta(t - t'_i) h(t - t'_i),\ t \in \text{projection}
```

where `h` is the kernel function, $t'_i$ is the spike time.
By default, it uses a Gaussian kernel function ([`SpikeTriggered.gaussian_kernel`](@ref)).

```julia
# Gaussian kernel with default sigma = 0.005
spike_histogram_smoothed(spike_or_raster, projection)

# change the Gaussian kernel width
spike_histogram_smoothed(spike_or_raster, projection; σ=0.010)
```

!!! note
    The kernel function shoudl be defined as `(dt; kwargs...) -> ...`:
    ```julia
    my_kernel(dt; thresh) = dt < thresh ? 0 : 1
    psth = spike_histogram_smoothed(spike_or_raster, projection, my_kernel; thresh=0.5)
    ```

!!! note
    There is also a `norm=true` flag, to normalize the PSTH by its maximum value,
    so the PSTH will be between `[0, 1]`.
