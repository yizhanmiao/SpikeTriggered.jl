# SpikeTriggered.jl Documentation

!!! warning
    This documentation is still work in progress.

The [SpikeTriggered](https://github.com/yizhanmiao/SpikeTriggered.jl) package provides a collection of functions for analysis related to spike trains. Including:

- converting spike trains into rasters.
- getting PSTH from spike trains, by using histogram or (Gaussian) filter.
- estimating RFs using forward correlation and reserve correlation (e.g., spike-triggered averaging).
- thresholding waveforms to extract spike events (for simple online estimation).
- etc.

To keep things simple, PSTH and spike trains are represented as `AbstractVector{T}`,
rasters are `AbstractVector{<:AbstractVector{T}}`, and other statistics are using
Julia's built-in types whenever possible.

Functions for spike train statistics are in the `Stats` submodule and functions for waveform processing are in the `Waveforms` submodule.
Functions in the `Stimulus` submodule are only for internal uses in HirschLab.

PSTHs and rasters are usually the foundation of all the analysis; so please double check
your marker times are correct (e.g., there might be some jitters during stimulus generation,
or there might be an offset between the marker and spike times between the stimulus generator and data collection).
