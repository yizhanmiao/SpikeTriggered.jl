# SpikeTriggered.jl

[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/JuliaDiff/BlueStyle)
[![][docs-img]][docs-url]

A Julia package for spike-train based neural data analysis, biased to vision.

## Features

- Spike rasters and PSTHs (histogram and Gaussian-smoothed)
- Reverse correlation (spike-triggered average/covariance)
- Forward correlation
- Spatiotemporal receptive field (STRF) analysis: 2D Gaussian fitting, spatial footprints, temporal response functions
- Burst detection and classification
- Drifting grating tuning indices (OSI, DSI)
- Spectrum, correlograms, and reliability measures
- Waveform event detection and filtering (`SpikeTriggered.Waveforms`)

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/yizhanmiao/SpikeTriggered.jl.git")
```

or through REPL
```julia
julia> ] add https://github.com/yizhanmiao/SpikeTriggered.jl.git
```

## Quick Start

```julia
using SpikeTriggered

# Build a raster from spike times and event markers
raster = spike_raster(spike_times, markers; pre=0.5, post=1.0)

# Spike-triggered average
sta = spike_triggered_average(stimulus, spike_times, markers; n=10)
```

[docs-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-url]: https://yizhanmiao.github.io/SpikeTriggered.jl
