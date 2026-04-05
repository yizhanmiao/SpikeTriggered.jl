# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Guidelines

- Do not re-read files you have already read in this session.
- Correctness and clarity first. Only optimise for performance when explicitly asked or when the issue is egregious.
- Explicit `return` at every exit point of a `function ... end` definition.
- Docstrings on every exported and public symbol using `@doc raw"""..."""` (not plain `"""`). Cover: what, args, return value. No docstrings required on _-prefixed internals.
- Only "why" comments — never describe what the code does; the code should be self-evident.
- Functional style over imperative loops. Prefer `map`, `filter`, `reduce`, broadcasting (`.`), and `|>` pipelines; use `for` only when clearly more readable.
- Descriptive error handling. Always `throw(...)` with context; never silently discard errors or return `nothing` on failure.
- No creative additions. Implement exactly what was asked; suggest extras in comments.
- Prefer existing files over new files. Only create a new file for a distinct logical component.
- Correctness check: logic errors, Julia gotchas (1-based indexing, column-major, copy vs view, scoping).
- Performance check (only when relevant): type instabilities, eachindex over 1:length, @code_warntype, BenchmarkTools.jl, pre-allocation, views, sizehint!, cautious @inbounds.
- Structure your review as: a brief summary, then a numbered list of specific issues (category, location, explanation, fix), then the corrected code.

## About

SpikeTriggered.jl is a Julia package for spike-train based neural data analysis, including spike statistics, reverse/forward correlation, spatiotemporal receptive field (STRF) analysis, and waveform processing. It follows the [Blue code style](https://github.com/JuliaDiff/BlueStyle). Requires Julia 1.11+.

## Architecture

### Core type aliases (defined in `src/SpikeTriggered.jl`)
- `SpikeRaster{T}` — `AbstractVector{<:AbstractVector{T}}` (collection of spike trains across trials)
- `AbstractSpikeTrain{T}` — `AbstractVector{T}` (single trial spike times)
- `AbstractPSTH{T}` — `AbstractVecOrMat{T}` (timepoints, or timepoints × repeats)
- `AbstractStimulus{T}` — `AbstractMatrix{T}` (dimensions × timepoints)
- `STRF{T}` — `Array{T,3}` (Height × Width × Time), defined in `src/strf/strf.jl`

These type aliases are used throughout all analysis functions for dispatch.

### Module organization
- **Top-level `src/`** — spike-train analysis functions (statistics, raster, PSTH, burst detection, spectrum, correlograms, reliability, drifting grating tuning, spike-triggered average/covariance)
- **`src/strf/`** — spatiotemporal receptive field analysis: 2D Gaussian fitting, spatial footprint extraction, pixel maps, temporal response functions, bright/dark polarity, ellipse overlap (Schiller method)
- **`Waveforms` submodule** (`src/Waveforms/`) — event detection and signal filtering (currently not re-exported; access via `SpikeTriggered.Waveforms`)

### Test structure
Tests in `test/runtests.jl` are organized into testsets: "receptive fields" and "spike statistics". Test files are included directly rather than being separate modules.

## Common Commands
```bash
# Run all tests
make test

# Run tests with coverage
make coverage-report

# Clean coverage artifacts
make clean

# start julia REPL with custom environment variables
pixi run julia
```
