# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Aqua.jl test suite for package quality checks (#31)
- Documentation build GitHub Actions workflow
- `Makefile` targets for running tests and launching the Julia REPL with the correct project environment

### Changed
- Functions renamed to a consistent snake\_case style; old names deprecated via `@deprecate` (#32)
- Fixed REPL prompt in README (`>>>` → `julia>`)

### Removed
- Deprecated `cdf` call style on matrix input (#30)

## [1.3.1] - 2026-03-19

### Changed
- Bumped Optim.jl compatibility to v2

## [1.3.0] - 2026-03-16

### Added
- Bright/dark polarity score for STRF analysis
- Biphasic temporal response function (TRF) segmentation

### Changed
- `get_footprint_map` now returns the spatially collapsed map by default
- Replaced GSL.jl dependency with Interpolations.jl
- Moved footprint analysis into the `strf/` submodule

### Fixed
- PSTH computation from spike rasters
- Smoothed PSTH computation from spike rasters

### Removed
- Packages no longer required after GSL → Interpolations migration

## [1.2.1] - 2026-01-28

### Changed
- Renamed `get_filtered_pixel_map` to `get_pixel_cloud`
- `spike_histogram` now accepts empty vectors and returns a zero array
- Raster construction now accepts trials with no spikes

### Fixed
- Burst split no longer errors when the input contains no bursts

## [1.2.0] - 2025-10-25

### Added
- Pixel threshold-based spatial masking (`get_filtered_pixel_map`) following the method in Fiona's paper (#21)

### Fixed
- Interpolations.jl call style updated for upstream API change (#23)

## [1.1.1] - 2025-10-16

### Fixed
- `floodfill!` was not exported/loaded (#19)

## [1.1.0] - 2025-10-13

### Added
- `floodfill` / `floodfill!` for connected-component filling on pixel maps
- `strf_vstack` for stacking STRFs along the time axis

### Changed
- Removed overly restrictive type annotations from STRF functions
- `benjamini_hochberg_constant` used by default for FDR correction

### Removed
- Stimulus submodule (functionality moved to parent scope)

## [1.0.1] - 2025-06-13

### Changed
- Refactored 2D Gaussian fit for spatial receptive fields (#13)

## [1.0.0] - 2025-05-18

### Added
- Forward correlation analysis
- `gOSI` and `gDSI` global orientation and direction selectivity indices
- Spatial receptive field 2D Gaussian fitting
- GitHub Actions documentation deployment

### Changed
- STRF array dimensions reordered to `[Height × Width × Time]`
- Refactored spike statistics module (#9)
- Adopted [Blue style](https://github.com/JuliaDiff/BlueStyle) code formatting

[Unreleased]: https://github.com/yizhanmiao/SpikeTriggered.jl/compare/v1.3.1...HEAD
[1.3.1]: https://github.com/yizhanmiao/SpikeTriggered.jl/compare/v1.3.0...v1.3.1
[1.3.0]: https://github.com/yizhanmiao/SpikeTriggered.jl/compare/v1.2.1...v1.3.0
[1.2.1]: https://github.com/yizhanmiao/SpikeTriggered.jl/compare/v1.2.0...v1.2.1
[1.2.0]: https://github.com/yizhanmiao/SpikeTriggered.jl/compare/v1.1.1...v1.2.0
[1.1.1]: https://github.com/yizhanmiao/SpikeTriggered.jl/compare/v1.1.0...v1.1.1
[1.1.0]: https://github.com/yizhanmiao/SpikeTriggered.jl/compare/v1.0.1...v1.1.0
[1.0.1]: https://github.com/yizhanmiao/SpikeTriggered.jl/compare/v1.0.0...v1.0.1
[1.0.0]: https://github.com/yizhanmiao/SpikeTriggered.jl/releases/tag/v1.0.0
