# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Next:
### Added
- Code to generate the ppp source with CECA
- Wave function for ppp from E. Garrido et al., PLB 868 (2025) 139731'
- `WaveFunction` class; fitting with Argonne v18 and pΛ wave functions, `.wf` format support
- Momentum smearing of CFs via resolution matrix or phase space (`Smear.py`)
- `SuperFitter`: support for `TGraph`/`TGraphMultiErrors`, x-axis rescaling and unit conversion
- 3-body source functions with resonances
- QA script (`RunQA.py`) for collisions, tracks and triplets with automatic particle-system detection
- Z test panel in `CompareGraphs.py`
- `utils.analysis.Divide`: division of `TGraph`s by histograms, graphs and functions, and of histograms by histograms and graphs

### Changed
- Massive refactoring of the code: merging several branches and renaming many files
- Started consistent use of tests in pre-commit hooks
- Renamed source functions consistently
- Units switched from GeV to MeV
- `ComputeRawCF.py` rewritten for the Run 3 femto framework (3-body)

### Fixed
- Ratio panel of `CompareGraphs.py`: works with `TGraph`s and no longer crashes on 2D objects
- application of smearing matrix as transposed
- Normalization of the ppr source function
- Wrong axes/binning in `ChangeUnits2D`

## 0.1.0
### Added
- Changelog
- [WIP] Script to simulate the source. Currently still
- root files with pT shapes or protons and antiprotons
- Script to compute the source size
- [WIP] Script to make the diff of two root files
- Header with source function
- Header with ROOT function wrappers for TF1s
- Logger macro for debugging and error handling
- Utility for projecting 2D histograms
- tests with different CECA configurations and refejkrence files
