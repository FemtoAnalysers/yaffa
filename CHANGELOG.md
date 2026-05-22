# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Next:
### Added
- Code to generate the ppp source with CECA
- Wave function for ppp from E. Garrido et al., PLB 868 (2025) 139731'

### Changed
- Massive refactoring of the code: merging several branches and renaming many files
- Started consistent use of tests in pre-commit hooks

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
