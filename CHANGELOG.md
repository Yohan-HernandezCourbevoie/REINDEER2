# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This project does **not** adheres to Semantic Versioning.

## [Unreleased]

## [1.1.0] - 2025-03-20

### Added

- User:
  - RD1-like output
  - Stranded mode
  - Prints statistics afer the indexation
  - Selection of the capacity of the index (in number of dataset)
  - Normalization option
  - Sampling mode
  - Infos mode
  - Add findere/fimpera support
- Technical:
  - This CHANGELOG file
  - Github CICD
  - System tests

### Changed

- Indexes have a .bin extension
- Indexes use 1024 files
- The input file is not read twice when querying using the colored output format
- The abundance chosen during a collision is the maximum possible abundance, instead of the minimum 
- Default filename does not include PACA anymore
- Indexation now crashes instead of producing a wrong index, when a file cannot be merged
- Replace integer by filename in median output
- Compilation now uses target native by default
- Uses simd minimizer
- New CLI design: most incompatible options are now incompatible in the CLI by design
- CLI now requires `reindeer2 query` instead of `reindeer2 --mode` 
- bin `Reindeer2` was renamed `reindeer2`

### Fixed

- Merge of chunks

## [1.0.0] - TODO date

### Added

- TODO