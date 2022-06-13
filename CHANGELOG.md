# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- CHANGELOG.md to track changes between versions.
- Header to display pipeline info like title and name on STDOUT and in log.
- Stub commands for easier testing of pipeline logic.
- Set seed for HISAT2 build and align pseudo-random processes.
- Add splice-aware option for building HISAT2 index.

### Changed
- Help message usage statement.
- Help message parameter help message display.
- Parameter specifications such as descriptions stored in `parameter_specifications.json`.
- Combined raw and trimmed reads FastQC steps to avoid duplication of modules and subworkflows.

### Deprecated

### Removed

### Fixed
- Made capitalization of SRAlign consistent (i.e. sralign to SRAlign).

### Security
