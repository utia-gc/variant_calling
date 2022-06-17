# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.2] - 2022-06-16
### Fixed
- Git submodule structure changed to allow for downloading pipeline with `nextflow pull`

## [1.0.1] - 2022-06-16
### Fixed
- HISAT2 build output only takes genome fasta and index files, no exon and splice site files.

## [1.0.0] - 2022-06-13
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
- Only MultiQC reports use the unique prefix name to keep other processes from having to be repeated upon `-resume`.

### Deprecated

### Removed

### Fixed
- Made capitalization of SRAlign consistent (i.e. sralign to SRAlign).
- README quick start message.
- README pipeline overview.

### Security
