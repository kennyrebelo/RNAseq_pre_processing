# Changelog

All notable changes to this pipeline are recorded here. Versions follow
[semantic versioning](https://semver.org): given `MAJOR.MINOR.PATCH`, a change
that breaks an existing `config_RNAseq_pipeline.yaml` or `samples.tsv` bumps
MAJOR, a new capability bumps MINOR, and a fix that leaves the interface
unchanged bumps PATCH.

## [1.0.1] - 2026-11-25

### Fixed
- `docs/INSTALLATION.md` told users to install Snakemake unpinned, which
  resolves to 8.x and is incompatible with this workflow. It now builds the
  environment from `environment_snakemake.yaml` and documents each pin.
- Documented that `environment_RNAseq.yaml` uses `defaults` rather than
  `nodefaults`, unlike the other three environment files.

### Added
- Troubleshooting entry for `ModuleNotFoundError: No module named 'pandas'`
  during DAG building, on the login node and on compute nodes.

## [1.0.0] - 2026-11-25

First public release.

### Added
- Pre-processing: SRA download, Trim Galore trimming with FastQC, STAR
  alignment, unique-read filtering and sorting, Salmon quantification,
  aggregate MultiQC report.
- Optional genome-browser branch: RSeQC strandedness inference, strand
  splitting, normalised bigwig coverage and a UCSC track hub, per sample via
  `generate_bw`.
- Optional differential gene expression with tximport and DESeq2, via
  `dge.enabled`.
- Optional differential splicing with StringTie and rMATS, including
  novel-isoform gene assignment and cryptic-exon annotation, via
  `splicing.enabled`.
- Mixed single-end/paired-end and mixed human/mouse samples within one run.
- PBS/Torque cluster profile, and documentation covering installation,
  configuration, usage, outputs and troubleshooting.
