# Livestock FASTQ to BAM with Nextflow and SLURM

A Nextflow DSL2 workflow for paired-end Illumina whole-genome sequencing of modern livestock. It performs read quality control, trimming, alignment, run merging, duplicate handling, and BAM/depth reporting. It does not call variants or process ancient DNA.

## Workflow

FastQC → Trimmomatic → FastQC → BWA-MEM/SAMtools sort → merge runs by animal → Picard MarkDuplicates → SAMtools statistics and depth summaries.

The workflow also provides optional read-header normalization, manually reviewed adapter additions, trimming overrides, and per-run/per-animal reports. Read groups retain animal, run and library identifiers. Inspect the exact filtering and duplicate-handling options in `main.nf` before applying them to another study.

## Status and validation

This is a publication-preparation copy of the livestock workflow. The related `nextflow_modern_batch` workflow has been run on the author's HPC. Local tests here cover Python QC/reporting helpers; they are not an end-to-end validation of this livestock workflow on other clusters.

Container integration is not implemented in this version. A separate udocker/FastQC 0.13.0 version test succeeded under SLURM, but that does not establish container support for the complete workflow.

## Requirements

Linux; Nextflow and compatible Java; SLURM; Python 3; FastQC; Trimmomatic; BWA; SAMtools; Picard. A reference FASTA and matching BWA indices are required. Inputs and reference data are supplied separately.

`nextflow.config` and `index_reference.sh` currently contain example site-specific module/Conda activation commands. Adapt these to your HPC; the workflow is not claimed to work unchanged on every cluster. No software binaries or research datasets are distributed here.

Versions observed in the related batch environment: FastQC 0.13.0, Trimmomatic 0.39, BWA 0.7.19-r1273, SAMtools/HTSlib 1.22.1, Picard 3.4.0 and Temurin Java 21.0.8+9. These are observed versions, not an automatically installed or container-locked environment.

## Configure and run

1. Copy `settings.env.example` to `settings.env` and configure the reference, account, partition, adapter FASTA and Nextflow activation script.
2. Replace tool-path placeholders and adapt module/Conda commands in `nextflow.config` and `index_reference.sh`.
3. Prepare `samples.tsv` using `samples.example.tsv` (one row per sequencing run).
4. Prepare matching reference indices if absent; submit `index_reference.sh` with appropriate SLURM resources.
5. Run the Nextflow controller in a persistent allocation/session allowed by your HPC. Child tasks use the SLURM executor.

```bash
bash run.sh /absolute/path/to/samples.tsv --stage raw_qc
# After reviewing QC and configuring any required overrides:
bash run.sh /absolute/path/to/samples.tsv --stage full -resume
```

Review the Nextflow version pinned in `run.sh` for your installation. Keep the work directory: most intermediate/published files are symbolic links. Use a separate execution directory for each reference/species and retain provenance alongside final outputs.

Detailed setup and review instructions: [Portuguese guide](GUIDE_PT.md).

## Local checks

```bash
python3 -m unittest discover -s tests -v
bash -n run.sh
bash -n index_reference.sh
```

## Reproducibility

The launcher records configuration snapshots, the validated manifest, adapter inputs, tool-controller information and Nextflow trace/report/timeline outputs in `provenance/`. Data, local settings and run outputs are excluded from Git by default. Store large outputs separately and cite their location when reporting an analysis.

## License

A distribution license has not yet been selected. Third-party tools retain their own licenses.

## Sample sheet

The input is a tab-separated file with one row per sequencing run:

| Column | Meaning |
| --- | --- |
| animal_id | Biological sample identifier; shared across runs from the same animal |
| run | Unique sequencing run identifier |
| library_id | Library identifier; shared only across runs from the same physical library |
| header_format | `standard`, or `hash_mate` for the supported special read-header format |
| R1 | Absolute path to the forward FASTQ file |
| R2 | Absolute path to the reverse FASTQ file |

See `samples.example.tsv` for an example.
