#!/bin/bash
set -euo pipefail
[[ $# -eq 1 && -s "$1" ]] || { echo 'Usage: sbatch ... index_reference.sh reference.fa' >&2; exit 1; }
REFERENCE=$(realpath "$1")
module purge
module load GCC/14.3.0 SAMtools/1.22.1
SAMTOOLS=$(command -v samtools)
module load site/langs Anaconda3/2025.12-1
CONDA_BASE=$(conda info --base)
set +u
source "${CONDA_BASE}/etc/profile.d/conda.sh"
conda activate /absolute/path/to/bwa_env
set -u
bwa index "$REFERENCE"
"$SAMTOOLS" faidx "$REFERENCE"
