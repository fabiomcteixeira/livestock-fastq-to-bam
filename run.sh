#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")"
[[ $# -ge 1 ]] || { echo 'Usage: bash run.sh /absolute/samples_runs.tsv [Nextflow options]' >&2; exit 1; }
MANIFEST=$(readlink -f "$1")
shift
source settings.env
[[ -s "$REFERENCE" && "$SLURM_ACCOUNT" != YOUR_ACCOUNT ]] || { echo 'Edit REFERENCE and SLURM_ACCOUNT in settings.env' >&2; exit 1; }
for EXT in amb ann bwt pac sa; do
    [[ -s "${REFERENCE}.${EXT}" ]] || { echo "Missing BWA index: ${REFERENCE}.${EXT}" >&2; exit 1; }
done
source "$NEXTFLOW_ACTIVATE"
export NXF_VER=26.10.5
mkdir -p provenance
python3 bin/validate_manifest.py "$MANIFEST" provenance/samples_runs.checked.tsv
# Set ADAPTERS_BASE when using a different standard adapter FASTA.
ADAPTERS_BASE=${ADAPTERS_BASE:-/absolute/path/to/Trimmomatic-0.39/adapters/TruSeq3-PE.fa}
COMBINED=$(python3 bin/prepare_review.py "$ADAPTERS_BASE" provenance/adapters)
STAMP=$(date +%Y%m%d_%H%M%S)_$$
SNAP="$PWD/provenance/$STAMP"
mkdir -p "$SNAP"
cp settings.env main.nf nextflow.config additional_adapters.fa approved_raw_adapter_runs.txt trim_overrides.tsv "$SNAP/"
cp provenance/samples_runs.checked.tsv "$SNAP/samples_runs.tsv"
cp "$ADAPTERS_BASE" "$SNAP/base_adapters.fa"
cp -R bin "$SNAP/"
module list > "$SNAP/controller_modules.txt" 2>&1
nextflow -version > "$SNAP/nextflow_version.txt"
printf '%s\n' "$COMBINED" > "$SNAP/combined_adapters_path.txt"
printf '%q ' "$@" > "$SNAP/arguments.txt"
NEXTFLOW_EXIT=0

nextflow run main.nf -c nextflow.config -profile slurm -ansi-log false \
    -with-trace "$SNAP/trace.tsv" \
    -with-report "$SNAP/report.html" \
    -with-timeline "$SNAP/timeline.html" \
    --reference "$REFERENCE" \
    --account "$SLURM_ACCOUNT" \
    --partition "$SLURM_PARTITION" \
    --manifest "$SNAP/samples_runs.tsv" \
    --adapters "$ADAPTERS_BASE" \
    --combined_adapters "$COMBINED" \
    --approvals "$SNAP/approved_raw_adapter_runs.txt" \
    --trim_overrides "$SNAP/trim_overrides.tsv" \
    "$@" || NEXTFLOW_EXIT=$?

printf '%s\n' "$NEXTFLOW_EXIT" > "$SNAP/nextflow_exit_code.txt"

printf '\nNextflow exit code: %s\n' "$NEXTFLOW_EXIT"
printf 'Execution records: %s\n' "$SNAP"
SUMMARY_EXIT=0

python3 "$SNAP/bin/summarize_execution.py" \
    "$SNAP/samples_runs.tsv" \
    "$SNAP/trace.tsv" \
    "$SNAP/execution_summary" || SUMMARY_EXIT=$?

if [[ "$SUMMARY_EXIT" -ne 0 ]]; then
    printf 'ERROR: execution summary failed (exit %s)\n' \
        "$SUMMARY_EXIT" >&2
    if [[ "$NEXTFLOW_EXIT" -eq 0 ]]; then
        NEXTFLOW_EXIT="$SUMMARY_EXIT"
    fi
fi
exit "$NEXTFLOW_EXIT"
