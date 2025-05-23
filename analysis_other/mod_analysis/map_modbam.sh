#!/bin/bash

# Exit on error
set -e

# Script configuration
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
LOG_DIR="${SCRIPT_DIR}/logs"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="${LOG_DIR}/map_modbam_${TIMESTAMP}.log"

# Input parameters
SAMPLE="TUE_02_03UL_old"
INPUT_BAM="data/basecalled/sup/24070LRa002_04503/24070LRa002_04503.sup.unmapped.bam"
REFERENCE="assembly/output/verkko/${SAMPLE}/assembly.fasta"
OUTPUT_BAM="${SCRIPT_DIR}/${SAMPLE}.mod.bam"
THREADS=40
SORT_THREADS=4
SORT_MEM="4G"

# Create log directory if it doesn't exist
mkdir -p "${LOG_DIR}"

# Logging function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "${LOG_FILE}"
}

# Error handling function
handle_error() {
    log "ERROR: An error occurred in the script at line $1"
    exit 1
}

# Set up error handling
trap 'handle_error $LINENO' ERR

# Start script
log "Starting nucflag analysis for sample ${SAMPLE}"
log "Working directory: $(pwd)"
log "Script directory: ${SCRIPT_DIR}"
log "Project root: ${PROJECT_ROOT}"

# Activate conda environment
log "Activating conda environment: minimap2"
source "$(conda info --base)/etc/profile.d/conda.sh" || { log "Failed to source conda"; exit 1; }
conda activate minimap2 || { log "Failed to activate conda environment"; exit 1; }

# Check if input files exist
if [ ! -f "${INPUT_BAM}" ]; then
    log "ERROR: Input BAM file not found: ${INPUT_BAM}"
    exit 1
fi

if [ ! -f "${REFERENCE}" ]; then
    log "ERROR: Reference FASTA file not found: ${REFERENCE}"
    exit 1
fi

# Run minimap2 alignment
log "Starting minimap2 alignment"
start_time=$(date +%s)

# Remove unmapped, non primary, secondary alignments
samtools fastq -TMM,ML --reference "${REFERENCE}" "${INPUT_BAM}" \
| minimap2 --MD -ax lr:hq --eqx -y \
    -t "${THREADS}" \
    "${REFERENCE}" - 2>>"${LOG_FILE}" \
| samtools view -h -F 2308 \
| samtools sort -m "${SORT_MEM}" -@ "${SORT_THREADS}" -o "${OUTPUT_BAM}" -O BAM - >>"${LOG_FILE}" 2>&1

end_time=$(date +%s)
duration=$((end_time - start_time))
log "Minimap2 alignment completed in ${duration} seconds"

# Index the BAM file
log "Indexing BAM file"
samtools index "${OUTPUT_BAM}" || { log "Failed to index BAM file"; exit 1; }

log "Analysis completed successfully"