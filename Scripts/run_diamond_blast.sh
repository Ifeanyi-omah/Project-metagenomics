#!/bin/bash

# Define directories and database
INPUT_DIR="/home/s2059633/ACEZAD_AGR_24072025_CONTIGS"
OUTPUT_DIR="/home/s2059633/Blast_nr"
DB_PATH="/home/s2059633/NR_DB/diamond_nr.dmnd"
TMP_DIR="/tmp"
EVALUE="1e-5"
THREADS=40

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

# Loop through each file ending with _modified.fa in the input directory
for QUERY_FILE in "${INPUT_DIR}"/*.fa; do
    # Extract the base name of the file (without directory and extension)
    BASENAME=$(basename "${QUERY_FILE}" .fa)
    
    # Define output file
    OUTPUT_FILE="${OUTPUT_DIR}/${BASENAME}_diamond_blast_nr.m9"
    
    # Run DIAMOND
    echo "Processing ${QUERY_FILE}..."
    diamond blastx \
        -d "${DB_PATH}" \
        -q "${QUERY_FILE}" \
        -o "${OUTPUT_FILE}" \
        --evalue "${EVALUE}" \
        --threads "${THREADS}" \
        --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames stitle \
        --tmpdir "${TMP_DIR}"
    
    echo "Output written to ${OUTPUT_FILE}"
done

echo "All files processed."
