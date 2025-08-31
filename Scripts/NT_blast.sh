#!/usr/bin/env bash
set -euo pipefail

# 1) point BLASTDB at the folder containing all nt.* volumes
export BLASTDB=/home/s2059633/NT_DB/NT_DB
DATABASE=nt

# 2) input/output dirs
INPUT_DIR=/home/s2059633/ACEZCAD_contigs
OUTPUT_DIR=/home/s2059633/ACECAD_virus_NT
mkdir -p "$OUTPUT_DIR"

# 3) verify the DB once
echo "Checking BLAST database '$DATABASE' in $BLASTDB..."
if ! blastdbcmd -info -db "$DATABASE" &>/dev/null; then
  echo "Error: BLAST database '$DATABASE' not found or invalid. Exiting." >&2
  exit 1
fi
echo "BLAST DB OK."

# 4) loop over each virus contig FASTA
shopt -s nullglob
for QUERY in "$INPUT_DIR"/*_virus_contigs.fa; do
  # derive a clean sample prefix
  # e.g. 1_01_23_0574_S41_virus_contigs.fa → 1_01_23_0574_S41
  SAMPLE=$(basename "$QUERY" "_virus_contigs.fa")
  OUT="${OUTPUT_DIR}/${SAMPLE}_virus_contigs_blast_nt.m9"

  echo "[$(date +'%Y-%m-%d %H:%M:%S')] Blasting $SAMPLE..."
  if ! blastn \
       -task dc-megablast \
       -db "$DATABASE" \
       -num_threads 4 \
       -query "$QUERY" \
       -outfmt "7 std staxid ssciname scomname stitle" \
       -evalue 1e-5 \
       -out "$OUT"; then
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] blastn FAILED for $SAMPLE" >&2
    continue
  fi
  echo "[$(date +'%Y-%m-%d %H:%M:%S')] Completed $SAMPLE → $(basename "$OUT")"
done

echo "All files processed."
