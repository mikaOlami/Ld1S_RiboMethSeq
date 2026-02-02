#!/usr/bin/env bash
set -euo pipefail

# L.donovani (Ld1S) RiboMethSeq analysis pipeline
# Starting from paired FASTQ (.fastq.gz or .fq.gz) in FASTQ/
#
# Repo layout assumed:
#   scripts/LD_RMS_Analysis_Pipeline.sh   (this script)
#   scripts/count_initiating.py
#   scripts/count3p.py
#   scripts/LD_RMS_Scores.py
#   DB/rRNA/LD_rRNA.genome
#   DB/rRNA/smalt_index/LD_rRNA_smalt_index.{smi,sma}
#
# Run from a working directory that contains FASTQ/ (recommended):
#   bash /path/to/repo/scripts/LD_RMS_Analysis_Pipeline.sh
#
# Optionally set REPO_DIR to the repo root if auto-detection fails:
#   REPO_DIR=/path/to/repo bash scripts/LD_RMS_Analysis_Pipeline.sh

# ---- Config ----
threads=8
max_jobs=25   # how many samples to process concurrently

# ---- Resolve repo root (so the script works from any working directory) ----
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="${REPO_DIR:-"$(cd -- "${SCRIPT_DIR}/.." && pwd)"}"

# ---- Reference paths inside repo ----
rRNA_DB="${REPO_DIR}/DB/rRNA/LD_rRNA_smalt_index"
genome="${REPO_DIR}/DB/rRNA/LD_rRNA.genome"

# ---- Script paths inside repo ----
COUNT_INIT="${REPO_DIR}/scripts/count_initiating.py"
COUNT_3P="${REPO_DIR}/scripts/count3p.py"
RMS_SCORES="${REPO_DIR}/scripts/LD_RMS_Scores.py"

# ---- Tools (assumed to be on PATH unless overridden) ----
SMALT="${SMALT:-smalt}"
SAMTOOLS="${SAMTOOLS:-samtools}"
BEDTOOLS="${BEDTOOLS:-bedtools}"
PYTHON="${PYTHON:-python3}"

# ---- Sanity checks ----
for tool in "$SMALT" "$SAMTOOLS" "$BEDTOOLS" "$PYTHON"; do
  command -v "$tool" >/dev/null 2>&1 || {
    echo "Error: required tool not found on PATH: $tool" >&2
    echo "Tip: set SMALT/SAMTOOLS/BEDTOOLS/PYTHON env vars if needed." >&2
    exit 1
  }
done

[[ -x "$(command -v "$SMALT")" ]] || { echo "Error: smalt not executable: $SMALT" >&2; exit 1; }
[[ -f "${rRNA_DB}.smi" && -f "${rRNA_DB}.sma" ]] || {
  echo "Error: SMALT index not found at:" >&2
  echo "  ${rRNA_DB}.smi / ${rRNA_DB}.sma" >&2
  exit 1
}
[[ -f "$genome" ]] || { echo "Error: genome file not found: $genome" >&2; exit 1; }
[[ -f "$COUNT_INIT" ]] || { echo "Error: missing script: $COUNT_INIT" >&2; exit 1; }
[[ -f "$COUNT_3P"   ]] || { echo "Error: missing script: $COUNT_3P" >&2; exit 1; }
[[ -f "$RMS_SCORES" ]] || { echo "Error: missing script: $RMS_SCORES" >&2; exit 1; }

[[ -d FASTQ ]] || {
  echo -e "Error: FASTQ directory does not exist in the current working directory.\n" \
          "Please create a directory named \"FASTQ\" and move the FASTQ files into it." >&2
  exit 1
}

mkdir -p Bams Beds Logs

# ---- simple job limiter ----
joblim_wait() {
  while (( $(jobs -rp | wc -l) >= max_jobs )); do
    wait -n || true
  done
}

# ---- per-sample pipeline (strictly sequential inside) ----
process_sample() {
  local r1="$1"
  local threads="$2"

  local r1_base
  r1_base="$(basename "$r1")"

  # derive R2 (supports .fastq.gz and .fq.gz; preserves extension)
  local r2_base
  r2_base="$(
    echo "$r1_base" | sed -E 's/(.*)_R1(_[0-9]+)?\.(fastq|fq)\.gz$/\1_R2\2.\3.gz/'
  )"
  local r2="FASTQ/$r2_base"

  if [[ ! -e "$r2" ]]; then
    echo "[WARN] Skipping: paired file not found for $r1 -> expected $r2" >&2
    return 0
  fi

  local base
  base="$(
    echo "$r1_base" | sed -E 's/_R1(_[0-9]+)?\.(fastq|fq)\.gz$//'
  )"

  local ubam="Bams/${base}_vs_rRNA.bam"
  local ubam_tmp="${ubam}.tmp"
  local bed="Beds/${base}_vs_rRNA.sorted.bed"
  local bed_tmp="${bed}.tmp"
  local init="${base}_vs_rRNA.sorted.init"
  local init_tmp="${init}.tmp"
  local threep="${base}_vs_rRNA.sorted.3p"
  local threep_tmp="${threep}.tmp"

  local log_map="Logs/${base}_vs_rRNA.smalt.log"
  local log_bed="Logs/${base}_vs_rRNA.bamtobed.log"
  local log_init="Logs/${base}_vs_rRNA.init.log"
  local log_3p="Logs/${base}_vs_rRNA.3p.log"

  echo "[START] $base"

  # --- MAP (blocking) ---
  if [[ ! -s "$ubam" || "$ubam" -ot "$r1" || "$ubam" -ot "$r2" ]]; then
    echo "[MAP: LD RMS vs rRNA] $base"
    "$SMALT" map -n "${threads}" "${rRNA_DB}" "$r1" "$r2" 2> "$log_map" \
      | "$SAMTOOLS" view -@ "${threads}" -b -f 0x02 -F 4 -o "$ubam_tmp" -
    "$SAMTOOLS" quickcheck -v "$ubam_tmp"
    mv -f "$ubam_tmp" "$ubam"
    [[ -s "$log_map" ]] || rm -f "$log_map"
  else
    echo "[MAP] SKIP (up-to-date) $base"
  fi

  if [[ ! -s "$ubam" ]]; then
    echo "[ERROR] BAM missing for $base" >&2
    return 1
  fi

  # --- BAM -> BED (blocking) ---
  if [[ ! -s "$bed" || "$bed" -ot "$ubam" ]]; then
    echo "[BED] $base"
    "$BEDTOOLS" bamtobed -bedpe -i "$ubam" 2> "$log_bed" \
      | cut -f1,2,6-9 \
      | sort -k1,1 -k2,2n > "$bed_tmp"
    [[ -s "$bed_tmp" ]] || { echo "[ERROR] BED empty for $base" >&2; rm -f "$bed_tmp"; return 1; }
    mv -f "$bed_tmp" "$bed"
    [[ -s "$log_bed" ]] || rm -f "$log_bed"
  else
    echo "[BED] SKIP (up-to-date) $base"
  fi

  if [[ ! -s "$bed" ]]; then
    echo "[ERROR] BED missing for $base" >&2
    return 1
  fi

  # --- INIT (blocking) ---
  if [[ ! -s "$init" || "$init" -ot "$bed" ]]; then
    echo "[INIT] $base"
    "$PYTHON" "$COUNT_INIT" "${genome}" "$bed" > "$init_tmp" 2> "$log_init" || true
    if [[ -s "$init_tmp" ]]; then mv -f "$init_tmp" "$init"; else rm -f "$init_tmp"; fi
    [[ -s "$log_init" ]] || rm -f "$log_init"
  else
    echo "[INIT] SKIP (up-to-date) $base"
  fi

  # --- 3P (blocking) ---
  if [[ ! -s "$threep" || "$threep" -ot "$bed" ]]; then
    echo "[3P] $base"
    "$PYTHON" "$COUNT_3P" "${genome}" "$bed" > "$threep_tmp" 2> "$log_3p" || true
    if [[ -s "$threep_tmp" ]]; then mv -f "$threep_tmp" "$threep"; else rm -f "$threep_tmp"; fi
    [[ -s "$log_3p" ]] || rm -f "$log_3p"
  else
    echo "[3P] SKIP (up-to-date) $base"
  fi

  echo "[DONE] $base"
}

# ---- discover R1s & launch each sample as ONE background job ----
shopt -s nullglob
found_any=false

for r1 in \
  FASTQ/*_R1_001.fastq.gz FASTQ/*_R1.fastq.gz \
  FASTQ/*_R1_001.fq.gz    FASTQ/*_R1.fq.gz
do
  [[ -e "$r1" ]] || continue
  found_any=true
  joblim_wait
  (
    set -euo pipefail
    process_sample "$r1" "$threads"
  ) &
done

if ! $found_any; then
  echo "No FASTQ/*_R1*.fastq.gz or FASTQ/*_R1*.fq.gz files found" >&2
  exit 0
fi

# wait for all samples to finish their full chain
wait
echo "[STAGE] All per-sample chains finished."

# ---- Global RMS step (depends on all .init/.3p) ----
echo "[RMS] Running LD_RMS_Scores.py"
"$PYTHON" "$RMS_SCORES" > Logs/RMS_Scores.log 2>&1 || true
[[ -s "Logs/RMS_Scores.log" ]] || rm -f "Logs/RMS_Scores.log"

# cleanup empty logs
find Logs -type f -name "*.log" -size 0 -delete 2>/dev/null || true

echo "SCRIPT FINISHED."
