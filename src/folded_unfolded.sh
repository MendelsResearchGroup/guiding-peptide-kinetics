#!/bin/bash
#PBS -N folded_unfolded
#PBS -q  mendels_comb_q
#PBS -o output.log
#PBS -l select=2:ncpus=4:mpiprocs=4
#PBS  -M  alexander.z@technion.ac.il

set -e

# Usage: ./folded_unfolded.sh <base_filename> [--force]
# Produces a folded unbiased trajectory, then stretches and continues unbiased for UF.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT_DEFAULT="$(cd "$SCRIPT_DIR/.." && pwd)"

if [ -n "${PBS_JOBID:-}" ]; then
  REPO_ROOT=${REPO_ROOT:-${PBS_O_WORKDIR:-$HOME/work/protein-toolkit}}
else
  REPO_ROOT=${REPO_ROOT:-$REPO_ROOT_DEFAULT}
fi

cd "$REPO_ROOT"

if [ -n "${1:-}" ]; then
  BASE=$1
elif [ -n "${BASE:-}" ]; then
  BASE=$BASE
else
  echo "Usage: $0 <base_filename> [--force] or export BASE=..."
  exit 1
fi

FORCE=false
for arg in "$@"; do
  if [[ "$arg" == "--force" || "$arg" == "-f" ]]; then
    FORCE=true
  fi
done

source "$REPO_ROOT/src/common/config.sh"

if [ ! -f "$OUTPUT_DIR/npt.gro" ] || [ ! -f "$OUTPUT_DIR/npt.cpt" ]; then
  echo "Missing npt.gro or npt.cpt in $OUTPUT_DIR; run nvt_npt.sh first."
  exit 1
fi

printf "${CYAN}\n========== [Folded: unbiased base run] ==========${NC}\n"
if [ "$FORCE" = true ]; then
  "$REPO_ROOT/src/base_run.sh" "$BASE" --force
else
  "$REPO_ROOT/src/base_run.sh" "$BASE"
fi

printf "${YELLOW}\n========== [Unfolded: stretch then unbiased] ==========${NC}\n"
export STRETCH_DEFFNM=${STRETCH_DEFFNM:-stretch}
export UF_DEFFNM=${UF_DEFFNM:-uf}
export STRETCH_TPR=${STRETCH_TPR:-md.tpr}
export STRETCH_CPI=${STRETCH_CPI:-md.cpt}
export UF_TPR=${UF_TPR:-$STRETCH_TPR}
export UF_CPI=${UF_CPI:-${STRETCH_DEFFNM}.cpt}

"$REPO_ROOT/src/stretch_run.sh" "$BASE"

UNBIASED_DIR="$OUTPUT_DIR/unbiased"
mkdir -p "$UNBIASED_DIR"

if [ -f "$OUTPUT_DIR/md.xtc" ]; then
  cp "$OUTPUT_DIR/md.xtc" "$UNBIASED_DIR/f.xtc"
fi
if [ -f "$OUTPUT_DIR/COLVAR_PIN" ]; then
  cp "$OUTPUT_DIR/COLVAR_PIN" "$UNBIASED_DIR/f.colvar"
fi
if [ -f "$OUTPUT_DIR/md_center.xtc" ]; then
  cp "$OUTPUT_DIR/md_center.xtc" "$UNBIASED_DIR/f_center.xtc"
fi

if [ -f "$OUTPUT_DIR/${UF_DEFFNM}.xtc" ]; then
  cp "$OUTPUT_DIR/${UF_DEFFNM}.xtc" "$UNBIASED_DIR/uf.xtc"
fi
if [ -f "$OUTPUT_DIR/COLVAR_FLAT" ]; then
  cp "$OUTPUT_DIR/COLVAR_FLAT" "$UNBIASED_DIR/uf.colvar"
fi
if [ -f "$OUTPUT_DIR/${UF_DEFFNM}_center.xtc" ]; then
  cp "$OUTPUT_DIR/${UF_DEFFNM}_center.xtc" "$UNBIASED_DIR/uf_center.xtc"
fi
