#!/usr/bin/env bash
#
# Usage:
#   ./fpt_run.sh <ID> <PROTEIN> [--force]
# Example:
#   ./fpt_run.sh 7 chignolin           # writes to data/output/chignolin/run_007/
#   ./fpt_run.sh 7 chignolin --force   # overwrite if already exists
#
set -euo pipefail

# ------------------ constants ------------------
MDP="md-charmm.mdp"
GRO="npt.gro"
CPT="npt.cpt"
TOP="topol.top"
REF="reference.pdb"
PLUMED_TEMPLATE="src/fpt_plumed.dat"
NSTEPS=3000000

# ------------------ parse args ------------------
if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <ID> <PROTEIN> [--force]" >&2
    exit 1
fi

RUN_ID=$(printf "%03d" "$1")
PROTEIN="$2"
BASE="$PROTEIN"

source "$(dirname "$0")/common/config.sh"

FORCE=false && [[ ${3:-} =~ ^(-f|--force)$ ]] && FORCE=true

DEFFNM="run_${RUN_ID}"
RUN_DIR="run_${RUN_ID}"
TPR="${DEFFNM}.tpr"
PLUMED="${DEFFNM}_plumed.dat"

(
    cd "$OUTPUT_DIR" || {
        echo "Output directory $OUTPUT_DIR not found! Exiting."
        exit 1
    }

    # ------------------ Setup replica folder ------------------
    if [[ -d "$RUN_DIR" && $FORCE == false ]]; then
        echo "Run folder $RUN_DIR exists – use --force to overwrite" >&2
        exit 1
    fi

    mkdir -p "$RUN_DIR"
    cp "../../$MDP" "$GRO" "$TOP" "$CPT" "$REF" "$RUN_DIR/"

    sed "s/__ID__/${RUN_ID}/g" "../../../$PLUMED_TEMPLATE" >"$RUN_DIR/$PLUMED"

    # ------------------ GROMACS run ------------------
    cd "$RUN_DIR"

    printf "\n${CYAN}>> [grompp] Generating TPR: $TPR${NC}\n"
    gmx grompp -f "$MDP" -c "$GRO" -t "$CPT" -p "$TOP" -o "$TPR"

    printf "\n${YELLOW}>> [mdrun] Starting MD with HLDA bias${NC}\n"
    if [[ "$FORCE" == true || ! -f ${DEFFNM}.edr ]]; then
        gmx mdrun -ntmpi 1 -ntomp 6 -v \
            -deffnm "$DEFFNM" -nsteps $NSTEPS \
            -nb gpu -pme gpu --plumed "$PLUMED"
    else
        echo "${DEFFNM}.edr exists – skipping mdrun"
    fi

    # ------------------ Center the protein ------------------
    if [[ -f ${DEFFNM}.xtc && -f ${DEFFNM}.tpr ]]; then
        printf "\n${CYAN}---------- [Center Protein in Box] ----------${NC}\n"
        printf "1\n1\n" | gmx trjconv -s "${DEFFNM}.tpr" -f "${DEFFNM}.xtc" \
            -o "${DEFFNM}_center.xtc" -center -pbc mol
    else
        echo "Warning: Cannot center – missing ${DEFFNM}.xtc or ${DEFFNM}.tpr"
    fi

    printf "\n${CYAN}Replica $RUN_ID finished; outputs in $OUTPUT_DIR/$RUN_DIR/${NC}\n"
)
