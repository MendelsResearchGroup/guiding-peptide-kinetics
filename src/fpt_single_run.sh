#!/usr/bin/env bash
#PBS -N fpt_run
#PBS -q  mendels_comb_q
#PBS -o output.log
#PBS -l select=1:ncpus=4:mpiprocs=4
#PBS  -M  alexander.z@technion.ac.il
#
# Usage:
#   ./fpt_single_run.sh <ID> <PROTEIN> [--force]
# Example:
#   ./fpt_single_run.sh 7 chignolin           # writes to data/chignolin/output/run_007/
#   ./fpt_single_run.sh 7 chignolin --force   # overwrite if already exists
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT_DEFAULT="$(cd "$SCRIPT_DIR/.." && pwd)"

if [ -n "${PBS_JOBID:-}" ]; then
    source ~/.bashrc
    conda activate gmx-plumed
    REPO_ROOT=${REPO_ROOT:-${PBS_O_WORKDIR:-$HOME/work/protein-toolkit}}
    GMX_CMD=${GMX_CMD:-gmx_mpi}
    FORCE_DEFAULT=true
    NSTEPS_DEFAULT=50000000
else
    REPO_ROOT=${REPO_ROOT:-$REPO_ROOT_DEFAULT}
    GMX_CMD=${GMX_CMD:-gmx}
    FORCE_DEFAULT=false
    NSTEPS_DEFAULT=3000000
fi

cd "$REPO_ROOT"

# ------------------ constants ------------------
MDP="md-charmm.mdp"
GRO="npt.gro"
CPT="npt.cpt"
TOP="topol.top"
REF="reference.pdb"

# ------------------ parse args ------------------
if [[ $# -ge 2 ]]; then
    RUN_ID_RAW="$1"
    PROTEIN="$2"
elif [[ -n "${ID:-}" && ( -n "${PROTEIN:-}" || -n "${BASE:-}" ) ]]; then
    RUN_ID_RAW="$ID"
    PROTEIN="${PROTEIN:-$BASE}"
else
    echo "Usage: $0 <ID> <PROTEIN> [--force] or export ID/PROTEIN (or BASE)" >&2
    exit 1
fi

RUN_ID=$(printf "%03d" "$RUN_ID_RAW")
BASE="$PROTEIN"

source "$REPO_ROOT/src/common/config.sh"

FORCE=${FORCE:-$FORCE_DEFAULT}
if [[ ${3:-} =~ ^(-f|--force)$ ]]; then
    FORCE=true
fi

NSTEPS=${NSTEPS:-$NSTEPS_DEFAULT}
PLUMED_BASE=${PLUMED_BASE:-"$REPO_ROOT/src/fpt_plumed/base.dat"}
PLUMED_TEMPLATE=${PLUMED_TEMPLATE:-"$REPO_ROOT/src/fpt_plumed/${PROTEIN}.dat"}
if [ -z "${MDRUN_FLAGS:-}" ]; then
    if [ -n "${PBS_JOBID:-}" ]; then
        MDRUN_FLAGS=""
    else
        MDRUN_FLAGS="-ntmpi 1 -ntomp 6 -nb gpu -pme gpu"
    fi
fi

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
    cp "$REPO_ROOT/data/$MDP" "$OUTPUT_DIR/$GRO" "$OUTPUT_DIR/$TOP" \
        "$OUTPUT_DIR/$CPT" "$REPO_ROOT/data/$BASE/$REF" "$RUN_DIR/"

    if [ -f "$PLUMED_BASE" ]; then
        cp "$PLUMED_BASE" "$RUN_DIR/"
    fi

    if [ ! -f "$PLUMED_TEMPLATE" ]; then
        echo "Plumed template not found: $PLUMED_TEMPLATE" >&2
        exit 1
    fi
    sed "s/__ID__/${RUN_ID}/g" "$PLUMED_TEMPLATE" >"$RUN_DIR/$PLUMED"

    # ------------------ GROMACS run ------------------
    cd "$RUN_DIR"

    printf "\n${CYAN}>> [grompp] Generating TPR: $TPR${NC}\n"
    $GMX_CMD grompp -f "$MDP" -c "$GRO" -t "$CPT" -p "$TOP" -o "$TPR"

    printf "\n${YELLOW}>> [mdrun] Starting MD with HLDA bias${NC}\n"
    if [[ "$FORCE" == true || ! -f ${DEFFNM}.edr ]]; then
        $GMX_CMD mdrun $MDRUN_FLAGS -v \
            -deffnm "$DEFFNM" -nsteps "$NSTEPS" \
            --plumed "$PLUMED"
    else
        echo "${DEFFNM}.edr exists – skipping mdrun"
    fi

    # ------------------ Center the protein ------------------
    if [[ -f ${DEFFNM}.xtc && -f ${DEFFNM}.tpr ]]; then
        printf "\n${CYAN}---------- [Center Protein in Box] ----------${NC}\n"
        printf "1\n1\n" | $GMX_CMD trjconv -s "${DEFFNM}.tpr" -f "${DEFFNM}.xtc" \
            -o "${DEFFNM}_center.xtc" -center -pbc mol
    else
        echo "Warning: Cannot center – missing ${DEFFNM}.xtc or ${DEFFNM}.tpr"
    fi

    printf "\n${CYAN}Replica $RUN_ID finished; outputs in $OUTPUT_DIR/$RUN_DIR/${NC}\n"
)
