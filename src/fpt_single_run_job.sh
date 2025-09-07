#!/bin/bash
#PBS -N fpt_run
#PBS -q  mendels_comb_q
#PBS -o output.log
#PBS -l select=1:ncpus=4:mpiprocs=4
#PBS  -M  alexander.z@technion.ac.il

PBS_O_WORKDIR=$HOME/work/protein-toolkit

source ~/.bashrc
conda activate gmx-plumed


cd $PBS_O_WORKDIR

set -euo pipefail
PROTEIN="$BASE"

# ------------------ constants ------------------
MDP="md-charmm.mdp"
GRO="npt.gro"
CPT="npt.cpt"
TOP="topol.top"
REF="reference.pdb"
PLUMED_TEMPLATE="src/fpt_plumed/$PROTEIN.dat"
NSTEPS=40000000

RUN_ID=$(printf "%03d" "$ID")
BASE="$PROTEIN"

source "$PBS_O_WORKDIR/src/common/config.sh"

# FORCE=false && [[ ${3:-} =~ ^(-f|--force)$ ]] && FORCE=true
FORCE=true

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
    cp "../../$MDP" "$GRO" "$TOP" "$CPT" "../$REF" "$RUN_DIR/"

    sed "s/__ID__/${RUN_ID}/g" "../../../$PLUMED_TEMPLATE" >"$RUN_DIR/$PLUMED"

    # ------------------ GROMACS run ------------------
    cd "$RUN_DIR"

    printf "\n${CYAN}>> [grompp] Generating TPR: $TPR${NC}\n"
    gmx_mpi grompp -f "$MDP" -c "$GRO" -t "$CPT" -p "$TOP" -o "$TPR"

    printf "\n${YELLOW}>> [mdrun] Starting MD with HLDA bias${NC}\n"
    if [[ "$FORCE" == true || ! -f ${DEFFNM}.edr ]]; then
        gmx_mpi mdrun -deffnm "$DEFFNM" -nsteps $NSTEPS --plumed "$PLUMED"
    else
        echo "${DEFFNM}.edr exists – skipping mdrun"
    fi

    # ------------------ Center the protein ------------------
    if [[ -f ${DEFFNM}.xtc && -f ${DEFFNM}.tpr ]]; then
        printf "\n${CYAN}---------- [Center Protein in Box] ----------${NC}\n"
        printf "1\n1\n" | gmx_mpi trjconv -s "${DEFFNM}.tpr" -f "${DEFFNM}.xtc" \
            -o "${DEFFNM}_center.xtc" -center -pbc mol
    else
        echo "Warning: Cannot center – missing ${DEFFNM}.xtc or ${DEFFNM}.tpr"
    fi

    printf "\n${CYAN}Replica $RUN_ID finished; outputs in $OUTPUT_DIR/$RUN_DIR/${NC}\n"
)
