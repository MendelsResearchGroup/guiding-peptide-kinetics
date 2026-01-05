#!/bin/bash
#PBS -N stretch_run
#PBS -q  mendels_q
#PBS -o output.log
#PBS -l select=2:ncpus=8:mpiprocs=2
#PBS  -M  alexander.z@technion.ac.il

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT_DEFAULT="$(cd "$SCRIPT_DIR/.." && pwd)"

if [ -n "${PBS_JOBID:-}" ]; then
    source ~/.bashrc
    conda activate gmx-plumed
    export OMP_NUM_THREADS=16
    REPO_ROOT=${REPO_ROOT:-${PBS_O_WORKDIR:-$HOME/work/protein-toolkit}}
    GMX_CMD=${GMX_CMD:-gmx_mpi}
else
    REPO_ROOT=${REPO_ROOT:-$REPO_ROOT_DEFAULT}
    GMX_CMD=${GMX_CMD:-gmx}
fi

cd "$REPO_ROOT"

if [ -n "${1:-}" ]; then
    BASE=$1
elif [ -n "${BASE:-}" ]; then
    BASE=$BASE
else
    echo "Usage: $0 <base_filename> or export BASE=..."
    exit 1
fi

source "$REPO_ROOT/src/common/config.sh"

PLUMED_STRETCH=${PLUMED_STRETCH:-"$REPO_ROOT/src/plumed/stretch.dat"}
PLUMED_MD=${PLUMED_MD:-"$REPO_ROOT/src/plumed/base.dat"}

if [ -z "${STRETCH_NSTEPS:-}" ]; then
    STRETCH_NSTEPS=400000
fi
if [ -z "${MD_NSTEPS:-}" ]; then
    if [ -n "${PBS_JOBID:-}" ]; then
        MD_NSTEPS=20000000
    else
        MD_NSTEPS=5000000
    fi
fi

if [ -z "${STRETCH_MDRUN_FLAGS:-}" ]; then
    if [ -n "${PBS_JOBID:-}" ]; then
        STRETCH_MDRUN_FLAGS="-pin on"
    else
        STRETCH_MDRUN_FLAGS="-ntmpi 1 -ntomp 12 -nb gpu -pme gpu -pin on"
    fi
fi
if [ -z "${MD_MDRUN_FLAGS:-}" ]; then
    if [ -n "${PBS_JOBID:-}" ]; then
        MD_MDRUN_FLAGS="-pin on"
    else
        MD_MDRUN_FLAGS="-ntmpi 1 -ntomp 12 -nb gpu -pme gpu -pin on"
    fi
fi

(
    cd "$OUTPUT_DIR"

    printf "${CYAN}\n---------- [Extract reference - disabled for now] ----------${NC}\n"
    # printf "1\n1\n" | $GMX_CMD trjconv -s md.tpr -f md_center.xtc -o reference.pdb -pbc mol -center -dump 200

    # Based on new GROMACS versions, OT1 and OT2 are not recognized as valid atom names.
    # awk '{name=substr($0,13,4); if (name ~ /OT1/) $0=substr($0,1,12) " O  " substr($0,17); else if (name ~ /OT2/) $0=substr($0,1,12) "OXT " substr($0,17); print }' reference.pdb > tmp.pdb && mv tmp.pdb reference.pdb

    printf "${YELLOW}\n---------- [Stretch protein] ----------${NC}\n"
    $GMX_CMD mdrun $STRETCH_MDRUN_FLAGS -v -deffnm md -nsteps "$STRETCH_NSTEPS" -plumed "$PLUMED_STRETCH"

    # printf "${CYAN}\n---------- [MD simulation] ----------${NC}\n"
    $GMX_CMD mdrun $MD_MDRUN_FLAGS -v -deffnm md -nsteps "$MD_NSTEPS" -plumed "$PLUMED_MD" -cpi md.cpt
    mv COLVAR COLVAR_FLAT

    # printf "\n${CYAN}---------- [Center Protein in Box] ----------${NC}\n"
    printf "1\n1\n" | $GMX_CMD trjconv -s md.tpr -f md.xtc -o md_center.xtc -center -pbc mol
)
