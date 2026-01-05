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
    MD_NSTEPS=100000
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

STRETCH_DEFFNM=${STRETCH_DEFFNM:-md}
UF_DEFFNM=${UF_DEFFNM:-$STRETCH_DEFFNM}
STRETCH_TPR=${STRETCH_TPR:-md.tpr}
UF_TPR=${UF_TPR:-$STRETCH_TPR}
STRETCH_CPI=${STRETCH_CPI:-}
UF_CPI=${UF_CPI:-${STRETCH_DEFFNM}.cpt}
UF_PLUMED=${UF_PLUMED:-$PLUMED_MD}

stretch_cpi_args=()
if [ -n "$STRETCH_CPI" ]; then
    stretch_cpi_args=(-cpi "$STRETCH_CPI")
fi

stretch_noappend_args=()
if [ "${STRETCH_NOAPPEND:-auto}" = "auto" ] && [ -n "$STRETCH_CPI" ] && [ "$STRETCH_DEFFNM" != "md" ]; then
    stretch_noappend_args=(-noappend)
elif [ "${STRETCH_NOAPPEND:-}" = "true" ]; then
    stretch_noappend_args=(-noappend)
fi

uf_cpi_args=()
if [ -n "$UF_CPI" ]; then
    uf_cpi_args=(-cpi "$UF_CPI")
fi

uf_plumed_args=()
if [ -n "$UF_PLUMED" ]; then
    uf_plumed_args=(-plumed "$UF_PLUMED")
fi

(
    cd "$OUTPUT_DIR"

    printf "${CYAN}\n---------- [Extract reference - disabled for now] ----------${NC}\n"
    # printf "1\n1\n" | $GMX_CMD trjconv -s md.tpr -f md_center.xtc -o reference.pdb -pbc mol -center -dump 200

    # Based on new GROMACS versions, OT1 and OT2 are not recognized as valid atom names.
    # awk '{name=substr($0,13,4); if (name ~ /OT1/) $0=substr($0,1,12) " O  " substr($0,17); else if (name ~ /OT2/) $0=substr($0,1,12) "OXT " substr($0,17); print }' reference.pdb > tmp.pdb && mv tmp.pdb reference.pdb

    printf "${YELLOW}\n---------- [Stretch protein] ----------${NC}\n"
    $GMX_CMD mdrun $STRETCH_MDRUN_FLAGS -v \
        -s "$STRETCH_TPR" -deffnm "$STRETCH_DEFFNM" -nsteps "$STRETCH_NSTEPS" \
        -plumed "$PLUMED_STRETCH" "${stretch_cpi_args[@]}" "${stretch_noappend_args[@]}"

    # printf "${CYAN}\n---------- [MD simulation] ----------${NC}\n"
    $GMX_CMD mdrun $MD_MDRUN_FLAGS -v \
        -s "$UF_TPR" -deffnm "$UF_DEFFNM" -nsteps "$MD_NSTEPS" \
        "${uf_plumed_args[@]}" "${uf_cpi_args[@]}"
    if [ -f COLVAR ]; then
        mv COLVAR COLVAR_FLAT
    fi

    # printf "\n${CYAN}---------- [Center Protein in Box] ----------${NC}\n"
    printf "1\n1\n" | $GMX_CMD trjconv -s "$UF_TPR" -f "${UF_DEFFNM}.xtc" -o "${UF_DEFFNM}_center.xtc" -center -pbc mol
)
