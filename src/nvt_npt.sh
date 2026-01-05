#!/bin/bash
#PBS -N nvt_npt
#PBS -q  mendels_comb_q
#PBS -o output.log
#PBS -l select=2:ncpus=2:mpiprocs=2
#PBS  -M  alexander.z@technion.ac.il

set -e

# Usage: ./nvt_npt.sh <base_filename> [-f|--force]
# Example: ./nvt_npt.sh 1fjs --force

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT_DEFAULT="$(cd "$SCRIPT_DIR/.." && pwd)"

if [ -n "${PBS_JOBID:-}" ]; then
    source ~/.bashrc
    conda activate gmx-plumed
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
    echo "Usage: $0 <base_filename> [-f|--force] or export BASE=..."
    exit 1
fi

source "$REPO_ROOT/src/common/config.sh"

FORCE=${FORCE:-false}
for arg in "$@"; do
    if [[ "$arg" == "--force" || "$arg" == "-f" ]]; then
        FORCE=true
    fi
done

(
    cd "$OUTPUT_DIR" || {
        echo "Output directory $OUTPUT_DIR not found! Exiting."
        exit 1
    }

    
    echo "Checking GROMACS version on HPC:"
    $GMX_CMD --version

    # Step 1: NVT grompp
    printf "${CYAN}\n---------- [Step 1: NVT grompp] ----------${NC}\n"
    if [[ "$FORCE" == true || ! -f nvt.tpr ]]; then
        $GMX_CMD grompp -f "../../nvt-charmm.mdp" -c em.gro -r em.gro -p topol.top -o nvt.tpr
    else
        echo "nvt.tpr already exists. Skipping grompp for NVT."
    fi

    # Step 2: NVT mdrun
    printf "${YELLOW}\n---------- [Step 2: NVT mdrun] ----------${NC}\n"
    if [[ "$FORCE" == true || ! -f nvt.edr ]]; then
        $GMX_CMD mdrun -v -deffnm nvt
    else
        echo "nvt.edr already exists. Skipping mdrun for NVT."
    fi

    # Step 3: NPT grompp
    printf "${CYAN}\n---------- [Step 3: NPT grompp] ----------${NC}\n"
    if [[ "$FORCE" == true || ! -f npt.tpr ]]; then
        $GMX_CMD grompp -f "../../npt-charmm.mdp" -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
    else
        echo "npt.tpr already exists. Skipping grompp for NPT."
    fi

    # Step 4: NPT mdrun
    printf "${YELLOW}\n---------- [Step 4: NPT mdrun] ----------${NC}\n"
    if [[ "$FORCE" == true || ! -f npt.edr ]]; then
        $GMX_CMD mdrun -v -deffnm npt
    else
        echo "npt.edr already exists. Skipping mdrun for NPT."
    fi

    # Step 5: Extract thermodynamic properties
    printf "${CYAN}\n---------- [Step 5: Extract Properties] ----------${NC}\n"
    echo "Temperature" | $GMX_CMD energy -f nvt.edr -o temperature.xvg -xvg none -b 20
    echo "Pressure" | $GMX_CMD energy -f npt.edr -o pressure.xvg -xvg none
    echo "Density" | $GMX_CMD energy -f npt.edr -o density.xvg -xvg none

    # # Step 6: Plot
    # printf "${YELLOW}\n---------- [Step 6: Plot Results] ----------${NC}\n"
    # python3 ../../../src/plots/property.py temperature.xvg temperature pressure.xvg pressure density.xvg density

    printf "\n${CYAN}✅ NVT and NPT simulations complete. Properties plotted.${NC}\n"
)
