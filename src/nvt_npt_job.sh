#!/bin/bash
#PBS -N gmx_test
#PBS -q  mendels_q
#PBS -o output.log
#PBS -l select=2:ncpus=16:mpiprocs=2
#PBS  -M  alexander.z@technion.ac.il

PBS_O_WORKDIR=$HOME/work/protein-toolkit

source ~/.bashrc
conda activate gmx-plumed

export OMP_NUM_THREADS=16

cd $PBS_O_WORKDIR

# Usage: ./run_nvt_npt.sh <base_filename> [-f|--force]
# Example: ./run_nvt_npt.sh 1fjs --force

BASE=$BASE

source "$PBS_O_WORKDIR/src/common/config.sh"

FORCE=false

for arg in "$@"; do
    if [[ "$arg" == "--force" || "$arg" == "-f" ]]; then
        FORCE=true
    fi
done

(
    cd "$PBS_O_WORKDIR/$OUTPUT_DIR" || exit

    
    echo "Checking GROMACS version on HPC:"
    gmx_mpi --version

    # Step 1: NVT grompp
    printf "${CYAN}\n---------- [Step 1: NVT grompp] ----------${NC}\n"
    if [[ "$FORCE" == true || ! -f nvt.tpr ]]; then
        gmx_mpi grompp -f "../../nvt-charmm.mdp" -c em.gro -r em.gro -p topol.top -o nvt.tpr
    else
        echo "nvt.tpr already exists. Skipping grompp for NVT."
    fi

    # Step 2: NVT mdrun
    printf "${YELLOW}\n---------- [Step 2: NVT mdrun] ----------${NC}\n"
    if [[ "$FORCE" == true || ! -f nvt.edr ]]; then
        gmx_mpi mdrun -v -deffnm nvt
    else
        echo "nvt.edr already exists. Skipping mdrun for NVT."
    fi

    # Step 3: NPT grompp
    printf "${CYAN}\n---------- [Step 3: NPT grompp] ----------${NC}\n"
    if [[ "$FORCE" == true || ! -f npt.tpr ]]; then
        gmx_mpi grompp -f "../../npt-charmm.mdp" -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
    else
        echo "npt.tpr already exists. Skipping grompp for NPT."
    fi

    # Step 4: NPT mdrun
    printf "${YELLOW}\n---------- [Step 4: NPT mdrun] ----------${NC}\n"
    if [[ "$FORCE" == true || ! -f npt.edr ]]; then
        gmx_mpi mdrun -v -deffnm npt
    else
        echo "npt.edr already exists. Skipping mdrun for NPT."
    fi

    # Step 5: Extract thermodynamic properties
    printf "${CYAN}\n---------- [Step 5: Extract Properties] ----------${NC}\n"
    echo "Temperature" | gmx_mpi energy -f nvt.edr -o temperature.xvg -xvg none -b 20
    echo "Pressure" | gmx_mpi energy -f npt.edr -o pressure.xvg -xvg none
    echo "Density" | gmx_mpi energy -f npt.edr -o density.xvg -xvg none

    # Step 6: Plot
    printf "${YELLOW}\n---------- [Step 6: Plot Results] ----------${NC}\n"
    python3 ../../../src/plots/property.py temperature.xvg temperature pressure.xvg pressure density.xvg density

    printf "\n${CYAN}✅ NVT and NPT simulations complete. Properties plotted.${NC}\n"
)
