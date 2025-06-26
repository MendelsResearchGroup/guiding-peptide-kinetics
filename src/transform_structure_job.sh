#!/bin/bash
#PBS -N transform_structure
#PBS -q  mendels_q
#PBS -o output.log
#PBS -l select=2:ncpus=16:mpiprocs=2
#PBS  -M  alexander.z@technion.ac.il

PBS_O_WORKDIR=$HOME/work/protein-toolkit

source ~/.bashrc
conda activate gmx-plumed

export OMP_NUM_THREADS=16

cd $PBS_O_WORKDIR

set -e

CYAN='\033[1;36m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

if [ -z "$BASE" ]; then
    echo -e "${YELLOW}Usage: $0 <base_filename>${NC}"
    exit 1
fi

BASE=$BASE
source "$PBS_O_WORKDIR/src/common/config.sh"

PDB_INPUT="$INPUT_DIR/${BASE}.pdb"

printf "${CYAN}\n---------- [Step 1: Create Output Directory] ----------${NC}\n"
mkdir -p "$OUTPUT_DIR"

printf "${YELLOW}\n---------- [Step 2: Clean PDB File] ----------${NC}\n"
grep -v "HETATM" "$PDB_INPUT" | grep -v "CONECT" >"$OUTPUT_DIR/${BASE}_protein.pdb"

printf "${CYAN}\n---------- [Step 3: Check for Missing Residues] ----------${NC}\n"
echo "Checking for missing residues:"
if grep "MISSING" "$PDB_INPUT"; then
    echo "Warning: Missing residues detected!"
else
    echo "No missing residues found."
fi

printf "${YELLOW}\n---------- [Step 4: GROMACS Processing] ----------${NC}\n"
(
    cd "$OUTPUT_DIR"


    printf "${CYAN}Generating topology and structure${NC}\n"
    gmx_mpi pdb2gmx -f "${BASE}_protein.pdb" -o "${BASE}_processed.gro" -water tip3p -ff charmm22star -ignh

    printf "${YELLOW}Defining simulation box${NC}\n"
    gmx_mpi editconf -f "${BASE}_processed.gro" -o "${BASE}_newbox.gro" -c -d 1.0 -bt dodecahedron

    printf "${CYAN}Adding solvent${NC}\n"
    gmx_mpi solvate -cp "${BASE}_newbox.gro" -cs spc216.gro -o "${BASE}_solv.gro" -p topol.top

    # Create empty ions.mdp file
    printf "${YELLOW}Preparing for ion addition...${NC}\n"
    touch ions.mdp

    gmx_mpi grompp -f ions.mdp -c "${BASE}_solv.gro" -p topol.top -o ions.tpr

    printf "${CYAN}Adding ions${NC}\n"
    printf "SOL\n" | gmx_mpi genion -s ions.tpr -o "${BASE}_solv_ions.gro" -conc 0.15 -p topol.top -pname NA -nname CL -neutral

    printf "${YELLOW}Running energy minimization${NC}\n"
    gmx_mpi grompp -f "../../emin-charmm.mdp" -c "${BASE}_solv_ions.gro" -p topol.top -o em.tpr
    gmx_mpi mdrun -v -deffnm em 

    # Extract potential energy
    printf "${CYAN}Extracting potential energy${NC}\n"
    printf "Potential\n0\n" | gmx_mpi energy -f em.edr -o potential.xvg -xvg none

    # Plot using external Python script
    printf "${YELLOW}Plotting potential${NC}\n"
    python3 ../../../src/plots/potential.py
)

printf "\n${CYAN}✅ Structure preparation and energy minimization complete! Check Epot and Fmax for success.${NC}\n"
