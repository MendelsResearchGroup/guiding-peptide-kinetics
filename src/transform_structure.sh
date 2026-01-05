#!/bin/bash
#PBS -N transform_structure
#PBS -q  mendels_comb_q
#PBS -o output.log
#PBS -l select=2:ncpus=2:mpiprocs=2
#PBS  -M  alexander.z@technion.ac.il

set -e

# Usage: ./transform_structure.sh <base_filename>
# Example: ./transform_structure.sh 1fjs

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT_DEFAULT="$(cd "$SCRIPT_DIR/.." && pwd)"

if [ -n "${PBS_JOBID:-}" ]; then
    source ~/.bashrc
    conda activate gmx-plumed
    export OMP_NUM_THREADS=16
    REPO_ROOT=${REPO_ROOT:-${PBS_O_WORKDIR:-$HOME/work/protein-toolkit}}
else
    REPO_ROOT=${REPO_ROOT:-$REPO_ROOT_DEFAULT}
fi

cd "$REPO_ROOT"

CYAN='\033[1;36m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

if [ -n "${1:-}" ]; then
    BASE=$1
elif [ -n "${BASE:-}" ]; then
    BASE=$BASE
else
    echo -e "${YELLOW}Usage: $0 <base_filename>${NC}"
    echo -e "${YELLOW}Or set BASE env var (e.g., BASE=1fjs $0)${NC}"
    exit 1
fi

source "$REPO_ROOT/src/common/config.sh"

PDB_INPUT="$INPUT_DIR/${BASE}.pdb"
if [ -z "${GMX_CMD:-}" ] && [ -n "${PBS_JOBID:-}" ]; then
    GMX_CMD=gmx_mpi
else
    GMX_CMD=${GMX_CMD:-gmx}
fi
if [ -z "${MDRUN_FLAGS:-}" ] && [ "$GMX_CMD" = "gmx" ]; then
    MDRUN_FLAGS="-ntmpi 1 -ntomp 1"
fi

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
    $GMX_CMD pdb2gmx -f "${BASE}_protein.pdb" -o "${BASE}_processed.gro" -water tip3p -ff charmm22star -ignh

    printf "${YELLOW}Defining simulation box${NC}\n"
    $GMX_CMD editconf -f "${BASE}_processed.gro" -o "${BASE}_newbox.gro" -c -d 1.0 -bt dodecahedron

    printf "${CYAN}Adding solvent${NC}\n"
    $GMX_CMD solvate -cp "${BASE}_newbox.gro" -cs spc216.gro -o "${BASE}_solv.gro" -p topol.top

    # Create empty ions.mdp file
    printf "${YELLOW}Preparing for ion addition...${NC}\n"
    touch ions.mdp

    $GMX_CMD grompp -f ions.mdp -c "${BASE}_solv.gro" -p topol.top -o ions.tpr

    printf "${CYAN}Adding ions${NC}\n"
    printf "SOL\n" | $GMX_CMD genion -s ions.tpr -o "${BASE}_solv_ions.gro" -conc 0.15 -p topol.top -pname NA -nname CL -neutral

    printf "${YELLOW}Running energy minimization${NC}\n"
    $GMX_CMD grompp -f "../../emin-charmm.mdp" -c "${BASE}_solv_ions.gro" -p topol.top -o em.tpr
    if [ -n "${MDRUN_FLAGS:-}" ]; then
        $GMX_CMD mdrun -v -deffnm em $MDRUN_FLAGS
    else
        $GMX_CMD mdrun -v -deffnm em
    fi

    # Extract potential energy
    printf "${CYAN}Extracting potential energy${NC}\n"
    printf "Potential\n0\n" | $GMX_CMD energy -f em.edr -o potential.xvg -xvg none

    # Plot using external Python script
    printf "${YELLOW}Plotting potential${NC}\n"
    python3 ../../../src/plots/potential.py
)

printf "\n${CYAN}✅ Structure preparation and energy minimization complete! Check Epot and Fmax for success.${NC}\n"
