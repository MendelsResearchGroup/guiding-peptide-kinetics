#!/bin/bash
set -e 

# Usage: ./prepare_structure.sh <base_filename>
# Example: ./prepare_structure.sh 1fjs

CYAN='\033[1;36m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

if [ -z "$1" ]; then
    echo -e "${YELLOW}Usage: $0 <base_filename>${NC}"
    exit 1
fi

BASE=$1
INPUT_DIR=data
OUTPUT_DIR=data/output

PDB_INPUT="$INPUT_DIR/${BASE}.pdb"

printf "${CYAN}\n---------- [Step 1: Create Output Directory] ----------${NC}\n"
mkdir -p "$OUTPUT_DIR"

printf "${YELLOW}\n---------- [Step 2: Clean PDB File] ----------${NC}\n"
grep -v "HETATM" "$PDB_INPUT" | grep -v "CONECT" > "$OUTPUT_DIR/${BASE}_protein.pdb"

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
  gmx pdb2gmx -f "${BASE}_protein.pdb" -o "${BASE}_processed.gro" -water tip3p -ff charmm22star

  printf "${YELLOW}Defining simulation box${NC}\n"
  gmx editconf -f "${BASE}_processed.gro" -o "${BASE}_newbox.gro" -c -d 1.0 -bt dodecahedron

  printf "${CYAN}Adding solvent${NC}\n"
  gmx solvate -cp "${BASE}_newbox.gro" -cs spc216.gro -o "${BASE}_solv.gro" -p topol.top

  # Create empty ions.mdp file
  printf "${YELLOW}Preparing for ion addition...${NC}\n"
  touch ions.mdp

  gmx grompp -f ions.mdp -c "${BASE}_solv.gro" -p topol.top -o ions.tpr

  printf "${CYAN}Adding ions${NC}\n"
  printf "SOL\n" | gmx genion -s ions.tpr -o "${BASE}_solv_ions.gro" -conc 0.15 -p topol.top -pname NA -nname CL -neutral

  printf "${YELLOW}Running energy minimization${NC}\n"
  gmx grompp -f "../emin-charmm.mdp" -c "${BASE}_solv_ions.gro" -p topol.top -o em.tpr
  gmx mdrun -v -deffnm em -ntmpi 1 -ntomp 1

  # Extract potential energy
  printf "${CYAN}Extracting potential energy${NC}\n"
  printf "Potential\n0\n" | gmx energy -f em.edr -o potential.xvg -xvg none

  # Plot using external Python script
  printf "${YELLOW}Plotting potential${NC}\n"
  python3 ../../src/plot_potential.py
)

printf "\n${CYAN}✅ Structure preparation and energy minimization complete! Check Epot and Fmax for success.${NC}\n"
