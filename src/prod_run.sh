#!/bin/bash

# Example: ./prod_run.sh [--force]

OUTPUT_DIR="data/output"
FORCE=false

# Parse arguments
for arg in "$@"; do
  if [[ "$arg" == "--force" || "$arg" == "-f" ]]; then
    FORCE=true
  fi
done

CYAN='\033[1;36m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

printf "${CYAN}===============================================\n"
echo "Starting GROMACS Simulation Script"
echo "Output Directory: $OUTPUT_DIR"
echo "Force Run: $FORCE"
printf "===============================================${NC}\n"

(
cd "$OUTPUT_DIR" || { echo "Output directory not found! Exiting."; exit 1; }

# Step 1: Preparation
printf "\n${CYAN}---------- [Step 1: Preparation] ----------${NC}\n"
gmx grompp -f ../md-charmm.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr

# Step 2: MD Simulation
printf "\n${YELLOW}---------- [Step 2: MD Simulation] ----------${NC}\n"
if [[ "$FORCE" == true || ! -f md.edr ]]; then
  echo "Running MD Simulation with mdrun..."
  gmx mdrun -ntmpi 1 -ntomp 16 -pin on -v -deffnm md -nsteps 30000
else
  echo "MD output (md.edr) already exists. Skipping mdrun."
fi

# Step 3: Center Protein in Box
printf "\n${CYAN}---------- [Step 3: Center Protein in Box] ----------${NC}\n"
printf "1\n1\n" | gmx trjconv -s md.tpr -f md.xtc -o md_center.xtc -center -pbc mol

# Step 4: Validate Periodic Distance
printf "\n${YELLOW}---------- [Step 4: Validate Periodic Distance] ----------${NC}\n"
printf "1\n" | gmx mindist -s md.tpr -f md_center.xtc -pi -od mindist.xvg

# Step 5: RMSD Stability Check
printf "\n${CYAN}---------- [Step 5: RMSD Stability Check] ----------${NC}\n"
printf "4\n1\n" | gmx rms -s em.tpr -f md_center.xtc -o rmsd_xray.xvg -tu ns -xvg none
python3 ../../src/plot_property.py rmsd_xray.xvg RMSD

# Final: Report Methods
printf "\n${YELLOW}---------- [Final: Report Methods] ----------${NC}\n"
gmx report-methods -s md.tpr

printf "\n${CYAN}===============================================\n"
echo "GROMACS Simulation Script Complete"
printf "===============================================${NC} \n"

)
