#!/bin/bash
set -e 

# Example: ./prod_run.sh [--force]

source "$(dirname "$0")/config.sh"
FORCE=false

# Parse arguments
for arg in "$@"; do
  if [[ "$arg" == "--force" || "$arg" == "-f" ]]; then
    FORCE=true
  fi
done



printf "${CYAN}===============================================\n"
echo "Starting GROMACS Simulation Script"
echo "Output Directory: $OUTPUT_DIR"
echo "Force Run: $FORCE"
printf "===============================================${NC}\n"

(
cd "$OUTPUT_DIR" || { echo "Output directory not found! Exiting."; exit 1; }

printf "\n${CYAN}---------- [Preparation] ----------${NC}\n"
gmx grompp -f ../md-charmm.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr

printf "\n${YELLOW}---------- [MD Simulation] ----------${NC}\n"
if [[ "$FORCE" == true || ! -f md.edr ]]; then
  echo "Running MD Simulation with mdrun..."
  gmx mdrun -ntmpi 1 -ntomp 12 -pin on -v -deffnm md -nb gpu -pme gpu -nsteps 50000000 -plumed ../../src/plumed/plumed.dat 
  cp COLVAR COLVAR_PIN 
else
  echo "MD output (md.edr) already exists. Skipping mdrun."
fi

printf "\n${CYAN}---------- [Center Protein in Box] ----------${NC}\n"
printf "1\n1\n" | gmx trjconv -s md.tpr -f md.xtc -o md_center.xtc -center -pbc mol

printf "\n${YELLOW}---------- [Validate Periodic Distance] ----------${NC}\n"
printf "1\n" | gmx mindist -s md.tpr -f md_center.xtc -pi -od mindist.xvg

printf "\n${CYAN}---------- [RMSD Stability Check] ----------${NC}\n"
printf "4\n1\n" | gmx rms -s em.tpr -f md_center.xtc -o rmsd_xray.xvg -tu ns -xvg none
python3 ../../src/plot_property.py rmsd_xray.xvg RMSD

printf "\n${YELLOW}---------- [Report Methods] ----------${NC}\n"
gmx report-methods -s md.tpr

printf "\n${CYAN}===============================================\n"
echo "GROMACS Simulation Script Complete"
printf "===============================================${NC} \n"

)
