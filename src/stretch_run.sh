#!/bin/bash
set -e 

source "$(dirname "$0")/config.sh"

if [ -z "$1" ]; then
    echo -e "${YELLOW}Usage: $0 <base_filename>${NC}"
    exit 1
fi

BASE=$1
(cd "$OUTPUT_DIR" 

printf "${CYAN}\n---------- [Extract reference] ----------${NC}\n"
printf "1\n1\n" | gmx trjconv -s md.tpr -f md_center.xtc -o reference.pdb -pbc mol -center -dump 200

# Based on new GROMACS versions, OT1 and OT2 are not recognized as valid atom names.
awk '{name=substr($0,13,4); if (name ~ /OT1/) $0=substr($0,1,12) " O  " substr($0,17); else if (name ~ /OT2/) $0=substr($0,1,12) "OXT " substr($0,17); print }' reference.pdb > tmp.pdb && mv tmp.pdb reference.pdb


printf "${YELLOW}\n---------- [Stretch protein] ----------${NC}\n"
gmx mdrun -ntmpi 1 -ntomp 16 -nb gpu -pme gpu -pin on -v -deffnm md -nsteps 400000 -plumed ../../src/plumed/stretch.dat


printf "${CYAN}\n---------- [MD simulation] ----------${NC}\n"
gmx mdrun -ntmpi 1 -ntomp 16 -nb gpu -pme gpu -pin on -v -deffnm md -nsteps 50000 -plumed ../../src/plumed/plumed.dat -cpi md.cpt

printf "\n${CYAN}---------- [Center Protein in Box] ----------${NC}\n"
printf "1\n1\n" | gmx trjconv -s md.tpr -f md.xtc -o md_center.xtc -center -pbc mol

)