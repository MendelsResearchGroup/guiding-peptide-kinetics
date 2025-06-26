#!/bin/bash
#PBS -N stretch_run
#PBS -q  mendels_q
#PBS -o output.log
#PBS -l select=2:ncpus=8:mpiprocs=2
#PBS  -M  alexander.z@technion.ac.il

PBS_O_WORKDIR=$HOME/work/protein-toolkit

source ~/.bashrc
conda activate gmx-plumed

export OMP_NUM_THREADS=16

cd $PBS_O_WORKDIR

set -e
set -x 

BASE=$BASE

source "$PBS_O_WORKDIR/src/common/config.sh"


BASE=$1
(
    cd "$OUTPUT_DIR"

    printf "${CYAN}\n---------- [Extract reference - disabled for now] ----------${NC}\n"
    # printf "1\n1\n" | gmx trjconv -s md.tpr -f md_center.xtc -o reference.pdb -pbc mol -center -dump 200

    # Based on new GROMACS versions, OT1 and OT2 are not recognized as valid atom names.
    # awk '{name=substr($0,13,4); if (name ~ /OT1/) $0=substr($0,1,12) " O  " substr($0,17); else if (name ~ /OT2/) $0=substr($0,1,12) "OXT " substr($0,17); print }' reference.pdb > tmp.pdb && mv tmp.pdb reference.pdb

    printf "${YELLOW}\n---------- [Stretch protein] ----------${NC}\n"
    gmx_mpi mdrun -pin on -deffnm md -nsteps 400000 -plumed ../../../src/plumed/stretch.dat

    # printf "${CYAN}\n---------- [MD simulation] ----------${NC}\n"
    gmx_mpi mdrun  -pin on -deffnm md -nsteps 20000000 -plumed ../../../src/plumed/base.dat -cpi md.cpt
    mv COLVAR COLVAR_FLAT

    # printf "\n${CYAN}---------- [Center Protein in Box] ----------${NC}\n"
    printf "1\n1\n" | gmx_mpi trjconv -s md.tpr -f md.xtc -o md_center.xtc -center -pbc mol

)
