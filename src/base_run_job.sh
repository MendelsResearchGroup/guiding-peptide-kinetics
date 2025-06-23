#!/bin/bash
#PBS -N base_run
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

# Example: ./prod_run.sh chignolin [--force]
BASE=$BASE

source "$PBS_O_WORKDIR/src/common/config.sh"
FORCE=false

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
  cd "$OUTPUT_DIR" || {
    echo "Output directory not found! Exiting."
    exit 1
  }

  printf "\n${CYAN}---------- [Preparation] ----------${NC}\n"
  gmx_mpi grompp -f ../../md-charmm.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr

  printf "\n${YELLOW}---------- [MD Simulation] ----------${NC}\n"
  if [[ "$FORCE" == true || ! -f md.edr ]]; then
    echo "Running MD Simulation with mdrun..."
    gmx_mpi mdrun -pin on -deffnm md -nsteps 200000 -plumed ../../../src/plumed/plumed.dat
    cp COLVAR COLVAR_PIN
  else
    echo "MD output (md.edr) already exists. Skipping mdrun."
  fi

  printf "\n${CYAN}---------- [Center Protein in Box] ----------${NC}\n"
  printf "1\n1\n" | gmx_mpi trjconv -s md.tpr -f md.xtc -o md_center.xtc -center -pbc mol

  printf "\n${YELLOW}---------- [Validate Periodic Distance] ----------${NC}\n"
  printf "1\n" | gmx_mpi mindist -s md.tpr -f md_center.xtc -pi -od mindist.xvg


  printf "\n${YELLOW}---------- [Report Methods] ----------${NC}\n"
  gmx_mpi report-methods -s md.tpr

  printf "\n${CYAN}===============================================\n"
  echo "GROMACS Simulation Script Complete"
  printf "===============================================${NC} \n"

)
