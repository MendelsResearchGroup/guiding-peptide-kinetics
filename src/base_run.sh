#!/bin/bash
#PBS -N base_run
#PBS -q  mendels_q
#PBS -o output.log
#PBS -l select=2:ncpus=4:mpiprocs=4:host=n151
#PBS  -M  alexander.z@technion.ac.il

set -e

# Example: ./base_run.sh WT [--force]

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT_DEFAULT="$(cd "$SCRIPT_DIR/.." && pwd)"

if [ -n "${PBS_JOBID:-}" ]; then
  source ~/.bashrc
  conda activate gmx-plumed
  export OMP_NUM_THREADS=16
  REPO_ROOT=${REPO_ROOT:-${PBS_O_WORKDIR:-$HOME/work/protein-toolkit}}
  GMX_CMD=${GMX_CMD:-gmx_mpi}
  FORCE=${FORCE:-true}
else
  REPO_ROOT=${REPO_ROOT:-$REPO_ROOT_DEFAULT}
  GMX_CMD=${GMX_CMD:-gmx}
  FORCE=${FORCE:-false}
fi

cd "$REPO_ROOT"

if [ -n "${1:-}" ]; then
  BASE=$1
elif [ -n "${BASE:-}" ]; then
  BASE=$BASE
else
  echo "Usage: $0 <base_filename> [--force] or export BASE=..."
  exit 1
fi

source "$REPO_ROOT/src/common/config.sh"

for arg in "$@"; do
  if [[ "$arg" == "--force" || "$arg" == "-f" ]]; then
    FORCE=true
  fi
done

if [ -z "${MD_NSTEPS:-}" ]; then
  if [ -n "${PBS_JOBID:-}" ]; then
    MD_NSTEPS=20000000
  else
    MD_NSTEPS=50000000
  fi
fi

PLUMED_FILE=${PLUMED_FILE:-"$REPO_ROOT/src/plumed/base.dat"}
if [ -z "${MDRUN_FLAGS:-}" ]; then
  if [ -n "${PBS_JOBID:-}" ]; then
    MDRUN_FLAGS="-pin on"
  else
    MDRUN_FLAGS="-ntmpi 1 -ntomp 12 -pin on -nb gpu -pme gpu"
  fi
fi

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
  $GMX_CMD grompp -f ../../md-charmm.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr

  printf "\n${YELLOW}---------- [MD Simulation] ----------${NC}\n"
  if [[ "$FORCE" == true || ! -f md.edr ]]; then
    echo "Running MD Simulation with mdrun..."
    $GMX_CMD mdrun $MDRUN_FLAGS -v -deffnm md -nsteps "$MD_NSTEPS" -plumed "$PLUMED_FILE"
    cp COLVAR COLVAR_PIN
  else
    echo "MD output (md.edr) already exists. Skipping mdrun."
  fi

  printf "\n${CYAN}---------- [Center Protein in Box] ----------${NC}\n"
  printf "1\n1\n" | $GMX_CMD trjconv -s md.tpr -f md.xtc -o md_center.xtc -center -pbc mol

  printf "\n${YELLOW}---------- [Validate Periodic Distance] ----------${NC}\n"
  printf "1\n" | $GMX_CMD mindist -s md.tpr -f md_center.xtc -pi -od mindist.xvg

  # printf "\n${CYAN}---------- [RMSD Stability Check] ----------${NC}\n"
  # printf "4\n1\n" | gmx rms -s em.tpr -f md_center.xtc -o rmsd_xray.xvg -tu ns -xvg none
  # python3 ../../../src/plots/property.py rmsd_xray.xvg RMSD

  printf "\n${YELLOW}---------- [Report Methods] ----------${NC}\n"
  $GMX_CMD report-methods -s md.tpr

  printf "\n${CYAN}===============================================\n"
  echo "GROMACS Simulation Script Complete"
  printf "===============================================${NC} \n"

)
