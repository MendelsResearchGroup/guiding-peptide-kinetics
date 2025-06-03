#!/bin/bash

# Usage: ./run_nvt_npt.sh <base_filename>
# Example: ./run_nvt_npt.sh 1fjs

BASE=$1
INPUT_DIR=data
OUTPUT_DIR=data/output

(
    cd "$OUTPUT_DIR" || exit

    # NVT
    if [ ! -f nvt.tpr ]; then
        gmx grompp -f "../nvt-charmm.mdp" -c em.gro -r em.gro -p topol.top -o nvt.tpr
    else
        echo "nvt.tpr already exists. Skipping grompp for NVT."
    fi

    if [ ! -f nvt.edr ]; then
        gmx mdrun -ntmpi 1 -ntomp 8 -v -deffnm nvt
    else
        echo "nvt.edr already exists. Skipping mdrun for NVT."
    fi


    # NPT
    if [ ! -f npt.tpr ]; then
        gmx grompp -f "../npt-charmm.mdp" -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
    else
        echo "npt.tpr already exists. Skipping grompp for NPT."
    fi

    if [ ! -f npt.edr ]; then
        gmx mdrun -ntmpi 1 -ntomp 8 -v -deffnm npt
    else
        echo "npt.edr already exists. Skipping mdrun for NPT."
    fi

    echo "Temperature" | gmx energy -f nvt.edr -o temperature.xvg -xvg none -b 20
    echo "Pressure" | gmx energy -f npt.edr -o pressure.xvg -xvg none
    echo "Density" | gmx energy -f npt.edr -o density.xvg -xvg none

    python3 ../../src/plot_property.py temperature.xvg temperature pressure.xvg pressure density.xvg density
)
