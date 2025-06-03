#!/bin/bash

# Usage: ./prepare_structure.sh <base_filename>
# Example: ./prepare_structure.sh 1fjs

if [ -z "$1" ]; then
    echo "Usage: $0 <base_filename>"
    exit 1
fi

BASE=$1
INPUT_DIR=data
OUTPUT_DIR=data/output

PDB_INPUT="$INPUT_DIR/${BASE}.pdb"

# Create output dir
mkdir -p "$OUTPUT_DIR"

# Clean PDB file
grep -v "HETATM" "$PDB_INPUT" | grep -v "CONECT" > "$OUTPUT_DIR/${BASE}_protein.pdb"

# Check for missing residues
echo "Checking for missing residues:"
if grep "MISSING" "$PDB_INPUT"; then
    echo "Warning: Missing residues detected!"
else
    echo "No missing residues found."
fi

# Run all subsequent steps inside output directory
(
cd "$OUTPUT_DIR"

# Generate topology and structure
gmx pdb2gmx -f "${BASE}_protein.pdb" -o "${BASE}_processed.gro" -water tip3p -ff charmm27

# Define box
gmx editconf -f "${BASE}_processed.gro" -o "${BASE}_newbox.gro" -c -d 1.0 -bt dodecahedron

# Solvate with water
gmx solvate -cp "${BASE}_newbox.gro" -cs spc216.gro -o "${BASE}_solv.gro" -p topol.top

# Create empty ions.mdp file
touch ions.mdp

# Assemble tpr file for ions
gmx grompp -f ions.mdp -c "${BASE}_solv.gro" -p topol.top -o ions.tpr

# Add ions
printf "SOL\n" | gmx genion -s ions.tpr -o "${BASE}_solv_ions.gro" -conc 0.15 -p topol.top -pname NA -nname CL -neutral

# Energy minimization
gmx grompp -f "../emin-charmm.mdp" -c "${BASE}_solv_ions.gro" -p topol.top -o em.tpr
gmx mdrun -v -deffnm em -ntmpi 1 -ntomp 1

# Extract potential energy
printf "Potential\n0\n" | gmx energy -f em.edr -o potential.xvg -xvg none

# Plot using external Python script
python3 ../../src/plot_potential.py
)

echo -e "\n✅ Structure preparation and energy minimization complete! Check Epot and Fmax for success."
