# protein-toolkit

To generate a metadynamics run biased by HLDA values, follow these steps:

1. **Prepare the protein structure**  
   Transform the input protein structure, solvate it, add ions, etc., by running:  
   ```bash
   ./src/transform_structure.sh 5awl
   ```

2. **Equilibrate temperature and pressure**  
   Run NVT and NPT equilibration to stabilize the system:  
   ```bash
   ./src/nvt_npt.sh 5awl
   ```

3. **Run a production MD simulation**  
   Perform an unbiased MD simulation starting from the equilibrated structure:  
   ```bash
   ./src/prod_run.sh
   ```

4. **Stretch the protein to its unfolded state**  
   Run a biased simulation that gradually stretches the protein and then follows with an unbiased run starting from the stretched conformation:  
   ```bash
   ./src/stretch.sh 5awl
   ```  
   This will generate a `COLVAR_STRETCH` file containing descriptors from the unfolded trajectory.

5. **Compute HLDA weights**  
   Once you have the COLVAR files for both the folded and unfolded states, calculate the HLDA projection:  
   ```bash
   python src/hlda.py
   ```  
   This script will output the HLDA weights.

6. **Update the PLUMED configuration**  
   Copy the resulting HLDA coefficients into your PLUMED input file (typically `plumed/hlda.dat`) to be used for biasing future runs.
