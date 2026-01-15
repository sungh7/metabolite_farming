#!/bin/bash
#
# GROMACS MD Simulation Workflow
# For isoflavonoid-enzyme binding validation
#
# Novel candidates:
#   1. 6-O-Malonylgenistin + enzyme (★★★)
#   2. 6-O-Acetyldaidzin + enzyme (★★★)
#   3. Daidzin + enzyme (★★)
#
# Usage: ./run_md_workflow.sh <complex.pdb> <output_dir>
#

set -e

# Check arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <complex.pdb> <output_dir>"
    echo "  complex.pdb: Protein-ligand complex from docking"
    echo "  output_dir:  Output directory for MD results"
    exit 1
fi

COMPLEX=$1
OUTDIR=$2
MDP_DIR="$(dirname $0)/mdp"

# Check GROMACS
if ! command -v gmx &> /dev/null; then
    echo "Error: GROMACS (gmx) not found. Install with:"
    echo "  conda install -c conda-forge gromacs"
    exit 1
fi

echo "========================================================"
echo "GROMACS MD Simulation Workflow"
echo "========================================================"
echo "Input: $COMPLEX"
echo "Output: $OUTDIR"
echo ""

# Create output directory
mkdir -p $OUTDIR
cd $OUTDIR

#----------------------------------------------
# Step 1: Generate topology
#----------------------------------------------
echo "[1/6] Generating topology..."
gmx pdb2gmx -f $COMPLEX -o complex.gro -water tip3p -ff amber99sb-ildn << EOF
1
EOF

#----------------------------------------------
# Step 2: Define box and solvate
#----------------------------------------------
echo "[2/6] Creating simulation box..."
gmx editconf -f complex.gro -o box.gro -c -d 1.2 -bt cubic

echo "[3/6] Adding solvent..."
gmx solvate -cp box.gro -cs spc216.gro -o solvated.gro -p topol.top

#----------------------------------------------
# Step 3: Add ions
#----------------------------------------------
echo "[4/6] Adding ions (neutralization)..."
gmx grompp -f $MDP_DIR/em.mdp -c solvated.gro -p topol.top -o ions.tpr -maxwarn 2
echo "SOL" | gmx genion -s ions.tpr -o system.gro -p topol.top -pname NA -nname CL -neutral

#----------------------------------------------
# Step 4: Energy minimization
#----------------------------------------------
echo "[5/6] Energy minimization..."
gmx grompp -f $MDP_DIR/em.mdp -c system.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em

# Check minimization
echo "Potential energy after EM:"
echo "Potential" | gmx energy -f em.edr -o potential.xvg

#----------------------------------------------
# Step 5: NVT equilibration
#----------------------------------------------
echo "[6a/6] NVT equilibration (100 ps)..."
gmx grompp -f $MDP_DIR/nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt

# Check temperature
echo "Temperature" | gmx energy -f nvt.edr -o temperature.xvg

#----------------------------------------------
# Step 6: NPT equilibration
#----------------------------------------------
echo "[6b/6] NPT equilibration (100 ps)..."
gmx grompp -f $MDP_DIR/npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt

# Check pressure and density
echo "Pressure" | gmx energy -f npt.edr -o pressure.xvg
echo "Density" | gmx energy -f npt.edr -o density.xvg

#----------------------------------------------
# Step 7: Production MD (100 ns)
#----------------------------------------------
echo "[7/7] Production MD (100 ns)..."
echo "  This will take several hours to days depending on hardware."
echo "  For GPU acceleration, use: gmx mdrun -deffnm md -nb gpu"

gmx grompp -f $MDP_DIR/md_100ns.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr

# Option A: Run directly
# gmx mdrun -deffnm md -nb gpu

# Option B: Submit to cluster (uncomment if using SLURM)
# sbatch << EOF
# #!/bin/bash
# #SBATCH --job-name=md_isoflavonoid
# #SBATCH --ntasks=1
# #SBATCH --cpus-per-task=8
# #SBATCH --gres=gpu:1
# #SBATCH --time=48:00:00
# gmx mdrun -deffnm md -nb gpu -ntmpi 1 -ntomp 8
# EOF

echo ""
echo "========================================================"
echo "WORKFLOW COMPLETE"
echo "========================================================"
echo ""
echo "Output files:"
echo "  em.gro     - Energy minimized structure"
echo "  nvt.gro    - NVT equilibrated"
echo "  npt.gro    - NPT equilibrated"
echo "  md.tpr     - Production run input (ready to run)"
echo ""
echo "To run production MD:"
echo "  gmx mdrun -deffnm md -nb gpu"
echo ""
echo "Analysis after MD:"
echo "  # RMSD"
echo "  gmx rms -s md.tpr -f md.xtc -o rmsd.xvg"
echo ""
echo "  # RMSF"
echo "  gmx rmsf -s md.tpr -f md.xtc -o rmsf.xvg"
echo ""
echo "  # Hydrogen bonds"
echo "  gmx hbond -s md.tpr -f md.xtc -num hbond.xvg"
echo ""
echo "  # Binding free energy (MM-PBSA)"
echo "  gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 13 -ct md.xtc"
