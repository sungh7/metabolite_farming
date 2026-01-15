# MD Simulation Protocol for Isoflavonoid-Enzyme Binding Validation

## Overview

This protocol validates novel isoflavonoid-enzyme binding interactions identified through computational analysis of ethylene-treated soybean metabolomics and proteomics data.

## Novel Candidates (Unreported in Literature)

| Rank | Ligand | Enzyme | Novelty | Rationale |
|------|--------|--------|---------|-----------|
| 1 | 6''-O-Malonylgenistin | 2-HIS (IFS homolog) | ★★★ | Malonylated form never docked |
| 2 | 6''-O-Acetyldaidzin | 2-HIS | ★★★ | Acetylated conjugate unstudied |
| 3 | 6''-O-Malonyldaidzin | 2-HIS | ★★★ | Major soybean isoflavone |
| 4 | Daidzin | GmIMaT | ★★ | Substrate with known Km |

## Directory Structure

```
md_protocol/
├── README.md                 # This file
├── run_md_workflow.sh        # Main execution script
└── mdp/
    ├── em.mdp               # Energy minimization
    ├── nvt.mdp              # NVT equilibration (100 ps)
    ├── npt.mdp              # NPT equilibration (100 ps)
    └── md_100ns.mdp         # Production MD (100 ns)
```

## Prerequisites

```bash
# Install GROMACS
conda create -n md gromacs -c conda-forge
conda activate md

# Verify installation
gmx --version
```

## Usage

### Step 1: Prepare Complex

After docking, combine protein and ligand:
```bash
# Convert docked pose to PDB
obabel docked.pdbqt -O ligand.pdb

# Merge with receptor
cat receptor.pdb ligand.pdb > complex.pdb
```

### Step 2: Run MD Workflow

```bash
./run_md_workflow.sh complex.pdb output_dir/
```

### Step 3: Production MD

```bash
cd output_dir/
gmx mdrun -deffnm md -nb gpu  # GPU acceleration
```

## Analysis Pipeline

### RMSD (Binding Stability)
```bash
echo "Backbone Backbone" | gmx rms -s md.tpr -f md.xtc -o rmsd.xvg
```

### RMSF (Flexibility)
```bash
echo "Protein" | gmx rmsf -s md.tpr -f md.xtc -o rmsf.xvg -res
```

### Hydrogen Bonds
```bash
gmx hbond -s md.tpr -f md.xtc -num hbond.xvg
```

### Binding Free Energy (MM-PBSA)

Using gmx_MMPBSA:
```bash
# Create input file
cat > mmpbsa.in << EOF
&general
startframe=1, endframe=1000, interval=10,
/
&pb
istrng=0.15, fillratio=4.0,
/
EOF

# Run MM-PBSA
gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 13 -ct md.xtc -o FINAL_RESULTS_MMPBSA.dat
```

## Expected Results

### Good Binding Indicators
- RMSD < 3 Å (stable binding)
- ΔG_bind < -7 kcal/mol (strong affinity)
- Persistent hydrogen bonds
- Low RMSF at binding site

### Comparison with Literature

| Parameter | Naringenin-CHI (Known) | Target for Novel |
|-----------|------------------------|------------------|
| Binding affinity | -7.5 kcal/mol | < -7 kcal/mol |
| RMSD | 1.5-2.5 Å | < 3 Å |
| H-bonds | 3-4 persistent | ≥ 2 persistent |

## References

1. [CHI-Naringenin structure (PDB 1EYQ)](https://www.rcsb.org/structure/1EYQ)
2. [IFS structure (PDB 8E83)](https://www.rcsb.org/structure/8E83)
3. [GmIMaT kinetics - Ahmad et al. 2017](https://www.frontiersin.org/articles/10.3389/fpls.2017.00735)

## Output Files

| File | Description |
|------|-------------|
| `md.xtc` | Trajectory (100 ns) |
| `md.edr` | Energy file |
| `rmsd.xvg` | RMSD over time |
| `rmsf.xvg` | Residue flexibility |
| `hbond.xvg` | H-bond count |
| `FINAL_RESULTS_MMPBSA.dat` | Binding free energy |
