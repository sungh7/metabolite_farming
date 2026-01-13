# Docking Pipeline Guide

## 1. Environment Setup
The docking script requires `openbabel` and `vina` which are not installed in the standard environment. Please install them using Conda:

```bash
# Create a new environment (highly recommended)
conda create -n docking -c conda-forge openbabel vina rdkit pandas requests tqdm

# Activate
conda activate docking
```

## 2. Run Docking
The script `src/run_docking.py` automates the following:
1.  Downloads AlphaFold protein structures (UniProt ID).
2.  Downloads Metabolite structures (KEGG ID).
3.  Converts PDB/MOL to PDBQT.
4.  Runs AutoDock Vina (Blind docking).

```bash
# Run from the project root
python src/run_docking.py \
    --input results/gnn/top1_enzymes_for_docking.csv \
    --outdir results/docking/output
```

## 3. Output
- `results/docking/docking_summary.csv`: Contains binding affinities (kcal/mol) for all pairs.
- `results/docking/output/`: Individual directories for each pair with PDBQT logs.

## 4. Expected Results
- **Binding Affinity**: Lower is better (e.g., -9.0 vs -6.0).
- **Target**: Pairs with **Strong Affinity (< -7.0 kcal/mol)** and **High GNN Score**.
