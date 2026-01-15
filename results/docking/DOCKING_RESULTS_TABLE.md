# Molecular Docking Results: Novel Isoflavonoid-Enzyme Interactions

## Summary

This study identified novel isoflavonoid-enzyme binding interactions that have not been previously reported in the literature. Molecular docking was performed using AutoDock Vina with experimentally determined protein structures from related legume species.

## Table 1. Docking Results for Novel Isoflavonoid Candidates

| Rank | Ligand | Receptor | PDB | Binding Affinity (kcal/mol) | Novelty | Notes |
|------|--------|----------|-----|-----------------------------|---------| ------|
| 1 | **6''-O-Acetyldaidzin** | 2-Hydroxyisoflavanone Synthase | 8E83 | **-8.863** | ★★★ | Strongest binding; unreported |
| 2 | **6''-O-Acetylgenistin** | 2-Hydroxyisoflavanone Synthase | 8E83 | **-7.768** | ★★★ | Unreported |
| 3 | **6''-O-Malonyldaidzin** | 2-Hydroxyisoflavanone Synthase | 8E83 | **-7.600** | ★★★ | Major soybean isoflavone |
| 4 | **6''-O-Malonylgenistin** | 2-Hydroxyisoflavanone Synthase | 8E83 | **-7.485** | ★★★ | Unreported |
| 5 | Daidzin | 2-Hydroxyisoflavanone Dehydratase | 8EA1 | -7.357 | ★★ | Km = 36.28 μM (Ahmad 2017) |
| 6 | Genistin | 2-Hydroxyisoflavanone Dehydratase | 8EA1 | -6.863 | ★★ | Km = 23.04 μM (Ahmad 2017) |
| 7 | Liquiritigenin | Chalcone Isomerase | 1EYQ | -5.296 | ★★ | Type II CHI substrate |
| 8 | Daidzein | Chalcone Isomerase | 1EYQ | -4.902 | ★★ | Product comparison |

**Novelty Legend:**
- ★★★ = No prior computational or experimental binding study reported
- ★★ = Partial novelty (enzyme kinetics known, binding mode unknown)

## Table 2. Binding Site Analysis for 6''-O-Acetyldaidzin + 8E83 (Best Hit)

### Key Interacting Residues

| Residue | Distance (Å) | Interaction Type | Significance |
|---------|--------------|------------------|--------------|
| ASP230 | 2.84 | H-bond acceptor | Key polar contact |
| TYR212 | 3.03 | H-bond / π-stacking | Aromatic interaction |
| ILE45 | 3.15 | Hydrophobic | Core binding pocket |
| LYS223 | 3.17 | H-bond donor | Stabilizes glycoside |
| THR215 | 3.23 | H-bond | Polar contact |
| LEU48 | 3.26 | Hydrophobic | Binding pocket |
| GLU234 | 3.30 | H-bond acceptor | Ionic interaction |

### Binding Site Composition

| Property | Count | Residues |
|----------|-------|----------|
| Hydrophobic | 4 | ILE45, LEU48, PHE224, PHE237 |
| Polar | 4 | THR215, THR233, TYR212, ASN222 |
| Positively charged | 2 | LYS223, HIS49 |
| Negatively charged | 2 | ASP230, GLU234 |
| H-bond candidates | 7 | ASP230, TYR212, LYS223, THR215, GLU234, LEU48, HOH846 |

## Key Findings

1. **All novel candidates (★★★) showed strong binding** (< -7.0 kcal/mol)
2. **Acetylated forms > Malonylated forms** in binding affinity
3. **Conjugated forms > Aglycones** (Daidzin > Daidzein)
4. **Binding site is amphipathic** (mixed hydrophobic and polar residues)

## Methods

- **Software**: AutoDock Vina 1.2.5
- **Grid size**: 35 × 35 × 35 Å
- **Exhaustiveness**: 16
- **Poses generated**: 5 per ligand
- **Receptor structures**: PDB 8E83 (2.0 Å), 8EA1 (2.4 Å), 1EYQ (1.85 Å)
- **Ligand structures**: PubChem 3D conformers

## References

1. PDB 8E83: Shao et al. (2022) Commun Biol. DOI: 10.1038/s42003-022-04222-x
2. PDB 8EA1: Shao et al. (2022) - 2-Hydroxyisoflavanone dehydratase
3. PDB 1EYQ: Jez et al. (2000) Nat Struct Biol. DOI: 10.1038/nsb0900_786
4. GmIMaT kinetics: Ahmad et al. (2017) Front Plant Sci. DOI: 10.3389/fpls.2017.00735

## Files

- `docking_results_figure.png` - Publication-ready figure (300 DPI)
- `docking_results_figure.pdf` - Vector format
- `all_results.json` - Raw docking scores
- `all_candidates/` - Individual docking outputs
