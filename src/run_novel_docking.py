#!/usr/bin/env python3
"""
AutoDock Vina docking script for novel isoflavonoid-enzyme binding validation.

Novel MD Candidates (UNREPORTED in literature):
1. 6-O-Malonylgenistin + enzyme (malonyltransferase product)
2. 6-O-Acetyldaidzin + enzyme (acetyltransferase product)
3. Daidzin + GmIMaT (substrate, Km known but binding mode unknown)

Uses homologous enzyme structures from related legumes:
- 8E83: 2-HIS from Medicago truncatula (IFS homolog)
- 1EYQ: CHI from Medicago sativa (CHI homolog)
- 8EA1: 2-HID from Pueraria lobata

Author: Generated for ethylene-isoflavonoid project
Date: January 2026
"""

import os
import subprocess
from pathlib import Path
import json

# Directories
BASE_DIR = Path("/data/ethylene")
STRUCT_DIR = BASE_DIR / "data/structures"
DOCKING_DIR = BASE_DIR / "results/docking/novel_candidates"
DOCKING_DIR.mkdir(parents=True, exist_ok=True)

# ============================================
# NOVEL DOCKING COMBINATIONS
# ============================================

NOVEL_DOCKING_PAIRS = [
    # Priority 1: Malonylated isoflavones (completely novel)
    {
        "receptor": "8E83",  # IFS homolog
        "receptor_name": "2HIS_MedicagoTruncatula",
        "ligand": "6-O-Malonylgenistin_CID5318568",
        "rationale": "IFS product analog - malonylated form never docked",
        "novelty": "★★★",
        "grid_center": None,  # Blind docking
        "grid_size": [40, 40, 40]
    },
    {
        "receptor": "8E83",
        "receptor_name": "2HIS_MedicagoTruncatula",
        "ligand": "6-O-Acetyldaidzin_CID14034712",
        "rationale": "Acetylated conjugate - never studied computationally",
        "novelty": "★★★",
        "grid_center": None,
        "grid_size": [40, 40, 40]
    },
    {
        "receptor": "8E83",
        "receptor_name": "2HIS_MedicagoTruncatula",
        "ligand": "6-O-Malonyldaidzin_CID5318574",
        "rationale": "Major isoflavone form in soybean - binding unknown",
        "novelty": "★★★",
        "grid_center": None,
        "grid_size": [40, 40, 40]
    },
    # Priority 2: CHI with novel substrates
    {
        "receptor": "1EYQ",
        "receptor_name": "CHI_MedicagoSativa",
        "ligand": "Liquiritigenin_CID114829",
        "rationale": "5-deoxy flavanone substrate comparison",
        "novelty": "★★",
        "grid_center": [12.5, 4.0, -7.5],  # Active site from crystal
        "grid_size": [25, 25, 25]
    },
    {
        "receptor": "1EYQ",
        "receptor_name": "CHI_MedicagoSativa",
        "ligand": "Daidzein_CID5281708",
        "rationale": "Product comparison - reverse binding",
        "novelty": "★★",
        "grid_center": [12.5, 4.0, -7.5],
        "grid_size": [25, 25, 25]
    },
    # Priority 3: 2-HID with isoflavone substrates
    {
        "receptor": "8EA1",
        "receptor_name": "2HID_PuerariaLobata",
        "ligand": "Daidzin_CID107971",
        "rationale": "Glycosylated isoflavone - dehydratase substrate analog",
        "novelty": "★★",
        "grid_center": None,
        "grid_size": [35, 35, 35]
    },
    {
        "receptor": "8EA1",
        "receptor_name": "2HID_PuerariaLobata",
        "ligand": "Genistin_CID5281377",
        "rationale": "Alternative glycoside substrate",
        "novelty": "★★",
        "grid_center": None,
        "grid_size": [35, 35, 35]
    },
]


def generate_vina_config(pair: dict, output_dir: Path) -> str:
    """Generate AutoDock Vina configuration file."""
    receptor_pdb = STRUCT_DIR / "pdb" / f"{pair['receptor']}.pdb"
    ligand_sdf = list((STRUCT_DIR / "ligands").glob(f"{pair['ligand']}*.sdf"))

    if not ligand_sdf:
        return None

    ligand_file = ligand_sdf[0]
    config_file = output_dir / f"config_{pair['receptor']}_{pair['ligand'][:20]}.txt"

    # Calculate grid center from PDB if not specified
    if pair['grid_center'] is None:
        center = calculate_protein_center(receptor_pdb)
    else:
        center = pair['grid_center']

    config_content = f"""# AutoDock Vina Configuration
# Receptor: {pair['receptor_name']}
# Ligand: {pair['ligand']}
# Novelty: {pair['novelty']}
# Rationale: {pair['rationale']}

receptor = {receptor_pdb}
ligand = {ligand_file}

center_x = {center[0]:.2f}
center_y = {center[1]:.2f}
center_z = {center[2]:.2f}

size_x = {pair['grid_size'][0]}
size_y = {pair['grid_size'][1]}
size_z = {pair['grid_size'][2]}

exhaustiveness = 32
num_modes = 20
energy_range = 3

out = {output_dir / f"docked_{pair['receptor']}_{pair['ligand'][:20]}.pdbqt"}
log = {output_dir / f"log_{pair['receptor']}_{pair['ligand'][:20]}.txt"}
"""

    config_file.write_text(config_content)
    return str(config_file)


def calculate_protein_center(pdb_file: Path) -> list:
    """Calculate geometric center of protein from PDB file."""
    coords = []

    if not pdb_file.exists():
        return [0.0, 0.0, 0.0]

    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                except ValueError:
                    continue

    if not coords:
        return [0.0, 0.0, 0.0]

    n = len(coords)
    center = [
        sum(c[0] for c in coords) / n,
        sum(c[1] for c in coords) / n,
        sum(c[2] for c in coords) / n
    ]
    return center


def check_vina_installed() -> bool:
    """Check if AutoDock Vina is installed."""
    try:
        result = subprocess.run(["vina", "--version"],
                                capture_output=True, text=True, timeout=10)
        return result.returncode == 0
    except (FileNotFoundError, subprocess.TimeoutExpired):
        return False


def main():
    print("=" * 60)
    print("Novel Isoflavonoid-Enzyme Docking Preparation")
    print("=" * 60)

    # Check Vina installation
    vina_installed = check_vina_installed()
    if not vina_installed:
        print("\n⚠ AutoDock Vina not found. Generating config files only.")
        print("  Install with: conda install -c conda-forge vina")

    # Generate configs for all pairs
    print(f"\nGenerating {len(NOVEL_DOCKING_PAIRS)} docking configurations...")

    configs = []
    for i, pair in enumerate(NOVEL_DOCKING_PAIRS, 1):
        print(f"\n{i}. {pair['ligand'][:30]} → {pair['receptor_name']}")
        print(f"   Novelty: {pair['novelty']} | {pair['rationale']}")

        config = generate_vina_config(pair, DOCKING_DIR)
        if config:
            configs.append(config)
            print(f"   ✓ Config: {Path(config).name}")
        else:
            print(f"   ✗ Ligand file not found")

    # Save summary
    summary = {
        "docking_pairs": NOVEL_DOCKING_PAIRS,
        "config_files": configs,
        "vina_installed": vina_installed,
        "output_directory": str(DOCKING_DIR)
    }

    summary_file = DOCKING_DIR / "docking_summary.json"
    with open(summary_file, "w") as f:
        json.dump(summary, f, indent=2)

    # Generate run script
    run_script = DOCKING_DIR / "run_all_docking.sh"
    with open(run_script, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("# Run all docking jobs\n\n")
        for config in configs:
            f.write(f"vina --config {config}\n")
    run_script.chmod(0o755)

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"\n총 {len(configs)}개 도킹 설정 생성 완료")
    print(f"\n신규성 분석:")
    novel_count = sum(1 for p in NOVEL_DOCKING_PAIRS if "★★★" in p['novelty'])
    partial_novel = sum(1 for p in NOVEL_DOCKING_PAIRS if "★★" in p['novelty'] and "★★★" not in p['novelty'])
    print(f"  ★★★ 완전 신규: {novel_count}개")
    print(f"  ★★  부분 신규: {partial_novel}개")

    print(f"\n출력 디렉토리: {DOCKING_DIR}")
    print(f"실행 스크립트: {run_script}")

    if not vina_installed:
        print("\n도킹 실행:")
        print("  conda activate docking")
        print(f"  bash {run_script}")


if __name__ == "__main__":
    main()
