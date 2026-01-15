#!/usr/bin/env python3
"""
Run docking for all novel isoflavonoid-enzyme candidates.
"""

import subprocess
import os
from pathlib import Path

BASE_DIR = Path("/data/ethylene")
STRUCT_DIR = BASE_DIR / "data/structures"
OUTPUT_DIR = BASE_DIR / "results/docking/all_candidates"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Docking pairs to run
DOCKING_PAIRS = [
    # Novel candidates (★★★)
    {"receptor": "8E83", "ligand": "6-O-Acetyldaidzin_CID14034712", "novelty": "★★★"},
    {"receptor": "8E83", "ligand": "6-O-Malonyldaidzin_CID5318574", "novelty": "★★★"},
    {"receptor": "8E83", "ligand": "6-O-Acetylgenistin_CID5320413", "novelty": "★★★"},
    # Partial novel (★★)
    {"receptor": "1EYQ", "ligand": "Liquiritigenin_CID114829", "novelty": "★★"},
    {"receptor": "1EYQ", "ligand": "Daidzein_CID5281708", "novelty": "★★"},
    {"receptor": "8EA1", "ligand": "Daidzin_CID107971", "novelty": "★★"},
    {"receptor": "8EA1", "ligand": "Genistin_CID5281377", "novelty": "★★"},
]

def get_protein_center(pdb_file):
    """Calculate geometric center of protein."""
    coords = []
    with open(pdb_file) as f:
        for line in f:
            if line.startswith("ATOM"):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                except:
                    pass
    if not coords:
        return [0, 0, 0]
    n = len(coords)
    return [sum(c[i] for c in coords) / n for i in range(3)]

def run_docking(pair):
    """Run docking for a single pair."""
    receptor_pdb = STRUCT_DIR / "pdb" / f"{pair['receptor']}.pdb"
    ligand_files = list((STRUCT_DIR / "ligands").glob(f"{pair['ligand']}*.sdf"))

    if not ligand_files:
        print(f"  ✗ Ligand not found: {pair['ligand']}")
        return None

    ligand_sdf = ligand_files[0]
    pair_name = f"{pair['receptor']}_{pair['ligand'][:25]}"
    work_dir = OUTPUT_DIR / pair_name
    work_dir.mkdir(exist_ok=True)

    # Prepare receptor
    receptor_pdbqt = work_dir / "receptor.pdbqt"
    subprocess.run(
        f"obabel {receptor_pdb} -O {receptor_pdbqt} -xr",
        shell=True, capture_output=True
    )

    # Prepare ligand
    ligand_pdbqt = work_dir / "ligand.pdbqt"
    subprocess.run(
        f"obabel {ligand_sdf} -O {ligand_pdbqt} --gen3d",
        shell=True, capture_output=True
    )

    # Get center
    center = get_protein_center(receptor_pdb)

    # Create config
    config_file = work_dir / "config.txt"
    config_file.write_text(f"""receptor = {receptor_pdbqt}
ligand = {ligand_pdbqt}
center_x = {center[0]:.2f}
center_y = {center[1]:.2f}
center_z = {center[2]:.2f}
size_x = 35
size_y = 35
size_z = 35
exhaustiveness = 16
num_modes = 5
energy_range = 3
out = {work_dir}/docked.pdbqt
""")

    # Run Vina
    result = subprocess.run(
        f"vina --config {config_file}",
        shell=True, capture_output=True, text=True
    )

    # Parse result
    output = result.stdout + result.stderr
    best_affinity = None
    for line in output.split('\n'):
        if line.strip().startswith('1 '):
            parts = line.split()
            if len(parts) >= 2:
                try:
                    best_affinity = float(parts[1])
                except:
                    pass

    # Save log
    (work_dir / "vina.log").write_text(output)

    return best_affinity

def main():
    print("=" * 60)
    print("Running All Docking Jobs")
    print("=" * 60)

    results = []

    for i, pair in enumerate(DOCKING_PAIRS, 1):
        print(f"\n[{i}/{len(DOCKING_PAIRS)}] {pair['ligand'][:30]}")
        print(f"    Receptor: {pair['receptor']} | Novelty: {pair['novelty']}")

        affinity = run_docking(pair)

        if affinity:
            status = "✓ Strong" if affinity < -7.0 else "○ Moderate" if affinity < -5.0 else "△ Weak"
            print(f"    Result: {affinity:.3f} kcal/mol {status}")
            results.append({
                "receptor": pair['receptor'],
                "ligand": pair['ligand'],
                "novelty": pair['novelty'],
                "affinity": affinity
            })
        else:
            print(f"    Result: Failed")

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"\n{'Ligand':<35} {'Receptor':<8} {'Affinity':>10} {'Novelty'}")
    print("-" * 60)

    for r in sorted(results, key=lambda x: x['affinity']):
        lig = r['ligand'][:33]
        print(f"{lig:<35} {r['receptor']:<8} {r['affinity']:>10.3f} {r['novelty']}")

    # Save results
    import json
    with open(OUTPUT_DIR / "all_results.json", "w") as f:
        json.dump(results, f, indent=2)

    print(f"\nResults saved to: {OUTPUT_DIR}")

if __name__ == "__main__":
    main()
