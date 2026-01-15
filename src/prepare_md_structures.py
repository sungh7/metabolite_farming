#!/usr/bin/env python3
"""
Prepare structures for MD simulation binding validation.
Downloads enzyme structures from PDB/AlphaFold and ligands from PubChem.

Novel MD candidates (unreported in literature):
1. 6''-O-Malonylgenistin + GmIMaT3
2. 6''-O-Acetyldaidzin + GmIMaT1
3. Daidzein + GmIFS1

Author: Generated for ethylene-isoflavonoid project
Date: January 2026
"""

import os
import requests
import json
from pathlib import Path

# Output directories
BASE_DIR = Path("/data/ethylene/data/structures")
PDB_DIR = BASE_DIR / "pdb"
ALPHAFOLD_DIR = BASE_DIR / "alphafold"
LIGAND_DIR = BASE_DIR / "ligands"

# Create directories
for d in [PDB_DIR, ALPHAFOLD_DIR, LIGAND_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# ============================================
# ENZYME STRUCTURES
# ============================================

# Known PDB structures (experimentally determined)
PDB_STRUCTURES = {
    "1EYQ": {
        "name": "CHI_MedicagoSativa_Naringenin",
        "description": "Chalcone isomerase with naringenin (1.85 Å)",
        "organism": "Medicago sativa",
        "ligand": "Naringenin"
    },
    "8E83": {
        "name": "2HIS_MedicagoTruncatula",
        "description": "2-Hydroxyisoflavanone synthase (IFS homolog, 2.0 Å)",
        "organism": "Medicago truncatula",
        "ligand": None
    },
    "8EA1": {
        "name": "2HID_PuerariaLobata",
        "description": "2-Hydroxyisoflavanone dehydratase (2.4 Å)",
        "organism": "Pueraria lobata",
        "ligand": "p-nitrophenol"
    }
}

# AlphaFold models for soybean enzymes (may not all be available)
ALPHAFOLD_MODELS = {
    # Soybean isoflavone pathway enzymes - UniProt IDs to try
    "GmIFR": ["I1LHU6", "Q9SDZ1"],  # Isoflavone reductase
    "GmIFS": ["Q9FRV8"],  # Isoflavone synthase (from legume)
    "GmCHS": ["P24826"],  # Chalcone synthase
    "GmCHI": ["O24520"],  # Chalcone isomerase
}

# ============================================
# LIGAND STRUCTURES (PubChem)
# ============================================

LIGANDS = {
    # Novel candidates - malonylated isoflavones (UNREPORTED in MD studies)
    "6-O-Malonyldaidzin": {
        "cid": 5318574,
        "chebi": "CHEBI:80371",
        "kegg": None,
        "novelty": "★★★ NOVEL - No MD simulation reported"
    },
    "6-O-Malonylgenistin": {
        "cid": 5318568,
        "chebi": "CHEBI:80372",
        "kegg": None,
        "novelty": "★★★ NOVEL - No MD simulation reported"
    },
    "6-O-Acetyldaidzin": {
        "cid": 14034712,
        "chebi": "CHEBI:133395",
        "kegg": None,
        "novelty": "★★★ NOVEL - No MD simulation reported"
    },
    "6-O-Acetylgenistin": {
        "cid": 5320413,
        "chebi": None,
        "kegg": None,
        "novelty": "★★★ NOVEL - No MD simulation reported"
    },
    # Standard isoflavones (for comparison)
    "Daidzein": {
        "cid": 5281708,
        "chebi": "CHEBI:28197",
        "kegg": "C10208",
        "novelty": "★ Known - binding to IFS studied"
    },
    "Genistein": {
        "cid": 5280961,
        "chebi": "CHEBI:74224",
        "kegg": "C06563",
        "novelty": "★ Known - binding to IFS studied"
    },
    "Daidzin": {
        "cid": 107971,
        "chebi": "CHEBI:4307",
        "kegg": "C00159",
        "novelty": "★★ Km reported, binding mode unknown"
    },
    "Genistin": {
        "cid": 5281377,
        "chebi": "CHEBI:28159",
        "kegg": "C09126",
        "novelty": "★★ Km reported, binding mode unknown"
    },
    "Naringenin": {
        "cid": 932,
        "chebi": "CHEBI:17846",
        "kegg": "C00509",
        "novelty": "★ Known - CHI crystal structure"
    },
    "Liquiritigenin": {
        "cid": 114829,
        "chebi": "CHEBI:28159",
        "kegg": "C09762",
        "novelty": "★★ IFS substrate"
    }
}


def download_pdb(pdb_id: str, output_dir: Path) -> bool:
    """Download PDB structure from RCSB."""
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    output_file = output_dir / f"{pdb_id}.pdb"

    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            output_file.write_text(response.text)
            print(f"  ✓ Downloaded {pdb_id}.pdb ({len(response.text)} bytes)")
            return True
        else:
            print(f"  ✗ Failed to download {pdb_id}: HTTP {response.status_code}")
            return False
    except Exception as e:
        print(f"  ✗ Error downloading {pdb_id}: {e}")
        return False


def download_alphafold(uniprot_id: str, name: str, output_dir: Path) -> bool:
    """Download AlphaFold model from EBI."""
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    output_file = output_dir / f"{name}_{uniprot_id}.pdb"

    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 200 and not response.text.startswith("<?xml"):
            output_file.write_text(response.text)
            print(f"  ✓ Downloaded AlphaFold {uniprot_id} ({len(response.text)} bytes)")
            return True
        else:
            print(f"  ✗ AlphaFold model not available for {uniprot_id}")
            return False
    except Exception as e:
        print(f"  ✗ Error downloading {uniprot_id}: {e}")
        return False


def download_ligand_sdf(cid: int, name: str, output_dir: Path) -> bool:
    """Download 3D ligand structure from PubChem in SDF format."""
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF?record_type=3d"
    output_file = output_dir / f"{name.replace(' ', '_')}_CID{cid}.sdf"

    try:
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            output_file.write_text(response.text)
            print(f"  ✓ Downloaded {name} (CID:{cid})")
            return True
        else:
            # Try 2D if 3D not available
            url_2d = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/SDF"
            response = requests.get(url_2d, timeout=30)
            if response.status_code == 200:
                output_file.write_text(response.text)
                print(f"  ✓ Downloaded {name} (CID:{cid}) [2D - needs 3D conversion]")
                return True
            print(f"  ✗ Failed to download {name}: HTTP {response.status_code}")
            return False
    except Exception as e:
        print(f"  ✗ Error downloading {name}: {e}")
        return False


def main():
    print("=" * 60)
    print("MD Simulation Structure Preparation")
    print("=" * 60)

    # 1. Download PDB structures
    print("\n[1/3] Downloading PDB structures (experimentally determined)...")
    for pdb_id, info in PDB_STRUCTURES.items():
        print(f"\n  {pdb_id}: {info['description']}")
        download_pdb(pdb_id, PDB_DIR)

    # 2. Download AlphaFold models
    print("\n[2/3] Downloading AlphaFold models (computational predictions)...")
    for enzyme, uniprot_ids in ALPHAFOLD_MODELS.items():
        for uniprot_id in uniprot_ids:
            download_alphafold(uniprot_id, enzyme, ALPHAFOLD_DIR)

    # 3. Download ligand structures
    print("\n[3/3] Downloading ligand structures from PubChem...")
    novel_count = 0
    for name, info in LIGANDS.items():
        print(f"\n  {name} ({info['novelty']})")
        if download_ligand_sdf(info["cid"], name, LIGAND_DIR):
            if "NOVEL" in info["novelty"]:
                novel_count += 1

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"\nStructures saved to: {BASE_DIR}")
    print(f"  - PDB structures: {PDB_DIR}")
    print(f"  - AlphaFold models: {ALPHAFOLD_DIR}")
    print(f"  - Ligand structures: {LIGAND_DIR}")

    print(f"\n신규 MD 후보 (미보고): {novel_count}개")
    print("\n추천 도킹 조합 (완전 신규):")
    print("  1. 6-O-Malonylgenistin + GmIMaT (말로닐전이효소)")
    print("  2. 6-O-Acetyldaidzin + GmIMaT (말로닐전이효소)")
    print("  3. Daidzin/Genistin + GmIMaT (기질, Km 보고됨)")

    # Save metadata
    metadata = {
        "pdb_structures": PDB_STRUCTURES,
        "alphafold_models": ALPHAFOLD_MODELS,
        "ligands": LIGANDS,
        "novel_candidates": [
            {
                "ligand": "6-O-Malonylgenistin",
                "enzyme": "GmIMaT3",
                "rationale": "Product of malonylation reaction, binding mode unknown"
            },
            {
                "ligand": "6-O-Acetyldaidzin",
                "enzyme": "GmIMaT1",
                "rationale": "Acetylated conjugate, never studied computationally"
            },
            {
                "ligand": "Daidzin",
                "enzyme": "GmIMaT1",
                "rationale": "Substrate with known Km (36.28 μM), binding mode unknown"
            }
        ]
    }

    with open(BASE_DIR / "structure_metadata.json", "w") as f:
        json.dump(metadata, f, indent=2)

    print(f"\nMetadata saved to: {BASE_DIR / 'structure_metadata.json'}")


if __name__ == "__main__":
    main()
