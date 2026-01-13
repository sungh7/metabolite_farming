"""
Automated Docking Pipeline for GNN Predictions
Requires: conda install -c conda-forge openbabel vina rdkit

Usage:
    python src/run_docking.py --input results/gnn/top1_enzymes_for_docking.csv --outdir results/docking
"""

import os
import argparse
import pandas as pd
import requests
import subprocess
from pathlib import Path
from tqdm import tqdm
import time

# Absolute paths for binaries
VINA_PATH = "/home/csh/mambaforge/bin/vina"
OBABEL_PATH = "/home/csh/mambaforge/bin/obabel"

def download_alphafold_pdb(uniprot_id, out_dir):
    """Download AlphaFold structure via API."""
    out_path = out_dir / f"{uniprot_id}.pdb"
    if out_path.exists():
        return out_path
    
    # 1. Get Metadata from API
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
    try:
        resp = requests.get(api_url, timeout=10)
        if resp.status_code != 200:
            return None
        
        data = resp.json()
        if not data or not isinstance(data, list):
            return None
            
        # Get the first entry (usually the correct one)
        pdb_url = data[0].get('pdbUrl')
        if not pdb_url:
            return None
            
        # 2. Download PDB
        pdb_resp = requests.get(pdb_url, timeout=30)
        if pdb_resp.status_code == 200:
            with open(out_path, 'wb') as f:
                f.write(pdb_resp.content)
            return out_path
            
    except Exception as e:
        print(f"Error downloading {uniprot_id}: {e}")
    
    return None

def download_kegg_mol(kegg_id, out_dir):
    """Download MOL file from KEGG."""
    url = f"https://rest.kegg.jp/get/{kegg_id}/mol"
    out_path = out_dir / f"{kegg_id}.mol"
    
    if out_path.exists():
        return out_path
        
    try:
        resp = requests.get(url, timeout=10)
        if resp.status_code == 200 and "NO DATA" not in resp.text:
            with open(out_path, 'wb') as f:
                f.write(resp.content)
            return out_path
    except Exception as e:
        print(f"Error downloading {kegg_id}: {e}")
    return None

def prepare_receptor(pdb_file, out_pdbqt):
    """Convert PDB to PDBQT using OpenBabel."""
    cmd = [OBABEL_PATH, str(pdb_file), "-xr", "-O", str(out_pdbqt), "-h", "--partialcharge", "gasteiger"]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        print(f"Error preparing receptor {pdb_file} (check if obabel is installed)")
        return False

def prepare_ligand(mol_file, out_pdbqt):
    """Convert MOL to PDBQT using OpenBabel."""
    cmd = [OBABEL_PATH, str(mol_file), "-O", str(out_pdbqt), "-h", "--gen3d", "--partialcharge", "gasteiger"]
    try:
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        print(f"Error preparing ligand {mol_file}")
        return False

def run_vina(receptor_pdbqt, ligand_pdbqt, out_pdbqt, log_file):
    """Run AutoDock Vina (Blind Docking)."""
    cmd = [
        VINA_PATH, 
        "--receptor", str(receptor_pdbqt), 
        "--ligand", str(ligand_pdbqt),
        "--out", str(out_pdbqt),
        "--center_x", "0", "--center_y", "0", "--center_z", "0",
        "--size_x", "80", "--size_y", "80", "--size_z", "80",
        "--cpu", "4", "--exhaustiveness", "8"
    ]
    
    try:
        with open(log_file, 'w') as f:
            subprocess.run(cmd, check=True, stdout=f, stderr=subprocess.PIPE)
        
        # Parse affinity from log
        best_affinity = None
        if os.path.exists(log_file):
            with open(log_file, 'r') as f:
                for line in f:
                    if line.strip().startswith("1"):
                        parts = line.split()
                        if len(parts) >= 2:
                            try:
                                best_affinity = float(parts[1])
                                break
                            except ValueError: pass
        return best_affinity
        
    except subprocess.CalledProcessError as e:
        print(f"Error running Vina: {e.stderr.decode()}")
        return None
    except FileNotFoundError:
        print(f"Error: Vina binary not found at {VINA_PATH}")
        return None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='CSV with metabolite_kegg, uniprot_id')
    parser.add_argument('--outdir', required=True)
    args = parser.parse_args()
    
    # Check dependencies
    if not os.path.exists(VINA_PATH) or not os.path.exists(OBABEL_PATH):
        print(f"Error: Tools not found at {VINA_PATH} or {OBABEL_PATH}")
        return

    df = pd.read_csv(args.input)
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)
    
    results = []
    
    print(f"Starting docking for {len(df)} pairs...")
    
    for _, row in tqdm(df.iterrows(), total=len(df)):
        met_id = row['metabolite_kegg']
        uniprot = row['uniprot_id']
        
        if uniprot == "Unknown": continue
        
        # Dirs
        pair_dir = out_dir / f"{met_id}_{uniprot}"
        pair_dir.mkdir(exist_ok=True)
        
        # 1. Structures
        receptor_pdb = download_alphafold_pdb(uniprot, pair_dir)
        ligand_mol = download_kegg_mol(met_id, pair_dir)
        
        if not receptor_pdb or not ligand_mol:
            continue
            
        # 2. Prepare
        receptor_pdbqt = pair_dir / "receptor.pdbqt"
        ligand_pdbqt = pair_dir / "ligand.pdbqt"
        
        if not prepare_receptor(receptor_pdb, receptor_pdbqt): continue
        if not prepare_ligand(ligand_mol, ligand_pdbqt): continue
        
        # 3. Dock
        dock_out = pair_dir / "docked.pdbqt"
        dock_log = pair_dir / "docking.log"
        
        affinity = run_vina(receptor_pdbqt, ligand_pdbqt, dock_out, dock_log)
        
        if affinity:
            results.append({
                'metabolite': met_id,
                'enzyme': uniprot,
                'affinity': affinity,
                'gnn_score': row.get('score', 0)
            })
            
    # Save results
    res_df = pd.DataFrame(results)
    res_df.to_csv(out_dir / "docking_summary.csv", index=False)
    print(f"Docking complete. Results saved to {out_dir / 'docking_summary.csv'}")

if __name__ == "__main__":
    main()
