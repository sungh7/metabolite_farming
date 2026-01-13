import torch
from src.model import HGT, LinkPredictor
from src.dataloader import StringDBLoader
import gzip

def find_valid_tf_pair():
    # 1. Load mappings
    loader = StringDBLoader()
    protein_map = loader.load_protein_info() # string_id -> name
    
    # 2. Identify Top TFs (ERF, WRKY, MYB)
    valid_tfs = []
    for sid, name in protein_map.items():
        n = name.lower()
        if 'erf' in n or 'wrky' in n or 'myb' in n:
            valid_tfs.append((sid, name))
            
    print(f"Found {len(valid_tfs)} valid TFs by name.")
    
    # 3. Load Model & Data
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    data = torch.load('data/processed/strict_bipartite.pt', map_location=device)
    
    # Reconstruct ID map
    sorted_string_ids = [None] * (data['Enzyme'].num_nodes + data['TF'].num_nodes + data['Signaling'].num_nodes + data['Protein'].num_nodes) 
    # This is hard because graph_builder didn't save the mapping.
    # But we can infer valid TFs if we just replicate the builder's logic slightly or trust the 'TF' node indices if we knew which one is which.
    # Actually, inference.py logic handles this. Let's just use the saved model and predict for *known* TF indices.
    
    # ... Wait, identifying which row in x_dict['TF'] corresponds to which string_id is tricky without the mapping list from builder.
    # graph_builder.py prints the mapping logic but doesn't save it.
    # CRITICAL FLASHLIGHT: graph_builder.py *re-sorts* IDs based on type!
    # "global_to_local[global_idx] = (ntype, local_idx)"
    
    # I need to rebuild the exact same mapping to be safe.
    # Fastest way: modify inference.py to filter candidates based on NAME *before* ranking.
    pass

if __name__ == "__main__":
    find_valid_tf_pair()
