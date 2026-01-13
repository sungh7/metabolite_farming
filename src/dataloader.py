import gzip
import pandas as pd
import os
from tqdm import tqdm

class StringDBLoader:
    def __init__(self, raw_dir='data/raw'):
        self.raw_dir = raw_dir
        self.protein_info_path = os.path.join(raw_dir, '3847.protein.info.v12.0.txt.gz')
        self.links_path = os.path.join(raw_dir, '3847.protein.links.full.v12.0.txt.gz')
        self.protein_map = {} # string_id -> preferred_name
        self.node_to_idx = {} # string_id -> int index
        self.idx_to_node = {} # int index -> string_id

    def load_protein_info(self):
        """
        Parses the protein info file.
        Returns:
            dict: string_id -> preferred_name
        """
        print(f"Loading protein info from {self.protein_info_path}...")
        self.protein_map = {}
        
        with gzip.open(self.protein_info_path, 'rt') as f:
            next(f) # Skip header
            for line in tqdm(f, desc="Parsing Info"):
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    string_id = parts[0]
                    preferred_name = parts[1]
                    self.protein_map[string_id] = preferred_name
        
        # Create an index mapping for all found proteins
        for idx, string_id in enumerate(self.protein_map.keys()):
            self.node_to_idx[string_id] = idx
            self.idx_to_node[idx] = string_id
            
        print(f"Loaded {len(self.protein_map)} proteins.")
        return self.protein_map

    def load_interactions(self, threshold=700, strict=False):
        """
        Parses the links file and filters by combined_score.
        Args:
            threshold (int): Minimum combined_score (0-1000). Default 700 (High confidence).
            strict (bool): If True, exclude 'textmining' channel and recalculate combined score.
        Returns:
            list: List of (source_idx, target_idx) tuples
        """
        print(f"Loading interactions from {self.links_path} with threshold={threshold}, strict={strict}...")
        edges = []
        
        if not self.node_to_idx:
            print("Warning: Protein info not loaded. Loading now to ensure index mapping.")
            self.load_protein_info()
            
        with gzip.open(self.links_path, 'rt') as f:
            header_line = next(f).strip()
            # Handle space-separated header
            header = header_line.split(' ')
            
            try:
                p1_idx = header.index('protein1')
                p2_idx = header.index('protein2')
                
                # Channels for strict mode
                if strict:
                    # Exclude textmining. Include others.
                    channels = ['neighborhood', 'fusion', 'cooccurence', 'coexpression', 'experimental', 'database']
                    channel_indices = []
                    for ch in channels:
                        if ch in header:
                            channel_indices.append(header.index(ch))
                        else:
                            # If channel missing, treat as 0
                            pass
                else:
                    # Normal mode: use combined_score
                    score_idx = header.index('combined_score')
                    
            except ValueError as e:
                print(f"Error parsing header: {e}")
                return []

            count = 0
            for line in tqdm(f, desc="Parsing Links"):
                parts = line.strip().split(' ')
                try:
                    if strict:
                        # Recalculate score
                        # Formula: 1 - Prod(1 - p_i)
                        # p_i = score / 1000.0
                        combined_prob = 1.0
                        
                        # Use Prior Correction? STRING scores in file are already prior-corrected and scaled.
                        # Standard combination is straightforward independence assumption.
                        
                        not_combined = 1.0
                        has_evidence = False
                        
                        for ch_idx in channel_indices:
                            # Check if part exists (sometimes lines are malformed?)
                            if ch_idx < len(parts):
                                val = int(parts[ch_idx])
                                if val > 0:
                                    prob = val / 1000.0
                                    not_combined *= (1.0 - prob)
                                    has_evidence = True
                        
                        if has_evidence:
                            final_prob = 1.0 - not_combined
                            score = int(final_prob * 1000)
                        else:
                            score = 0
                    else:
                        score = int(parts[score_idx])

                    if score >= threshold:
                        p1 = parts[p1_idx]
                        p2 = parts[p2_idx]
                        
                        if p1 in self.node_to_idx and p2 in self.node_to_idx:
                            idx1 = self.node_to_idx[p1]
                            idx2 = self.node_to_idx[p2]
                            edges.append((idx1, idx2))
                            count += 1
                            
                except (ValueError, IndexError):
                    continue
                    
        print(f"Loaded {count} interactions (score >= {threshold}, strict={strict}).")
        return edges

    def get_node_name(self, idx):
        string_id = self.idx_to_node.get(idx)
        return self.protein_map.get(string_id, "Unknown")

if __name__ == "__main__":
    loader = StringDBLoader()
    info = loader.load_protein_info()
    links = loader.load_interactions(threshold=900) # Testing with very high threshold
    
    # Print some stats
    print(f"Total Proteins: {len(info)}")
    print(f"Total Edges (Conf >= 900): {len(links)}")
    
    # Show first 5 edges
    print("Sample Edges:")
    for u, v in links[:5]:
        print(f"{loader.get_node_name(u)} <--> {loader.get_node_name(v)}")
