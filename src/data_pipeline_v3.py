"""
Data Pipeline v3 for Ethylene-Isoflavonoid GNN Project

학술적으로 방어 가능한 그래프 설계 (Academically Defensible Design)

Key Changes from v2:
1. Tier-P → Metabolite↔Metabolite (rxn_neighbor) edges based on reaction distance
2. Supervision uses Tier-R only (Enzyme→Metabolite)
3. Omics features stored in data, learnable embeddings in model
4. TF domain one-hot features
5. Currency metabolite validation (key + neighbor)
6. tier_r_lookup for filtered evaluation

Usage:
    python src/data_pipeline_v3.py --output data/processed/graph_v3.pt
"""

import torch
from torch_geometric.data import HeteroData
import pandas as pd
import numpy as np
import gzip
from pathlib import Path
from tqdm import tqdm
from collections import defaultdict
from typing import Optional, Dict, Set, Tuple, List

import sys
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.config import DATA_PATHS, GRAPH_CONFIG, V3_CONFIG
from src.utils.seed import set_seed


# =============================================================================
# Currency Metabolites (ATP, GTP, H2O, etc.)
# =============================================================================
CURRENCY_METABOLITES = {
    'C00001',  # H2O
    'C00002',  # ATP
    'C00003',  # NAD+
    'C00004',  # NADH
    'C00005',  # NADPH
    'C00006',  # NADP+
    'C00008',  # ADP
    'C00009',  # Orthophosphate
    'C00010',  # CoA
    'C00011',  # CO2
    'C00013',  # Diphosphate
    'C00014',  # Ammonia
    'C00020',  # AMP
    'C00027',  # H2O2
    'C00044',  # GTP
    'C00080',  # H+
}


def load_string_ncbi_mapping(aliases_path: Path) -> Dict[str, str]:
    """Load NCBI Gene ID to STRING protein mapping from aliases file."""
    if not aliases_path.exists():
        print(f"Warning: {aliases_path} not found. Falling back to UniProt only.")
        return {}

    ncbi_to_string = {}
    print(f"Loading STRING aliases from {aliases_path}...")

    with gzip.open(aliases_path, 'rt') as f:
        next(f)  # Skip header
        for line in tqdm(f, desc="Parsing aliases"):
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                string_id = parts[0]
                alias = parts[1]
                source = parts[2]

                if 'GeneID' in source or 'NCBI' in source or 'Entrez' in source:
                    ncbi_to_string[alias] = string_id

    print(f"NCBI to STRING mappings: {len(ncbi_to_string)}")
    return ncbi_to_string


def build_reaction_graph(
    edges_df: pd.DataFrame,
    valid_mets: Set[str],
    currency_mets: Set[str] = CURRENCY_METABOLITES,
) -> Dict[str, Set[str]]:
    """
    Build reaction adjacency graph for Metabolite↔Metabolite edges.

    Two metabolites are connected if they participate in the same reaction.
    Currency metabolites are excluded from both keys and neighbors.

    Args:
        edges_df: DataFrame with columns [reaction_id, metabolite_id]
        valid_mets: Set of valid metabolite IDs to include
        currency_mets: Set of currency metabolites to exclude

    Returns:
        Dictionary mapping metabolite_id → Set of neighbor metabolite_ids
    """
    # Group metabolites by reaction
    reaction_to_mets: Dict[str, Set[str]] = defaultdict(set)

    for _, row in edges_df.iterrows():
        reaction_id = str(row['reaction_id'])
        met_id = str(row['metabolite_id'])

        # Exclude currency metabolites
        if met_id in currency_mets:
            continue

        # Only include valid metabolites
        if met_id in valid_mets:
            reaction_to_mets[reaction_id].add(met_id)

    # Build adjacency: metabolites in same reaction are neighbors
    reaction_adj: Dict[str, Set[str]] = defaultdict(set)

    for reaction_id, mets in reaction_to_mets.items():
        mets_list = list(mets)
        for i, met_i in enumerate(mets_list):
            for met_j in mets_list[i + 1:]:
                reaction_adj[met_i].add(met_j)
                reaction_adj[met_j].add(met_i)

    print(f"Reaction graph: {len(reaction_adj)} metabolites with neighbors")
    return dict(reaction_adj)


def build_rxn_neighbor_edges(
    reaction_adj: Dict[str, Set[str]],
    met_to_idx: Dict[str, int],
    max_dist: int = 2,
) -> torch.Tensor:
    """
    Build Metabolite↔Metabolite edges based on reaction distance.

    Distance 1: Same reaction
    Distance 2: Share a common neighbor (via BFS)

    rxn_neighbor represents reaction neighborhood proximity (undirected).
    This is NOT a directional flux model.

    Args:
        reaction_adj: Adjacency dict from build_reaction_graph
        met_to_idx: Metabolite ID to index mapping
        max_dist: Maximum reaction distance (1 or 2)

    Returns:
        Edge index tensor [2, num_edges]
    """
    edges = set()

    for met_id, neighbors in reaction_adj.items():
        if met_id not in met_to_idx:
            continue

        met_idx = met_to_idx[met_id]

        # Distance 1: Direct neighbors
        for neighbor_id in neighbors:
            if neighbor_id in met_to_idx:
                neighbor_idx = met_to_idx[neighbor_id]
                # Add both directions for undirected edge
                if met_idx != neighbor_idx:  # No self-loops
                    edges.add((min(met_idx, neighbor_idx), max(met_idx, neighbor_idx)))

        # Distance 2: Neighbors of neighbors (if max_dist >= 2)
        if max_dist >= 2:
            for neighbor_id in neighbors:
                if neighbor_id not in reaction_adj:
                    continue
                for hop2_id in reaction_adj[neighbor_id]:
                    if hop2_id in met_to_idx and hop2_id != met_id:
                        hop2_idx = met_to_idx[hop2_id]
                        if met_idx != hop2_idx:
                            edges.add((min(met_idx, hop2_idx), max(met_idx, hop2_idx)))

    # Convert to bidirectional edge index
    edge_list = list(edges)
    if not edge_list:
        return torch.zeros((2, 0), dtype=torch.long)

    # Create bidirectional edges (undirected graph)
    src = [e[0] for e in edge_list] + [e[1] for e in edge_list]
    dst = [e[1] for e in edge_list] + [e[0] for e in edge_list]

    return torch.tensor([src, dst], dtype=torch.long)


def validate_no_currency(
    reaction_adj: Dict[str, Set[str]],
    met_to_idx: Dict[str, int],
    rxn_neighbor_edge_index: torch.Tensor,
    currency_mets: Set[str] = CURRENCY_METABOLITES,
):
    """
    Validate that currency metabolites are excluded from both keys and neighbors.

    Args:
        reaction_adj: Reaction adjacency dictionary
        met_to_idx: Metabolite ID to index mapping
        rxn_neighbor_edge_index: rxn_neighbor edge index
        currency_mets: Set of currency metabolite IDs
    """
    idx_to_met = {v: k for k, v in met_to_idx.items()}

    # Check reaction_adj keys and neighbors
    for met_id, neighbors in reaction_adj.items():
        assert met_id not in currency_mets, f"Currency {met_id} in reaction_adj keys"
        for neighbor_id in neighbors:
            assert neighbor_id not in currency_mets, \
                f"Currency {neighbor_id} in neighbors of {met_id}"

    # Check rxn_neighbor edges
    for idx in rxn_neighbor_edge_index.flatten().tolist():
        if idx in idx_to_met:
            met_id = idx_to_met[idx]
            assert met_id not in currency_mets, \
                f"Currency {met_id} (idx={idx}) in rxn_neighbor edges"

    print("Currency validation passed: no currency metabolites in graph")


def zscore_normalize(features: np.ndarray) -> torch.Tensor:
    """
    Z-score normalize features, handling NaN and constant columns.

    Args:
        features: NumPy array of shape [num_nodes, num_features]

    Returns:
        Normalized tensor
    """
    # Replace NaN with 0
    features = np.nan_to_num(features, nan=0.0)

    # Z-score normalize each column
    result = np.zeros_like(features)
    for i in range(features.shape[1]):
        col = features[:, i]
        std = np.std(col)
        if std > 1e-8:
            result[:, i] = (col - np.mean(col)) / std
        else:
            result[:, i] = 0.0  # Constant column

    return torch.tensor(result, dtype=torch.float32)


def build_metabolite_omics_features(
    met_list: List[str],
    metabolomics_path: Path,
    kegg_metabolites_path: Path,
) -> torch.Tensor:
    """
    Build metabolite omics features: [log2fc, -log10(pvalue), n_pathways].

    Args:
        met_list: List of metabolite KEGG IDs
        metabolomics_path: Path to metabolomics differential data
        kegg_metabolites_path: Path to KEGG metabolites with pathway info

    Returns:
        Tensor of shape [num_mets, 3]
    """
    n_mets = len(met_list)
    features = np.zeros((n_mets, 3))

    # Load KEGG metabolites (has pathway info)
    kegg_df = pd.read_csv(kegg_metabolites_path) if kegg_metabolites_path.exists() else pd.DataFrame()

    # Create lookup from KEGG
    kegg_lookup = {}
    if len(kegg_df) > 0 and 'compound_id' in kegg_df.columns:
        for _, row in kegg_df.iterrows():
            cid = row['compound_id']
            kegg_lookup[cid] = {
                'log2fc': row.get('log2fc', 0.0),
                'pvalue': row.get('pvalue', 1.0),
                'n_pathways': row.get('n_pathways', 0),
            }

    # Load metabolomics if available
    met_diff_lookup = {}
    if metabolomics_path.exists():
        met_df = pd.read_csv(metabolomics_path)
        if 'KEGG' in met_df.columns:
            for _, row in met_df.iterrows():
                kegg_id = row.get('KEGG')
                if pd.notna(kegg_id):
                    met_diff_lookup[kegg_id] = {
                        'log2fc': row.get('Log2FC', 0.0),
                        'pvalue': row.get('P_Value', 1.0),
                    }

    # Fill features
    for i, met_id in enumerate(met_list):
        log2fc = 0.0
        pvalue = 1.0
        n_pathways = 0

        # Try KEGG lookup first
        if met_id in kegg_lookup:
            info = kegg_lookup[met_id]
            log2fc = info.get('log2fc', 0.0) or 0.0
            pvalue = info.get('pvalue', 1.0) or 1.0
            n_pathways = info.get('n_pathways', 0) or 0

        # Override with metabolomics if available
        if met_id in met_diff_lookup:
            info = met_diff_lookup[met_id]
            log2fc = info.get('log2fc', log2fc) or log2fc
            pvalue = info.get('pvalue', pvalue) or pvalue

        # Convert pvalue to -log10(pvalue) for better scaling
        neg_log_pvalue = -np.log10(max(pvalue, 1e-300))

        features[i] = [log2fc, neg_log_pvalue, n_pathways]

    return zscore_normalize(features)


def build_enzyme_omics_features(
    num_enzymes: int,
    enzyme_df: pd.DataFrame,
    proteomics_path: Path,
) -> torch.Tensor:
    """
    Build enzyme omics features: [log2fc, -log10(pvalue), abundance].

    Args:
        num_enzymes: Number of enzyme nodes
        enzyme_df: DataFrame with enzyme string mapping
        proteomics_path: Path to proteomics differential data

    Returns:
        Tensor of shape [num_enzymes, 3]
    """
    features = np.zeros((num_enzymes, 3))

    if not proteomics_path.exists():
        print(f"Warning: Proteomics data not found at {proteomics_path}")
        return zscore_normalize(features)

    prot_df = pd.read_csv(proteomics_path)

    # Create lookup by protein ID
    prot_lookup = {}
    for _, row in prot_df.iterrows():
        protein_ids = str(row.get('Protein IDs', ''))
        for pid in protein_ids.split(';'):
            pid = pid.strip()
            if pid:
                prot_lookup[pid] = {
                    'log2fc': row.get('Log2FC', 0.0),
                    'pvalue': row.get('P_Value', 1.0),
                    'abundance': row.get('Mean_Ethylene', 0.0),
                }

    # Map enzymes to proteomics
    if 'uniprot_id' in enzyme_df.columns:
        for _, row in enzyme_df.iterrows():
            enz_idx = row['enzyme_idx']
            uniprot = row.get('uniprot_id', '')

            if enz_idx < num_enzymes and uniprot in prot_lookup:
                info = prot_lookup[uniprot]
                log2fc = info.get('log2fc', 0.0) or 0.0
                pvalue = info.get('pvalue', 1.0) or 1.0
                abundance = info.get('abundance', 0.0) or 0.0

                neg_log_pvalue = -np.log10(max(pvalue, 1e-300))
                features[enz_idx] = [log2fc, neg_log_pvalue, np.log1p(abundance)]

    return zscore_normalize(features)


def build_tf_domain_onehot(
    tf_names: List[str],
    domains: Optional[List[str]] = None,
) -> torch.Tensor:
    """
    Build TF domain one-hot features.

    Args:
        tf_names: List of TF names/descriptions
        domains: List of domain categories (default from V3_CONFIG)

    Returns:
        Tensor of shape [num_tfs, num_domains]
    """
    if domains is None:
        domains = list(V3_CONFIG['tf_domains'])

    n_tfs = len(tf_names)
    n_domains = len(domains)
    features = torch.zeros(n_tfs, n_domains)

    for i, tf_name in enumerate(tf_names):
        tf_lower = tf_name.lower()
        matched = False

        for j, domain in enumerate(domains[:-1]):  # Exclude 'other'
            if domain.lower() in tf_lower:
                features[i, j] = 1.0
                matched = True
                break

        if not matched:
            # Assign to 'other' category
            features[i, -1] = 1.0

    return features


def build_graph_v3(
    base_graph_path: Optional[Path] = None,
    output_path: Optional[Path] = None,
    seed: int = 42,
) -> HeteroData:
    """
    Build the v3 heterogeneous graph with academically defensible design.

    Key Features:
    - Tier-R: Enzyme→Metabolite (catalyzes_R) for supervision
    - rxn_neighbor: Metabolite↔Metabolite based on reaction distance
    - Omics features stored in data (learnable embeddings in model)
    - TF domain one-hot features
    - Currency metabolite validation

    Args:
        base_graph_path: Path to base PPI graph
        output_path: Path to save output graph
        seed: Random seed

    Returns:
        HeteroData graph
    """
    set_seed(seed)

    # Use config defaults
    base_graph_path = base_graph_path or DATA_PATHS['base_graph']
    output_path = output_path or DATA_PATHS['graph_v3']
    kegg_dir = DATA_PATHS['kegg_dir']
    processed_dir = DATA_PATHS['processed_dir']

    max_context_mets = GRAPH_CONFIG['max_context_metabolites']
    rxn_max_dist = V3_CONFIG['rxn_neighbor_max_dist']

    print("=" * 70)
    print("Building v3 Graph (Academically Defensible Design)")
    print("=" * 70)
    print(f"- Tier-R: Enzyme→Metabolite (supervision)")
    print(f"- rxn_neighbor: Metabolite↔Metabolite (reaction distance ≤{rxn_max_dist})")
    print(f"- Omics features: z-score normalized")
    print("=" * 70)

    # 1. Load Base PPI Graph
    if not base_graph_path.exists():
        raise FileNotFoundError(f"Base graph not found: {base_graph_path}")

    base_graph = torch.load(base_graph_path)
    data = base_graph.clone()
    num_enzymes = data['Enzyme'].num_nodes
    print(f"Base graph: {num_enzymes} enzymes")

    # Check for TF nodes
    num_tfs = data['TF'].num_nodes if 'TF' in data.node_types else 0
    print(f"TF nodes: {num_tfs}")

    # 2. Load enzyme mapping
    enzyme_df = pd.read_csv(processed_dir / "enzyme_string_mapping.csv")
    uniprot_to_enzyme_idx = dict(zip(enzyme_df['uniprot_id'], enzyme_df['enzyme_idx']))
    string_to_enzyme_idx = dict(zip(enzyme_df['string_id'], enzyme_df['enzyme_idx']))

    # 3. Load NCBI to STRING mapping
    aliases_path = DATA_PATHS['raw_dir'] / "3847.protein.aliases.v12.0.txt.gz"
    ncbi_to_string = load_string_ncbi_mapping(aliases_path)

    ncbi_to_enzyme_idx = {}
    for ncbi, string_id in ncbi_to_string.items():
        if string_id in string_to_enzyme_idx:
            ncbi_to_enzyme_idx[ncbi] = string_to_enzyme_idx[string_id]

    # 4. Load KEGG gene-EC mapping
    gene_ec_df = pd.read_csv(kegg_dir / "gene_ec_mapping.tsv", sep='\t')

    # 5. Load KEGG -> UniProt mapping
    kegg_uniprot_path = kegg_dir / "kegg_uniprot_mapping.csv"
    kegg_to_uniprot = {}
    if kegg_uniprot_path.exists():
        kegg_uniprot_df = pd.read_csv(kegg_uniprot_path)
        kegg_to_uniprot = dict(zip(
            kegg_uniprot_df['kegg_gene'].astype(str),
            kegg_uniprot_df['uniprot']
        ))

    # 6. Create EC -> Enzyme indices
    ec_to_enzyme_indices = defaultdict(set)
    for _, row in gene_ec_df.iterrows():
        gene_id = str(row['gene_id'])
        ec = row['ec']

        enzyme_idx = None
        uniprot = kegg_to_uniprot.get(gene_id)
        if uniprot and uniprot in uniprot_to_enzyme_idx:
            enzyme_idx = uniprot_to_enzyme_idx[uniprot]
        elif gene_id in ncbi_to_enzyme_idx:
            enzyme_idx = ncbi_to_enzyme_idx[gene_id]

        if enzyme_idx is not None:
            ec_to_enzyme_indices[ec].add(enzyme_idx)

    print(f"ECs with enzyme mapping: {len(ec_to_enzyme_indices)}")

    # 7. Load KEGG edges
    edges_df = pd.read_csv(kegg_dir / "full_enzyme_metabolite_edges.tsv", sep='\t')
    print(f"Full KEGG edges: {len(edges_df)}")

    # 8. Create metabolite index (excluding currency)
    exp_met_path = kegg_dir / "metabolites.csv"
    exp_mets: Set[str] = set()
    if exp_met_path.exists():
        exp_df = pd.read_csv(exp_met_path)
        exp_mets = set(exp_df['compound_id'].unique())

    all_kegg_mets_raw = set(str(m) for m in edges_df['metabolite_id'].unique())
    all_kegg_mets = all_kegg_mets_raw - CURRENCY_METABOLITES
    valid_exp_mets = exp_mets.intersection(all_kegg_mets_raw)
    valid_exp_mets = valid_exp_mets - CURRENCY_METABOLITES

    print(f"Experimental mets (valid): {len(valid_exp_mets)}")
    print(f"Total KEGG metabolites (non-currency): {len(all_kegg_mets)}")

    met_counts = edges_df['metabolite_id'].value_counts()
    top_context_mets = set(str(m) for m in met_counts.index[:max_context_mets]) - CURRENCY_METABOLITES

    final_met_set = valid_exp_mets.union(top_context_mets)
    met_list = sorted(final_met_set)
    met_to_idx = {m: i for i, m in enumerate(met_list)}
    n_metabolites = len(met_list)
    print(f"Final metabolites: {n_metabolites}")

    # 9. Build Tier-R edges (Enzyme→Metabolite)
    tier_r_edges: Set[Tuple[int, int]] = set()

    for _, row in tqdm(edges_df.iterrows(), total=len(edges_df), desc="Building Tier-R edges"):
        met_id = row['metabolite_id']
        ec = row['enzyme_ec']

        if met_id not in met_to_idx:
            continue
        if ec not in ec_to_enzyme_indices:
            continue

        met_idx = met_to_idx[met_id]
        for enz_idx in ec_to_enzyme_indices[ec]:
            tier_r_edges.add((enz_idx, met_idx))

    print(f"Tier-R edges: {len(tier_r_edges)}")

    # 10. Build rxn_neighbor edges (Metabolite↔Metabolite)
    reaction_adj = build_reaction_graph(edges_df, set(met_list), CURRENCY_METABOLITES)
    rxn_neighbor_edge_index = build_rxn_neighbor_edges(
        reaction_adj, met_to_idx, max_dist=rxn_max_dist
    )
    print(f"rxn_neighbor edges: {rxn_neighbor_edge_index.shape[1]}")

    # 11. Validate currency exclusion
    validate_no_currency(reaction_adj, met_to_idx, rxn_neighbor_edge_index)

    # Validate no self-loops
    assert (rxn_neighbor_edge_index[0] != rxn_neighbor_edge_index[1]).all(), \
        "Self-loops detected in rxn_neighbor edges"

    # 12. Build tier_r_lookup for filtered evaluation
    tier_r_lookup: Dict[int, Set[int]] = defaultdict(set)
    for enz_idx, met_idx in tier_r_edges:
        tier_r_lookup[enz_idx].add(met_idx)
    tier_r_lookup = dict(tier_r_lookup)

    # 13. Build omics features
    # Metabolite omics
    met_omics = build_metabolite_omics_features(
        met_list,
        processed_dir / "mtbls531_differential.csv",
        kegg_dir / "metabolites.csv",
    )
    print(f"Metabolite omics shape: {met_omics.shape}")

    # Enzyme omics
    enz_omics = build_enzyme_omics_features(
        num_enzymes,
        enzyme_df,
        processed_dir / "pxd006989_differential.csv",
    )
    print(f"Enzyme omics shape: {enz_omics.shape}")

    # TF domain features
    tf_domain_features = None
    if num_tfs > 0:
        # Get TF names from mapping or node attributes
        tf_mapping_path = processed_dir / "tf_string_mapping.csv"
        if tf_mapping_path.exists():
            tf_df = pd.read_csv(tf_mapping_path)
            tf_names = tf_df['name'].fillna('unknown').tolist() if 'name' in tf_df.columns else ['unknown'] * num_tfs
        else:
            tf_names = ['unknown'] * num_tfs

        tf_domain_features = build_tf_domain_onehot(tf_names)
        print(f"TF domain features shape: {tf_domain_features.shape}")

    # 14. Create edge tensors
    tier_r_src = [e[0] for e in tier_r_edges]
    tier_r_dst = [e[1] for e in tier_r_edges]
    tier_r_edge_index = torch.tensor([tier_r_src, tier_r_dst], dtype=torch.long)

    # 15. Update HeteroData
    # Metabolite nodes
    data['Metabolite'].num_nodes = n_metabolites
    data['Metabolite'].omics_x = met_omics
    data['Metabolite'].node_ids = torch.arange(n_metabolites)
    data['Metabolite'].compound_ids = met_list

    # Mark experimental metabolites
    is_experimental = torch.zeros(n_metabolites, dtype=torch.bool)
    for i, met_id in enumerate(met_list):
        if met_id in valid_exp_mets:
            is_experimental[i] = True
    data['Metabolite'].is_experimental = is_experimental
    print(f"Experimental metabolites in graph: {is_experimental.sum().item()}")

    # EC mappings for negative sampling
    ec_to_met_indices = defaultdict(list)
    for _, row in edges_df.iterrows():
        met_id = row['metabolite_id']
        ec = row['enzyme_ec']
        if met_id in met_to_idx:
            ec_to_met_indices[ec].append(met_to_idx[met_id])
    data['Metabolite'].ec_to_indices = dict(ec_to_met_indices)

    met_to_ecs = defaultdict(set)
    for ec, met_indices in ec_to_met_indices.items():
        for met_idx in met_indices:
            met_to_ecs[met_idx].add(ec)
    data['Metabolite'].met_to_ecs = dict(met_to_ecs)

    # Enzyme nodes
    data['Enzyme'].omics_x = enz_omics
    data['Enzyme'].node_ids = torch.arange(num_enzymes)

    # TF nodes
    if num_tfs > 0 and tf_domain_features is not None:
        data['TF'].domain_x = tf_domain_features
        data['TF'].node_ids = torch.arange(num_tfs)

    # 16. Set edge indices
    # Tier-R: Enzyme→Metabolite (supervision)
    data['Enzyme', 'catalyzes_R', 'Metabolite'].edge_index = tier_r_edge_index

    # Reverse edges
    data['Metabolite', 'rev_catalyzes_R', 'Enzyme'].edge_index = torch.stack([
        tier_r_edge_index[1], tier_r_edge_index[0]
    ])

    # rxn_neighbor: Metabolite↔Metabolite (message passing only)
    data['Metabolite', 'rxn_neighbor', 'Metabolite'].edge_index = rxn_neighbor_edge_index

    # Store tier_r_lookup for filtered evaluation
    data.tier_r_lookup = tier_r_lookup

    # Rename TF edges if present (interacts → associates)
    if ('TF', 'interacts', 'Enzyme') in data.edge_types:
        tf_enz_edges = data['TF', 'interacts', 'Enzyme'].edge_index
        del data['TF', 'interacts', 'Enzyme']
        data['TF', 'associates', 'Enzyme'].edge_index = tf_enz_edges
        data['Enzyme', 'rev_associates', 'TF'].edge_index = torch.stack([
            tf_enz_edges[1], tf_enz_edges[0]
        ])

    # 17. Verify indices
    max_enz_idx = tier_r_edge_index[0].max().item() if tier_r_edge_index.shape[1] > 0 else 0
    max_met_idx = tier_r_edge_index[1].max().item() if tier_r_edge_index.shape[1] > 0 else 0
    assert max_enz_idx < num_enzymes, f"Enzyme index {max_enz_idx} >= {num_enzymes}"
    assert max_met_idx < n_metabolites, f"Metabolite index {max_met_idx} >= {n_metabolites}"

    # 18. Summary
    print("\n" + "=" * 70)
    print("v3 Graph Summary")
    print("=" * 70)
    print(f"Nodes:")
    print(f"  - Enzymes: {num_enzymes}")
    print(f"  - Metabolites: {n_metabolites} (experimental: {is_experimental.sum().item()})")
    print(f"  - TFs: {num_tfs}")
    print(f"Edges:")
    print(f"  - catalyzes_R (Tier-R): {tier_r_edge_index.shape[1]}")
    print(f"  - rxn_neighbor: {rxn_neighbor_edge_index.shape[1]}")
    if 'ppi' in str(data.edge_types):
        ppi_edges = data['Enzyme', 'ppi', 'Enzyme'].edge_index.shape[1]
        print(f"  - PPI: {ppi_edges}")
    print(f"Features:")
    print(f"  - Metabolite omics: {met_omics.shape}")
    print(f"  - Enzyme omics: {enz_omics.shape}")
    if tf_domain_features is not None:
        print(f"  - TF domain: {tf_domain_features.shape}")
    print(f"Evaluation:")
    print(f"  - tier_r_lookup: {len(tier_r_lookup)} enzymes")

    # Save
    torch.save(data, output_path)
    print(f"\nSaved to: {output_path}")

    return data


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Build v3 graph with academically defensible design")
    parser.add_argument(
        '--base-graph', type=str,
        default=str(DATA_PATHS['base_graph']),
        help='Path to base PPI graph'
    )
    parser.add_argument(
        '--output', type=str,
        default=str(DATA_PATHS['graph_v3']),
        help='Output path for v3 graph'
    )
    parser.add_argument(
        '--seed', type=int, default=42,
        help='Random seed'
    )
    args = parser.parse_args()

    build_graph_v3(
        base_graph_path=Path(args.base_graph),
        output_path=Path(args.output),
        seed=args.seed,
    )
