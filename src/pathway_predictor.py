"""
Hybrid Pathway Prediction System

Combines:
1. Known pathway templates (reaction order from literature)
2. GNN predictions (specific enzyme isoform ranking)

Input: Target metabolite (e.g., "Daidzein")
Output: Ordered pathway with predicted enzyme isoforms
"""

import torch
import torch.nn as nn
import numpy as np
import pandas as pd
import gzip
import os
import sys
sys.path.append(os.getcwd())

from torch_geometric.nn import HGTConv, Linear


# ==============================================================================
# ISOFLAVONOID BIOSYNTHESIS PATHWAY TEMPLATE
# ==============================================================================
# Based on established literature (Dixon & Paiva 1995, Liu et al. 2006)

PATHWAY_TEMPLATES = {
    'isoflavonoid': {
        'name': 'Isoflavonoid Biosynthesis',
        'description': 'Phenylalanine → Isoflavonoids (Daidzein, Genistein, etc.)',
        'steps': [
            {
                'step': 1,
                'reaction': 'Phenylalanine → Cinnamic acid',
                'enzyme_class': 'PAL',
                'enzyme_name': 'Phenylalanine ammonia-lyase',
                'ec': '4.3.1.24',
                'substrate': 'L-Phenylalanine',
                'product': 'trans-Cinnamic acid',
                'kegg_substrate': 'C00079',
                'kegg_product': 'C00423',
            },
            {
                'step': 2,
                'reaction': 'Cinnamic acid → p-Coumaric acid',
                'enzyme_class': 'C4H',
                'enzyme_name': 'Cinnamate 4-hydroxylase',
                'ec': '1.14.14.91',
                'substrate': 'trans-Cinnamic acid',
                'product': 'p-Coumaric acid',
                'kegg_substrate': 'C00423',
                'kegg_product': 'C00811',
            },
            {
                'step': 3,
                'reaction': 'p-Coumaric acid → p-Coumaroyl-CoA',
                'enzyme_class': '4CL',
                'enzyme_name': '4-Coumarate-CoA ligase',
                'ec': '6.2.1.12',
                'substrate': 'p-Coumaric acid',
                'product': 'p-Coumaroyl-CoA',
                'kegg_substrate': 'C00811',
                'kegg_product': 'C00223',
            },
            {
                'step': 4,
                'reaction': 'p-Coumaroyl-CoA + 3 Malonyl-CoA → Naringenin chalcone',
                'enzyme_class': 'CHS',
                'enzyme_name': 'Chalcone synthase',
                'ec': '2.3.1.74',
                'substrate': 'p-Coumaroyl-CoA',
                'product': 'Naringenin chalcone',
                'kegg_substrate': 'C00223',
                'kegg_product': 'C06561',
            },
            {
                'step': 5,
                'reaction': 'Naringenin chalcone → Naringenin / Liquiritigenin',
                'enzyme_class': 'CHI',
                'enzyme_name': 'Chalcone isomerase',
                'ec': '5.5.1.6',
                'substrate': 'Naringenin chalcone',
                'product': 'Naringenin',
                'kegg_substrate': 'C06561',
                'kegg_product': 'C00509',
                'branch': {
                    'product_alt': 'Liquiritigenin',
                    'kegg_product_alt': 'C09762',
                }
            },
            {
                'step': 6,
                'reaction': 'Naringenin/Liquiritigenin → Genistein/Daidzein',
                'enzyme_class': 'IFS',
                'enzyme_name': 'Isoflavone synthase',
                'ec': '1.14.14.87',
                'substrate': 'Liquiritigenin',
                'product': 'Daidzein',
                'kegg_substrate': 'C09762',
                'kegg_product': 'C02495',
                'branch': {
                    'substrate_alt': 'Naringenin',
                    'product_alt': 'Genistein',
                    'kegg_substrate_alt': 'C00509',
                    'kegg_product_alt': 'C06563',
                }
            },
            {
                'step': 7,
                'reaction': '2-Hydroxyisoflavanone → Isoflavone (dehydration)',
                'enzyme_class': 'HID',
                'enzyme_name': '2-Hydroxyisoflavanone dehydratase',
                'ec': '4.2.1.105',
                'substrate': '2-Hydroxyisoflavanone',
                'product': 'Isoflavone',
                'kegg_substrate': None,
                'kegg_product': None,
                'note': 'Often coupled with IFS'
            },
        ],
        'final_products': ['Daidzein', 'Genistein', 'Formononetin', 'Glycitein'],
        'kegg_final': ['C02495', 'C06563', 'C00858'],
    },
    'phaseollin': {
        'name': 'Phaseollin (Phytoalexin) Biosynthesis',
        'description': 'Daidzein → Phaseollin',
        'prerequisite': 'isoflavonoid',
        'steps': [
            {
                'step': 1,
                'reaction': 'Daidzein → 2\'-Hydroxydaidzein',
                'enzyme_class': 'I2\'H',
                'enzyme_name': 'Isoflavone 2\'-hydroxylase',
                'ec': '1.14.14.-',
            },
            {
                'step': 2,
                'reaction': '2\'-Hydroxydaidzein → Phaseollin',
                'enzyme_class': 'PTS',
                'enzyme_name': 'Pterocarpan synthase',
                'ec': '1.1.1.-',
            },
        ],
        'final_products': ['Phaseollin'],
        'kegg_final': ['C10514'],
    }
}


# ==============================================================================
# MODEL DEFINITION (Same as ablation scripts)
# ==============================================================================

class HGTWithFeatures(nn.Module):
    def __init__(self, metadata, in_channels, hidden_channels, out_channels,
                 num_heads=4, num_layers=3, feature_dims=None):
        super().__init__()
        self.feature_dims = feature_dims or {}

        self.lin_dict = nn.ModuleDict()
        for node_type in metadata[0]:
            feat_dim = self.feature_dims.get(node_type, in_channels)
            self.lin_dict[node_type] = Linear(feat_dim, hidden_channels)

        self.convs = nn.ModuleList()
        for _ in range(num_layers):
            conv = HGTConv(hidden_channels, hidden_channels, metadata, heads=num_heads)
            self.convs.append(conv)

    def forward(self, x_dict, edge_index_dict):
        x_dict = {
            node_type: self.lin_dict[node_type](x).relu_()
            for node_type, x in x_dict.items()
        }
        for conv in self.convs:
            x_dict = conv(x_dict, edge_index_dict)
        return x_dict


class LinkPredictor(nn.Module):
    def __init__(self, in_channels):
        super().__init__()

    def forward(self, x_src, x_dst, edge_label_index):
        row, col = edge_label_index
        return (x_src[row] * x_dst[col]).sum(dim=-1)


# ==============================================================================
# ENZYME DATABASE
# ==============================================================================

def load_enzyme_database(raw_dir='data/raw'):
    """Load enzyme annotations from STRING protein info."""
    protein_info_path = os.path.join(raw_dir, '3847.protein.info.v12.0.txt.gz')

    enzymes = {}
    with gzip.open(protein_info_path, 'rt') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                string_id = parts[0]
                name = parts[1]
                annotation = parts[3]

                uniprot = string_id.split('.')[-1] if '.' in string_id else string_id
                enzymes[uniprot] = {
                    'string_id': string_id,
                    'name': name,
                    'annotation': annotation
                }

    return enzymes


def find_enzymes_by_class(enzyme_db, enzyme_class, annotation_keywords):
    """Find enzymes matching a class by name/annotation keywords."""
    matches = []

    class_lower = enzyme_class.lower()
    keywords_lower = [k.lower() for k in annotation_keywords]

    for uniprot, info in enzyme_db.items():
        name = info['name'].lower()
        ann = info['annotation'].lower()

        # Check name match
        if class_lower in name:
            matches.append((uniprot, info, 'name_match'))
            continue

        # Check annotation match
        for kw in keywords_lower:
            if kw in ann:
                matches.append((uniprot, info, 'annotation_match'))
                break

    return matches


# ==============================================================================
# GNN ENZYME RANKING
# ==============================================================================

def load_proteomics_features(proteomics_path):
    """Load proteomics data for combined features."""
    prot_df = pd.read_csv(proteomics_path)

    proteomics_features = {}
    for _, row in prot_df.iterrows():
        string_id = row['STRING_ID']
        if pd.isna(string_id):
            continue

        uniprot_id = string_id.split('.')[-1] if '.' in string_id else string_id

        features = [
            row['Log2FC'] if not pd.isna(row['Log2FC']) else 0.0,
            row['P_Value'] if not pd.isna(row['P_Value']) else 1.0,
            row['Mean_Control'] if not pd.isna(row['Mean_Control']) else 0.0,
            row['Mean_Ethylene'] if not pd.isna(row['Mean_Ethylene']) else 0.0,
            abs(row['Log2FC']) if not pd.isna(row['Log2FC']) else 0.0,
            -np.log10(row['P_Value'] + 1e-10) if not pd.isna(row['P_Value']) else 0.0,
            1 if row['P_Value'] < 0.05 else 0,
        ]

        proteomics_features[uniprot_id] = features

    return proteomics_features


def add_proteomics_to_graph(data, proteomics_features, enzyme_mapping_df):
    """Add proteomics features to enzyme nodes."""
    device = data['Enzyme'].x.device
    num_enzymes = data['Enzyme'].num_nodes
    n_prot_features = 7

    enzyme_idx_to_uniprot = {}
    for _, row in enzyme_mapping_df.iterrows():
        enzyme_idx_to_uniprot[int(row['enzyme_idx'])] = row['uniprot_id']

    random_x = data['Enzyme'].x
    prot_x = torch.zeros(num_enzymes, n_prot_features, device=device)

    for idx in range(num_enzymes):
        if idx in enzyme_idx_to_uniprot:
            uniprot_id = enzyme_idx_to_uniprot[idx]
            if uniprot_id in proteomics_features:
                feats = proteomics_features[uniprot_id]
                prot_x[idx, :len(feats)] = torch.tensor(feats, dtype=torch.float32)

    # Normalize
    prot_mean = prot_x.mean(dim=0)
    prot_std = prot_x.std(dim=0) + 1e-8
    prot_x = (prot_x - prot_mean) / prot_std

    data['Enzyme'].x = torch.cat([random_x, prot_x], dim=1)

    return data


def train_gnn_model(data, device, epochs=30, num_layers=3):
    """Train GNN model and return trained model."""
    feature_dims = {}
    for node_type in data.node_types:
        feature_dims[node_type] = data[node_type].x.shape[1]

    model = HGTWithFeatures(
        data.metadata(),
        in_channels=64,
        hidden_channels=64,
        out_channels=64,
        num_heads=4,
        num_layers=num_layers,
        feature_dims=feature_dims
    ).to(device)

    predictor = LinkPredictor(64).to(device)

    optimizer = torch.optim.Adam(
        list(model.parameters()) + list(predictor.parameters()),
        lr=0.01, weight_decay=1e-5
    )

    edge_index = data['Enzyme', 'catalyzes', 'Metabolite'].edge_index

    model.train()
    for epoch in range(1, epochs + 1):
        optimizer.zero_grad()

        x_dict = model(data.x_dict, data.edge_index_dict)

        num_pos = edge_index.size(1)
        num_metabolites = data['Metabolite'].num_nodes

        neg_src = edge_index[0]
        offset = torch.randint(1, 6, (num_pos,), device=device) * \
                 (2 * torch.randint(0, 2, (num_pos,), device=device) - 1)
        neg_dst = (edge_index[1] + offset) % num_metabolites

        pos_out = predictor(x_dict['Enzyme'], x_dict['Metabolite'], edge_index)
        neg_out = predictor(x_dict['Enzyme'], x_dict['Metabolite'],
                           torch.stack([neg_src, neg_dst]))

        loss = -torch.log(torch.sigmoid(pos_out) + 1e-15).mean() - \
               torch.log(1 - torch.sigmoid(neg_out) + 1e-15).mean()

        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()

    return model, predictor


def rank_enzymes_for_metabolite(model, predictor, data, metabolite_idx, candidate_indices, device):
    """Rank enzyme candidates for a specific metabolite."""
    model.eval()

    with torch.no_grad():
        x_dict = model(data.x_dict, data.edge_index_dict)
        enz_emb = x_dict['Enzyme']
        met_emb = x_dict['Metabolite']

        # Create edges for scoring
        candidate_tensor = torch.tensor(candidate_indices, device=device)
        met_tensor = torch.tensor([metabolite_idx] * len(candidate_indices), device=device)
        eval_edges = torch.stack([candidate_tensor, met_tensor])

        scores = predictor(enz_emb, met_emb, eval_edges).sigmoid()

    return scores.cpu().numpy()


# ==============================================================================
# PATHWAY PREDICTION
# ==============================================================================

ENZYME_KEYWORDS = {
    'PAL': ['phenylalanine ammonia-lyase', 'PAL'],
    'C4H': ['cinnamate 4-hydroxylase', 'cinnamate-4-hydroxylase', 'C4H'],
    '4CL': ['4-coumarate:CoA ligase', '4-coumarate-CoA ligase', '4CL'],
    'CHS': ['chalcone synthase', 'CHS'],
    'CHI': ['chalcone isomerase', 'chalcone-flavonone isomerase', 'CHI'],
    'IFS': ['isoflavone synthase', 'IFS', '2-hydroxyisoflavanone synthase'],
    'HID': ['2-hydroxyisoflavanone dehydratase', 'HID', 'HIDH'],
    'IFR': ['isoflavone reductase', 'IFR'],
    'I2\'H': ['isoflavone 2\'-hydroxylase', 'isoflavone hydroxylase'],
    'IOMT': ['isoflavone O-methyltransferase', 'IOMT'],
    'IF7GT': ['isoflavone 7-O-glucosyltransferase', 'IF7GT'],
}


def predict_pathway(target_metabolite, pathway_type='isoflavonoid', top_k=3, device=None):
    """
    Predict complete pathway to target metabolite.

    Args:
        target_metabolite: Target metabolite name (e.g., 'Daidzein')
        pathway_type: Pathway template to use
        top_k: Number of top enzyme isoforms to show per step
        device: Torch device

    Returns:
        Ordered pathway with predicted enzymes
    """
    if device is None:
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    print(f"\n{'='*70}")
    print(f"PATHWAY PREDICTION: {target_metabolite}")
    print(f"{'='*70}")

    # Load data
    print("\n[1/4] Loading data...")
    data = torch.load('data/processed/strict_bipartite_v2.pt')
    data = data.to(device)

    enzyme_mapping_df = pd.read_csv('data/processed/enzyme_string_mapping.csv')
    uniprot_to_idx = dict(zip(enzyme_mapping_df['uniprot_id'], enzyme_mapping_df['enzyme_idx']))
    idx_to_uniprot = dict(zip(enzyme_mapping_df['enzyme_idx'], enzyme_mapping_df['uniprot_id']))

    enzyme_db = load_enzyme_database()

    # Add proteomics features
    print("[2/4] Adding proteomics features...")
    proteomics_features = load_proteomics_features('data/processed/pxd006989_mapped.csv')
    data = add_proteomics_to_graph(data, proteomics_features, enzyme_mapping_df)

    # Train model
    print("[3/4] Training GNN model...")
    model, predictor = train_gnn_model(data, device)

    # Get pathway template
    template = PATHWAY_TEMPLATES.get(pathway_type)
    if not template:
        print(f"Unknown pathway type: {pathway_type}")
        return None

    print(f"\n[4/4] Predicting pathway: {template['name']}")
    print(f"Description: {template['description']}")

    # Predict each step
    pathway_result = []

    print(f"\n{'='*70}")
    print("PREDICTED PATHWAY")
    print(f"{'='*70}")

    for step_info in template['steps']:
        step_num = step_info['step']
        enzyme_class = step_info['enzyme_class']
        enzyme_name = step_info['enzyme_name']
        reaction = step_info['reaction']

        print(f"\n--- Step {step_num}: {reaction} ---")
        print(f"Enzyme class: {enzyme_class} ({enzyme_name})")

        # Find candidate enzymes
        keywords = ENZYME_KEYWORDS.get(enzyme_class, [enzyme_name])
        candidates = find_enzymes_by_class(enzyme_db, enzyme_class, keywords)

        if not candidates:
            print(f"  No candidates found for {enzyme_class}")
            pathway_result.append({
                'step': step_num,
                'reaction': reaction,
                'enzyme_class': enzyme_class,
                'candidates': [],
                'status': 'no_candidates'
            })
            continue

        # Get indices for candidates that are in our graph
        candidate_data = []
        for uniprot, info, match_type in candidates:
            if uniprot in uniprot_to_idx:
                idx = uniprot_to_idx[uniprot]

                # Get proteomics info
                prot_info = proteomics_features.get(uniprot, None)
                log2fc = prot_info[0] if prot_info else None
                pvalue = prot_info[1] if prot_info else None

                candidate_data.append({
                    'uniprot': uniprot,
                    'idx': idx,
                    'name': info['name'],
                    'annotation': info['annotation'][:80] + '...' if len(info['annotation']) > 80 else info['annotation'],
                    'match_type': match_type,
                    'log2fc': log2fc,
                    'pvalue': pvalue
                })

        if not candidate_data:
            print(f"  Candidates found but none in graph")
            pathway_result.append({
                'step': step_num,
                'reaction': reaction,
                'enzyme_class': enzyme_class,
                'candidates': candidates,
                'status': 'not_in_graph'
            })
            continue

        # Rank by GNN (use first metabolite as proxy, or use proteomics)
        candidate_indices = [c['idx'] for c in candidate_data]

        # Use metabolite 0 as proxy for ranking (or could average across related metabolites)
        gnn_scores = rank_enzymes_for_metabolite(model, predictor, data, 0, candidate_indices, device)

        for i, cand in enumerate(candidate_data):
            cand['gnn_score'] = gnn_scores[i]

            # Combined score: GNN + proteomics significance
            prot_boost = 0
            if cand['log2fc'] is not None and cand['pvalue'] is not None:
                if cand['pvalue'] < 0.05:
                    prot_boost = 0.2 * abs(cand['log2fc'])
            cand['combined_score'] = cand['gnn_score'] + prot_boost

        # Sort by combined score
        candidate_data.sort(key=lambda x: x['combined_score'], reverse=True)

        # Display top candidates
        print(f"  Top-{top_k} predicted isoforms:")
        for i, cand in enumerate(candidate_data[:top_k]):
            prot_str = ""
            if cand['log2fc'] is not None:
                sig = "*" if cand['pvalue'] < 0.05 else ""
                prot_str = f" [Log2FC: {cand['log2fc']:.2f}{sig}]"

            print(f"    {i+1}. {cand['name']} ({cand['uniprot']})")
            print(f"       Score: {cand['combined_score']:.4f} (GNN: {cand['gnn_score']:.4f}){prot_str}")

        pathway_result.append({
            'step': step_num,
            'reaction': reaction,
            'enzyme_class': enzyme_class,
            'top_candidates': candidate_data[:top_k],
            'all_candidates': len(candidate_data),
            'status': 'predicted'
        })

    # Summary
    print(f"\n{'='*70}")
    print("PATHWAY SUMMARY")
    print(f"{'='*70}")

    print(f"\nTarget: {target_metabolite}")
    print(f"Pathway: {template['name']}\n")

    print("Precursor → ", end="")
    for i, step in enumerate(pathway_result):
        if step['status'] == 'predicted' and step['top_candidates']:
            best = step['top_candidates'][0]
            print(f"[{best['name']}]", end="")
        else:
            print(f"[{step['enzyme_class']}?]", end="")

        if i < len(pathway_result) - 1:
            print(" → ", end="")
    print(f" → {target_metabolite}")

    return pathway_result


def find_associated_proteins(data, enzyme_indices, enzyme_mapping_df, tf_mapping_df,
                              enzyme_db, proteomics_features, top_k=5):
    """
    Find TFs and proteins that interact with pathway enzymes.

    Returns:
        dict with 'tfs' and 'proteins' lists
    """
    # Build reverse mappings
    idx_to_enzyme_uniprot = dict(zip(enzyme_mapping_df['enzyme_idx'], enzyme_mapping_df['uniprot_id']))
    idx_to_tf_uniprot = dict(zip(tf_mapping_df['tf_idx'], tf_mapping_df['uniprot_id']))
    tf_idx_to_name = dict(zip(tf_mapping_df['tf_idx'], tf_mapping_df['name']))

    # Get TF -> Enzyme edges
    tf_enz_edges = data['TF', 'interacts', 'Enzyme'].edge_index

    # Get Protein -> Enzyme edges
    prot_enz_edges = data['Protein', 'interacts', 'Enzyme'].edge_index

    # Find TFs interacting with pathway enzymes
    interacting_tfs = {}
    for enz_idx in enzyme_indices:
        # Find TFs connected to this enzyme
        mask = tf_enz_edges[1] == enz_idx
        tf_indices = tf_enz_edges[0][mask].tolist()

        for tf_idx in tf_indices:
            if tf_idx not in interacting_tfs:
                interacting_tfs[tf_idx] = {
                    'idx': tf_idx,
                    'uniprot': idx_to_tf_uniprot.get(tf_idx, f'TF_{tf_idx}'),
                    'name': tf_idx_to_name.get(tf_idx, f'TF_{tf_idx}'),
                    'connected_enzymes': [],
                    'connection_count': 0
                }
            interacting_tfs[tf_idx]['connected_enzymes'].append(enz_idx)
            interacting_tfs[tf_idx]['connection_count'] += 1

    # Add proteomics info to TFs
    for tf_idx, tf_info in interacting_tfs.items():
        uniprot = tf_info['uniprot']
        if uniprot in proteomics_features:
            prot = proteomics_features[uniprot]
            tf_info['log2fc'] = prot[0]
            tf_info['pvalue'] = prot[1]
            tf_info['significant'] = prot[1] < 0.05 if prot[1] else False
        else:
            tf_info['log2fc'] = None
            tf_info['pvalue'] = None
            tf_info['significant'] = False

    # Sort TFs by connection count and significance
    tf_list = list(interacting_tfs.values())
    tf_list.sort(key=lambda x: (x['significant'], x['connection_count'],
                                 abs(x['log2fc']) if x['log2fc'] else 0), reverse=True)

    # Find signaling/other proteins (sample from Protein nodes)
    # Note: We don't have direct Protein mapping, so we report TFs mainly

    return {
        'tfs': tf_list[:top_k*2],
        'tf_count': len(tf_list),
    }


def predict_pathway_with_context(target_metabolite, pathway_type='isoflavonoid',
                                  top_k=3, device=None):
    """
    Predict pathway with associated regulatory proteins.
    """
    if device is None:
        device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    print(f"\n{'='*70}")
    print(f"PATHWAY PREDICTION WITH REGULATORY CONTEXT: {target_metabolite}")
    print(f"{'='*70}")

    # Load data
    print("\n[1/5] Loading data...")
    data = torch.load('data/processed/strict_bipartite_v2.pt')
    data = data.to(device)

    enzyme_mapping_df = pd.read_csv('data/processed/enzyme_string_mapping.csv')
    tf_mapping_df = pd.read_csv('data/processed/tf_string_mapping.csv')

    uniprot_to_idx = dict(zip(enzyme_mapping_df['uniprot_id'], enzyme_mapping_df['enzyme_idx']))
    idx_to_uniprot = dict(zip(enzyme_mapping_df['enzyme_idx'], enzyme_mapping_df['uniprot_id']))

    enzyme_db = load_enzyme_database()

    # Add proteomics features
    print("[2/5] Adding proteomics features...")
    proteomics_features = load_proteomics_features('data/processed/pxd006989_mapped.csv')
    data = add_proteomics_to_graph(data, proteomics_features, enzyme_mapping_df)

    # Train model
    print("[3/5] Training GNN model...")
    model, predictor = train_gnn_model(data, device)

    # Get pathway template
    template = PATHWAY_TEMPLATES.get(pathway_type)
    if not template:
        print(f"Unknown pathway type: {pathway_type}")
        return None

    print(f"\n[4/5] Predicting pathway: {template['name']}")

    # Collect pathway enzymes
    pathway_enzyme_indices = []
    pathway_result = []

    print(f"\n{'='*70}")
    print("PREDICTED BIOSYNTHETIC ENZYMES")
    print(f"{'='*70}")

    for step_info in template['steps']:
        step_num = step_info['step']
        enzyme_class = step_info['enzyme_class']
        enzyme_name = step_info['enzyme_name']
        reaction = step_info['reaction']

        print(f"\n--- Step {step_num}: {reaction} ---")
        print(f"Enzyme class: {enzyme_class} ({enzyme_name})")

        # Find candidate enzymes
        keywords = ENZYME_KEYWORDS.get(enzyme_class, [enzyme_name])
        candidates = find_enzymes_by_class(enzyme_db, enzyme_class, keywords)

        if not candidates:
            print(f"  No candidates found for {enzyme_class}")
            continue

        # Get indices for candidates in graph
        candidate_data = []
        for uniprot, info, match_type in candidates:
            if uniprot in uniprot_to_idx:
                idx = uniprot_to_idx[uniprot]
                prot_info = proteomics_features.get(uniprot, None)
                log2fc = prot_info[0] if prot_info else None
                pvalue = prot_info[1] if prot_info else None

                candidate_data.append({
                    'uniprot': uniprot,
                    'idx': idx,
                    'name': info['name'],
                    'log2fc': log2fc,
                    'pvalue': pvalue
                })

        if not candidate_data:
            print(f"  Candidates found but none in graph")
            continue

        # Rank by GNN
        candidate_indices = [c['idx'] for c in candidate_data]
        gnn_scores = rank_enzymes_for_metabolite(model, predictor, data, 0, candidate_indices, device)

        for i, cand in enumerate(candidate_data):
            cand['gnn_score'] = gnn_scores[i]
            prot_boost = 0
            if cand['log2fc'] is not None and cand['pvalue'] is not None:
                if cand['pvalue'] < 0.05:
                    prot_boost = 0.2 * abs(cand['log2fc'])
            cand['combined_score'] = cand['gnn_score'] + prot_boost

        candidate_data.sort(key=lambda x: x['combined_score'], reverse=True)

        # Display and collect top candidates
        print(f"  Top-{top_k} predicted isoforms:")
        for i, cand in enumerate(candidate_data[:top_k]):
            prot_str = ""
            if cand['log2fc'] is not None:
                sig = "*" if cand['pvalue'] < 0.05 else ""
                prot_str = f" [Log2FC: {cand['log2fc']:.2f}{sig}]"
            print(f"    {i+1}. {cand['name']} ({cand['uniprot']})")
            print(f"       Score: {cand['combined_score']:.4f}{prot_str}")

            pathway_enzyme_indices.append(cand['idx'])

        pathway_result.append({
            'step': step_num,
            'enzyme_class': enzyme_class,
            'top_candidates': candidate_data[:top_k]
        })

    # Find associated regulatory proteins
    print(f"\n[5/5] Finding regulatory proteins...")

    associated = find_associated_proteins(
        data, pathway_enzyme_indices, enzyme_mapping_df, tf_mapping_df,
        enzyme_db, proteomics_features, top_k=top_k
    )

    # Display associated TFs
    print(f"\n{'='*70}")
    print("ASSOCIATED TRANSCRIPTION FACTORS")
    print(f"{'='*70}")
    print(f"Total TFs interacting with pathway: {associated['tf_count']}")
    print(f"\nTop regulatory TFs:")

    for i, tf in enumerate(associated['tfs'][:10]):
        prot_str = ""
        if tf['log2fc'] is not None:
            sig = "**" if tf['significant'] else ""
            prot_str = f" [Log2FC: {tf['log2fc']:.2f}{sig}]"

        print(f"  {i+1}. {tf['name']} ({tf['uniprot']})")
        print(f"     Connected to {tf['connection_count']} pathway enzyme(s){prot_str}")

    # Summary visualization
    print(f"\n{'='*70}")
    print("PATHWAY SUMMARY")
    print(f"{'='*70}")

    print(f"\nTarget: {target_metabolite}")
    print(f"Pathway: {template['name']}\n")

    # Pathway flow
    print("BIOSYNTHETIC FLOW:")
    print("Precursor → ", end="")
    for i, step in enumerate(pathway_result):
        if step['top_candidates']:
            best = step['top_candidates'][0]
            print(f"[{best['name']}]", end="")
        else:
            print(f"[{step['enzyme_class']}?]", end="")
        if i < len(pathway_result) - 1:
            print(" → ", end="")
    print(f" → {target_metabolite}")

    # Regulatory context
    print("\nREGULATORY CONTEXT:")
    sig_tfs = [tf for tf in associated['tfs'] if tf['significant']]
    if sig_tfs:
        print(f"  Ethylene-responsive TFs: {len(sig_tfs)}")
        for tf in sig_tfs[:5]:
            direction = "↑" if tf['log2fc'] > 0 else "↓"
            print(f"    - {tf['name']}: {direction} ({tf['log2fc']:.1f}x, p<0.05)")
    else:
        print("  No significantly changed TFs with proteomics data")

    return {
        'pathway': pathway_result,
        'associated_tfs': associated['tfs'],
        'pathway_enzymes': pathway_enzyme_indices
    }


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(description='Predict biosynthetic pathway')
    parser.add_argument('--target', type=str, default='Daidzein',
                        help='Target metabolite name')
    parser.add_argument('--pathway', type=str, default='isoflavonoid',
                        choices=['isoflavonoid', 'phaseollin'],
                        help='Pathway type')
    parser.add_argument('--top_k', type=int, default=3,
                        help='Number of top candidates per step')
    parser.add_argument('--with-context', action='store_true',
                        help='Include regulatory proteins')

    args = parser.parse_args()

    if args.with_context:
        predict_pathway_with_context(
            target_metabolite=args.target,
            pathway_type=args.pathway,
            top_k=args.top_k
        )
    else:
        predict_pathway(
            target_metabolite=args.target,
            pathway_type=args.pathway,
            top_k=args.top_k
        )


if __name__ == "__main__":
    main()
