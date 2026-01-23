"""
Pathway Activation Scorer Module

Computes pathway activation scores based on integrated metabolomics and proteomics data.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple
from pathlib import Path

from .pathway_mapper import PathwayMapper


class PathwayActivationScorer:
    """
    Computes pathway activation scores from metabolomics and proteomics data.

    The score reflects how well the observed metabolite changes are explained
    by enzyme activity changes in each pathway.
    """

    def __init__(
        self,
        metabolomics_df: pd.DataFrame,
        proteomics_df: pd.DataFrame,
        mapper: PathwayMapper,
        significance_threshold: float = 0.05
    ):
        """
        Initialize the scorer with omics data and pathway mapper.

        Args:
            metabolomics_df: DataFrame with columns [KEGG, Log2FC, P_Value]
            proteomics_df: DataFrame with columns [STRING_ID, Log2FC, P_Value]
            mapper: PathwayMapper instance
            significance_threshold: P-value threshold for significance weighting
        """
        self.metabolomics_df = metabolomics_df.copy()
        self.proteomics_df = proteomics_df.copy()
        self.mapper = mapper
        self.significance_threshold = significance_threshold

        # Preprocess data
        self._preprocess_metabolomics()
        self._preprocess_proteomics()

    def _preprocess_metabolomics(self):
        """Preprocess metabolomics data for analysis."""
        # Create KEGG -> data mapping
        self.metabolite_data = {}

        for _, row in self.metabolomics_df.iterrows():
            kegg_id = row.get('KEGG')
            if pd.isna(kegg_id) or str(kegg_id).strip() == '':
                continue

            kegg_id = str(kegg_id).strip()
            log2fc = float(row['Log2FC']) if pd.notna(row['Log2FC']) else 0.0
            p_value = float(row['P_Value']) if pd.notna(row['P_Value']) else 1.0

            # Store or update (keep the one with lower p-value if duplicate)
            if kegg_id not in self.metabolite_data:
                self.metabolite_data[kegg_id] = {
                    'log2fc': log2fc,
                    'p_value': p_value,
                    'sig_weight': self._compute_significance_weight(p_value)
                }
            else:
                if p_value < self.metabolite_data[kegg_id]['p_value']:
                    self.metabolite_data[kegg_id] = {
                        'log2fc': log2fc,
                        'p_value': p_value,
                        'sig_weight': self._compute_significance_weight(p_value)
                    }

    def _preprocess_proteomics(self):
        """Preprocess proteomics data for analysis."""
        # Create UniProt -> data mapping
        self.protein_data = {}

        for _, row in self.proteomics_df.iterrows():
            string_id = row.get('STRING_ID')
            if pd.isna(string_id):
                continue

            string_id = str(string_id).strip()
            uniprot_id = self.mapper.string_to_uniprot_id(string_id)

            if uniprot_id is None:
                # Try extracting from string_id (format: 3847.XXXX)
                if '.' in string_id:
                    uniprot_id = string_id.split('.')[-1]
                else:
                    continue

            log2fc = float(row['Log2FC']) if pd.notna(row['Log2FC']) else 0.0
            p_value = float(row['P_Value']) if pd.notna(row['P_Value']) else 1.0

            self.protein_data[uniprot_id] = {
                'log2fc': log2fc,
                'p_value': p_value,
                'sig_weight': self._compute_significance_weight(p_value),
                'string_id': string_id
            }

    def _compute_significance_weight(self, p_value: float) -> float:
        """
        Compute significance weight from p-value.

        Args:
            p_value: P-value

        Returns:
            Weight based on -log10(p_value), normalized
        """
        if p_value <= 0:
            p_value = 1e-300
        return -np.log10(p_value)

    def compute_metabolite_signal(self, pathway_id: str) -> Dict:
        """
        Compute metabolite signal for a pathway.

        Args:
            pathway_id: KEGG pathway ID

        Returns:
            Dictionary with mean_fc, weighted_fc, n_detected, compounds
        """
        pathway_compounds = self.mapper.get_pathway_metabolites(pathway_id)

        if not pathway_compounds:
            return {
                'mean_fc': 0.0,
                'weighted_fc': 0.0,
                'n_detected': 0,
                'n_total': 0,
                'compounds': [],
                'coverage': 0.0
            }

        detected = []
        fc_values = []
        weights = []

        for compound_id in pathway_compounds:
            if compound_id in self.metabolite_data:
                data = self.metabolite_data[compound_id]
                detected.append({
                    'compound_id': compound_id,
                    'log2fc': data['log2fc'],
                    'p_value': data['p_value'],
                    'sig_weight': data['sig_weight']
                })
                fc_values.append(data['log2fc'])
                weights.append(data['sig_weight'])

        n_detected = len(detected)
        n_total = len(pathway_compounds)
        coverage = n_detected / n_total if n_total > 0 else 0.0

        if n_detected == 0:
            return {
                'mean_fc': 0.0,
                'weighted_fc': 0.0,
                'n_detected': 0,
                'n_total': n_total,
                'compounds': [],
                'coverage': 0.0
            }

        mean_fc = np.mean(fc_values)
        total_weight = sum(weights)
        weighted_fc = sum(fc * w for fc, w in zip(fc_values, weights)) / total_weight if total_weight > 0 else mean_fc

        return {
            'mean_fc': mean_fc,
            'weighted_fc': weighted_fc,
            'n_detected': n_detected,
            'n_total': n_total,
            'compounds': detected,
            'coverage': coverage
        }

    def compute_enzyme_signal(self, pathway_id: str) -> Dict:
        """
        Compute enzyme signal for a pathway.

        Args:
            pathway_id: KEGG pathway ID

        Returns:
            Dictionary with mean_fc, weighted_fc, n_detected, proteins
        """
        pathway_proteins = self.mapper.get_pathway_proteins(pathway_id)

        if not pathway_proteins:
            return {
                'mean_fc': 0.0,
                'weighted_fc': 0.0,
                'n_detected': 0,
                'n_total': 0,
                'proteins': [],
                'coverage': 0.0
            }

        detected = []
        fc_values = []
        weights = []

        for uniprot_id in pathway_proteins:
            if uniprot_id in self.protein_data:
                data = self.protein_data[uniprot_id]
                detected.append({
                    'uniprot_id': uniprot_id,
                    'log2fc': data['log2fc'],
                    'p_value': data['p_value'],
                    'sig_weight': data['sig_weight'],
                    'string_id': data.get('string_id', '')
                })
                fc_values.append(data['log2fc'])
                weights.append(data['sig_weight'])

        n_detected = len(detected)
        n_total = len(pathway_proteins)
        coverage = n_detected / n_total if n_total > 0 else 0.0

        if n_detected == 0:
            return {
                'mean_fc': 0.0,
                'weighted_fc': 0.0,
                'n_detected': 0,
                'n_total': n_total,
                'proteins': [],
                'coverage': 0.0
            }

        mean_fc = np.mean(fc_values)
        total_weight = sum(weights)
        weighted_fc = sum(fc * w for fc, w in zip(fc_values, weights)) / total_weight if total_weight > 0 else mean_fc

        return {
            'mean_fc': mean_fc,
            'weighted_fc': weighted_fc,
            'n_detected': n_detected,
            'n_total': n_total,
            'proteins': detected,
            'coverage': coverage
        }

    def compute_activation_score(self, pathway_id: str) -> Dict:
        """
        Compute the final pathway activation score.

        Score = combined evidence from metabolite and enzyme changes,
        weighted by significance and coverage.

        Args:
            pathway_id: KEGG pathway ID

        Returns:
            Dictionary with activation_score and component scores
        """
        metabolite_signal = self.compute_metabolite_signal(pathway_id)
        enzyme_signal = self.compute_enzyme_signal(pathway_id)

        # Base scores from fold changes (using absolute values for upregulation)
        metabolite_fc = metabolite_signal['weighted_fc']
        enzyme_fc = enzyme_signal['weighted_fc']

        # Coverage factors
        metabolite_coverage = metabolite_signal['coverage']
        enzyme_coverage = enzyme_signal['coverage']

        # Evidence consistency score
        # Higher if both metabolites and enzymes show coordinated upregulation
        if metabolite_fc > 0 and enzyme_fc > 0:
            # Both upregulated - strong positive evidence
            consistency_bonus = 1.0 + min(metabolite_fc, enzyme_fc) * 0.1
        elif metabolite_fc > 0 and enzyme_fc < 0:
            # Metabolites up but enzymes down - weaker evidence
            consistency_bonus = 0.7
        elif metabolite_fc < 0:
            # Metabolites downregulated
            consistency_bonus = 0.5
        else:
            consistency_bonus = 0.8

        # Compute individual component scores
        metabolite_score = abs(metabolite_fc) * metabolite_signal.get('n_detected', 0) * consistency_bonus
        enzyme_score = abs(enzyme_fc) * enzyme_signal.get('n_detected', 0)

        # Metabolite detection bonus: reward pathways with detected metabolites
        n_metabolites = metabolite_signal['n_detected']
        metabolite_bonus = 1.0 + n_metabolites * 0.3

        # Combined score
        # Prioritize metabolite evidence - model goal is to explain metabolite changes
        if metabolite_signal['n_detected'] > 0 and enzyme_signal['n_detected'] > 0:
            # Both metabolites and enzymes detected: highest score, no coverage penalty
            activation_score = (
                0.6 * metabolite_score +
                0.4 * enzyme_score
            ) * metabolite_bonus
        elif metabolite_signal['n_detected'] > 0:
            # Only metabolites detected
            activation_score = metabolite_score * 0.8
        elif enzyme_signal['n_detected'] > 0:
            # Only enzymes detected: heavy penalty (0.5 -> 0.1)
            # Pathways without metabolite evidence should not rank high
            activation_score = enzyme_score * 0.1
        else:
            activation_score = 0.0

        return {
            'pathway_id': pathway_id,
            'pathway_name': self.mapper.get_pathway_name(pathway_id),
            'activation_score': activation_score,
            'normalized_score': 0.0,  # Will be normalized across all pathways
            'metabolite_fc': metabolite_fc,
            'enzyme_fc': enzyme_fc,
            'metabolite_coverage': metabolite_coverage,
            'enzyme_coverage': enzyme_coverage,
            'n_metabolites_detected': metabolite_signal['n_detected'],
            'n_metabolites_total': metabolite_signal['n_total'],
            'n_enzymes_detected': enzyme_signal['n_detected'],
            'n_enzymes_total': enzyme_signal['n_total'],
            'consistency_bonus': consistency_bonus,
            'metabolite_details': metabolite_signal['compounds'],
            'enzyme_details': enzyme_signal['proteins']
        }

    def rank_all_pathways(self) -> pd.DataFrame:
        """
        Compute and rank activation scores for all pathways.

        Returns:
            DataFrame with pathway rankings sorted by activation score
        """
        all_pathways = self.mapper.get_all_pathways()
        results = []

        for pathway_id in all_pathways:
            score_info = self.compute_activation_score(pathway_id)
            results.append(score_info)

        # Normalize scores
        if results:
            max_score = max(r['activation_score'] for r in results)
            if max_score > 0:
                for r in results:
                    r['normalized_score'] = r['activation_score'] / max_score

        # Convert to DataFrame
        df = pd.DataFrame(results)

        # Sort by activation score
        df = df.sort_values('activation_score', ascending=False).reset_index(drop=True)

        return df

    def get_pathway_details(self, pathway_id: str) -> Dict:
        """
        Get detailed information about a pathway's activation.

        Args:
            pathway_id: KEGG pathway ID

        Returns:
            Detailed dictionary including all molecules and their contributions
        """
        score_info = self.compute_activation_score(pathway_id)

        # Add additional context
        pathway_info = self.mapper.get_pathway_info(pathway_id)
        score_info['pathway_info'] = pathway_info

        return score_info

    def summary(self) -> Dict:
        """Get summary statistics."""
        return {
            'n_metabolites_with_kegg': len(self.metabolite_data),
            'n_proteins_mapped': len(self.protein_data),
            'n_pathways_analyzed': len(self.mapper.get_all_pathways())
        }


def compute_embedding_scores(
    graph_path: str,
    model_path: str,
    mapper: 'PathwayMapper',
    output_path: Optional[str] = None
) -> pd.DataFrame:
    """
    Compute pathway activation scores using GNN embeddings.

    Uses pre-trained HGT model to compute enzyme-metabolite similarity scores
    for pathway ranking.

    Args:
        graph_path: Path to PyTorch Geometric graph data (.pt)
        model_path: Path to trained HGT model weights (.pth)
        mapper: PathwayMapper instance
        output_path: Optional path to save results CSV

    Returns:
        DataFrame with pathway embedding scores
    """
    try:
        import torch
        from src.model import HGT
    except ImportError as e:
        print(f"Warning: Could not import torch/model: {e}")
        return pd.DataFrame()

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    # Load graph data
    try:
        data = torch.load(graph_path, map_location=device)
    except FileNotFoundError:
        print(f"Warning: Graph file not found: {graph_path}")
        return pd.DataFrame()

    # Load model
    try:
        model = HGT(data.metadata(), 64, 64, 64, 4, 2).to(device)
        state_dict = torch.load(model_path, map_location=device)

        # Filter out incompatible keys (size mismatches)
        model_state = model.state_dict()
        filtered_state = {}
        for k, v in state_dict.items():
            if k in model_state and model_state[k].shape == v.shape:
                filtered_state[k] = v

        if len(filtered_state) > 0:
            model.load_state_dict(filtered_state, strict=False)
            print(f"  Loaded {len(filtered_state)}/{len(state_dict)} model parameters")
        else:
            print(f"  Warning: No compatible parameters found in model checkpoint")
            # Continue with randomly initialized model for basic similarity computation

        model.eval()
    except FileNotFoundError:
        print(f"Warning: Model file not found: {model_path}")
        return pd.DataFrame()
    except Exception as e:
        print(f"Warning: Could not load model: {e}")
        return pd.DataFrame()

    # Compute embeddings
    with torch.no_grad():
        x_dict = model(data.x_dict, data.edge_index_dict)
        enzyme_emb = x_dict['Enzyme'].cpu().numpy()
        metabolite_emb = x_dict['Metabolite'].cpu().numpy()

    # Load ID mappings from files
    graph_dir = Path(graph_path).parent
    project_root = graph_dir.parent.parent

    # Enzyme mapping: enzyme_idx -> uniprot_id
    enzyme_mapping_path = graph_dir / 'enzyme_string_mapping.csv'
    enzyme_id_to_idx = {}
    if enzyme_mapping_path.exists():
        enzyme_map_df = pd.read_csv(enzyme_mapping_path)
        for _, row in enzyme_map_df.iterrows():
            enzyme_id_to_idx[row['uniprot_id']] = int(row['enzyme_idx'])
        print(f"  Loaded enzyme mapping: {len(enzyme_id_to_idx)} entries")

    # Metabolite mapping: idx -> compound_id
    metabolite_mapping_path = project_root / 'data' / 'kegg' / 'metabolite_index_mapping.csv'
    metabolite_id_to_idx = {}
    if metabolite_mapping_path.exists():
        met_map_df = pd.read_csv(metabolite_mapping_path)
        for _, row in met_map_df.iterrows():
            metabolite_id_to_idx[row['compound_id']] = int(row['idx'])
        print(f"  Loaded metabolite mapping: {len(metabolite_id_to_idx)} entries")

    # Compute pathway scores based on embedding similarity
    results = []
    all_pathways = mapper.get_all_pathways()

    for pathway_id in all_pathways:
        pathway_name = mapper.get_pathway_name(pathway_id)
        pathway_proteins = mapper.get_pathway_proteins(pathway_id)
        pathway_metabolites = mapper.get_pathway_metabolites(pathway_id)

        # Compute average embedding similarity for pathway
        n_matched_enzymes = 0
        n_matched_metabolites = 0
        enzyme_indices = []
        metabolite_indices = []

        # Map pathway proteins to graph indices using mapping file
        for uniprot_id in pathway_proteins:
            if uniprot_id in enzyme_id_to_idx:
                idx = enzyme_id_to_idx[uniprot_id]
                if idx < enzyme_emb.shape[0]:
                    enzyme_indices.append(idx)
                    n_matched_enzymes += 1

        # Map pathway metabolites to graph indices using mapping file
        for kegg_id in pathway_metabolites:
            if kegg_id in metabolite_id_to_idx:
                idx = metabolite_id_to_idx[kegg_id]
                if idx < metabolite_emb.shape[0]:
                    metabolite_indices.append(idx)
                    n_matched_metabolites += 1

        # Compute embedding score
        embedding_score = 0.0
        if len(enzyme_indices) > 0 and len(metabolite_indices) > 0:
            # Compute cosine similarity between pathway enzymes and metabolites
            pathway_enz_emb = enzyme_emb[enzyme_indices]
            pathway_met_emb = metabolite_emb[metabolite_indices]

            # Normalize embeddings
            pathway_enz_emb_norm = pathway_enz_emb / (np.linalg.norm(pathway_enz_emb, axis=1, keepdims=True) + 1e-8)
            pathway_met_emb_norm = pathway_met_emb / (np.linalg.norm(pathway_met_emb, axis=1, keepdims=True) + 1e-8)

            # Average pairwise similarity
            similarity_matrix = pathway_enz_emb_norm @ pathway_met_emb_norm.T
            embedding_score = similarity_matrix.mean()
        elif len(enzyme_indices) > 0:
            # Only enzymes - use intra-pathway enzyme similarity
            pathway_enz_emb = enzyme_emb[enzyme_indices]
            pathway_enz_emb_norm = pathway_enz_emb / (np.linalg.norm(pathway_enz_emb, axis=1, keepdims=True) + 1e-8)
            if len(enzyme_indices) > 1:
                similarity_matrix = pathway_enz_emb_norm @ pathway_enz_emb_norm.T
                # Exclude diagonal
                mask = ~np.eye(similarity_matrix.shape[0], dtype=bool)
                embedding_score = similarity_matrix[mask].mean() * 0.5  # Lower weight
            else:
                embedding_score = 0.0

        results.append({
            'pathway_id': pathway_id,
            'pathway_name': pathway_name,
            'embedding_score': float(embedding_score),
            'n_enzymes_in_graph': n_matched_enzymes,
            'n_metabolites_in_graph': n_matched_metabolites,
            'n_enzymes_total': len(pathway_proteins),
            'n_metabolites_total': len(pathway_metabolites)
        })

    # Create DataFrame and normalize
    df = pd.DataFrame(results)
    if len(df) > 0 and df['embedding_score'].max() > 0:
        df['normalized_embedding_score'] = df['embedding_score'] / df['embedding_score'].max()
    else:
        df['normalized_embedding_score'] = 0.0

    df = df.sort_values('embedding_score', ascending=False).reset_index(drop=True)

    # Save if output path provided
    if output_path:
        df.to_csv(output_path, index=False)
        print(f"Saved embedding scores to: {output_path}")

    return df
