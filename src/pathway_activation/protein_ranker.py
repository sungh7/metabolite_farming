"""
Protein Contribution Ranker Module

Ranks proteins by their contribution to pathway activation.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Optional
from pathlib import Path

from .pathway_mapper import PathwayMapper


# Known enzyme names for isoflavonoid biosynthesis pathway
ENZYME_NAMES = {
    # Isoflavonoid biosynthesis key enzymes
    'Q9M6D6': {'name': 'IFS', 'full_name': 'Isoflavone synthase', 'ec': '1.14.14.87'},
    'Q53B70': {'name': 'CHI', 'full_name': 'Chalcone isomerase', 'ec': '5.5.1.6'},
    'Q9LLQ8': {'name': 'HIDH', 'full_name': '2-hydroxyisoflavanone dehydratase', 'ec': '4.2.1.105'},
    'Q9XHC6': {'name': 'PAL', 'full_name': 'Phenylalanine ammonia-lyase', 'ec': '4.3.1.24'},
    'Q8L7J4': {'name': '4CL', 'full_name': '4-coumarate-CoA ligase', 'ec': '6.2.1.12'},
    'B0M1A4': {'name': 'C4H', 'full_name': 'Cinnamate 4-hydroxylase', 'ec': '1.14.14.91'},
    'Q5ECI5': {'name': 'CHS', 'full_name': 'Chalcone synthase', 'ec': '2.3.1.74'},
    'I1MQD2': {'name': 'CHR', 'full_name': 'Chalcone reductase', 'ec': '1.1.1.34'},
    'A0A0R0GNE3': {'name': 'HI4OMT', 'full_name': 'Hydroxyisoflavanone 4-O-methyltransferase', 'ec': '2.1.1.212'},
    'A0A0R0E8N5': {'name': 'IF7GT', 'full_name': 'Isoflavone 7-O-glucosyltransferase', 'ec': '2.4.1.170'},

    # Extended soybean isoflavonoid/flavonoid pathway enzymes
    'I1M622': {'name': 'IFS2', 'full_name': 'Isoflavone synthase 2', 'ec': '1.14.14.87'},
    'A7BIC9': {'name': 'CYP81E', 'full_name': 'Cytochrome P450 81E subfamily', 'ec': '1.14.-.-'},
    'I1JNN8': {'name': 'UGT', 'full_name': 'UDP-glycosyltransferase', 'ec': '2.4.1.-'},
    'C6TDP1': {'name': 'OMT', 'full_name': 'O-methyltransferase', 'ec': '2.1.1.-'},
    'C6TNS8': {'name': 'F3H', 'full_name': 'Flavanone 3-hydroxylase', 'ec': '1.14.11.9'},
    'I1KM16': {'name': 'DFR', 'full_name': 'Dihydroflavonol 4-reductase', 'ec': '1.1.1.219'},
    'A0A0R0KSM2': {'name': 'CAD', 'full_name': 'Cinnamyl alcohol dehydrogenase', 'ec': '1.1.1.195'},
    'I1MNR6': {'name': 'CCR', 'full_name': 'Cinnamoyl-CoA reductase', 'ec': '1.2.1.44'},
    'A0A0R4J4D2': {'name': 'HCT', 'full_name': 'Hydroxycinnamoyl-CoA transferase', 'ec': '2.3.1.133'},
    'C6SVC4': {'name': 'COMT', 'full_name': 'Caffeic acid O-methyltransferase', 'ec': '2.1.1.68'},
    'I1MG66': {'name': 'CCoAOMT', 'full_name': 'Caffeoyl-CoA O-methyltransferase', 'ec': '2.1.1.104'},

    # Additional phenylpropanoid pathway enzymes
    'Q43153': {'name': 'PAL1', 'full_name': 'Phenylalanine ammonia-lyase 1', 'ec': '4.3.1.24'},
    'P27991': {'name': 'CHS1', 'full_name': 'Chalcone synthase 1', 'ec': '2.3.1.74'},
    'P48408': {'name': 'CHS7', 'full_name': 'Chalcone synthase 7', 'ec': '2.3.1.74'},
    'O04058': {'name': 'CHI1A', 'full_name': 'Chalcone isomerase 1A', 'ec': '5.5.1.6'},
    'P93510': {'name': 'CHI1B', 'full_name': 'Chalcone isomerase 1B', 'ec': '5.5.1.6'},
    'Q9SWR5': {'name': 'IFS1', 'full_name': 'Isoflavone synthase 1', 'ec': '1.14.14.87'},
    'Q9XHG7': {'name': 'IOMT', 'full_name': 'Isoflavone O-methyltransferase', 'ec': '2.1.1.46'},
    'Q84QV3': {'name': 'I2H', 'full_name': 'Isoflavone 2-hydroxylase', 'ec': '1.14.14.89'},
    'Q9FVE0': {'name': 'VR', 'full_name': 'Vestitone reductase', 'ec': '1.1.1.348'},
    'Q5NTD4': {'name': 'DMID', 'full_name': "7,2'-dihydroxy-4'-methoxyisoflavanol dehydratase", 'ec': '4.2.1.139'},

    # Glucosyltransferases
    'Q9XJ28': {'name': 'UGT73C', 'full_name': 'UDP-glucosyltransferase 73C', 'ec': '2.4.1.-'},
    'Q9FVF4': {'name': 'IF7MaT', 'full_name': 'Isoflavone 7-O-glucoside 6-O-malonyltransferase', 'ec': '2.3.1.159'},
}


class ProteinContributionRanker:
    """
    Ranks proteins by their contribution to pathway activation.
    """

    def __init__(
        self,
        proteomics_df: pd.DataFrame,
        mapper: PathwayMapper,
        centrality_weights: Optional[Dict[str, float]] = None
    ):
        """
        Initialize the ranker.

        Args:
            proteomics_df: DataFrame with columns [STRING_ID, Log2FC, P_Value]
            mapper: PathwayMapper instance
            centrality_weights: Optional dictionary of protein centrality scores
        """
        self.proteomics_df = proteomics_df.copy()
        self.mapper = mapper
        self.centrality_weights = centrality_weights or {}

        self._preprocess_proteomics()

    def _preprocess_proteomics(self):
        """Preprocess proteomics data."""
        self.protein_data = {}

        for _, row in self.proteomics_df.iterrows():
            string_id = row.get('STRING_ID')
            if pd.isna(string_id):
                continue

            string_id = str(string_id).strip()
            uniprot_id = self.mapper.string_to_uniprot_id(string_id)

            if uniprot_id is None:
                if '.' in string_id:
                    uniprot_id = string_id.split('.')[-1]
                else:
                    continue

            log2fc = float(row['Log2FC']) if pd.notna(row['Log2FC']) else 0.0
            p_value = float(row['P_Value']) if pd.notna(row['P_Value']) else 1.0

            # Get gene names and protein names if available
            gene_names = row.get('Gene names', '')
            protein_names = row.get('Protein names', '')

            self.protein_data[uniprot_id] = {
                'log2fc': log2fc,
                'p_value': p_value,
                'string_id': string_id,
                'gene_names': gene_names if pd.notna(gene_names) else '',
                'protein_names': protein_names if pd.notna(protein_names) else ''
            }

    def compute_contribution(
        self,
        pathway_id: str,
        protein_id: str
    ) -> float:
        """
        Compute protein contribution score.

        Contribution = |Log2FC| * -log10(p_value) * centrality

        Args:
            pathway_id: KEGG pathway ID
            protein_id: UniProt ID

        Returns:
            Contribution score
        """
        if protein_id not in self.protein_data:
            return 0.0

        data = self.protein_data[protein_id]

        # Base components
        fc_component = abs(data['log2fc'])
        p_value = max(data['p_value'], 1e-300)
        sig_component = -np.log10(p_value)

        # Centrality (default to 1.0 if not available)
        centrality = self.centrality_weights.get(protein_id, 1.0)

        # Contribution score
        contribution = fc_component * sig_component * centrality

        return contribution

    def rank_proteins(self, pathway_id: str) -> pd.DataFrame:
        """
        Rank proteins by contribution within a pathway.

        Args:
            pathway_id: KEGG pathway ID

        Returns:
            DataFrame with protein rankings
        """
        pathway_proteins = self.mapper.get_pathway_proteins(pathway_id)
        results = []

        for uniprot_id in pathway_proteins:
            if uniprot_id not in self.protein_data:
                continue

            data = self.protein_data[uniprot_id]
            contribution = self.compute_contribution(pathway_id, uniprot_id)

            # Get enzyme name if known
            enzyme_info = ENZYME_NAMES.get(uniprot_id, {})

            results.append({
                'uniprot_id': uniprot_id,
                'string_id': data['string_id'],
                'log2fc': data['log2fc'],
                'p_value': data['p_value'],
                'contribution_score': contribution,
                'enzyme_name': enzyme_info.get('name', ''),
                'enzyme_full_name': enzyme_info.get('full_name', data['protein_names']),
                'ec_number': enzyme_info.get('ec', ''),
                'gene_names': data['gene_names'],
                'protein_names': data['protein_names'],
                'centrality': self.centrality_weights.get(uniprot_id, 1.0)
            })

        # Convert to DataFrame and sort
        if not results:
            return pd.DataFrame(columns=[
                'uniprot_id', 'string_id', 'log2fc', 'p_value',
                'contribution_score', 'enzyme_name', 'enzyme_full_name',
                'ec_number', 'gene_names', 'protein_names', 'centrality'
            ])

        df = pd.DataFrame(results)
        df = df.sort_values('contribution_score', ascending=False).reset_index(drop=True)

        # Normalize contribution scores
        max_contribution = df['contribution_score'].max()
        if max_contribution > 0:
            df['normalized_contribution'] = df['contribution_score'] / max_contribution
        else:
            df['normalized_contribution'] = 0.0

        return df

    def rank_all_pathway_proteins(self) -> Dict[str, pd.DataFrame]:
        """
        Rank proteins for all pathways.

        Returns:
            Dictionary mapping pathway_id to protein ranking DataFrame
        """
        all_pathways = self.mapper.get_all_pathways()
        results = {}

        for pathway_id in all_pathways:
            results[pathway_id] = self.rank_proteins(pathway_id)

        return results

    def get_top_contributors(
        self,
        pathway_id: str,
        n: int = 10
    ) -> List[Dict]:
        """
        Get top N contributing proteins for a pathway.

        Args:
            pathway_id: KEGG pathway ID
            n: Number of top proteins to return

        Returns:
            List of dictionaries with protein info
        """
        df = self.rank_proteins(pathway_id)
        if len(df) == 0:
            return []

        top_df = df.head(n)
        return top_df.to_dict('records')

    def compute_pathway_protein_summary(self, pathway_id: str) -> Dict:
        """
        Compute summary statistics for proteins in a pathway.

        Args:
            pathway_id: KEGG pathway ID

        Returns:
            Summary dictionary
        """
        df = self.rank_proteins(pathway_id)

        if len(df) == 0:
            return {
                'pathway_id': pathway_id,
                'pathway_name': self.mapper.get_pathway_name(pathway_id),
                'n_detected': 0,
                'n_significant': 0,
                'mean_log2fc': 0.0,
                'mean_contribution': 0.0,
                'top_contributors': []
            }

        significant_mask = df['p_value'] < 0.05

        return {
            'pathway_id': pathway_id,
            'pathway_name': self.mapper.get_pathway_name(pathway_id),
            'n_detected': len(df),
            'n_significant': significant_mask.sum(),
            'mean_log2fc': df['log2fc'].mean(),
            'mean_contribution': df['contribution_score'].mean(),
            'top_contributors': self.get_top_contributors(pathway_id, 5)
        }

    def export_contributions(
        self,
        output_path: str,
        top_n_per_pathway: int = 20
    ):
        """
        Export protein contributions to CSV.

        Args:
            output_path: Output file path
            top_n_per_pathway: Number of top proteins per pathway
        """
        all_results = []

        for pathway_id in self.mapper.get_all_pathways():
            df = self.rank_proteins(pathway_id)
            if len(df) == 0:
                continue

            df = df.head(top_n_per_pathway)
            df['pathway_id'] = pathway_id
            df['pathway_name'] = self.mapper.get_pathway_name(pathway_id)
            all_results.append(df)

        if all_results:
            combined_df = pd.concat(all_results, ignore_index=True)
            combined_df.to_csv(output_path, index=False)
        else:
            pd.DataFrame().to_csv(output_path, index=False)
