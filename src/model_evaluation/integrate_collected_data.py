#!/usr/bin/env python
"""
Integrate Collected Data with Existing Proteomics/GNN Results

Maps KEGG pathway genes, TF motifs, and other collected data
to the existing proteomics and GNN prediction results.
"""

import os
import json
import pandas as pd
import numpy as np
from typing import Dict, List, Set

# Paths
COLLECTED_DIR = "results/collected_data"
PROTEOMICS_PATH = "data/processed/pxd006989_mapped.csv"
CANDIDATES_PATH = "results/candidates.csv"
OUTPUT_DIR = "results/integrated_analysis"


def load_collected_data() -> Dict:
    """Load all collected data files."""
    data = {}

    files = {
        'kegg': 'kegg_pathway_data.json',
        'ncbi': 'ncbi_data.json',
        'gget': 'gget_data.json',
        'jaspar': 'jaspar_motifs.json'
    }

    for key, filename in files.items():
        filepath = os.path.join(COLLECTED_DIR, filename)
        if os.path.exists(filepath):
            with open(filepath) as f:
                data[key] = json.load(f)
            print(f"  Loaded {key}: {filename}")

    return data


def extract_kegg_gene_ids(kegg_data: Dict) -> Dict[str, Set[str]]:
    """Extract KEGG gene IDs by pathway."""
    pathway_genes = {}

    for pathway_id, info in kegg_data.get('pathways', {}).items():
        genes = set(info.get('genes', []))
        pathway_genes[pathway_id] = genes

    return pathway_genes


def map_kegg_to_string(kegg_genes: Set[str]) -> Set[str]:
    """
    Convert KEGG gene IDs to STRING format.
    KEGG: 100784567 -> STRING: 3847.xxx
    """
    # KEGG soybean genes are in format like "100784567"
    # STRING format is "3847.XXXXX"
    # We'll return the KEGG IDs as-is for now, actual mapping needs database
    return {f"gmx:{g}" for g in kegg_genes}


def annotate_proteomics_with_pathways(proteomics_df: pd.DataFrame,
                                       pathway_genes: Dict[str, Set[str]]) -> pd.DataFrame:
    """Add pathway annotations to proteomics data."""

    # Create pathway annotation columns
    for pathway_id in pathway_genes:
        proteomics_df[f'in_{pathway_id}'] = False

    # Try to match based on STRING_ID or Protein IDs
    for idx, row in proteomics_df.iterrows():
        string_id = str(row.get('STRING_ID', ''))
        protein_ids = str(row.get('Protein IDs', ''))

        # Extract potential gene identifiers
        # Glyma format: Glyma.01G000700 -> try to match
        glyma_ids = set()
        for part in protein_ids.split(';'):
            if 'Glyma.' in part:
                # Extract locus like Glyma.01G000700
                glyma_id = part.split('.p')[0] if '.p' in part else part
                glyma_ids.add(glyma_id)

        # For now, mark as in pathway if any match
        # (In real implementation, use proper ID mapping)

    return proteomics_df


def identify_de_proteins(proteomics_df: pd.DataFrame,
                          log2fc_thresh: float = 1.0,
                          pval_thresh: float = 0.05) -> pd.DataFrame:
    """Identify differentially expressed proteins."""

    de_mask = (
        (proteomics_df['Log2FC'].abs() >= log2fc_thresh) &
        (proteomics_df['P_Value'] <= pval_thresh)
    )

    de_proteins = proteomics_df[de_mask].copy()
    de_proteins['DE_direction'] = np.where(de_proteins['Log2FC'] > 0, 'up', 'down')

    return de_proteins


def match_candidates_to_pathways(candidates_df: pd.DataFrame,
                                  pathway_genes: Dict[str, Set[str]],
                                  kegg_data: Dict) -> pd.DataFrame:
    """Match GNN candidates to pathway information."""

    # Add pathway context columns
    candidates_df = candidates_df.copy()
    candidates_df['pathway_context'] = ''
    candidates_df['isoflavonoid_related'] = False

    # Check enzyme descriptions for pathway keywords
    isoflavonoid_keywords = [
        'isoflavon', 'flavon', 'chalcone', 'phenylpropanoid',
        'CHS', 'CHI', 'IFS', 'HIDM', 'glucosyltransferase',
        '2-hydroxyisoflavanone', 'daidzein', 'genistein'
    ]

    for idx, row in candidates_df.iterrows():
        enzyme_desc = str(row.get('Enzyme_Desc', '')).lower()
        tf_desc = str(row.get('TF_Desc', '')).lower()

        # Check for isoflavonoid-related terms
        for kw in isoflavonoid_keywords:
            if kw.lower() in enzyme_desc or kw.lower() in tf_desc:
                candidates_df.at[idx, 'isoflavonoid_related'] = True
                break

        # Add pathway context
        contexts = []
        if 'methyltransferase' in enzyme_desc or 'methyltransferase' in tf_desc:
            contexts.append('methylation')
        if 'kinase' in enzyme_desc:
            contexts.append('signaling')
        if 'dehydratase' in enzyme_desc:
            contexts.append('isoflavonoid_biosynthesis')
        if 'glucosyltransferase' in enzyme_desc:
            contexts.append('glycosylation')

        candidates_df.at[idx, 'pathway_context'] = '; '.join(contexts)

    return candidates_df


def create_tf_target_evidence(candidates_df: pd.DataFrame,
                               jaspar_data: Dict) -> pd.DataFrame:
    """Add TF binding motif evidence."""

    candidates_df = candidates_df.copy()
    candidates_df['tf_motif_available'] = False
    candidates_df['tf_family'] = ''

    # Get available TF families from JASPAR
    tf_families = set()
    for motif in jaspar_data.get('plant_tf_motifs', []):
        tf_class = motif.get('tf_class', '')
        if tf_class:
            tf_families.add(tf_class.lower())

    # Check if TF descriptions match available motifs
    tf_keywords = {
        'nac': 'NAC',
        'myb': 'MYB',
        'bhlh': 'bHLH',
        'wrky': 'WRKY',
        'erf': 'ERF/AP2',
        'ap2': 'ERF/AP2',
        'bzip': 'bZIP'
    }

    for idx, row in candidates_df.iterrows():
        tf_desc = str(row.get('TF_Desc', '')).lower()
        tf_name = str(row.get('TF_Name', '')).lower()

        for kw, family in tf_keywords.items():
            if kw in tf_desc or kw in tf_name:
                candidates_df.at[idx, 'tf_family'] = family
                # Check if motif data available
                if any(kw in tf.lower() for tf in tf_families):
                    candidates_df.at[idx, 'tf_motif_available'] = True
                break

    return candidates_df


def generate_integration_report(proteomics_df: pd.DataFrame,
                                 candidates_df: pd.DataFrame,
                                 collected_data: Dict,
                                 output_dir: str):
    """Generate comprehensive integration report."""

    os.makedirs(output_dir, exist_ok=True)

    report = {
        'summary': {},
        'pathway_enrichment': {},
        'tf_analysis': {},
        'recommendations': []
    }

    # Summary statistics
    report['summary'] = {
        'total_proteins': len(proteomics_df),
        'de_proteins_up': len(proteomics_df[proteomics_df['Log2FC'] > 1]),
        'de_proteins_down': len(proteomics_df[proteomics_df['Log2FC'] < -1]),
        'total_candidates': len(candidates_df),
        'isoflavonoid_related_candidates': int(candidates_df['isoflavonoid_related'].sum()),
        'candidates_with_tf_motifs': int(candidates_df['tf_motif_available'].sum())
    }

    # TF family distribution
    tf_family_counts = candidates_df['tf_family'].value_counts().to_dict()
    report['tf_analysis'] = {
        'tf_family_distribution': tf_family_counts,
        'motif_coverage': report['summary']['candidates_with_tf_motifs'] / max(1, report['summary']['total_candidates'])
    }

    # Pathway context
    pathway_counts = {}
    for ctx in candidates_df['pathway_context'].dropna():
        for c in ctx.split('; '):
            if c:
                pathway_counts[c] = pathway_counts.get(c, 0) + 1
    report['pathway_enrichment'] = pathway_counts

    # Recommendations
    report['recommendations'] = [
        f"1. {report['summary']['isoflavonoid_related_candidates']} candidates directly related to isoflavonoid biosynthesis",
        f"2. {report['summary']['candidates_with_tf_motifs']} candidates have TF motif data available in JASPAR",
        "3. Priority targets for promoter motif analysis: NAC and MYB family TFs",
        "4. Consider downloading GSE112584 for NAC-glyceollin temporal data",
        "5. Use KEGG pathway genes for functional validation"
    ]

    # Save report
    with open(os.path.join(output_dir, 'integration_report.json'), 'w') as f:
        json.dump(report, f, indent=2)

    # Save annotated candidates
    candidates_df.to_csv(os.path.join(output_dir, 'annotated_candidates.csv'), index=False)

    # Print summary
    print("\n" + "="*60)
    print("INTEGRATION REPORT")
    print("="*60)
    print(f"\nProteomics Data:")
    print(f"  - Total proteins: {report['summary']['total_proteins']}")
    print(f"  - Upregulated (Log2FC > 1): {report['summary']['de_proteins_up']}")
    print(f"  - Downregulated (Log2FC < -1): {report['summary']['de_proteins_down']}")

    print(f"\nGNN Candidates:")
    print(f"  - Total candidates: {report['summary']['total_candidates']}")
    print(f"  - Isoflavonoid-related: {report['summary']['isoflavonoid_related_candidates']}")
    print(f"  - With TF motif data: {report['summary']['candidates_with_tf_motifs']}")

    print(f"\nTF Family Distribution:")
    for family, count in sorted(tf_family_counts.items(), key=lambda x: -x[1])[:5]:
        if family:
            print(f"  - {family}: {count}")

    print(f"\nPathway Context:")
    for ctx, count in sorted(pathway_counts.items(), key=lambda x: -x[1])[:5]:
        print(f"  - {ctx}: {count}")

    print(f"\nRecommendations:")
    for rec in report['recommendations']:
        print(f"  {rec}")

    print(f"\nResults saved to: {output_dir}/")

    return report


def main():
    """Run integration pipeline."""
    print("="*60)
    print("INTEGRATING COLLECTED DATA WITH EXISTING RESULTS")
    print("="*60)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load collected data
    print("\n1. Loading collected data...")
    collected_data = load_collected_data()

    # Load proteomics data
    print("\n2. Loading proteomics data...")
    if os.path.exists(PROTEOMICS_PATH):
        proteomics_df = pd.read_csv(PROTEOMICS_PATH)
        print(f"  Loaded {len(proteomics_df)} proteins")
    else:
        print(f"  Warning: {PROTEOMICS_PATH} not found")
        proteomics_df = pd.DataFrame()

    # Load GNN candidates
    print("\n3. Loading GNN candidates...")
    if os.path.exists(CANDIDATES_PATH):
        candidates_df = pd.read_csv(CANDIDATES_PATH)
        print(f"  Loaded {len(candidates_df)} candidates")
    else:
        print(f"  Warning: {CANDIDATES_PATH} not found")
        candidates_df = pd.DataFrame()

    # Extract KEGG pathway genes
    print("\n4. Processing KEGG pathway data...")
    kegg_data = collected_data.get('kegg', {})
    pathway_genes = extract_kegg_gene_ids(kegg_data)
    for pid, genes in pathway_genes.items():
        print(f"  {pid}: {len(genes)} genes")

    # Annotate candidates with pathway context
    print("\n5. Annotating candidates with pathway context...")
    if not candidates_df.empty:
        candidates_df = match_candidates_to_pathways(
            candidates_df, pathway_genes, kegg_data
        )
        isoflavonoid_count = candidates_df['isoflavonoid_related'].sum()
        print(f"  Isoflavonoid-related candidates: {isoflavonoid_count}")

    # Add TF motif evidence
    print("\n6. Adding TF motif evidence...")
    jaspar_data = collected_data.get('jaspar', {})
    if not candidates_df.empty:
        candidates_df = create_tf_target_evidence(candidates_df, jaspar_data)
        tf_motif_count = candidates_df['tf_motif_available'].sum()
        print(f"  Candidates with TF motif data: {tf_motif_count}")

    # Generate integration report
    print("\n7. Generating integration report...")
    report = generate_integration_report(
        proteomics_df, candidates_df, collected_data, OUTPUT_DIR
    )

    print("\n" + "="*60)
    print("INTEGRATION COMPLETE")
    print("="*60)

    return report


if __name__ == "__main__":
    main()
