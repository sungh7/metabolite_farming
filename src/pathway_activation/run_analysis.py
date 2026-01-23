"""
Run Pathway Activation Analysis

Main script to execute pathway activation scoring analysis.
"""

import argparse
import pandas as pd
from pathlib import Path
import sys

# Add project root to path
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

from src.pathway_activation.pathway_mapper import PathwayMapper
from src.pathway_activation.activation_scorer import PathwayActivationScorer, compute_embedding_scores
from src.pathway_activation.protein_ranker import ProteinContributionRanker
from src.pathway_activation.evidence_reporter import EvidenceReporter


def load_data(data_dir: Path) -> tuple:
    """
    Load required data files.

    Args:
        data_dir: Path to data directory

    Returns:
        Tuple of (metabolomics_df, proteomics_df, kegg_path, uniprot_mapping_path, string_mapping_path)
    """
    # Metabolomics data
    metabolomics_path = data_dir / 'processed' / 'mtbls531_differential.csv'
    metabolomics_df = pd.read_csv(metabolomics_path)
    print(f"Loaded metabolomics data: {len(metabolomics_df)} rows")

    # Proteomics data
    proteomics_path = data_dir / 'processed' / 'pxd006989_mapped.csv'
    proteomics_df = pd.read_csv(proteomics_path)
    print(f"Loaded proteomics data: {len(proteomics_df)} rows")

    # KEGG pathway data
    kegg_path = project_root / 'results' / 'collected_data' / 'kegg_pathway_data.json'

    # ID mappings
    uniprot_mapping_path = data_dir / 'kegg' / 'kegg_uniprot_mapping.csv'
    string_mapping_path = data_dir / 'processed' / 'enzyme_string_mapping.csv'

    return metabolomics_df, proteomics_df, kegg_path, uniprot_mapping_path, string_mapping_path


def run_analysis(
    data_dir: str = None,
    output_dir: str = None,
    generate_figures: bool = True
):
    """
    Run the complete pathway activation analysis.

    Args:
        data_dir: Path to data directory
        output_dir: Path to output directory
        generate_figures: Whether to generate figures
    """
    # Set default paths
    if data_dir is None:
        data_dir = project_root / 'data'
    else:
        data_dir = Path(data_dir)

    if output_dir is None:
        output_dir = project_root / 'results' / 'pathway_activation'
    else:
        output_dir = Path(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("Pathway Activation Analysis")
    print("=" * 60)
    print()

    # Load data
    print("Loading data...")
    metabolomics_df, proteomics_df, kegg_path, uniprot_mapping_path, string_mapping_path = load_data(data_dir)

    # Initialize PathwayMapper
    print("\nInitializing pathway mapper...")
    mapper = PathwayMapper(
        kegg_pathway_path=str(kegg_path),
        kegg_uniprot_mapping_path=str(uniprot_mapping_path),
        enzyme_string_mapping_path=str(string_mapping_path)
    )

    mapper_summary = mapper.summary()
    print(f"  Pathways: {mapper_summary['n_pathways']}")
    print(f"  KEGG genes mapped: {mapper_summary['n_kegg_genes']}")
    print(f"  STRING proteins mapped: {mapper_summary['n_string_proteins']}")
    print(f"  Compound mappings: {mapper_summary['n_compound_mappings']}")

    # Initialize PathwayActivationScorer
    print("\nInitializing activation scorer...")
    scorer = PathwayActivationScorer(
        metabolomics_df=metabolomics_df,
        proteomics_df=proteomics_df,
        mapper=mapper
    )

    scorer_summary = scorer.summary()
    print(f"  Metabolites with KEGG ID: {scorer_summary['n_metabolites_with_kegg']}")
    print(f"  Proteins mapped: {scorer_summary['n_proteins_mapped']}")

    # Compute pathway activation scores
    print("\nComputing pathway activation scores...")
    pathway_scores_df = scorer.rank_all_pathways()
    print(f"  Analyzed {len(pathway_scores_df)} pathways")

    # Display top pathways
    print("\n" + "=" * 60)
    print("Top Activated Pathways")
    print("=" * 60)
    for i, row in pathway_scores_df.head(5).iterrows():
        print(f"\n{i+1}. {row['pathway_name']} ({row['pathway_id']})")
        print(f"   Activation Score: {row['normalized_score']:.3f}")
        print(f"   Metabolite FC: {row['metabolite_fc']:.3f} (n={row['n_metabolites_detected']})")
        print(f"   Enzyme FC: {row['enzyme_fc']:.3f} (n={row['n_enzymes_detected']})")
        print(f"   Coverage: M={row['metabolite_coverage']:.0%}, E={row['enzyme_coverage']:.0%}")

    # Initialize ProteinContributionRanker
    print("\n" + "=" * 60)
    print("Computing protein contributions...")
    print("=" * 60)
    ranker = ProteinContributionRanker(
        proteomics_df=proteomics_df,
        mapper=mapper
    )

    protein_contributions = ranker.rank_all_pathway_proteins()

    # Display top proteins for top pathway
    if len(pathway_scores_df) > 0:
        top_pathway_id = pathway_scores_df.iloc[0]['pathway_id']
        top_proteins = ranker.get_top_contributors(top_pathway_id, 5)

        print(f"\nTop contributors to {top_pathway_id}:")
        for i, protein in enumerate(top_proteins):
            name = protein.get('enzyme_name') or protein.get('uniprot_id', 'Unknown')
            print(f"  {i+1}. {name}: score={protein.get('normalized_contribution', 0):.3f}, "
                  f"FC={protein.get('log2fc', 0):.2f}, p={protein.get('p_value', 1):.2e}")

    # Compute GNN embedding scores
    print("\n" + "=" * 60)
    print("Computing GNN embedding scores...")
    print("=" * 60)

    graph_path = project_root / 'data' / 'processed' / 'strict_bipartite_v2.pt'
    model_path = project_root / 'data' / 'models' / 'refined_hgt_strict.pth'
    embedding_output_path = output_dir / 'embedding_scores.csv'

    embedding_scores_df = compute_embedding_scores(
        graph_path=str(graph_path),
        model_path=str(model_path),
        mapper=mapper,
        output_path=str(embedding_output_path)
    )

    if len(embedding_scores_df) > 0:
        print(f"  Computed embedding scores for {len(embedding_scores_df)} pathways")
        print("\n  Top pathways by embedding score:")
        for i, row in embedding_scores_df.head(5).iterrows():
            print(f"    {i+1}. {row['pathway_name']} ({row['pathway_id']})")
            print(f"       Score: {row['normalized_embedding_score']:.3f} "
                  f"(E={row['n_enzymes_in_graph']}, M={row['n_metabolites_in_graph']})")
    else:
        print("  Warning: Could not compute embedding scores (missing model or data)")

    # Generate report
    print("\n" + "=" * 60)
    print("Generating reports...")
    print("=" * 60)
    reporter = EvidenceReporter(output_dir=str(output_dir))

    # Export tables
    exported_files = reporter.export_tables(pathway_scores_df, protein_contributions)
    for name, path in exported_files.items():
        print(f"  Exported: {path}")

    # Generate summary report
    results = {
        'pathway_scores_df': pathway_scores_df,
        'protein_contributions': protein_contributions
    }
    report_path = reporter.summary_report(results)
    print(f"  Generated: {report_path}")

    # Generate figures
    if generate_figures:
        print("\nGenerating figures...")
        figure_paths = reporter.generate_all_figures(
            pathway_scores_df,
            protein_contributions
        )
        for path in figure_paths:
            print(f"  Generated: {path}")

    print("\n" + "=" * 60)
    print("Analysis Complete!")
    print("=" * 60)
    print(f"\nOutput directory: {output_dir}")

    # Return results for programmatic use
    return {
        'pathway_scores': pathway_scores_df,
        'protein_contributions': protein_contributions,
        'embedding_scores': embedding_scores_df,
        'mapper': mapper,
        'scorer': scorer,
        'ranker': ranker
    }


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Run pathway activation analysis for ethylene treatment'
    )
    parser.add_argument(
        '--data-dir',
        type=str,
        default=None,
        help='Path to data directory'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default=None,
        help='Path to output directory'
    )
    parser.add_argument(
        '--no-figures',
        action='store_true',
        help='Skip figure generation'
    )

    args = parser.parse_args()

    run_analysis(
        data_dir=args.data_dir,
        output_dir=args.output_dir,
        generate_figures=not args.no_figures
    )


if __name__ == '__main__':
    main()
