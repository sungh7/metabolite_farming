"""
Evidence Reporter Module

Generates visualizations and reports for pathway activation analysis.
"""

import json
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional
from datetime import datetime


class EvidenceReporter:
    """
    Generates visualizations and reports for pathway activation analysis.
    """

    def __init__(self, output_dir: str):
        """
        Initialize the reporter.

        Args:
            output_dir: Directory for output files
        """
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.figures_dir = self.output_dir / 'figures'
        self.figures_dir.mkdir(exist_ok=True)

    def pathway_heatmap(
        self,
        scores_df: pd.DataFrame,
        output_path: Optional[str] = None
    ) -> str:
        """
        Generate pathway activation heatmap.

        Args:
            scores_df: DataFrame with pathway scores
            output_path: Optional output path

        Returns:
            Path to saved figure
        """
        try:
            import matplotlib.pyplot as plt
            import matplotlib.colors as mcolors
        except ImportError:
            print("Warning: matplotlib not available, skipping heatmap")
            return ""

        if output_path is None:
            output_path = str(self.figures_dir / 'pathway_activation_heatmap.png')

        # Prepare data
        df = scores_df.copy()
        if len(df) == 0:
            return ""

        # Sort by activation score
        df = df.sort_values('activation_score', ascending=True)

        # Create figure
        fig, ax = plt.subplots(figsize=(12, max(6, len(df) * 0.5)))

        # Create heatmap data
        data = df[['normalized_score', 'metabolite_coverage', 'enzyme_coverage']].values

        # Color mapping
        cmap = plt.cm.RdYlBu_r

        # Plot bars
        y_positions = np.arange(len(df))
        bar_height = 0.6

        # Activation score bars
        colors = [cmap(score) for score in df['normalized_score'].values]
        bars = ax.barh(y_positions, df['normalized_score'].values,
                       height=bar_height, color=colors, edgecolor='black', linewidth=0.5)

        # Add pathway names
        pathway_labels = [
            f"{row['pathway_name']} ({row['pathway_id']})"
            for _, row in df.iterrows()
        ]
        ax.set_yticks(y_positions)
        ax.set_yticklabels(pathway_labels, fontsize=10)

        # Add score annotations
        for i, (score, row) in enumerate(zip(df['normalized_score'].values, df.itertuples())):
            text_x = score + 0.02 if score < 0.8 else score - 0.15
            text_color = 'black' if score < 0.8 else 'white'
            ax.text(text_x, i, f'{score:.2f}', va='center', ha='left',
                    fontsize=9, color=text_color, fontweight='bold')

            # Add coverage info
            coverage_text = f'M:{row.metabolite_coverage:.0%} E:{row.enzyme_coverage:.0%}'
            ax.text(1.05, i, coverage_text, va='center', ha='left',
                    fontsize=8, color='gray', transform=ax.get_yaxis_transform())

        ax.set_xlabel('Normalized Activation Score', fontsize=12)
        ax.set_title('Pathway Activation Scores (Ethylene Treatment)', fontsize=14, fontweight='bold')
        ax.set_xlim(0, 1.0)

        # Add colorbar
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=mcolors.Normalize(vmin=0, vmax=1))
        sm.set_array([])
        cbar = plt.colorbar(sm, ax=ax, shrink=0.5, aspect=20)
        cbar.set_label('Activation Level', fontsize=10)

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

        return output_path

    def protein_contribution_plot(
        self,
        pathway_id: str,
        contributions_df: pd.DataFrame,
        output_path: Optional[str] = None
    ) -> str:
        """
        Generate protein contribution plot for a pathway.

        Args:
            pathway_id: KEGG pathway ID
            contributions_df: DataFrame with protein contributions
            output_path: Optional output path

        Returns:
            Path to saved figure
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("Warning: matplotlib not available, skipping contribution plot")
            return ""

        if output_path is None:
            output_path = str(self.figures_dir / f'protein_contributions_{pathway_id}.png')

        if len(contributions_df) == 0:
            return ""

        # Take top 15 proteins
        df = contributions_df.head(15).copy()
        df = df.iloc[::-1]  # Reverse for plotting

        fig, ax = plt.subplots(figsize=(10, max(6, len(df) * 0.4)))

        # Create labels
        labels = []
        for _, row in df.iterrows():
            if row['enzyme_name']:
                label = f"{row['enzyme_name']} ({row['uniprot_id']})"
            else:
                label = row['uniprot_id']
            labels.append(label)

        y_positions = np.arange(len(df))

        # Color by Log2FC direction
        colors = ['#e74c3c' if fc > 0 else '#3498db' for fc in df['log2fc'].values]

        bars = ax.barh(y_positions, df['normalized_contribution'].values,
                       color=colors, edgecolor='black', linewidth=0.5, height=0.6)

        ax.set_yticks(y_positions)
        ax.set_yticklabels(labels, fontsize=9)

        # Add Log2FC annotations
        for i, row in enumerate(df.itertuples()):
            ax.text(row.normalized_contribution + 0.02, i,
                    f'FC:{row.log2fc:.2f} p:{row.p_value:.2e}',
                    va='center', ha='left', fontsize=8, color='gray')

        ax.set_xlabel('Normalized Contribution Score', fontsize=11)
        ax.set_title(f'Protein Contributions to {pathway_id}', fontsize=12, fontweight='bold')
        ax.set_xlim(0, 1.3)

        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='#e74c3c', edgecolor='black', label='Upregulated'),
            Patch(facecolor='#3498db', edgecolor='black', label='Downregulated')
        ]
        ax.legend(handles=legend_elements, loc='lower right', fontsize=9)

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

        return output_path

    def fc_comparison_plot(
        self,
        scores_df: pd.DataFrame,
        output_path: Optional[str] = None
    ) -> str:
        """
        Generate metabolite vs enzyme FC comparison plot.

        Args:
            scores_df: DataFrame with pathway scores
            output_path: Optional output path

        Returns:
            Path to saved figure
        """
        try:
            import matplotlib.pyplot as plt
        except ImportError:
            print("Warning: matplotlib not available, skipping FC comparison")
            return ""

        if output_path is None:
            output_path = str(self.figures_dir / 'fc_comparison.png')

        if len(scores_df) == 0:
            return ""

        fig, ax = plt.subplots(figsize=(10, 8))

        df = scores_df.copy()

        # Size based on activation score
        sizes = 100 + df['normalized_score'] * 400

        # Color based on coverage
        colors = df['metabolite_coverage'] + df['enzyme_coverage']

        scatter = ax.scatter(
            df['metabolite_fc'],
            df['enzyme_fc'],
            s=sizes,
            c=colors,
            cmap='viridis',
            alpha=0.7,
            edgecolors='black',
            linewidth=0.5
        )

        # Add labels for top pathways
        top_pathways = df.nlargest(5, 'activation_score')
        for _, row in top_pathways.iterrows():
            ax.annotate(
                row['pathway_id'],
                (row['metabolite_fc'], row['enzyme_fc']),
                xytext=(5, 5),
                textcoords='offset points',
                fontsize=9,
                fontweight='bold'
            )

        # Add reference lines
        ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)
        ax.axvline(x=0, color='gray', linestyle='--', linewidth=0.5, alpha=0.5)

        ax.set_xlabel('Metabolite Log2FC (weighted)', fontsize=12)
        ax.set_ylabel('Enzyme Log2FC (weighted)', fontsize=12)
        ax.set_title('Metabolite vs Enzyme Changes by Pathway', fontsize=14, fontweight='bold')

        # Colorbar
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label('Combined Coverage', fontsize=10)

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()

        return output_path

    def summary_report(
        self,
        results: Dict,
        output_path: Optional[str] = None
    ) -> str:
        """
        Generate summary report as JSON.

        Args:
            results: Dictionary with analysis results
            output_path: Optional output path

        Returns:
            Path to saved report
        """
        if output_path is None:
            output_path = str(self.output_dir / 'pathway_report.json')

        report = {
            'analysis_date': datetime.now().isoformat(),
            'summary': {
                'n_pathways_analyzed': len(results.get('pathway_scores', [])),
                'top_pathway': None,
                'method': 'FC correlation-based pathway activation scoring'
            },
            'pathway_rankings': [],
            'top_proteins_by_pathway': {}
        }

        # Process pathway scores
        if 'pathway_scores_df' in results:
            scores_df = results['pathway_scores_df']
            for _, row in scores_df.iterrows():
                report['pathway_rankings'].append({
                    'pathway_id': row['pathway_id'],
                    'pathway_name': row['pathway_name'],
                    'activation_score': float(row['activation_score']),
                    'normalized_score': float(row['normalized_score']),
                    'metabolite_fc': float(row['metabolite_fc']),
                    'enzyme_fc': float(row['enzyme_fc']),
                    'metabolite_coverage': float(row['metabolite_coverage']),
                    'enzyme_coverage': float(row['enzyme_coverage']),
                    'n_metabolites': int(row['n_metabolites_detected']),
                    'n_enzymes': int(row['n_enzymes_detected'])
                })

            if len(scores_df) > 0:
                top_row = scores_df.iloc[0]
                report['summary']['top_pathway'] = {
                    'id': top_row['pathway_id'],
                    'name': top_row['pathway_name'],
                    'score': float(top_row['normalized_score'])
                }

        # Process protein contributions
        if 'protein_contributions' in results:
            for pathway_id, contrib_df in results['protein_contributions'].items():
                if len(contrib_df) > 0:
                    top_proteins = contrib_df.head(10).to_dict('records')
                    # Clean up for JSON serialization
                    for p in top_proteins:
                        for k, v in p.items():
                            if isinstance(v, (np.integer, np.floating)):
                                p[k] = float(v) if np.isfinite(v) else 0.0
                    report['top_proteins_by_pathway'][pathway_id] = top_proteins

        # Save report
        with open(output_path, 'w') as f:
            json.dump(report, f, indent=2, default=str)

        return output_path

    def generate_all_figures(
        self,
        pathway_scores_df: pd.DataFrame,
        protein_contributions: Dict[str, pd.DataFrame]
    ) -> List[str]:
        """
        Generate all figures for the analysis.

        Args:
            pathway_scores_df: DataFrame with pathway scores
            protein_contributions: Dict mapping pathway_id to protein contributions

        Returns:
            List of generated figure paths
        """
        generated_figures = []

        # Pathway heatmap
        heatmap_path = self.pathway_heatmap(pathway_scores_df)
        if heatmap_path:
            generated_figures.append(heatmap_path)
            print(f"Generated: {heatmap_path}")

        # FC comparison plot
        fc_path = self.fc_comparison_plot(pathway_scores_df)
        if fc_path:
            generated_figures.append(fc_path)
            print(f"Generated: {fc_path}")

        # Protein contribution plots for top pathways
        if len(pathway_scores_df) > 0:
            top_pathways = pathway_scores_df.head(4)['pathway_id'].tolist()
            for pathway_id in top_pathways:
                if pathway_id in protein_contributions:
                    contrib_path = self.protein_contribution_plot(
                        pathway_id,
                        protein_contributions[pathway_id]
                    )
                    if contrib_path:
                        generated_figures.append(contrib_path)
                        print(f"Generated: {contrib_path}")

        return generated_figures

    def export_tables(
        self,
        pathway_scores_df: pd.DataFrame,
        protein_contributions: Dict[str, pd.DataFrame]
    ) -> Dict[str, str]:
        """
        Export analysis results as CSV tables.

        Args:
            pathway_scores_df: DataFrame with pathway scores
            protein_contributions: Dict mapping pathway_id to protein contributions

        Returns:
            Dictionary mapping table name to file path
        """
        exported_files = {}

        # Pathway scores
        scores_path = str(self.output_dir / 'pathway_scores.csv')
        export_cols = [
            'pathway_id', 'pathway_name', 'activation_score', 'normalized_score',
            'metabolite_fc', 'enzyme_fc', 'metabolite_coverage', 'enzyme_coverage',
            'n_metabolites_detected', 'n_metabolites_total',
            'n_enzymes_detected', 'n_enzymes_total'
        ]
        available_cols = [c for c in export_cols if c in pathway_scores_df.columns]
        pathway_scores_df[available_cols].to_csv(scores_path, index=False)
        exported_files['pathway_scores'] = scores_path

        # Protein contributions (combined)
        all_contributions = []
        for pathway_id, contrib_df in protein_contributions.items():
            if len(contrib_df) > 0:
                df = contrib_df.copy()
                df['pathway_id'] = pathway_id
                all_contributions.append(df)

        if all_contributions:
            combined_df = pd.concat(all_contributions, ignore_index=True)
            contrib_path = str(self.output_dir / 'protein_contributions.csv')
            combined_df.to_csv(contrib_path, index=False)
            exported_files['protein_contributions'] = contrib_path

        return exported_files
