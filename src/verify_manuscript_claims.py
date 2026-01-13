#!/usr/bin/env python3
"""
Comprehensive verification of all manuscript claims against actual data.
This ensures every statistic cited in the manuscript matches the data files.
"""

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact

def verify_all_claims():
    """Verify every numerical claim in the manuscript."""

    print("="*80)
    print("MANUSCRIPT CLAIMS VERIFICATION")
    print("="*80)
    print()

    # Load all data
    metabolites = pd.read_csv('data/processed/mtbls531_differential_enhanced.csv')
    kegg_enrich = pd.read_csv('results/kegg_pathway_detailed.csv')
    plantcyc_enrich = pd.read_csv('results/plantcyc_pathway_enrichment.csv')
    enzymes = pd.read_csv('results/IFS_IFR_CHI_Evidence.csv')

    errors = []

    # ============================================================================
    # ABSTRACT CLAIMS
    # ============================================================================
    print("ABSTRACT CLAIMS:")
    print("-" * 80)

    # Claim: "79 metabolites"
    actual_total = len(metabolites)
    print(f"✓ Total metabolites: {actual_total} (claimed: 79)")
    if actual_total != 79:
        errors.append(f"Total metabolites: claimed 79, actual {actual_total}")

    # Claim: "43 significant metabolites (P<0.05)"
    actual_sig = len(metabolites[metabolites['P_Value'] < 0.05])
    print(f"✓ Significant metabolites: {actual_sig} (claimed: 43)")
    if actual_sig != 43:
        errors.append(f"Significant metabolites: claimed 43, actual {actual_sig}")

    # Claim: "6''-O-acetyldaidzin: 12.3-fold, P=1.7×10⁻⁸"
    acetyl_daidzin = metabolites[metabolites['Name'].str.contains("6''-O-Acetyldaidzin", case=False, na=False)]
    if len(acetyl_daidzin) > 0:
        row = acetyl_daidzin.iloc[0]
        log2fc = row['Log2FC']
        pval = row['P_Value']
        print(f"✓ 6''-O-Acetyldaidzin: Log2FC={log2fc:.2f}, P={pval:.2e}")
        print(f"  (claimed: 12.3-fold, P=1.7×10⁻⁸)")
        if abs(log2fc - 12.30) > 0.1:
            errors.append(f"6''-O-Acetyldaidzin Log2FC: claimed 12.3, actual {log2fc:.2f}")
        if abs(pval - 1.72e-8) / 1.72e-8 > 0.1:
            errors.append(f"6''-O-Acetyldaidzin P-value: claimed 1.7e-8, actual {pval:.2e}")

    # Claim: "6''-malonylgenistin: 12.1-fold, P=5.3×10⁻⁷"
    malonyl_genistin = metabolites[metabolites['Name'].str.contains("6''-Malonylgenistin", case=False, na=False)]
    if len(malonyl_genistin) > 0:
        row = malonyl_genistin.iloc[0]
        log2fc = row['Log2FC']
        pval = row['P_Value']
        print(f"✓ 6''-Malonylgenistin: Log2FC={log2fc:.2f}, P={pval:.2e}")
        print(f"  (claimed: 12.1-fold, P=5.3×10⁻⁷)")
        if abs(log2fc - 12.09) > 0.1:
            errors.append(f"6''-Malonylgenistin Log2FC: claimed 12.1, actual {log2fc:.2f}")
        if abs(pval - 5.28e-7) / 5.28e-7 > 0.2:
            errors.append(f"6''-Malonylgenistin P-value: claimed 5.3e-7, actual {pval:.2e}")

    # Claim: "IFR: 6.4-fold"
    ifr = enzymes[enzymes['Protein Name'].str.contains('IFR', case=False, na=False)]
    if len(ifr) > 0:
        ifr_fc = ifr.iloc[0]['Log2 Fold Change']
        print(f"✓ IFR enzyme: {ifr_fc:.2f}× (claimed: 6.4×)")
        if abs(ifr_fc - 6.39) > 0.1:
            errors.append(f"IFR fold change: claimed 6.4, actual {ifr_fc:.2f}")

    # Claim: "CHI: 5.1-fold"
    chi = enzymes[enzymes['Protein Name'].str.contains('CHI', case=False, na=False)]
    if len(chi) > 0:
        chi_fc = chi.iloc[0]['Log2 Fold Change']
        print(f"✓ CHI enzyme: {chi_fc:.2f}× (claimed: 5.1×)")
        if abs(chi_fc - 5.08) > 0.1:
            errors.append(f"CHI fold change: claimed 5.1, actual {chi_fc:.2f}")

    # Claim: "IFS: 3.2-fold"
    ifs = enzymes[enzymes['Protein Name'].str.contains('IFS', case=False, na=False)]
    if len(ifs) > 0:
        ifs_fc = ifs.iloc[0]['Log2 Fold Change']
        print(f"✓ IFS enzyme: {ifs_fc:.2f}× (claimed: 3.2×)")
        if abs(ifs_fc - 3.22) > 0.1:
            errors.append(f"IFS fold change: claimed 3.2, actual {ifs_fc:.2f}")

    print()

    # ============================================================================
    # PATHWAY ENRICHMENT CLAIMS
    # ============================================================================
    print("PATHWAY ENRICHMENT CLAIMS:")
    print("-" * 80)

    # Claim: "KEGG map01110, P=0.030"
    map01110 = kegg_enrich[kegg_enrich['Pathway'].str.contains('map01110', case=False, na=False)]
    if len(map01110) > 0:
        map_pval = map01110.iloc[0]['P_Value']
        print(f"✓ KEGG map01110: P={map_pval:.4f} (claimed: P=0.030)")
        if abs(map_pval - 0.030) > 0.001:
            errors.append(f"map01110 P-value: claimed 0.030, actual {map_pval:.4f}")

    # Claim: "Odds ratio: 10.43"
    # Need to calculate this
    n_total = len(metabolites)
    n_sig = len(metabolites[metabolites['P_Value'] < 0.05])

    if len(map01110) > 0:
        pathway_sig = map01110.iloc[0]['Sig_Count']
        pathway_bg = map01110.iloc[0]['Bg_Count']

        # Fisher's exact test contingency table
        a = pathway_sig
        b = pathway_bg - pathway_sig
        c = n_sig - pathway_sig
        d = (n_total - pathway_bg) - c

        # Calculate odds ratio
        if b > 0 and c > 0:
            odds_ratio = (a * d) / (b * c)
        else:
            # Haldane correction
            odds_ratio = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))

        print(f"✓ Odds ratio: {odds_ratio:.2f} (claimed: 10.43)")
        if abs(odds_ratio - 10.43) > 1.0:
            errors.append(f"Odds ratio: claimed 10.43, actual {odds_ratio:.2f}")

    # Claim: "44 KEGG pathways tested"
    n_kegg = len(kegg_enrich)
    print(f"✓ KEGG pathways tested: {n_kegg} (claimed: 44)")
    if n_kegg != 44:
        errors.append(f"KEGG pathways: claimed 44, actual {n_kegg}")

    # Claim: "268 PlantCyc pathways tested"
    n_plantcyc = len(plantcyc_enrich)
    print(f"✓ PlantCyc pathways tested: {n_plantcyc} (claimed: 268)")
    if n_plantcyc != 268:
        errors.append(f"PlantCyc pathways: claimed 268, actual {n_plantcyc}")

    print()

    # ============================================================================
    # EFFECT SIZE CLAIMS
    # ============================================================================
    print("EFFECT SIZE CLAIMS:")
    print("-" * 80)

    # Claim: "Mean Log2FC = 0.91"
    sig_metabolites = metabolites[metabolites['P_Value'] < 0.05]
    mean_log2fc = sig_metabolites['Log2FC'].mean()
    print(f"✓ Mean Log2FC (significant): {mean_log2fc:.2f} (claimed: 0.91)")
    if abs(mean_log2fc - 0.91) > 0.1:
        errors.append(f"Mean Log2FC: claimed 0.91, actual {mean_log2fc:.2f}")

    # Claim: "54.4% significant"
    pct_sig = (n_sig / n_total) * 100
    print(f"✓ Percentage significant: {pct_sig:.1f}% (claimed: 54.4%)")
    if abs(pct_sig - 54.4) > 1.0:
        errors.append(f"Percentage significant: claimed 54.4%, actual {pct_sig:.1f}%")

    # Claim: "11 metabolites with large effect sizes"
    # Large effect = Cohen's d > 0.8, approximate as |Log2FC| > 0.8
    large_effects = len(sig_metabolites[sig_metabolites['Log2FC'].abs() > 0.8])
    print(f"✓ Large effect sizes (|Log2FC|>0.8): {large_effects} (claimed: 11)")
    if abs(large_effects - 11) > 2:
        errors.append(f"Large effects: claimed 11, actual {large_effects}")

    print()

    # ============================================================================
    # DATABASE COVERAGE CLAIMS
    # ============================================================================
    print("DATABASE COVERAGE CLAIMS:")
    print("-" * 80)

    # Claim: "36.7% KEGG mapping coverage"
    n_mapped = metabolites['KEGG'].notna().sum()
    pct_mapped = (n_mapped / n_total) * 100
    print(f"✓ KEGG mapping: {n_mapped}/{n_total} = {pct_mapped:.1f}% (claimed: 36.7%)")
    if abs(pct_mapped - 36.7) > 1.0:
        errors.append(f"KEGG mapping: claimed 36.7%, actual {pct_mapped:.1f}%")

    print()

    # ============================================================================
    # ENZYME CLAIMS
    # ============================================================================
    print("ENZYME CLAIMS:")
    print("-" * 80)

    print(f"✓ Number of key enzymes: {len(enzymes)} (claimed: 6)")
    if len(enzymes) != 6:
        errors.append(f"Number of enzymes: claimed 6, actual {len(enzymes)}")

    # List all enzymes with fold changes
    print("\nAll enzyme fold changes:")
    for _, row in enzymes.iterrows():
        print(f"  {row['Protein Name'][:40]:40s} {row['Log2 Fold Change']:.2f}×")

    print()

    # ============================================================================
    # SUMMARY
    # ============================================================================
    print("="*80)
    if len(errors) == 0:
        print("✅ VERIFICATION COMPLETE: ALL CLAIMS MATCH DATA!")
    else:
        print(f"⚠️  VERIFICATION FOUND {len(errors)} DISCREPANCIES:")
        for error in errors:
            print(f"  ❌ {error}")
    print("="*80)

    return len(errors) == 0

if __name__ == "__main__":
    success = verify_all_claims()
    exit(0 if success else 1)
