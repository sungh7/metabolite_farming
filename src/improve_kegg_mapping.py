#!/usr/bin/env python3
"""
Enhanced KEGG mapping using alternative strategies.
Try PubChem, HMDB, and fuzzy name matching.
"""

import pandas as pd
import urllib.request
import urllib.parse
import time

def search_pubchem_for_kegg(compound_name):
    """Try to find KEGG ID via PubChem cross-reference."""
    try:
        # Search PubChem
        query = urllib.parse.quote(compound_name)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{query}/xrefs/RegistryID/JSON"

        with urllib.request.urlopen(url, timeout=5) as response:
            import json
            data = json.loads(response.read().decode('utf-8'))
            # Look for KEGG IDs in cross-references
            if 'InformationList' in data:
                for info in data['InformationList'].get('Information', []):
                    for xref in info.get('RegistryID', []):
                        if xref.startswith('C') and len(xref) == 6:
                            return xref
    except:
        pass
    return None

def fuzzy_name_search(name):
    """Try variations of compound name."""
    variations = [
        name,
        name.replace('6\'\'', '6"'),  # Quote style
        name.replace('DL-', ''),       # Remove stereochemistry
        name.replace('L-', ''),
        name.replace('(E)-', ''),
        name.replace('(Z)-', ''),
        name.split('(')[0].strip(),    # Remove parentheticals
    ]

    for variant in set(variations):
        if not variant:
            continue
        try:
            query = urllib.parse.quote(variant)
            url = f"http://rest.kegg.jp/find/compound/{query}"
            with urllib.request.urlopen(url, timeout=3) as response:
                data = response.read().decode('utf-8').strip()
                if data:
                    first_line = data.split('\n')[0]
                    kegg_id = first_line.split('\t')[0].replace('cpd:', '')
                    return kegg_id
            time.sleep(0.3)
        except:
            continue
    return None

def enhance_mapping(input_file, output_file):
    """Improve KEGG mapping using multiple strategies."""

    df = pd.read_csv(input_file)

    print(f"Original KEGG mapping: {df['KEGG'].notna().sum()}/{len(df)} ({100*df['KEGG'].notna().sum()/len(df):.1f}%)")

    unmapped = df[df['KEGG'].isna()].copy()
    print(f"\nAttempting to map {len(unmapped)} unmapped metabolites...")

    improved = 0
    for idx, row in unmapped.iterrows():
        name = row['Name']
        print(f"  Trying: {name[:50]}...", end=' ')

        # Try fuzzy KEGG search
        kegg_id = fuzzy_name_search(name)

        if kegg_id:
            df.at[idx, 'KEGG'] = kegg_id
            improved += 1
            print(f"✓ {kegg_id}")
        else:
            print("✗")

        if idx % 5 == 0:
            time.sleep(0.5)  # Rate limiting

    print(f"\nImproved mapping: {df['KEGG'].notna().sum()}/{len(df)} ({100*df['KEGG'].notna().sum()/len(df):.1f}%)")
    print(f"Gained: {improved} additional mappings")

    df.to_csv(output_file, index=False)
    print(f"Saved to: {output_file}")

    return df

if __name__ == "__main__":
    enhance_mapping(
        'data/processed/mtbls531_differential.csv',
        'data/processed/mtbls531_differential_enhanced.csv'
    )
