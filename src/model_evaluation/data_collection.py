#!/usr/bin/env python
"""
Data Collection Script using Scientific Skills

Collects additional data to strengthen mechanism claims:
1. KEGG pathway data for soybean isoflavonoid biosynthesis
2. Promoter sequences for TF binding motif analysis
3. TF binding motif data
4. Gene expression correlation data
5. Time-course RNA-seq dataset search

Based on GNN vs Transformer review recommendations.
"""

import os
import sys
import json
import time
from datetime import datetime

# Create output directory
OUTPUT_DIR = "results/collected_data"
os.makedirs(OUTPUT_DIR, exist_ok=True)


def collect_kegg_pathway_data():
    """
    Collect KEGG pathway data for soybean isoflavonoid biosynthesis.
    """
    print("\n" + "="*60)
    print("1. KEGG Pathway Data Collection")
    print("="*60)

    try:
        from bioservices import KEGG
        k = KEGG()

        results = {
            'organism': 'gmx',  # Glycine max
            'pathways': {},
            'genes': {},
            'compounds': {}
        }

        # Key pathways for isoflavonoid biosynthesis
        target_pathways = [
            'gmx00940',  # Phenylpropanoid biosynthesis
            'gmx00941',  # Flavonoid biosynthesis
            'gmx00943',  # Isoflavonoid biosynthesis
            'gmx04075',  # Plant hormone signal transduction (ethylene)
        ]

        for pathway_id in target_pathways:
            print(f"\n  Fetching {pathway_id}...")
            try:
                data = k.get(pathway_id)
                if data:
                    parsed = k.parse(data)
                    results['pathways'][pathway_id] = {
                        'name': parsed.get('NAME', ['Unknown'])[0] if 'NAME' in parsed else 'Unknown',
                        'genes': list(parsed.get('GENE', {}).keys()) if 'GENE' in parsed else [],
                        'compounds': list(parsed.get('COMPOUND', {}).keys()) if 'COMPOUND' in parsed else [],
                    }
                    print(f"    Found {len(results['pathways'][pathway_id]['genes'])} genes")
            except Exception as e:
                print(f"    Error: {e}")

        # Search for ethylene-related genes in soybean
        print("\n  Searching ethylene signaling genes...")
        ethylene_genes = ['EIN2', 'EIN3', 'ERF', 'ETR1', 'CTR1']
        for gene in ethylene_genes:
            try:
                search_result = k.find('gmx', gene)
                if search_result:
                    lines = search_result.strip().split('\n')[:5]
                    results['genes'][gene] = [line.split('\t')[0] for line in lines if line]
                    print(f"    {gene}: {len(results['genes'][gene])} hits")
            except Exception as e:
                print(f"    {gene} search error: {e}")

        # Key isoflavonoid compounds
        print("\n  Fetching isoflavonoid compound data...")
        isoflavonoids = ['C00814', 'C00858', 'C05623', 'C05631']  # Genistein, Daidzein, etc.
        for cpd_id in isoflavonoids:
            try:
                cpd_data = k.get(f'cpd:{cpd_id}')
                if cpd_data:
                    parsed = k.parse(cpd_data)
                    results['compounds'][cpd_id] = {
                        'name': parsed.get('NAME', ['Unknown'])[0] if 'NAME' in parsed else 'Unknown',
                        'formula': parsed.get('FORMULA', 'Unknown'),
                    }
                    print(f"    {cpd_id}: {results['compounds'][cpd_id]['name']}")
            except Exception as e:
                print(f"    {cpd_id} error: {e}")

        # Save results
        output_file = os.path.join(OUTPUT_DIR, 'kegg_pathway_data.json')
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\n  Saved to {output_file}")

        return results

    except ImportError:
        print("  bioservices not installed. Run: pip install bioservices")
        return None
    except Exception as e:
        print(f"  Error: {e}")
        return None


def collect_uniprot_tf_data():
    """
    Collect TF information from UniProt for soybean.
    """
    print("\n" + "="*60)
    print("2. UniProt TF Data Collection")
    print("="*60)

    try:
        from bioservices import UniProt
        u = UniProt(verbose=False)

        results = {
            'transcription_factors': [],
            'ethylene_responsive': [],
            'nac_family': [],
            'myb_family': []
        }

        # Search for soybean TFs
        tf_queries = [
            ('NAC domain', 'nac_family'),
            ('MYB', 'myb_family'),
            ('ethylene responsive', 'ethylene_responsive'),
        ]

        for query, key in tf_queries:
            print(f"\n  Searching: {query}...")
            try:
                search_query = f'organism_id:3847 AND ({query}) AND reviewed:false'
                result = u.search(search_query,
                                 frmt="tab",
                                 columns="accession,id,gene_names,protein_name",
                                 limit=50)
                if result:
                    lines = result.strip().split('\n')[1:]  # Skip header
                    for line in lines[:20]:
                        parts = line.split('\t')
                        if len(parts) >= 4:
                            results[key].append({
                                'accession': parts[0],
                                'entry_name': parts[1],
                                'gene_names': parts[2],
                                'protein_name': parts[3]
                            })
                    print(f"    Found {len(results[key])} entries")
            except Exception as e:
                print(f"    Error: {e}")

        # Save results
        output_file = os.path.join(OUTPUT_DIR, 'uniprot_tf_data.json')
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\n  Saved to {output_file}")

        return results

    except ImportError:
        print("  bioservices not installed")
        return None
    except Exception as e:
        print(f"  Error: {e}")
        return None


def search_geo_datasets():
    """
    Search GEO for soybean ethylene time-course datasets.
    """
    print("\n" + "="*60)
    print("3. GEO Dataset Search (Ethylene Time-course)")
    print("="*60)

    try:
        from Bio import Entrez
        Entrez.email = "researcher@example.com"  # Required

        results = {
            'ethylene_timecourse': [],
            'isoflavonoid_related': [],
            'soybean_rnaseq': []
        }

        search_terms = [
            ('"Glycine max"[Organism] AND ethylene AND time', 'ethylene_timecourse'),
            ('"Glycine max"[Organism] AND isoflavon*', 'isoflavonoid_related'),
            ('"Glycine max"[Organism] AND RNA-seq AND treatment', 'soybean_rnaseq'),
        ]

        for term, key in search_terms:
            print(f"\n  Searching: {term[:50]}...")
            try:
                handle = Entrez.esearch(db="gds", term=term, retmax=20)
                record = Entrez.read(handle)
                handle.close()

                count = int(record['Count'])
                ids = record['IdList']
                print(f"    Found {count} datasets, fetching top {len(ids)}...")

                if ids:
                    # Fetch details
                    handle = Entrez.esummary(db="gds", id=",".join(ids[:10]))
                    summaries = Entrez.read(handle)
                    handle.close()

                    for summary in summaries:
                        if isinstance(summary, dict):
                            results[key].append({
                                'id': summary.get('Id', ''),
                                'accession': summary.get('Accession', ''),
                                'title': summary.get('title', ''),
                                'summary': summary.get('summary', '')[:200],
                                'platform': summary.get('GPL', ''),
                                'samples': summary.get('n_samples', 0)
                            })

                    print(f"    Retrieved {len(results[key])} dataset details")

                time.sleep(0.5)  # Rate limiting

            except Exception as e:
                print(f"    Error: {e}")

        # Save results
        output_file = os.path.join(OUTPUT_DIR, 'geo_datasets.json')
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\n  Saved to {output_file}")

        return results

    except ImportError:
        print("  Biopython not installed")
        return None
    except Exception as e:
        print(f"  Error: {e}")
        return None


def collect_gene_info_gget():
    """
    Use gget to collect gene information and enrichment data.
    """
    print("\n" + "="*60)
    print("4. gget Gene Information Collection")
    print("="*60)

    try:
        import gget

        results = {
            'gene_info': {},
            'enrichment': None,
            'orthologs': {}
        }

        # Key isoflavonoid genes (Arabidopsis orthologs for reference)
        key_genes = [
            'AT5G13930',  # CHS (Chalcone synthase)
            'AT3G55120',  # CHI (Chalcone isomerase)
            'AT5G07990',  # F3H (Flavanone 3-hydroxylase)
        ]

        print("\n  Fetching gene information...")
        for gene_id in key_genes:
            try:
                info = gget.info([gene_id])
                if info is not None and not info.empty:
                    results['gene_info'][gene_id] = info.to_dict('records')[0]
                    print(f"    {gene_id}: OK")
            except Exception as e:
                print(f"    {gene_id}: {e}")

        # Enrichment analysis for phenylpropanoid genes
        print("\n  Running enrichment analysis...")
        gene_symbols = ['CHS', 'CHI', 'F3H', 'IFS', 'HIDM', 'I2H']
        try:
            enrichment = gget.enrichr(gene_symbols, database='pathway')
            if enrichment is not None and not enrichment.empty:
                results['enrichment'] = enrichment.head(20).to_dict('records')
                print(f"    Found {len(results['enrichment'])} enriched pathways")
        except Exception as e:
            print(f"    Enrichment error: {e}")

        # Save results
        output_file = os.path.join(OUTPUT_DIR, 'gget_gene_data.json')
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        print(f"\n  Saved to {output_file}")

        return results

    except ImportError:
        print("  gget not installed. Run: pip install gget")
        return None
    except Exception as e:
        print(f"  Error: {e}")
        return None


def search_sra_timecourse():
    """
    Search SRA for time-course RNA-seq experiments.
    """
    print("\n" + "="*60)
    print("5. SRA Time-course Dataset Search")
    print("="*60)

    try:
        from Bio import Entrez
        Entrez.email = "researcher@example.com"

        results = {
            'experiments': []
        }

        # Search for soybean ethylene RNA-seq
        search_term = '(Glycine max[Organism]) AND (RNA-Seq[Strategy]) AND (ethylene OR "hormone treatment")'

        print(f"\n  Searching SRA: {search_term[:60]}...")

        handle = Entrez.esearch(db="sra", term=search_term, retmax=30)
        record = Entrez.read(handle)
        handle.close()

        count = int(record['Count'])
        ids = record['IdList']
        print(f"    Found {count} experiments")

        if ids:
            handle = Entrez.efetch(db="sra", id=",".join(ids[:15]), rettype="full", retmode="xml")
            data = handle.read()
            handle.close()

            # Parse basic info from XML
            import re
            titles = re.findall(r'<TITLE>([^<]+)</TITLE>', data.decode() if isinstance(data, bytes) else data)
            accessions = re.findall(r'<PRIMARY_ID>([^<]+)</PRIMARY_ID>', data.decode() if isinstance(data, bytes) else data)

            for i, (acc, title) in enumerate(zip(accessions[:10], titles[:10])):
                results['experiments'].append({
                    'accession': acc,
                    'title': title[:200]
                })

            print(f"    Retrieved {len(results['experiments'])} experiment details")

        # Save results
        output_file = os.path.join(OUTPUT_DIR, 'sra_experiments.json')
        with open(output_file, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\n  Saved to {output_file}")

        return results

    except Exception as e:
        print(f"  Error: {e}")
        return None


def generate_summary_report():
    """
    Generate summary report of all collected data.
    """
    print("\n" + "="*60)
    print("GENERATING SUMMARY REPORT")
    print("="*60)

    # Load all collected data
    collected_files = [
        'kegg_pathway_data.json',
        'uniprot_tf_data.json',
        'geo_datasets.json',
        'gget_gene_data.json',
        'sra_experiments.json'
    ]

    summary = {
        'timestamp': datetime.now().isoformat(),
        'data_sources': {},
        'recommendations': []
    }

    for filename in collected_files:
        filepath = os.path.join(OUTPUT_DIR, filename)
        if os.path.exists(filepath):
            with open(filepath) as f:
                data = json.load(f)

            source_name = filename.replace('.json', '')

            if 'kegg' in filename:
                summary['data_sources']['KEGG'] = {
                    'pathways': len(data.get('pathways', {})),
                    'genes': sum(len(v) for v in data.get('genes', {}).values()),
                    'compounds': len(data.get('compounds', {}))
                }
            elif 'uniprot' in filename:
                summary['data_sources']['UniProt'] = {
                    'nac_tfs': len(data.get('nac_family', [])),
                    'myb_tfs': len(data.get('myb_family', [])),
                    'ethylene_responsive': len(data.get('ethylene_responsive', []))
                }
            elif 'geo' in filename:
                summary['data_sources']['GEO'] = {
                    'ethylene_timecourse': len(data.get('ethylene_timecourse', [])),
                    'isoflavonoid_related': len(data.get('isoflavonoid_related', [])),
                    'soybean_rnaseq': len(data.get('soybean_rnaseq', []))
                }
            elif 'sra' in filename:
                summary['data_sources']['SRA'] = {
                    'experiments': len(data.get('experiments', []))
                }

    # Add recommendations
    summary['recommendations'] = [
        "1. Use KEGG pathway genes for promoter motif analysis",
        "2. Cross-reference UniProt TFs with GNN predictions",
        "3. Download top GEO time-course datasets for temporal validation",
        "4. Use Arabidopsis orthologs for ChIP-seq data mapping",
        "5. Integrate SRA RNA-seq for pseudo-temporal ordering"
    ]

    # Save summary
    output_file = os.path.join(OUTPUT_DIR, 'collection_summary.json')
    with open(output_file, 'w') as f:
        json.dump(summary, f, indent=2)

    # Print summary
    print(f"\nData Collection Summary:")
    print(f"-" * 40)
    for source, stats in summary['data_sources'].items():
        print(f"\n  {source}:")
        for key, val in stats.items():
            print(f"    - {key}: {val}")

    print(f"\n  Results saved to: {OUTPUT_DIR}/")

    return summary


def main():
    """Run all data collection tasks."""
    print("="*60)
    print("DATA COLLECTION FOR MECHANISM STRENGTHENING")
    print("="*60)
    print(f"\nTimestamp: {datetime.now().isoformat()}")
    print(f"Output directory: {OUTPUT_DIR}")

    # Run collection tasks
    collect_kegg_pathway_data()
    collect_uniprot_tf_data()
    search_geo_datasets()
    collect_gene_info_gget()
    search_sra_timecourse()

    # Generate summary
    summary = generate_summary_report()

    print("\n" + "="*60)
    print("DATA COLLECTION COMPLETE")
    print("="*60)

    return summary


if __name__ == "__main__":
    main()
