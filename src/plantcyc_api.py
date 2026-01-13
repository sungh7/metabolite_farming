"""
PlantCyc/SoyCyc API Client for metabolite pathway mapping.

This module provides functions to query PlantCyc/SoyCyc databases via BioCyc web services
to map metabolites to pathways and retrieve pathway information.

References:
- BioCyc Web Services: https://biocyc.org/web-services.shtml
- SoyCyc Database: https://plantcyc.org/databases/soycyc/
- Plant Metabolic Network: https://plantcyc.org/
"""

import requests
import time
import os
import json
import pandas as pd
from typing import List, Dict, Optional, Tuple
from xml.etree import ElementTree as ET


class PlantCycClient:
    """Client for querying PlantCyc/SoyCyc databases via BioCyc web services."""

    # Common organism IDs in Plant Metabolic Network
    ORGIDS = {
        'plantcyc': 'PLANT',      # PlantCyc reference pathways (may not work via API)
        'soycyc': 'GMAX',         # Glycine max (soybean) - testing
        'metacyc': 'META',        # MetaCyc (all organisms) - VERIFIED WORKING
        'aracyc': 'ARA',          # Arabidopsis thaliana - VERIFIED WORKING
    }

    BASE_URL = "https://websvc.biocyc.org"

    def __init__(self, email: Optional[str] = None, password: Optional[str] = None, orgid: str = 'PLANT'):
        """
        Initialize PlantCyc API client.

        Args:
            email: BioCyc account email (optional, from env var BIOCYC_EMAIL)
            password: BioCyc account password (optional, from env var BIOCYC_PASSWORD)
            orgid: Organism database ID (default: PLANT for PlantCyc)
        """
        self.email = email or os.getenv('BIOCYC_EMAIL')
        self.password = password or os.getenv('BIOCYC_PASSWORD')
        self.orgid = orgid
        self.session = requests.Session()
        self.authenticated = False

        # Rate limiting: max 1 request per second
        self.last_request_time = 0
        self.min_request_interval = 1.0

    def authenticate(self) -> bool:
        """
        Authenticate with BioCyc web services.

        Returns:
            True if authentication successful, False otherwise
        """
        if not self.email or not self.password:
            print("Warning: No BioCyc credentials provided. Some features may be limited.")
            return False

        try:
            response = self.session.post(
                f"{self.BASE_URL}/credentials/login/",
                data={'email': self.email, 'password': self.password}
            )

            if response.status_code == 200:
                self.authenticated = True
                print(f"Successfully authenticated as {self.email}")
                return True
            else:
                print(f"Authentication failed: {response.status_code}")
                return False

        except Exception as e:
            print(f"Authentication error: {e}")
            return False

    def _rate_limit(self):
        """Enforce rate limiting (max 1 request per second)."""
        current_time = time.time()
        elapsed = current_time - self.last_request_time

        if elapsed < self.min_request_interval:
            time.sleep(self.min_request_interval - elapsed)

        self.last_request_time = time.time()

    def _make_request(self, endpoint: str, params: Dict = None) -> Optional[requests.Response]:
        """
        Make a rate-limited request to BioCyc API.

        Args:
            endpoint: API endpoint path
            params: Query parameters

        Returns:
            Response object or None if failed
        """
        self._rate_limit()

        url = f"{self.BASE_URL}/{endpoint}"

        try:
            response = self.session.get(url, params=params, timeout=30)

            if response.status_code == 200:
                return response
            else:
                print(f"Request failed: {response.status_code} - {url}")
                return None

        except Exception as e:
            print(f"Request error: {e}")
            return None

    def search_compound_by_name(self, name: str, fmt: str = 'json') -> Optional[List[Dict]]:
        """
        Search for compounds by name.

        Args:
            name: Compound name to search
            fmt: Response format ('json' or 'xml')

        Returns:
            List of compound objects or None
        """
        params = {
            'object': name,
            'class': 'Compounds',
            'fmt': fmt
        }

        response = self._make_request(f"{self.orgid}/name-search", params)

        if response and fmt == 'json':
            try:
                data = response.json()
                return data.get('RESULTS', [])
            except json.JSONDecodeError:
                print(f"Failed to parse JSON response for {name}")
                return None
        elif response:
            return response.text

        return None

    def get_pathways_of_compound(self, compound_id: str, detail: str = 'low') -> Optional[List[str]]:
        """
        Get pathways containing a specific compound.

        Args:
            compound_id: BioCyc compound ID (e.g., 'DAIDZEIN')
            detail: Detail level ('low', 'full')

        Returns:
            List of pathway IDs or None
        """
        # Ensure compound ID has orgid prefix
        if ':' not in compound_id:
            compound_id = f"{self.orgid}:{compound_id}"

        params = {
            'fn': 'pathways-of-compound',
            'id': compound_id,
            'detail': detail
        }

        response = self._make_request("apixml", params)

        if response:
            return self._parse_pathway_xml(response.text)

        return None

    def get_compounds_of_pathway(self, pathway_id: str, detail: str = 'low') -> Optional[List[str]]:
        """
        Get compounds in a specific pathway.

        Args:
            pathway_id: BioCyc pathway ID
            detail: Detail level ('low', 'full')

        Returns:
            List of compound IDs or None
        """
        if ':' not in pathway_id:
            pathway_id = f"{self.orgid}:{pathway_id}"

        params = {
            'fn': 'compounds-of-pathway',
            'id': pathway_id,
            'detail': detail
        }

        response = self._make_request("apixml", params)

        if response:
            return self._parse_compound_xml(response.text)

        return None

    def get_pathway_info(self, pathway_id: str) -> Optional[Dict]:
        """
        Get detailed information about a pathway.

        Args:
            pathway_id: BioCyc pathway ID

        Returns:
            Dictionary with pathway information or None
        """
        if ':' not in pathway_id:
            pathway_id = f"{self.orgid}:{pathway_id}"

        params = {
            'id': pathway_id,
            'detail': 'full'
        }

        response = self._make_request("getxml", params)

        if response:
            return self._parse_pathway_info_xml(response.text)

        return None

    def _parse_pathway_xml(self, xml_text: str) -> List[str]:
        """Parse pathway IDs from XML response."""
        try:
            root = ET.fromstring(xml_text)
            pathways = []

            # Look for Pathway elements
            for elem in root.iter():
                if 'Pathway' in elem.tag and 'frameid' in elem.attrib:
                    pathways.append(elem.attrib['frameid'])

            return pathways
        except ET.ParseError as e:
            print(f"XML parsing error: {e}")
            return []

    def _parse_compound_xml(self, xml_text: str) -> List[str]:
        """Parse compound IDs from XML response."""
        try:
            root = ET.fromstring(xml_text)
            compounds = []

            for elem in root.iter():
                if 'Compound' in elem.tag and 'frameid' in elem.attrib:
                    compounds.append(elem.attrib['frameid'])

            return compounds
        except ET.ParseError as e:
            print(f"XML parsing error: {e}")
            return []

    def _parse_pathway_info_xml(self, xml_text: str) -> Dict:
        """Parse pathway information from XML response."""
        try:
            root = ET.fromstring(xml_text)
            info = {}

            # Extract common-name
            for elem in root.iter():
                if 'common-name' in elem.tag:
                    info['name'] = elem.text
                elif 'comment' in elem.tag:
                    info['description'] = elem.text

            return info
        except ET.ParseError as e:
            print(f"XML parsing error: {e}")
            return {}


def test_plantcyc_api():
    """Test PlantCyc API connectivity and ORGID."""
    print("=" * 60)
    print("Testing PlantCyc/SoyCyc API Connectivity")
    print("=" * 60)

    # Test different ORGIDs
    test_orgids = ['GMAX', 'META', 'ARA']  # Updated: GMAX for soybean
    test_compounds = ['DAIDZEIN', 'FORMONONETIN', 'GENISTEIN']

    for orgid in test_orgids:
        print(f"\n--- Testing ORGID: {orgid} ---")
        client = PlantCycClient(orgid=orgid)

        # Test compound search
        for compound in test_compounds[:1]:  # Test only first compound
            print(f"\nSearching for compound: {compound}")
            results = client.search_compound_by_name(compound)

            if results:
                print(f"  Found {len(results)} results")
                for r in results[:2]:  # Show first 2 results
                    print(f"    - {r.get('OBJECT-ID')}: {r.get('COMMON-NAME')}")

                # Test pathway retrieval for first result
                if len(results) > 0:
                    compound_id = results[0].get('OBJECT-ID')
                    print(f"\n  Testing pathway retrieval for {compound_id}...")
                    pathways = client.get_pathways_of_compound(compound_id)
                    if pathways:
                        print(f"    Found {len(pathways)} pathways:")
                        for p in pathways[:3]:  # Show first 3 pathways
                            print(f"      - {p}")
                    else:
                        print(f"    No pathways found or request failed")
            else:
                print(f"  No results or request failed")

            # Small delay between compounds
            time.sleep(0.5)


def map_metabolites_to_plantcyc_pathways(
    metabolites_csv: str,
    output_csv: str,
    orgid: str = 'PLANT',
    name_column: str = 'Name',
    chebi_column: str = 'ChEBI'
) -> pd.DataFrame:
    """
    Map metabolites to PlantCyc pathways.

    Args:
        metabolites_csv: Path to metabolites CSV (from MTBLS531)
        output_csv: Path to output CSV with pathway mappings
        orgid: Organism database ID
        name_column: Column name containing metabolite names
        chebi_column: Column name containing ChEBI IDs

    Returns:
        DataFrame with pathway mappings
    """
    print(f"Loading metabolites from {metabolites_csv}...")
    df = pd.read_csv(metabolites_csv)

    print(f"Initializing PlantCyc client (ORGID: {orgid})...")
    client = PlantCycClient(orgid=orgid)

    # Try to authenticate (optional)
    client.authenticate()

    # Store results
    results = []

    print(f"\nMapping {len(df)} metabolites to PlantCyc pathways...")

    for idx, row in df.iterrows():
        metabolite_name = row.get(name_column, '')
        chebi_id = row.get(chebi_column, '')

        if pd.isna(metabolite_name) or metabolite_name.strip() == '':
            continue

        print(f"\n[{idx+1}/{len(df)}] Processing: {metabolite_name}")

        # Search by name
        compounds = client.search_compound_by_name(metabolite_name)

        if compounds and len(compounds) > 0:
            # Use first match
            compound_id = compounds[0].get('OBJECT-ID')
            print(f"  Found compound: {compound_id}")

            # Get pathways for this compound
            pathways = client.get_pathways_of_compound(compound_id)

            if pathways and len(pathways) > 0:
                print(f"  Found {len(pathways)} pathways")
                for pathway_id in pathways:
                    results.append({
                        'Metabolite_Name': metabolite_name,
                        'ChEBI': chebi_id,
                        'PlantCyc_Compound_ID': compound_id,
                        'PlantCyc_Pathway_ID': pathway_id,
                        'ORGID': orgid
                    })
            else:
                print(f"  No pathways found")
                results.append({
                    'Metabolite_Name': metabolite_name,
                    'ChEBI': chebi_id,
                    'PlantCyc_Compound_ID': compound_id,
                    'PlantCyc_Pathway_ID': None,
                    'ORGID': orgid
                })
        else:
            print(f"  Compound not found")
            results.append({
                'Metabolite_Name': metabolite_name,
                'ChEBI': chebi_id,
                'PlantCyc_Compound_ID': None,
                'PlantCyc_Pathway_ID': None,
                'ORGID': orgid
            })

    # Convert to DataFrame
    results_df = pd.DataFrame(results)

    # Save results
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    results_df.to_csv(output_csv, index=False)
    print(f"\n{'=' * 60}")
    print(f"Saved {len(results_df)} mappings to {output_csv}")
    print(f"{'=' * 60}")

    # Print summary
    total_metabolites = len(df)
    mapped_compounds = results_df['PlantCyc_Compound_ID'].notna().sum()
    mapped_pathways = results_df['PlantCyc_Pathway_ID'].notna().sum()

    print(f"\nSummary:")
    print(f"  Total metabolites: {total_metabolites}")
    print(f"  Mapped to PlantCyc compounds: {mapped_compounds} ({100*mapped_compounds/total_metabolites:.1f}%)")
    print(f"  Mapped to pathways: {mapped_pathways} ({100*mapped_pathways/len(results_df):.1f}% of compound matches)")

    return results_df


if __name__ == "__main__":
    import sys

    # Check command line arguments
    if len(sys.argv) > 1 and sys.argv[1] == 'map':
        # Run pathway mapping with MetaCyc (no authentication required)
        print("\nRunning PlantCyc/MetaCyc pathway mapping...")
        print("This will take approximately 3-5 minutes due to API rate limiting.\n")

        map_metabolites_to_plantcyc_pathways(
            metabolites_csv='data/processed/mtbls531_differential.csv',
            output_csv='data/processed/plantcyc_metabolite_pathways.csv',
            orgid='META'  # Use MetaCyc - verified working
        )
    else:
        # Run API connectivity test
        test_plantcyc_api()
