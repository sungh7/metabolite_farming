"""Debug PlantCyc API responses to diagnose JSON parsing issues."""

import requests
import json
import time

BASE_URL = "https://websvc.biocyc.org"
ORGID = "META"

def test_compound_search(compound_name):
    """Test compound search and print raw response."""
    print(f"\n{'='*60}")
    print(f"Testing compound search: {compound_name}")
    print(f"{'='*60}")

    url = f"{BASE_URL}/{ORGID}/name-search"
    params = {
        'object': compound_name,
        'class': 'Compounds',
        'fmt': 'json'
    }

    print(f"URL: {url}")
    print(f"Params: {params}")

    try:
        response = requests.get(url, params=params, timeout=30)
        print(f"\nStatus Code: {response.status_code}")
        print(f"Content-Type: {response.headers.get('Content-Type')}")
        print(f"\nRaw Response (first 500 chars):")
        print(response.text[:500])

        if response.status_code == 200:
            try:
                data = response.json()
                print(f"\nParsed JSON successfully!")
                print(f"Results: {json.dumps(data, indent=2)}")
            except json.JSONDecodeError as e:
                print(f"\nJSON Parsing Error: {e}")
                print(f"Full response text:")
                print(response.text)
        else:
            print(f"\nRequest failed with status {response.status_code}")

    except Exception as e:
        print(f"\nException: {e}")

    time.sleep(1.5)  # Rate limit

if __name__ == "__main__":
    # Test compounds that failed
    test_compounds = [
        "Daidzein",           # Should work (worked in test)
        "daidzein",           # Try lowercase
        "DAIDZEIN",           # Try uppercase
        "Formononetin",
        "Genistein",
        "L-Arginine",
        "DL-Phenylalanine",   # This one failed with JSON error
        "Phenylalanine",      # Try simpler name
        "6''-Malonylgenistin",  # Complex name - likely to fail
        "malonylgenistin",      # Simpler version
    ]

    for compound in test_compounds[:5]:  # Test first 5
        test_compound_search(compound)
