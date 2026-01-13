import urllib.request
import urllib.parse

def test_url(url):
    print(f"Testing: {url}")
    try:
        with urllib.request.urlopen(url) as response:
            data = response.read().decode('utf-8').strip()
            print(f"Result: {data[:100]}...") # Print first 100 chars
    except Exception as e:
        print(f"Error: {e}")

# Test 1: Simple Name
test_url("http://rest.kegg.jp/find/compound/phenylalanine")

# Test 2: ChEBI
test_url("http://rest.kegg.jp/conv/compound/chebi:17650") # Phenylalanine
test_url("http://rest.kegg.jp/conv/compound/chebi:28044") # DL-Phenylalanine

# Test 3: Complex Name
name = "(E)-3-Hexadecenoic acid"
query = urllib.parse.quote(name)
test_url(f"http://rest.kegg.jp/find/compound/{query}")
