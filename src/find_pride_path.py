import urllib.request
import sys

base_url = "ftp://ftp.pride.ebi.ac.uk/pride/data/archive"
years = [2017, 2018, 2019]
months = range(1, 13)

accession = "PXD006989"

print(f"Searching for {accession}...")

for y in years:
    for m in months:
        m_str = f"{m:02d}"
        url = f"{base_url}/{y}/{m_str}/{accession}/"
        try:
            # print(f"Checking {url}")
            with urllib.request.urlopen(url, timeout=2) as response:
                print(f"FOUND: {url}")
                # List files
                content = response.read().decode('utf-8')
                print(content)
                sys.exit(0)
        except Exception as e:
            # print(e)
            pass
            
print("Not found in 2017-2019.")
