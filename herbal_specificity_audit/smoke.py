#!/usr/bin/env python3
import json
import urllib.parse
import urllib.request
from pathlib import Path

query = 'OPEN_ACCESS:Y AND "network pharmacology" AND (herbal OR "traditional Chinese medicine" OR Kampo OR "Korean medicine")'
url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?' + urllib.parse.urlencode({
    'query': query,
    'format': 'json',
    'pageSize': 5,
    'resultType': 'core',
})
with urllib.request.urlopen(url, timeout=60) as r:
    data = json.load(r)
out = {
    'url': url,
    'hitCount': data.get('hitCount'),
    'titles': [x.get('title') for x in data.get('resultList', {}).get('result', [])],
}
Path('artifact').mkdir(exist_ok=True)
Path('artifact/smoke.json').write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding='utf-8')
print(json.dumps(out, indent=2, ensure_ascii=False))
