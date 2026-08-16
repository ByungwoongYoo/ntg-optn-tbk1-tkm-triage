#!/usr/bin/env python3
import json
import urllib.parse
import urllib.request
from pathlib import Path

queries = {
    'broad_fulltext': 'OPEN_ACCESS:Y AND "network pharmacology" AND (herbal OR "traditional Chinese medicine" OR Kampo OR "Korean medicine")',
    'title_abs': 'OPEN_ACCESS:Y AND TITLE_ABS:"network pharmacology" AND (herbal OR "traditional Chinese medicine" OR Kampo OR "Korean medicine")',
    'title_abs_year': 'OPEN_ACCESS:Y AND FIRST_PDATE:[2015-01-01 TO 2026-08-16] AND TITLE_ABS:"network pharmacology" AND (herbal OR "traditional Chinese medicine" OR Kampo OR "Korean medicine")',
    'title_abs_formula': 'OPEN_ACCESS:Y AND FIRST_PDATE:[2015-01-01 TO 2026-08-16] AND TITLE_ABS:"network pharmacology" AND (decoction OR formula OR herb OR herbal OR phytochemical OR "traditional Chinese medicine")',
    'title_only': 'OPEN_ACCESS:Y AND TITLE:"network pharmacology" AND (decoction OR formula OR herb OR herbal OR phytochemical OR "traditional Chinese medicine")',
}

def fetch(query):
    url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?' + urllib.parse.urlencode({
        'query': query,
        'format': 'json',
        'pageSize': 10,
        'resultType': 'core',
    })
    with urllib.request.urlopen(url, timeout=60) as r:
        data = json.load(r)
    return {
        'url': url,
        'hitCount': data.get('hitCount'),
        'titles': [x.get('title') for x in data.get('resultList', {}).get('result', [])],
    }

out = {name: fetch(q) for name, q in queries.items()}
Path('artifact').mkdir(exist_ok=True)
Path('artifact/smoke.json').write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding='utf-8')
print(json.dumps(out, indent=2, ensure_ascii=False))
