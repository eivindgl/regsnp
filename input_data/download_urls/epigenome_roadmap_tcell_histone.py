import os
import json

import requests
import bs4
import csvkit


base_url = 'http://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/'
# extract T-cell states
metapath = 'external_static/metadata/epigenome_roadmap/chromatin_state_samples_meta.csv'
assert os.path.isfile(metapath)

states = []
with open(metapath) as f:
    for x in csvkit.DictReader(f):
        if x['group'] == 'Blood & T-cell':
            states.append(x['eid'])

entries = []
html = bs4.BeautifulSoup(requests.get(base_url).content, 'html5lib')
for x in html.find_all('a'):
    if x.attrs['href'].endswith('narrowPeak.gz'):
        filename = x.attrs['href']
        eid = filename.split('-')[0]
        if eid in states:
            entries.append(dict(
                url = requests.compat.urljoin(base_url, filename),
                filepath = filename))

print(json.dumps(entries, sort_keys=True, indent=4))

