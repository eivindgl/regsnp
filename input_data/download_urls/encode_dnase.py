'''
Parses metadata for uniform encode DNase peaks
'''
import collections
import json

import requests
from requests.compat import urljoin


def parse_meta(raw_lines, **kv):
    for line in raw_lines:
        filename, rest = line.split('\t')
        d = dict(x.split('=', 1) for x in rest.split('; '))
        if d['type'] != 'narrowPeak':
            continue
        d['filename'] = filename
        for k in kv:
            d[k] = kv[k]
        yield d

def get_filename(d, seen):
    name = d['cell']
    if d['treatment'] != 'None':
        name = f'{name}_{d["treatment"]}'
    if d['track_type'] == 'uniform':
        name = 'u_' + name
    else:
        name = f's_{name}_r{d["replicate"]}'
    seen[name] += 1
    return '{}_s{}.narrowPeak.gz'.format(name, seen[name])

BASE_URLS = dict(uniform='http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgDnaseUniform/',
                 single='http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeUwDnase/')
entries = []
for URL_NAME, BASE_URL in BASE_URLS.items():
    meta_path = urljoin(BASE_URL, 'files.txt')
    lines = requests.get(meta_path).content.decode('UTF-8').splitlines()
    #lines = open('external_static/metadata/encode/Dnase_metadata_raw.txt').readlines()
    seen_names = collections.defaultdict(int)
    for d in parse_meta(lines, track_type=URL_NAME):
        local_filename = get_filename(d, seen_names)
        file_url = urljoin(BASE_URL, d['filename'])
        x = dict(filepath=local_filename, url=file_url)
        entries.append(x)

print(json.dumps(entries, sort_keys=True, indent=4))
