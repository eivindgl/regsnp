'''
Parses metadata for uniform encode DNase peaks
'''
import collections
import json

import requests
from requests.compat import urljoin


def parse_meta(raw_lines):
    dict_list = []
    for line in raw_lines:
        filename, rest = line.split('\t')
        d = dict(x.split('=', 1) for x in rest.split('; '))
        d['filename'] = filename
        dict_list.append(d)
    return dict_list

def get_name(name, seen):
    seen[name] += 1
    return '{}_{}.narrowPeak.gz'.format(name, seen[name])

BASE_URL = 'http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgDnaseUniform/'
meta_path = urljoin(BASE_URL, 'files.txt')
lines = requests.get(meta_path).content.decode('UTF-8').splitlines()
#lines = open('external_static/metadata/encode/Dnase_metadata_raw.txt').readlines()
meta = parse_meta(lines)
seen_names = collections.defaultdict(int)

names = (get_name(x['cell'], seen_names) for x in meta)
urls = (urljoin(BASE_URL, x['filename']) for x in meta)
xs = [dict(filepath=name, url=url) for name, url in zip(names, urls)]
print(json.dumps(xs, sort_keys=True, indent=4))
