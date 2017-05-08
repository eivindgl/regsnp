import requests
from requests.compat import urljoin
import collections
import json

def parse_meta(lines):
    dict_list = []
    for line in lines:
        filename, rest = line.split('\t')
        d = dict(x.split('=', 1) for x in rest.split('; '))
        d['filename'] = filename
        dict_list.append(d)
    return dict_list

def get_name(name, seen):
    seen[name] += 1
    return '{}_{}.narrowPeak.gz'.format(name, seen[name])

base_url='http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeAwgDnaseUniform/'
meta_path = urljoin(base_url, 'files.txt')
lines = requests.get(meta_path).content.decode('UTF-8').splitlines()
#lines = open('external_static/metadata/encode/Dnase_metadata_raw.txt').readlines()
meta = parse_meta(lines)
seen = collections.defaultdict(int)

xs = [dict(filepath = get_name(d['cell'], seen), url = urljoin(base_url, d['filename'])) for d in meta]
print(json.dumps(xs, sort_keys=True, indent=4))
