#!/usr/bin/env python
import json
import argparse
import requests
import os
import shutil
import sys

p = argparse.ArgumentParser()
p.add_argument('-i', '--input', type=argparse.FileType('r'), default='-')
p.add_argument('-d', '--dst-dir', default='out')
args = p.parse_args()

if not os.path.isdir(args.dst_dir):
    os.path.makedirs(args.dst_dir)
for d in json.load(args.input):
    dst_path = os.path.join(args.dst_dir, d['filepath'])
    if os.path.isfile(dst_path):
        print(f'File exists. Skipping {dst_path}')
        continue
    else:
        print(f'downloading {dst_path}')
    r = requests.get(d['url'], stream=True)
    if r.status_code == 200:
        try:
            with open(dst_path, 'wb') as f:
                shutil.copyfileobj(r.raw, f)
        except:
            print(f'Interrupted. Deleting file {dst_path}', file=sys.stderr)
            os.unlink(dst_path)
            raise
    else:
        print(f'unable to download {dst_path} from {url}')



