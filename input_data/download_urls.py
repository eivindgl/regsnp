#!/usr/bin/env python
import json
import argparse
import requests
import os
import shutil
import sys
import tempfile
import tarfile

def download_file(url, dst):
    r = requests.get(url, stream=True)
    if r.status_code == 200:
        try:
            with open(dst, 'wb') as f:
                shutil.copyfileobj(r.raw, f)
        except:
            print(f'Interrupted. Deleting file {dst}', file=sys.stderr)
            os.unlink(dst)
            raise
    else:
        raise Exception(f'unable to download {dst} from {url}')

def decompress(url, dst_folder):
    t = tempfile.NamedTemporaryFile()
    download_file(url, t.name)

    tb = tarfile.open(t.name)
    member_files = (x for x in tb.getmembers() if x.isfile())
    for mem in member_files:
        filename = os.path.basename(mem.name)
        new_dst = os.path.join(dst_folder, filename)
        print('unzipping', new_dst)
        mem.name = new_dst
        tb.extract(mem)

p = argparse.ArgumentParser()
p.add_argument('-i', '--input', type=argparse.FileType('r'), default='-')
p.add_argument('-d', '--dst-dir', default='out')
args = p.parse_args()

if not os.path.isdir(args.dst_dir):
    os.makedirs(args.dst_dir)
for d in json.load(args.input):
    if d.get('compressed_folder', False):
        decompress(d['url'], args.dst_dir)
    else:
        dst_path = os.path.join(args.dst_dir, d['filepath'])
        if os.path.isfile(dst_path):
            print(f'File exists. Skipping {dst_path}')
        else:
            print(f'downloading {dst_path}')
            try:
                download_file(d['url'], dst_path)
            except Exception as e:
                print(e, file=sys.stderr)

