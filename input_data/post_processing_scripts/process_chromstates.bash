#!/usr/bin/env bash
set -u
set -e

ROOT_DIR=$(git rev-parse --show-toplevel)
STATES_FILE="$ROOT_DIR/input_data/external/epigenome_roadmap/all.dense.browserFiles.tgz"
OUT_DIR="$ROOT_DIR/input_data/external/epigenome_roadmap/states"
[ -f "$STATES_FILE" ] || (echo "Missing $STATES_FILE. Run download scripts first" ; exit 1)
mkdir -p "$OUT_DIR"
TDIR=$(mktemp -d -p .)
tar zxvf "$STATES_FILE" -C "$TDIR"
for x in $TDIR/*.bed.gz ; do
    name=$(basename $x)
    zcat $x | \
        tail -n +2 | \
        cut -f 1-4 | \
        gzip > "$OUT_DIR/$name"
done
rm -fr $TDIR
