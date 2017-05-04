#!/usr/bin/env bash
set -u
set -e

base_dir=external

for url_file in download_urls/*.txt ; do
    dstd=$base_dir/$(basename $url_file .txt)
    mkdir -p "$dstd"
    while read url ; do
        filename=$( echo "$url" | perl -MURI -le 'chomp($url = <>); print URI->new($url)->path')
        filename=$(basename "$filename")
        storage_path="$dstd/$filename"
        if [ -f $storage_path ] ; then
            echo "File exists. Skipping $storage_path"
        else
            echo "Downloading $storage_path"
            curl -o "$storage_path" "$url"
        fi
    done < $url_file
done
