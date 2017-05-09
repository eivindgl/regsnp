#!/usr/bin/env bash
set -u
set -e
shopt -s nullglob

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

for url_gen in download_urls/*.py ; do
    dstd=$base_dir/$(basename $url_gen .py)
    mkdir -p "$dstd"
    python $url_gen | ./download_urls.py --dst-dir $dstd
done

for json in download_urls/*.json ; do
    dstd=$base_dir/$(basename $json .json)
  ./download_urls.py --input $json --dst-dir $dstd
done
