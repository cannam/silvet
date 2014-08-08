#!/bin/bash
tag=$(hg tags | grep -v ^tip | head -1 | awk '{ print $1 }')
out="/tmp/silvet-$tag.tar.bz2"
echo "Archiving from tag $tag to file $out"
if [ -z "$tag" ]; then 
    echo "ERROR: No tag found!"
    exit 1
fi
hg archive -r"$tag" -S --exclude notes --exclude papers --exclude mirex2012-matlab --exclude testdata --exclude yeti --exclude constant-q-cpp/misc "$out"
