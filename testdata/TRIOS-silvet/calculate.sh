#!/bin/bash

# Run this from the directory that contains it

trios_path="/home/cannam/Music/TRIOS_dataset"

if [ ! -d "$trios_path" ]; then
    echo "TRIOS dataset directory $trios_path not found, giving up"
    exit 1
fi

outbase="`pwd`"
echo "Will read TRIOS files from $trios_path"
echo "Will write output files below $outbase"
echo "If either of these is incorrect, hit ctrl-C now!"
sleep 8

if ! sonic-annotator -v ; then
    echo "Failed to run sonic-annotator (not in PATH?), giving up"
    exit 1
fi

VAMP_PATH=../.. sonic-annotator \
    -w csv \
    --csv-basedir "$outbase" \
    -d vamp:silvet:silvet:notes \
    -r "$trios_path"

