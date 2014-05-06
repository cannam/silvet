#!/bin/sh

trios_path="/home/cannam/Music/TRIOS_dataset"

if [ ! -d "$trios_path" ]; then
    echo "TRIOS dataset directory $trios_path not found, giving up"
    exit 1
fi

if ! sonic-annotator -v ; then
    echo "Failed to run sonic-annotator (not in PATH?), giving up"
    exit 1
fi

VAMP_PATH=../..
export VAMP_PATH

time sonic-annotator \
    --writer csv \
    --csv-one-file /dev/null \
    --csv-force \
    --default vamp:silvet:silvet:notes \
    "$trios_path/take_five/mix.wav"

