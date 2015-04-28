#!/bin/sh

trios_path="/home/cannam/Music/TRIOS_dataset"
yc="/home/cannam/code/may/bin/yc"

if [ ! -d "$trios_path" ]; then
    echo "TRIOS dataset directory $trios_path not found, giving up"
    exit 1
fi

if ! "$yc" -v ; then
    echo "Failed to run Yeti compiler yc at $yc_path, giving up";
fi

if ! sonic-annotator -v ; then
    echo "Failed to run sonic-annotator (not in PATH?), giving up"
    exit 1
fi

VAMP_PATH=../..
export VAMP_PATH

tfile=transform.ttl

outfile="/tmp/$$"

tmpwav="/tmp/$$norm.wav"

trap 'rm -f "$outfile" "$tmpwav" "$outfile.lab"' 0

sox "$trios_path/take_five/mix.wav" "$tmpwav" gain -n -6.020599913279624

time sonic-annotator \
    --writer csv \
    --csv-one-file "$outfile" \
    --csv-force \
    --transform "$tfile" \
    "$tmpwav"

cat "$outfile" | \
    sed 's/^[^,]*,//' | \
    while IFS=, read start duration frequency level label; do
    end=`echo "$start $duration + p" | dc`
    echo -e "$start\t$end\t$frequency"
    done > "$outfile.lab"

for ms in 50 100 150; do
    echo
    echo "Validating against ground truth at $ms ms:"
    "$yc" ../scripts/evaluate_lab.yeti "$ms" "../TRIOS-groundtruth/take_five/mix.lab" "$outfile.lab"
    echo
    echo "Validating against MIREX submission at $ms ms:"
    "$yc" ../scripts/evaluate_lab.yeti "$ms" "../TRIOS-mirex2012-matlab/take_five/mix.lab" "$outfile.lab"
    echo
    echo "Validating MIREX against ground truth at $ms ms":
    "$yc" ../scripts/evaluate_lab.yeti "$ms" "../TRIOS-groundtruth/take_five/mix.lab" "../TRIOS-mirex2012-matlab/take_five/mix.lab"
done;

echo

