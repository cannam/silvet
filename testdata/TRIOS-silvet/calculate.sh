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

for d in brahms lussier mozart schubert take_five; do
    dir="$trios_path/$d"
    outdir="$outbase/$d"
    if [ ! -d "$dir" ]; then
        echo "TRIOS subdir $dir not found, skipping it"
    else 
        mkdir -p "$outdir"
	VAMP_PATH=../.. sonic-annotator \
	    --writer csv \
	    --csv-basedir "$outdir" \
	    --csv-force \
	    --default vamp:silvet:silvet:notes \
	    --recursive \
	    "$dir"
    fi
    echo "Converting to lab files..."
    for csv in "$outdir"/*.csv; do
	cbase=`basename "$csv" .csv`
	cat "$csv" | while IFS=, read start duration frequency level label; do
	    end=`echo "$start $duration + p" | dc`
	    echo -e "$start\t$end\t$frequency"
	done > "$outdir/$cbase.lab"
    done
done

