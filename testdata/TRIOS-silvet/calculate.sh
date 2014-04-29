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

if ! sox --version ; then
    echo "Failed to run sox (not in PATH?), giving up"
    exit 1
fi

for d in brahms lussier mozart schubert take_five; do
    dir="$trios_path/$d"
    outdir="$outbase/$d"
    if [ ! -d "$dir" ]; then
        echo "TRIOS subdir $dir not found, skipping it"
    else 
        mkdir -p "$outdir"
        for w in "$dir"/*.wav; do
            wbase=`basename "$w" .wav`
            outlab="$outdir/$wbase.lab"
            echo "Processing wav file $w, writing to lab file $outlab"
	    # The MATLAB method starts by normalising to a peak 0.5
	    # (approx -3dBFS amplitude or -6dB power). We can't do
	    # that in the plugin, so must do it here
	    tmpwav="$outdir/$wbase.norm.wav"
	    sox "$w" "$tmpwav" gain -n -6.020599913279624
	    VAMP_PATH=../.. sonic-annotator \
		--writer csv \
		--csv-stdout \
		--csv-force \
		--default vamp:silvet:silvet:notes \
		"$tmpwav" | \
		while IFS=, read start duration frequency level label; do
		end=`echo "$start $duration + p" | dc`
		echo -e "$start\t$end\t$frequency"
	    done > "$outlab"
	    rm "$tmpwav"
	done
    fi
done

