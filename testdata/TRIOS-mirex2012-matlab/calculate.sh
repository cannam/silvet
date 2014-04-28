#!/bin/bash

# Run this from the directory that contains it

trios_path="/import/c4dm-datasets/TRIOS_dataset"

if [ ! -d "$trios_path" ]; then
    echo "TRIOS dataset directory $trios_path not found, giving up"
    exit 1
fi

matlab_path="../../mirex2012-matlab"

if [ ! -d "$matlab_path" ] || [ ! -f "$matlab_path/doMultiF0.m" ]; then
    echo "Required MATLAB code not found in $matlab_path, giving up"
    exit 1
fi

outbase="`pwd`"
echo "Will read TRIOS files from $trios_path"
echo "Will write output files below $outbase"
echo "If either of these is incorrect, hit ctrl-C now!"
sleep 8

if echo quit | matlab -nojvm ; then echo
else
    echo "Failed to start MATLAB to check that it works, giving up"
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
            out="$outdir/$wbase.lab"
            echo "Processing wav file $w, writing to lab file $out"
            time ( cd "$matlab_path" ; echo "doMultiF0('$w','$out')" | matlab -nojvm )
            echo "Done"
        done
    fi
    echo
done


