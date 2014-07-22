#!/bin/sh

# Assumption: instrument parameter value 0 of the Silvet plugin is
# always "multiple or unknown instruments" and 1 is always "piano".

piano_path="/home/cannam/Music/piano_small_dataset"
yc="/home/cannam/code/may/bin/yc"

if [ ! -d "$piano_path" ]; then
    echo "Piano dataset directory $piano_path not found, giving up"
    exit 1
fi

if ! "$yc" -v ; then
    echo "Failed to run Yeti compiler yc at $yc_path, giving up";
fi

if ! sonic-annotator -v ; then
    echo "Failed to run sonic-annotator (not in PATH?), giving up"
    exit 1
fi

rdffile="../../silvet.n3"
if [ ! -f "$rdffile" ] ; then
    echo "Failed to find plugin RDF file at $rdffile, giving up"
    exit 1
fi

case "$piano_path" in
*\ *) echo "Piano dataset path $piano_path has a space in it, this script won't handle that"; exit 1;;
esac

( cd ../.. ; make -f Makefile.linux ) || exit 1

VAMP_PATH=../..
export VAMP_PATH

outfile="/tmp/$$"
reference="/tmp/$$ref"

tmpwav="/tmp/$$.wav"

transfile="/tmp/$$transform.ttl"

trap 'rm -f "$outfile" "$tmpwav" "$instfile" "$transfile" "$outfile.lab"' 0

infiles=$(find "$piano_path" -name \*.wav | sort)

echo
echo "Input files are:"
echo $infiles | fmt -1

time for infile in $infiles; do

    echo
    echo "Evaluating for file $infile..."

    intended_instrument=1 ## assumption: 1 == piano

    # We run this twice, once using the default instrument
    # (i.e. "multiple or unknown") and once using the intended
    # instrument preset (piano).

    filename=$(basename "$infile" .wav)

    duration=30

    for instrument in $intended_instrument ; do

	for norm in no yes; do

	    echo
	    echo "For file $filename, instrument $instrument, norm $norm..."

	    if [ "$norm" = "no" ]; then
		# Don't normalise; plugin is now supposed to do it
		sox "$infile" "$tmpwav" trim 0 $duration
	    else
		# Normalise as reference
		sox "$infile" "$tmpwav" trim 0 $duration gain -n -6.020599913279624
	    fi

	    # generate the transform by interpolating the instrument parameter
	    cat transform.ttl | sed "s/INSTRUMENT_PARAMETER/$instrument/" > "$transfile"

	    sonic-annotator \
		--writer csv \
		--csv-one-file "$outfile" \
		--csv-force \
		--transform "$transfile" \
		"$tmpwav"

	    cat "$outfile" | \
		sed 's/^[^,]*,//' | \
		while IFS=, read start duration frequency level label; do
		end=`echo "$start $duration + p" | dc`
		echo -e "$start\t$end\t$frequency"
	    done > "$outfile.lab"

	    for ms in 50 100 150; do
		mark=""
		if [ "$ms" = "50" ]; then
		    if [ "$instrument" = "0" ]; then
			mark="  <-- main generic preset for $filename (norm = $norm)"; 
		    else
			mark="  <-- main piano preset for $filename (norm = $norm)";
		    fi
		fi;
		echo
		echo "Validating against ground truth at $ms ms:"
		egrep '(^[0-9]\.)|(^[012][0-9]\.)' "../piano-groundtruth/$filename.lab" > "$reference.lab"
		"$yc" ../scripts/evaluate_lab.yeti "$ms" "$reference.lab" "$outfile.lab" | sed 's,$,'"$mark"','
		cp "$reference.lab" /tmp/reference.lab
		cp "$outfile.lab" /tmp/detected.lab
	    done;
	    echo
	done
    done
done
