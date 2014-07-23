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

rdffile="../../silvet.n3"
if [ ! -f "$rdffile" ] ; then
    echo "Failed to find plugin RDF file at $rdffile, giving up"
    exit 1
fi

case "$trios_path" in
*\ *) echo "TRIOS dataset path $trios_path has a space in it, this script won't handle that"; exit 1;;
esac

( cd ../.. ; make -f Makefile.linux ) || exit 1

VAMP_PATH=../..
export VAMP_PATH

outfile="/tmp/$$"

tmpwav="/tmp/$$norm.wav"

instfile="/tmp/$$instruments.txt"

transfile="/tmp/$$transform.ttl"

trap 'rm -f "$outfile" "$tmpwav" "$instfile" "$transfile" "$outfile.lab"' 0

# Use the single-instrument monophonic non-synthetic files for a
# (varied) subset of the TRIOS dataset. We take only those instruments
# for which we have a preset available
infiles="$trios_path/lussier/bassoon.wav $trios_path/take_five/saxophone.wav $trios_path/schubert/violin.wav $trios_path/mozart/clarinet.wav"

grep Piano "$rdffile" | sed 's/^.*( *//' | sed 's/ *).*$//' | sed 's/ "/\n/g' | sed 's/"//g' | tr '[A-Z]' '[a-z]' | tail -n +2 | cat -n > "$instfile"

instrument_for() {
    filename="$1"
    base=`basename "$filename" .wav`
    if [ "$base" = "saxophone" ]; then base="tenorsax"; fi
    instrument_no=`grep "$base" "$instfile" | awk '{ print $1; }'`
    if [ -z "$instrument_no" ] || [ -z "$base" ]; 
    then echo 0
    else echo "$instrument_no"
    fi
}

echo
echo "Input files are:"
echo $infiles | fmt -1

time for infile in $infiles; do

    echo
    echo "Evaluating for file $infile..."

    intended_instrument=`instrument_for "$infile"`
    case "$intended_instrument" in
	[0-9]*) ;;
	*) echo "Instrument extraction failed for infile $infile -- not even default multi-instrument setting returned?"; exit 1;;
    esac
    
    piece=`basename \`dirname "$infile" \``
    arrangement=`basename "$infile" .wav`

    # We run this twice, once using the default instrument
    # (i.e. "multiple or unknown") and once using the intended
    # instrument preset (the solo one).

    for instrument in $intended_instrument 0; do

	echo
	echo "For piece $piece, arrangement $arrangement, using instrument $instrument..."

	sox "$infile" "$tmpwav" gain -n -6.020599913279624

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

	for ms in 50 100; do
	    mark=""
	    if [ "$ms" = "50" ]; then mark="  <-- main $piece/$arrangement"; fi;
	    echo
	    echo "Validating against ground truth at $ms ms:"
	    "$yc" ../scripts/evaluate_lab.yeti "$ms" "../TRIOS-groundtruth/$piece/$arrangement.lab" "$outfile.lab" | sed 's,$,'"$mark"','
	    echo
	    echo "Validating against MIREX submission at $ms ms:"
	    "$yc" ../scripts/evaluate_lab.yeti "$ms" "../TRIOS-mirex2012-matlab/$piece/$arrangement.lab" "$outfile.lab"
	    echo
	    echo "Validating MIREX against ground truth at $ms ms":
	    "$yc" ../scripts/evaluate_lab.yeti "$ms" "../TRIOS-groundtruth/$piece/$arrangement.lab" "../TRIOS-mirex2012-matlab/$piece/$arrangement.lab"
	done;

	echo
    done
done
