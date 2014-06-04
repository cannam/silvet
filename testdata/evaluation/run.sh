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

rdffile="../silvet.n3"
if [ ! -f "$rdffile" ] ; then
    

case "$trios_path" in
*\ *) echo "TRIOS dataset path $trios_path has a space in it, this script won't handle that!"; exit 1;;
esac

VAMP_PATH=../..
export VAMP_PATH

outfile="/tmp/$$"

tmpwav="/tmp/$$norm.wav"

instfile="/tmp/$$instruments.txt"

trap 'rm -f "$outfile" "$tmpwav" "$instfile" "$outfile.lab"' 0

# Use the mix and single-instrument non-synthetic files for a (varied)
# subset of the TRIOS dataset
infiles=`find "$trios_path" -name \*.wav -print | egrep '(mozart|lussier|take_five)' | grep -v _syn`

grep Piano "$rdffile" | sed 's/^.*( *//' | sed 's/ *).*$//' | sed 's/ "/\n/g' | sed 's/"//g' | tr '[A-Z]' '[a-z]' | tail -n +2 | cat -n > "$instfile"

instrument_for() {
    filename="$0"
    base=`basename "$filename" .wav`
    instrument_no=`grep "$base" "$instfile" | awk '{ print $1; }'`
    if [ -z "$instrument_no" ]; 
    then echo 0
    else echo "$instrument_no"
    fi
}

for infile in $infiles; do

    echo
    echo "Evaluating for file $infile..."

    instrument=`instrument_for "$infile"`
    case "$instrument" in
	[0-9]*) ;;
	*) echo "Instrument extraction failed for infile $infile -- not even default multi-instrument setting returned?"; exit 1;;
    esac

    echo
    echo "For file $infile, using instrument setting $instrument..."

    sox "$infile" "$tmpwav" gain -n -6.020599913279624

    ##!!! todo: actually apply the instrument setting!

    time sonic-annotator \
	--writer csv \
	--csv-one-file "$outfile" \
	--csv-force \
	--default vamp:silvet:silvet:notes \
	"$tmpwav"

    cat "$outfile" | \
	sed 's/^[^,]*,//' | \
	while IFS=, read start duration frequency level label; do
	    end=`echo "$start $duration + p" | dc`
	    echo -e "$start\t$end\t$frequency"
    done > "$outfile.lab"
    
    piece=`basename \`dirname "$infile" \``
    arrangement=`basename "$infile" .wav`

    for ms in 50 100; do
	echo
	echo "Validating against ground truth at $ms ms:"
	"$yc" ./evaluate_lab.yeti "$ms" "../TRIOS-groundtruth/$piece/$arrangement.lab" "$outfile.lab"
	echo
	echo "Validating against MIREX submission at $ms ms:"
	"$yc" ./evaluate_lab.yeti "$ms" "../TRIOS-mirex2012-matlab/$piece/$arrangement.lab" "$outfile.lab"
	echo
	echo "Validating MIREX against ground truth at $ms ms":
	"$yc" ./evaluate_lab.yeti "$ms" "../TRIOS-groundtruth/$piece/$arrangement.lab" "../TRIOS-mirex2012-matlab/$piece/$arrangement.lab"
    done;

    echo
done
