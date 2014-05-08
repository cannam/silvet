#!/bin/sh

revision="$1"

if [ -z "$revision" ]; then
    echo "Usage: run-for-revision.sh <hg-revision-id>"
    exit 2
fi

tmppath="../../../silvet-hg-$revision-$$"

echo
echo "Running timing and validation test script for revision $revision..."

echo
echo "Cloning into $tmppath..."

hg clone -r"$revision" ../.. "$tmppath" || exit 1

echo 
echo "Building..."

( cd "$tmppath" && make -f Makefile.linux ) || exit 1

echo
echo "Copying plugin to local directory..."

cp "$tmppath"/silvet.so ../.. || exit 1

rm -r "$tmppath" 

echo
echo "Running timing and evaluation script..."

./run.sh

echo
echo "Done, for revision $revision"
