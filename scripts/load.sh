#!/bin/bash

if [ $# -ne 2 ]
then
    echo "usage: ./load runid label"
    echo "        "
    exit
fi

runid=$1
label=$2

echo ""
echo "Copying model files in" $runid "with label" $label "to the model input files."
echo ""


# Load hdf5 files
for model in $runid/*_$label.h5
do
    input=${model/$runid\//}
    input=${input/$label/input}
    cp -v $model $input
done

# load fort.44
cp -v $runid/fort.44 fort.44

