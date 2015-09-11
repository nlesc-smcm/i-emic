#!/bin/bash
# This script transforms the profile output to a format that matlab can read
# Author: Erik -> t.e.mulder@uu.nl

# Specify in and output files
if [ $# -ne 2 ]
then
	iname=profile_16x16
	oname=gatheredprofile
else
	iname=$1
	oname=$2
fi

units=89 # Magic number

delim=" " # Delimiting character

echo -n > $oname # Initialize output file

# Write the different program component names to the output
for i in {1..89}
do
	grep -m 1 "($i)" $iname | sed 's/.*)\ *//' | sed 's/\ * :.*//' >> $oname
done

# Write the number of experiments to the output
echo $delim >> $oname
grep -c "(1)" $iname >> $oname
echo $delim >> $oname

# Write the total time results to the output
for i in {1..89}
do
	grep "($i)" $iname | sed 's/.*:\ *//' | sed 's/\ .*//' >> $oname
	echo $delim >> $oname
done
