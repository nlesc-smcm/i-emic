#!/bin/bash

iname=profile_16x16
oname=gatheredprofile

echo > $oname

for i in {1..20}
do
	grep -m 1 "($i)" $iname | sed 's/.*)\ *//' | sed 's/\ * :.*//' >> $oname
	echo >> $oname
	grep "($i)" $iname | sed 's/.*:\ *//' | sed 's/\ .*//' >> $oname
	echo >> $oname
done
