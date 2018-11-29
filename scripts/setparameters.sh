#!/bin/bash

if [ $# -ne 4 ]
then
    echo "usage: setparameters.sh <restart> <cpar> <dest> <initstep>"
    echo "       restart:   true/false"
    echo "       cpar:      \"Combined Forcing\"/\"Solar Forcing\"/etc"
    echo "       dest:      destination value, 1.0 for spinups"
    echo "       initstep:  initial step size, i.e., 1e-2"

else
	iname=$1
	oname=$2
fi


# start from existing state
restart=$1

# set continuation parameter
cpar=$2

for xml in coupledmodel_params.xml ocean_params.xml seaice_params.xml atmosphere_params.xml
do
    sed -i "s/Continuation parameter.*value.*/Continuation parameter\" type=\"string\" value=\"$cpar\"\/>/" $xml
    sed -i "s/Load state.*value.*/Load state\" type=\"bool\" value=\"$restart\"\/>/" $xml
done

# set destination parameter
dest=$3
sed -i "s/destination 0.*value.*/destination 0\" type=\"double\" value=\"$dest\"\/>/" continuation_params.xml

# set initial step size
initstep=$4
sed -i "s/initial step size.*value=.*/initial step size\"  type=\"double\"  value=\"$initstep\"\/>/" continuation_params.xml
