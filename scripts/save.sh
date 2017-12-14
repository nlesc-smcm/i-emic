#!/bin/bash

if [ $# -ne 2 ]
then
    echo "usage: ./save runid label"
    echo "          "
    exit
fi

runid=$1
label=$2

read -p "Please describe label '$label': " -e label_description

echo ""
echo "Storing model states with label" $label "in directory" $runid
echo ""

mkdir -vp $runid

for model in *_output.h5
do
    newfile=${model/output/$label}
    cp -v $model $runid/$newfile
done

# Copy fort.44
cp -v fort.44 $runid/fort.44

# Copy parameter files-------------------------------------------------

logdir=$runid/log_label_$label
echo ""
echo "Storing parameterlists in directory" $logdir
echo ""
mkdir -vp $logdir

for paramlist in *.xml
do
    newparamlist=${paramlist/.xml/_$label.xml}
    echo '<!-- LABEL' $label ': ' $label_description '-->' > $logdir/$newparamlist
    echo $paramlist '->' $logdir/$newparamlist
    cat $paramlist >> $logdir/$newparamlist
done

echo ""
echo "Copy outputfiles to" $logdir
echo ""

## Copy info files -----------------------------------------------------

for infofile in info_*.txt
do
    newinfofile=${infofile/.txt/_lbl$label.txt}
    echo 'LABEL' $label ': ' $label_description > $logdir/$newinfofile
    echo   $infofile '->' $logdir/$newinfofile
    cat    $infofile >> $logdir/$newinfofile
done

## Copy cdata file ----------------------------------------------------

cdatafile=cdata.txt
newcdatafile=${cdatafile/.txt/_lbl$label.txt}
echo   $cdatafile '->' $logdir/$newcdatafile
cat    $cdatafile > $logdir/$newcdatafile

echo ""
echo "Building log" $logdir
echo ""

contpar=`grep '\"Continuation parameter\"' continuation_params.xml | sed 's/.*value=\"//' | sed 's/\".*//'`
initstep=`grep '\"initial step size\"' continuation_params.xml | sed 's/.*value=\"//' | sed 's/\".*//'`
maxnumstep=`grep '\"maximum number of steps\"' continuation_params.xml | sed 's/.*value=\"//' | sed 's/\".*//'`
destvalue=`grep 'summary' -A20  info_0.txt | tail -n 21 | grep 'destination value' | sed 's/.*: //'`
parvalue=`grep 'summary' -A20  info_0.txt | tail -n 21 | grep 'parameter value' | sed 's/.*: //'`
startvalue=`grep 'summary' -A20  info_0.txt | tail -n 21 | grep 'starting value' | sed 's/.*: //'`
numsteps=`grep 'summary' -A20  info_0.txt | tail -n 21 | grep 'step:' | sed 's/.*: //'`
numresets=`grep 'summary' -A20  info_0.txt | tail -n 21 | grep 'resets:' | sed 's/.*: //'`

logfile=log
echo 'LABEL' $label ': ' $label_description > $logdir/$logfile
echo '      Continuation summary ' >> $logdir/$logfile
echo '      starting label: ' $startlabel >> $logdir/$logfile
echo '           parameter: ' $contpar >> $logdir/$logfile
echo '     parameter value: ' $parvalue >> $logdir/$logfile
echo '      starting value: ' $startvalue >> $logdir/$logfile
echo '   destination value: ' $destvalue >> $logdir/$logfile
echo '        initial step: ' $initstep >> $logdir/$logfile
echo '         total steps: ' $numsteps >> $logdir/$logfile
echo '   max allowed steps: ' $maxnumstep >> $logdir/$logfile
echo '              resets: ' $numresets >> $logdir/$logfile
echo "" >> $logdir/$logfile

all_logs=$runid/cont.log

echo "Collecting all logs in " $all_logs
echo ""                        

echo "Overview of labels; the order is based on time of label creation." > $all_logs
echo "" >> $all_logs

for logfile in `find $runid -name $logfile -printf "%T+%p\n" | sort`
do
    cat ${logfile/*$runid/$runid} >> $all_logs
done

