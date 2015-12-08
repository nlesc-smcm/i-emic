#!/bin/bash

if [ $# -eq 0 ]
then
    echo "  usage: <solver history> <newton history> <horizontal size> <vertical size>"
    echo "  using defaults..."
    solvehistory=200
    newtonhistory=100
    horizontal=70
    vertical=20
else
    solvehistory=$1
    newtonhistory=$2
    horizontal=$3
    vertical=$4
fi


date=`date +%m%d%y-%k%M`

# You might have to change this
less dump | grep 'impl res' | tail -n $solvehistory | sed 's/.*iteration://' | sed 's/impl res://' \
	| sed 's/expl.*//' > solveresiduals
less info_0.txt | grep '||R||' | tail -n $newtonhistory | sed 's/.*||R||://' > newtonresiduals


if [ -s "solveresiduals" ]
then
	gnuplot -e "set terminal dumb $horizontal $vertical;\
                set logscale y;\
                unset xtics;\
                set title 'linear solver';\
                plot 'solveresiduals' using 2 with lines linetype 1 notitle;" > solveresidual.plot
	cat solveresidual.plot
fi

if [ -s "newtonresiduals" ]
then
	gnuplot -e "set terminal dumb $horizontal $vertical;\
                set logscale y;\
                unset xtics;\
                set title 'Newton iterations';\
                plot 'newtonresiduals' using 1 with lines linetype 0 notitle" > newtonresidual.plot
	cat newtonresidual.plot
fi
