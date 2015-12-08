#!/bin/bash

# You might have to change this
solvehistory=$1
newtonhistory=$2
horizontal=100
vertical=30
date=`date +%m%d%y-%k%M`

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
