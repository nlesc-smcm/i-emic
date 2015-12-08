#!/bin/bash

# You might have to change this
solvehistory=50
newtonhistory=30
less dump | grep 'impl res' | tail -n $solvehistory | sed 's/.*iteration://' | sed 's/impl res://' \
	| sed 's/expl.*//' > solveresiduals
less info_0.txt | grep '||R||' | tail -n $newtonhistory | sed 's/.*||R||://' > newtonresiduals


if [ -s "solveresiduals" ]
then
	gnuplot -e "set terminal dumb 80 20;\
                set logscale y;\
                unset xtics;\
                set xrange [0:$solvehistory];\
                set title 'linear solver';\
                plot 'solveresiduals' using 2 with lines linetype 0 notitle;"
fi

if [ -s "newtonresiduals" ]
then
	gnuplot -e "set terminal dumb 80 20;\
                set logscale y;\
                unset xtics;\
                set title 'Newton iterations';\
                plot 'newtonresiduals' using 1 with lines linetype 0 notitle"		
fi
