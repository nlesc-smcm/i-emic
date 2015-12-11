#!/bin/bash

if [ $# -eq 0 ]
then
	measure=norm
else
	measure=$1
fi

horizontal=$((`tput cols` - 5))
vertical=$((`tput lines` * 4 / 5 ))
padding=$(( (`tput lines` - 2 * $vertical) / 2 ))

for i in `eval echo {1..$padding}`
do
	echo ""
done

less info_0.txt | grep '  parameter value:' | sed 's/.*parameter value://' > .lambda

testr=moc+
if [ $measure = $testr ]
then
	less info_0.txt | grep 'MOC+:' | sed 's/.*||x||\:// ' > .state
else
	less info_0.txt | grep '||x||\:' | sed 's/.*||x||\:// ' > .state
fi
	
paste .lambda .state > .bifdata

if [ -s ".state" ]
then
	gnuplot -e "set terminal dumb $horizontal $vertical;\
                set title 'bifurcation diagram';\
                plot '.bifdata' using 1:2 with lines linetype 1 notitle;" > bif.plot
	cat bif.plot
fi

