#!/bin/bash

# Check arguments
if [ $# -eq 0 ]
then
    echo "  plot 4 columns against first in datafile "
    echo "  usage: plot_data.sh  <file> <col1> <col2> <col3> <col4> "
    exit
fi

xlabel=$(head $1 -n 1 | awk '{print $1}')
title1=$(head $1 -n 1 | awk '{print $'$2'}')
title2=$(head $1 -n 1 | awk '{print $'$3'}')
title3=$(head $1 -n 1 | awk '{print $'$4'}')
title4=$(head $1 -n 1 | awk '{print $'$5'}')


while true;
do
    horizontal=$((`tput cols`))
    vertical=$((`tput lines`))

    gnuplot -e \
            "set terminal dumb $horizontal $vertical; \
             set autoscale;\
             set multiplot layout 2,2 columnsfirst margins 0.1,0.9,0.1,0.9 spacing 0.1;\
             set xlabel '$xlabel'; \
             set tics scale 0.1;\
             set label 1 '$title1' at graph 0.5,1.05 center;\
             plot [$6:$7] '$1' using 1:$2 notitle with lines linetype 0;\
             set xlabel '$xlabel'; \
             set label 1 '$title2' ;\
             plot [$6:$7] '$1' using 1:$3 notitle with lines linetype 0;\
             set xlabel '$xlabel'; \
             set label 1 '$title3' ;\
             plot [$6:$7] '$1' using 1:$4 notitle with lines linetype 0;\
             set xlabel '$xlabel'; \
             set label 1 '$title4' ;\
             set label 2 '${PWD}' at  screen 0.5,0.95 right;\
             plot [$6:$7] '$1' using 1:$5 notitle with lines linetype 0;"

    sleep 5;
done
