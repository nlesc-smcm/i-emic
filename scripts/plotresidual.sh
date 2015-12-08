#!/bin/bash

# You might have to change this
less dump | grep 'impl res' | tail -n 200 | sed 's/.*iteration://' | sed 's/impl res://' \
    | sed 's/expl.*//' > residuals

gnuplot -e "set terminal dumb; plot 'residuals' using 1:2"
