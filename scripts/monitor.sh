#!/bin/bash

if [ $# -eq 0 ]
	echo "  usage: <solver history> <newton history> <horizontal size> <vertical size>"
    echo "  using defaults..."
	server=cartesius
	dir=rundir_test
else
	server=$1
	dir=$2
fi

remote_command='"watch ls ~/Projects/I-EMIC/'$1'/"'
compl_command="'exec ssh -t "$server" $remote_command'"

tailhist=1000

echo $compl_command

tmux new-session -d -s status 'exec ssh -t '$server' \
"tail -f -n '$tailhist' ~/Projects/I-EMIC/'$dir'/dump" '
tmux rename-window 'Status'

tmux split-window -h 'exec ssh -t '$server' \
"cd ~/Projects/I-EMIC/'$dir'/ \
&& watch -n 5 -t ./plotresidual.sh 3000 1000" '

tmux split-window -v -t 0 'exec ssh -t '$server' \
"tail -f -n '$tailhist' ~/Projects/I-EMIC/'$dir'/info_0.txt" '

tmux resize-pane -L 10

#'exec ssh -t '$server' "watch ls ~/Projects/I-EMIC/'$dir'/" '

tmux -2 attach-session -t status
