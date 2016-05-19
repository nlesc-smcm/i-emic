#!/bin/bash

# check arguments
if [ $# -eq 0 ]
then
    echo "  using defaults..."
	server=cartesius
	dir=run/default/
else
	server=$1
	dir=$2
fi

if [ $# -eq 4 ]
then
	linhist=$3
	newthist=$4
else
	linhist=200
	newthist=100
fi

# Random number generator
rand()
{
	tr -dc 0-9 < /dev/urandom | head -c10;
	echo
}

# Session name needs a random number...
session_name=status_$(rand)

teststr=local
local=0
if [ $server = $teststr ]
then
	local=1
fi

tailhist=1000

if [ $local -eq 1 ]
then
	echo "LOCAL!" ${SHARED_DIR}
	sleep 1
	tmux new-session -d -s $session_name 'exec tail -f -n '$tailhist' '${SHARED_DIR}'/i-emic/'$dir'/dump '
#	tmux split-window -h 'cd '${SHARED_DIR}'/i-emic/'$dir'/ && watch -n 5 -t ./plotresidual.sh '$linhist' '$newthist' '
	tmux split-window -v -t 0 'exec  tail -f -n '$tailhist'  '${SHARED_DIR}'/i-emic/'$dir'/info_0.txt '
	tmux resize-pane -L 30
	
else

	tmux new-session -d -s $session_name 'exec ssh -t '$server' \
"tail -f -n '$tailhist' \${SHARED_DIR}/i-emic/'$dir'/dump" '
	
	tmux split-window -h 'exec ssh -t '$server' "cd \${SHARED_DIR}/i-emic/'$dir'/ \
&& watch -n 5 -t ./plotresidual.sh '$linhist' '$newthist'" '
	
	tmux split-window -v -t 0 'exec ssh -t '$server' \
"tail -f -n '$tailhist' \${SHARED_DIR}/i-emic/'$dir'/info_0.txt" '
	
	tmux resize-pane -L 30

	tmux split-window -v -t 1 'exec ssh -t '$server' "cd \${SHARED_DIR}/i-emic/'$dir'/; bash" '
	
	tmux resize-pane -D 20
fi

tmux rename-window 'Status'

tmux -2 attach-session -t $session_name
