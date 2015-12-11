#!/bin/bash

if [ $# -eq 0 ]
then
    echo "  using defaults..."
	server=cartesius
	dir=rundir_test
else
	server=$1
	dir=$2
fi

rand()
{
	tr -dc 0-9 < /dev/urandom | head -c10;
	echo
}

session_name=status_$(rand)

teststr=local
local=0
if [ $server = $teststr ]
then
	local=1
	echo "LOCAL!"
fi

tailhist=1000

if [ $local -eq 1 ]
then
	tmux new-session -d -s $session_name 'exec tail -f -n '$tailhist' '${SHARED_DIR}'/i-emic/'$dir'/dump '
	tmux split-window -h 'cd '${SHARED_DIR}'/i-emic/'$dir'/ && watch -n 5 -t ./plotresidual.sh 3000 1000'
	tmux split-window -v -t 0 'exec  tail -f -n '$tailhist'  '${SHARED_DIR}'/i-emic/'$dir'/info_0.txt '
	tmux resize-pane -L 30
# cd '${SHARED_DIR}'/i-emic/'$dir'/; watch -n 5 -t ls '
else

	tmux new-session -d -s $session_name 'exec ssh -t '$server' \
"tail -f -n '$tailhist' \${SHARED_DIR}/i-emic/'$dir'/dump" '

	tmux split-window -h 'exec ssh -t '$server' "cd \${SHARED_DIR}/i-emic/'$dir'/ \
&& watch -n 5 -t ./plotresidual.sh 3000 1000" '
	
	tmux split-window -v -t 0 'exec ssh -t '$server' \
"tail -f -n '$tailhist' \${SHARED_DIR}/i-emic/'$dir'/info_0.txt" '

	tmux resize-pane -L 30

	tmux split-window -v -t 1 'exec ssh -t '$server' "cd \${SHARED_DIR}/i-emic/'$dir'/; bash" '

	tmux resize-pane -D 20

fi

tmux rename-window 'Status'

tmux -2 attach-session -t $session_name
