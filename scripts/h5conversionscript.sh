parameterfile=.parameter

h5dump $1 | grep -n '(0):' | grep 8: | sed 's/[0-9].*)://' > $parameterfile
backup=new_$1
cp -v $1 $outfile

h5import $parameterfile -dims 1 -path Parameters/Combined\ Forcing -type TEXTFP -size 64 -o $outfile

