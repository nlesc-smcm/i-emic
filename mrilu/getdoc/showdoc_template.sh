#!/bin/sh
#
# File:  showdoc  (1997-08-13)
#
# Show the documententation files for the routines mentioned in the
# parameters:
#
#  Last update:  2002-09-23
#
if [ ${#} -gt 0 ]
then
  for name in $*
  do
    echo
    echo ""Documentation of: "$name"""
    echo
    more $NUMBASE/doc/$name.txt
  done
else
  echo
  echo "Usage:"
  echo "  showdoc <subroutine/module name>..."
  echo
fi
#
# End of  showdoc
