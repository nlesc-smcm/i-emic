.EXPORT_ALL_VARIABLES:

NUMBASE?=/iwi200/home/users/rudy/fortran
OSNAME ?=
PATH   ?=/iwi200/home/users/rudy/fortran/bin:/iwi200/home/users/rudy/fortran/usr/bin:.:/iwi200/home/users/rudy/scripts:/iwi200/home/users/rudy/bin:.:/iwi200/home/users/rudy/scripts:/iwi200/home/users/rudy/bin:/usr/local/bin:/usr/bin:/usr/X11R6/bin:/bin:/usr/lib/java/bin:/usr/games/bin:/usr/games:/opt/gnome/bin:/opt/kde2/bin:/opt/kde/bin:.:/opt/local/bin
CC     ?=icc
CFLAGS ?=-ansi -Wall
CXX    ?=icc
MPI_CC ?=/iwi200/home/users/rudy/fortran/usr/bin/icc
F77    ?=ifort
MPI_F77?=/usr/bin/f77
F90    ?=ifort
MPI_F90?=/iwi200/home/users/rudy/fortran/usr/bin/ifc
