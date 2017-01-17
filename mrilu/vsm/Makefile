#                                                           2002-10-07
#
# make
#   make executable program  vsm  in the current directory
# make tar
#   make a tar file  matrix.tar  from the executables
# make install
#   make executable program  vsm  and move it to  $NUMBASE/bin/
# make clean
#   remove all obsolete files and object files.
# make mailtje
#   make a shar file in the file  mailtje
#
###############################################################################
#
#                        Definitions:
#                        ============
# include userdefined constants
include ../makefile.inc
#
# Define the base names of the modules involved:
#
modules =
#
# Define the base names of the source files involved:
#
objects =  LXlib vissmat smatrix
#
# Define the base names of the main programs involved:
#
progs = vsm
#
# Define the name of the library to be made:
#
thislib = libvsm
#
# Define the loader and its environment
#
LD       = $(CC)
LDLIBDIR = $(MRILU_LIB) $(LIB_X11)
LDLIBES  += vsm X11 m
#
#
# Definition of search path variable:
#
VPATH = $(MRILU_LIB):$(MRILU_MOD):$(MRILU_OBJ)
#
#
#
# Delete target in case of errors during construction of target
#
.DELETE_ON_ERROR:
#
#
.SUFFIXES:                          # Delete the default suffixes
.SUFFIXES: .o .c .h .F90 .txt       # Define the used suffix list
#
###############################################################################

#                    Pattern rules:
#                    ==============
#
#
include $(SOURCEDIR)/mrilu/scripts/rules.mk
#
#
#                        Rules implementing Actions:
#                        ===========================
#
# Default goal: compile all modules
#
.PHONY: all
#
all: default_all
#
# Installation
#
.PHONY: install
#
install: default_install
#
# Making the documentation
#
.PHONY: install_doc
#
install_doc: default_install_doc
#
# Adding the object files to the library
#
.PHONY: install_obj
#
install_obj: default_install_obj
#
.PHONY: install_prg
#
install_prg: default_install_prg
#
#
.PHONY: install_others
install_others:
#
#
.PHONY: uninstall
#
uninstall: default_uninstall
#
.PHONY: clean
#
clean: default_clean
#
.PHONY: clean_bak
#
clean_bak: default_clean_bak
#
.PHONY: clean_obj
#
default_clean_obj: clean_obj 

.PHONY: clean_prog
#
clean_prog: default_clean_prog

.PHONY: default_clean_doc
#
clean_doc: default_clean_doc
#

#


# Overige afhankelijkheden:
#
globals.h: inttypes.h
LXlib.o:   globals.h
smatrix.o: globals.h
vissmat.o: globals.h
vsm.o:     globals.h LXlib.o smatrix.o vissmat.o
#
#
# End of Makefile
