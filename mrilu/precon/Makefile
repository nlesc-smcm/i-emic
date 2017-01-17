#                                                            2003-03-04
#
#                        Description:
#                        ============
#
# Makefile for the subroutines, which create or use a preconditioner.
# 
# Make can be used to:
# . produce or update the object files in  $(NUMBASE)/lib/libprecon.a
# . produce or update the documentation files in  $(NUMBASE)/doc/
# . update the include files in  $(NUMBASE)/include/
#
#
# Usage:
#
#   make [ OPT = [ debug | optim | prof ] ]
#        [ CC = <name C-Compiler> ] [ F77 = <name Fortran 77-Compiler ]
#        <action>
#
#   The contents of the C and Fortran 77 source files depend on the value of
#   the OPTion  variable definition, OPT.
#   The preprocessor(s) and compiler(s) depend on the value of the
#   the options CC or F77
#   The compilation options depend on the value of the  OPTion  variable
#   definition, OPT, and the compilers to be used.
#
#   For the different definitions we have:
#   OPT=  or omitted
#     The compiler optimisation option '-O' is used to produce the object
#     files.  The range checking and profiling are disabled.
#   OPT=debug
#     The preprocessor symbol DEBUG is defined.
#     The compiler optimisation option is disabled and the options for
#     range checking and profiling are enabled to produce the object
#     files.
#   OPT=optim
#     The preprocessor symbol OPTIM is defined.
#     The compiler optimisation options are enabled to produce the object
#     files.  Range checking and profiling are disabled.
#   OPT=prof
#     The profiling is enabled.  The compiler optimisation option is enabled,
#     if possible, and the range checking is disabled.
#   CC=<name C-compiler>
#     The name of C compiler is <name C-compiler>.
#   F77=<name Fortran 77-compiler>
#     The name of Fortran 77 compiler is <name Fortran 77-compiler>.
#
#   The following actions are defined:
#   all  (the default action)
#     produce the object files and the executable program(s).
#   clean
#     remove all redundant files.
#   clean_doc
#     remove all the document files (extension '.txt') from the current
#     directory.
#   clean_ftn
#     remove each '.f', with an existing pre-Fortran '.F' file, from the
#     current directory.
#   clean_obj
#     remove all the object files from the current directory.
#   doc
#     extract the documentation from the source files.
#   fortran_src
#     produce Fortran source files, with extension '.f', from the pre-Fortran
#     source files with extension '.F'.
#   install
#     update the document files in    $(NUMBASE)/doc/,
#            the include files in     $(NUMBASE)/include  and
#            the object files in      $(NUMBASE)/lib/libprecon.a
#   install_doc
#     place the document files of the changed source files in the directory
#     $(NUMBASE)/doc/
#   install_obj
#     place the object files of the changed source files in the library
#     $(NUMBASE)/lib/libprecon.a
#
#
# N.B.
# The make file uses the following environment variables:
#   NUMBASE  The base of the subdirectories in which binaries, libraries
#            etc. reside.
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
modules =	bepnum		bicgstab	bicgstabr	cg	\
		cmpprc		cmpsol 		copymt    	crbper	\
	    	csflpar		csslpar		diavec  	eblkdia	\
		filgrd    	gmres   	glbpars		hernumsch\
		ilduk     	incldup		iniglb	\
		iniprc		inisol		inivis		invbdia	\
		invblkcutmck	lumpdrop	lumpspace	matvec	\
		matvecp  	order		permbdi	   	permpr	\
		permrc    	possred		prcpars		prepro	\
		presred   	prpars		rdpblpars	rdprcpars 	rdsolpars\
		rdvispars 	readpars  	redblack	reordrb	\
		reordwrap 	scalmat   	scd2csr   	scd2fm	\
		schaak    	schurcmpl	solpars		solprc	\
		solve     	stomat    	stopldu   	vispars	\
		visscd    	applmlp		applprc         chmat
#
# Define the base names of the main programs involved:
#
progs = 
#
# Define the name of the library to be made:
#
thislib = libprecon
#
#
# Define the base names of the source files of the programs involved:
#
# Define the loader and its environment
#
LD       = $(F90)
LDLIBDIR = $(MRILU_LIB)
LDLIBES  = precon mtstor misc iosrc
#
ifeq "$(MKL)" "1"
mkl95_blas = mkl95_precision.mod mkl95_blas.mod
endif
#
ifdef LIBATLAS
FFLAGS   += -DWITH_ATLAS
LDLIBDIR += $(LIBATLAS)
LDLIBES  += cblas f77blas atlas g2c
else
dgetrf    = m_dgetrf.mod
dgetrs    = m_dgetrs.mod
endif
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
.SUFFIXES: .o .f .F .txt            # Define the used suffix list
#
###############################################################################
#
#                    Pattern rules:
#                    ==============
#
include $(SOURCEDIR)/mrilu/scripts/rules.mk
#
#
#                        Rules implementing Actions:
#                        ===========================
#
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
#
.PHONY: install_others
install_others:
#	cd ../getdoc;     $(MAKE) install
#	cd ../misc;       $(MAKE) install
	cd ../iosrclib;   $(MAKE) install
#	cd ../mtstor;     $(MAKE) install
#
#
.PHONY: uninstall
#
uninstall: default_uninstall

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
.PHONY: doc
#
doc: $(addsuffix .txt, $(modules) )
#
#
#                        Other Rules:
#                        ============
#
#
# Explicit dependencies for included files:
#
#
m_applmlp.mod:		m_build.mod	m_solve.mod	m_dperv.mod
m_applprc.mod:		m_build.mod	m_possred.mod	m_chkcnt.mod	m_applmlp.mod\
			m_glbpars.mod	m_presred.mod	m_dperv.mod
m_bepnum.mod:     	m_build.mod	m_csslpar.mod	m_csflpar.mod	m_fstrlen.mod\
			m_glbpars.mod	m_invbdia.mod   m_lumpspace.mod	m_permpr.mod\
			m_prcpars.mod	m_reordwrap.mod	m_schurcmpl.mod	m_stopldu.mod\
			m_vispars.mod	m_visscd.mod	m_wappnd.mod\
			m_wennz.mod	m_wfree.mod	m_xfnminr.mod
m_bicgstab.mod:     	m_build.mod	m_chkcnt.mod   	m_glbpars.mod  	m_matvec.mod\
			m_matvecp.mod	m_solpars.mod	m_solve.mod	m_wennz.mod	$(mkl95_blas)
m_bicgstabr.mod:    	m_build.mod	m_chkcnt.mod   	m_glbpars.mod  	m_matvec.mod\
			m_solpars.mod	m_solve.mod	m_wennz.mod	$(mkl95_blas)
m_cg.mod:           	m_chkcnt.mod   	m_glbpars.mod	m_matvecp.mod\
			m_solpars.mod	m_solve.mod	m_wennz.mod	$(mkl95_blas)
m_chmat.mod:	m_build.mod	
m_cmpprc.mod:       	m_bepnum.mod	m_build.mod	m_glbpars.mod	m_prepro.mod\
			m_prpars.mod	m_waprc.mod
m_cmpsol.mod:      	m_bicgstab.mod	m_bicgstabr.mod	m_build.mod	m_cg.mod\
			m_dperv.mod	m_gmres.mod	m_solpars.mod
m_copymt.mod:       	m_build.mod    	m_glbpars.mod	m_wascde.mod	m_wcompr.mod
m_crbper.mod:       	m_build.mod     m_invblkcutmck.mod		m_prcpars.mod\
			m_schaak.mod
m_csflpar.mod:      	m_build.mod    	m_prcpars.mod	m_scd2fm.mod\
	  		m_wapffp.mod	m_wcompr.mod	m_wfree.mod	$(dgetrf)
m_csslpar.mod:      	m_build.mod    	m_glbpars.mod  	m_ilduk.mod\
			m_incldup.mod	m_ioerrmsg.mod	m_prcpars.mod	m_scd2csr.mod\
			m_vispars.mod	m_wacsr.mod	m_wfree.mod	m_wrtbldu.mod
m_defvals.mod:    	m_glbpars.mod  	m_prcpars.mod  	m_solpars.mod   m_vispars.mod
m_eblkdia.mod:    	m_glbpars.mod
m_gmres.mod:      	m_build.mod	m_chkcnt.mod	m_glbpars.mod  	m_matvecp.mod\
			m_solpars.mod	m_solve.mod	m_wennz.mod	$(mkl95_blas)
m_ilduk.mod:      	m_glbpars.mod	m_csrresize.mod
m_incldup.mod:    	m_build.mod	m_glbpars.mod	m_csrresize.mod
m_iniglb.mod:     	m_glbpars.mod
m_iniprc.mod:     	m_prcpars.mod
m_inisol.mod:       	m_solpars.mod
m_inivis.mod:       	m_vispars.mod
m_invbdia.mod:      	m_dgeco.mod	m_dgedi.mod
m_invblkcutmck.mod: 	
m_lumpdrop.mod:     	m_prcpars.mod
m_lumpspace.mod:    	m_build.mod    	m_dgeco.mod	m_dgedi.mod	m_glbpars.mod\
		  	m_prcpars.mod	m_wadia.mod\
			m_wfree.mod	
m_matvec.mod:       	m_build.mod	m_cscvec.mod	m_csrvec.mod	m_diavec.mod
m_matvecp.mod:      	m_build.mod     m_dperv.mod	m_matvec.mod
m_order.mod:        	m_prcpars.mod	m_redblack.mod	m_wacsr.mod
m_permbdi.mod:		m_dperv.mod
m_permpr.mod:       	m_build.mod	m_mterrmsg.mod	m_permrc.mod
m_possred.mod:      	m_build.mod	m_diavec.mod 	m_csrvec.mod
m_prepro.mod:       	m_build.mod	m_copymt.mod	m_crbper.mod	m_eblkdia.mod\
		 	m_glbpars.mod	m_invbdia.mod\
			m_prcpars.mod	m_reordrb.mod	m_scalmat.mod	m_schurcmpl.mod\
			m_stomat.mod	m_wacsr.mod	m_wacsrd.mod	m_wamlp.mod\
			m_wfree.mod
m_presred.mod:      	m_build.mod	m_matvec.mod
m_prpars.mod:		m_prcpars.mod   m_solpars.mod	m_vispars.mod
m_rdpblpars.mod:    	m_glbpars.mod   m_dump.mod
m_rdprcpars.mod:    	m_glbpars.mod	m_dump.mod	m_iniprc.mod
m_rdsolpars.mod:    	m_dump.mod	m_glbpars.mod	m_inisol.mod
m_rdvispars.mod:	m_dump.mod	m_inivis.mod
m_readpars.mod:		m_getunit.mod	m_rdfilenm.mod	m_rdpblpars.mod	m_rdsolpars.mod\
			m_rdvispars.mod
m_redblack.mod:     	
m_reordrb.mod:    	m_fstrlen.mod	m_glbpars.mod  	m_hernumsch.mod	m_ioerrmsg.mod\
			m_permbdi.mod	m_vispars.mod	m_wrtmtd.mod	m_xfnminr.mod
m_reordwrap.mod:  	m_build.mod	m_fstrlen.mod	m_dperv.mod  	m_glbpars.mod\
			m_hernumsch.mod	m_ioerrmsg.mod	m_iperv.mod	m_lumpdrop.mod\
			m_lumpspace.mod	m_order.mod	m_permbdi.mod	m_prcpars.mod\
			m_vispars.mod	m_wacsr.mod	m_wacsrd.mod	m_wcompr.mod\
			m_wrtmtd.mod	m_xfnminr.mod
m_scalmat.mod:    	m_glbpars.mod  	m_dump.mod	m_ioerrmsg.mod	m_prcpars.mod\
		  	m_vispars.mod	m_wrtmt.mod
m_schaak.mod:     	m_dgeco.mod	m_redblack.mod
m_schurcmpl.mod:  	m_build.mod    	m_dump.mod	m_wascde.mod	m_wcompr.mod
m_solprc.mod:		m_dump.mod	m_build.mod	m_chkcnt.mod	m_cmpsol.mod	m_dperv.mod\
			m_glbpars.mod	m_possred.mod	m_presred.mod
m_solve.mod:      	m_build.mod	m_cscvec.mod	m_csrvec.mod	$(dgetrs)\
			m_diavec.mod	m_mterrmsg.mod	$(mkl95_blas)
m_stomat.mod:		m_build.mod	m_wacsr.mod	m_wascbm.mod
m_stopldu.mod:    	m_build.mod	m_wapldu.mod
m_visscd.mod:     	m_build.mod    	m_glbpars.mod	m_ioerrmsg.mod  m_wrtmtd.mod
#
# End of Makefile
