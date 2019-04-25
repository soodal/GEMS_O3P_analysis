##########################################################
#                                                         #
#   This is the Makefile for SOMIPROF, the SAO's          #
#                                                         #
#     BOas REtrieval for Atmospheric Spectra              #
#                                                         #
#                                                         #
#                                                         #
#   Author:  Thomas P. Kurosu, Kelly Chance               #
#                                                         #
#            Harvard-Smithsonian Center for Astrophysics  #
#            Atomic and Molecular Physics Division        #
#            60 Garden Street (MS 50)                     #
#            Cambridge, MA 02138, USA                     #
#            Tel/Fax: +1 -- 617 - 495 7213                #
#            EMail tkurosu@cfa.harvard.edu                #
#                                                         #
#                                                         #
#   Last modified: June 2003                              #
#                                                         #
###########################################################

# ================================================
# Define some environment variables that
# are unlikely to change between host-types
# ================================================
AR      = ar
ARFLAGS = ruv
MV      = mv
RM      = rm -f
SHELL   = bash

# HDF/HDF-EOS/Toolkit/Pspline libraries
PSPLINELIBS  = -lpspline -lezcdf
NETCDFLIBS   = -lnetcdf
HE4LIBS      = -lhdfeos -lmfhdf -lGctp -ldf -ljpeg -lsz -lz
HE5LIBS      = -lhdf5_fortran -lhe5_hdfeos -lhdf5 
PGSLIBS      = -lPGSTK
OMILIBS      = -lOMIutil_v2.1.0

#===============================================================


# =======================================================
# Target Directories for Modules, Library, and Executable
# =======================================================
PGETMPDIR  = $(PGEHOME)/tmp
PGEBINDIR  = $(PGEHOME)/bin
PGEOUTDIR  = $(PGEIODIR)
PGEOBJDIR  = $(PGEHOME)/obj

$(PGEHOME=$(PGEHOME))

OMSAOLIB   = SYNT_O3P
MAINDIR    = synt-pge
MAIN       = $(OMSAOLIB)_Main
$(MAIN=$(MAIN))

TARGET     = ../bin/$(SAOPGE)_exec
$(TARGET=$(TARGET))

EXECUTABLE = $(TARGET) 
$(EXECUTABLE=$(EXECUTABLE))

include ../make/make.flags

# =============================================================
#  Library for all compiled object files. Rather than keeping
#  the *.o in the current directory, they are archived in
#  LIBFILE. This way, we can maintain object files for various
#  compilation flags (opt,dbx,...) or even architectures
# =============================================================
LOCALLIB = $(PGETMPDIR)/lib$(OMSAOLIB).a
VPATH = ./$(PGEBINDIR):$(PGETMPDIR):$(PGEHOME)/src:$(PGEHOME)/src/mod:\
	$(PGEHOME)/src/omi-pge:$(PGEHOME)/src/o3prof:$(PGEHOME)/src/gems-pge:\
	$(PGEHOME)/src/synt-share:$(PGEHOME)/src/synt-pge

# ================================================
# All sources are "sourced out" to an include file
# ================================================
include ../make/make.synt_share_sources
include ../make/make.mod_sources
include ../make/make.synt_mod_sources
include ../make/make.sources
include ../make/make.o3prof_sources
include ../make/make.synt_sources
include ../make/make.gems_sources

GEMS_LIB += -I/usr/local/mpi/include -L/usr/local/mpi/lib -lmpi -lmpi_mpifh
GEMS_LIB += -L$(LAPACK_PATH_LIB) -llapack -lblas
# =================================================
#  Everything we want to make in this distribution
# =================================================
all : envcheck dirs messages mainlink $(EXECUTABLE)
	@for target in $^; do                          \
           $(MAKE) $$target;     \
        done


# ============================================
#  Here's how the111133.executable is made from the
#  archive (LOCALLIB) and the various Modules
# ============================================
$(EXECUTABLE) : $(LOCALLIB) $(MAIN).f90  $(F90MODSOURCES) $(F90SOURCES) $(F90GEMSSOURCES) $(F90SYNTSOURCES) $(OZSOURCES) \
                $(LIDORTPATH)/$(LIDORTLIBFILE) 
	$(FC) $(MAIN).f90 -o $(TARGET) $(LNFLAGS) $(FFLAGS) $(IFLAGS) $(LFLAGS) $(LIBS) $(GEMS_LIB) -I/usr/local/mpi/include -L/usr/local/mpi/lib -lmpi -lmpi_mpifh



# ======================================
#  The SAO_OMOCLO object library archive
# ======================================
$(LOCALLIB):$(MODOBJECTS) $(OBJECTS) $(GEMSOBJECTS) $(OZOBJECTS) $(SYNTOBJECTS)
$(SYNTOBJECTS=$(SYNTOBJECTS))

# =============================================================
# A list of directories that are empty at CVS checkin, and are
# thus not version tagged (files only). They must be created at
# the start of the PGE make process.
# =============================================================
emptydirs = $(PGETMPDIR) $(PGSMSG) $(PGEBINDIR) $(PGEOUTDIR) $(PGEOBJDIR)

# =====================
#  Various Make targets
# =====================
build: all

install: all
	$(MAKE) cleanlog cleanoutput
	cp -p $(EXECUTABLE) $(PGEOUTDIR)/
	mv 	*.mod $(PGETMPDIR)
$(PGEOUTDIR=$(PRGEOUTDIR))

test: all
	$(MAKE) install
	cd $(PGEOUTDIR); $(PGEOUTDIR)/$(EXECUTABLE) 2> $(SAOPGE).stderr > $(SAOPGE).stdout
	@grep -i "exit code" $(PGEOUTDIR)/LogStatus

vtest: all
	$(MAKE) install
	cd $(PGEOUTDIR); $(PGEOUTDIR)/$(EXECUTABLE)
	@grep -i "exit code" $(PGEOUTDIR)/LogStatus

run:
	cd $(PGEOUTDIR); $(PGEOUTDIR)/$(EXECUTABLE)

include ../make/make.targets
include ../make/make.synt_share_patterns
include ../make/make.mod_patterns
include ../make/make.synt_mod_patterns
include ../make/make.patterns
include ../make/make.o3prof_patterns
include ../make/make.gems_patterns
include ../make/make.synt_patterns
# =================
#  End of Makefile
# =================
