#
# Makefile for PM/Tree code pyul.exe - April 23, 2001 - JJD
#


RANLIB = ranlib



#####################################
# System specific definitions
#####################################

#########################################################
####      KIAS PG Compiler                          #####
#########################################################
AR = ar rcv
FC = mpif90
CC = mpicc
F90C = mpif90
#OPT = -DPGCC -mcmodel=medium -tp nehalem-64 -fast -mp -fastsse
OPT = -DINTEL -g -fcommon

INCLUDES = 

#LIBS = -L$(FFTWDIR)/lib/ -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw   $(CAMBLIBS) 
#####
#####
F90INCLUDES =
COMFLAGS = -DINDEX -DVarPM   -DXYZDBL 

FDFLAGS =  -DINCLUDE_TREE_FORCE  

CDFLAGS = -DWGROUPSIZE=8 -DNMEG=17000L -DINCLUDE_TREE_FORCE \
        -D_LARGE_FILES -DSAVESLICE  -DPMSEEDFORCE -DQUADHILBERT -DNENER=0 -DNPRE=8 -DREAD_SINK -DNCHEM=9 -DNDUST=4

RAMLIBS = -L../ -lmyram

LIBS = $(RAMLIBS) -lm

FFLAGS = $(FDFLAGS) $(OPT) $(COMFLAGS)  -std=legacy -fallow-argument-mismatch
CFLAGS = $(OPT) $(CDFLAGS)  $(COMFLAGS)  -fpermissive
LDFLAGS = $(OPT) $(RAMLIBS) 



#################################
# Compaq Alpha
#################################
#FC = f77
#CC = cc 
#INCLUDES = -I/home/kjhan/dolphin/fftw/include
#FDFLAGS = -DBIT64 -DINCLUDE_TREE_FORCE -DTREE -DCOMPACT
#DFLAGS = -DBIT64 -DNMEG=256 -DINCLUDE_TREE_FORCE -DTREE -DCOMPACT -DTREEFIX
#FFLAGS = $(FDFLAGS)  -fast -nofor_main
#CFLAGS = $(DFLAGS)  -fast
#LDFLAGS =  -fast -nofor_main
#LIBS = -L/home/kjhan/dolphin/fftw/lib -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw -lm -lz -lmpi -lm
##################################

#--- C Compiler information
#  Leave the rest untouched


#--- Suffix-based compilation rules
.SUFFIXES: .exe .o .c .f .F .f90

#rules to build binary from source


.c.o :
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

.f90.o :
	$(F90C) $(FFLAGS) $(INCLUDES) $(F90INCLUDES) -c $<

.f.o :
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

.for.o :
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

.F.o :
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

.cu.o :
	$(NVCC) $(CUFLAGS) $(INCLUDES) -c $<

.c.exe :
	$(CC) $(LDFLAGS) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS)

.f.exe :
	$(FC) $(LDFLAGS) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS) $(INCLUDES)

.f90.exe :
	$(F90C) $(LDFLAGS) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS) $(INCLUDES)

.for.exe :
	$(FC) $(LDFLAGS) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS) $(INCLUDES)

.F.exe :
	$(FC) $(LDFLAGS) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS) $(INCLUDES)

.cu.exe :
	$(FC) $(LDFLAGS) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS)

#--- Targets
