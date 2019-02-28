# user specific directories
#--------------------------
ROOT   = ../..
CODE   = $(ROOT)
DCODE  = $(ROOT)
ECODE  = $(ROOT)/runable

# general directories
#--------------------------
INCLUDES = $(DCODE)/include/pentium_nopar
EXE      = $(ECODE)/piny_md_linux
CMALLOC  =
FFTW     = /opt/ohpc/pub/libs/gnu8/openmpi3/fftw/3.3.8
# HP compiler
#--------------------------
FC = ifort -nofor_main -DLINUX -DFFTW3
CC = icc -DLINUX -DFFTW3
OPT = -Ofast
OPT_CARE = $(OPT)
OPT_GRP = $(OPT)
CFLAGS = -I/$(FFTW)/include
FFLAGS = -I/$(FFTW)/include
LIBS =   $(LIB_PATH) $(MALLOC) -L/$(FFTW)/lib -lfftw3
