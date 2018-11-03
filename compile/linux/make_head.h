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
CMALLOC =

# HP compiler
#--------------------------
FC = gfortran -DLINUX -DFFTW
CC = gcc -DLINUX -DFFTW
OPT = -Ofast
OPT_CARE = $(OPT)
OPT_GRP = $(OPT)
CFLAGS =
FFLAGS =
LIBS =   $(LIB_PATH) $(MALLOC) -lfftw
