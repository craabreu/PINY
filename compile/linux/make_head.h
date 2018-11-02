# user specific directories
#--------------------------
ROOT   = ../..
CODE   = $(ROOT)
DCODE  = $(ROOT)
ECODE  = $(ROOT)/runable

# general directories
#--------------------------
INCLUDES = $(DCODE)/include/pentium_nopar
#EXE      = $(ECODE)/piny_md_pentium_fftw3_opt
EXE      = $(ECODE)/piny_md_pentium_test
CMALLOC = 

# HP compiler
#--------------------------
FC = gfortran -DLINUX -DFFTW3
CC = gcc -DLINUX -DFFTW3
OPT = -Ofast
OPT_CARE = -Ofast
OPT_GRP = -Ofast
CFLAGS =  
FFLAGS = 
LIBS =   $(LIB_PATH) $(MALLOC) -lfftw3

