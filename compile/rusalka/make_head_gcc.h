# user specific directories
#--------------------------
ROOT   = /home/d/daniel/RESPA2_DANIEL
CODE   = $(ROOT)
DCODE  = $(ROOT)
ECODE  = $(ROOT)/runable

# general directories
#--------------------------
INCLUDES = $(DCODE)/include/pentium_nopar
EXE      = $(ECODE)/piny_md_pentium_gcc
CMALLOC = 

# HP compiler
#--------------------------
FC = gfortran -DLINUX -DFFTW3 -fno-second-underscore
CC = gcc -DLINUX -DFFTW3
OPT = -O2 
OPT_CARE = -O 
OPT_GRP = -O
CFLAGS =  -I /home/d/daniel/fftw332/include
FFLAGS = 
LIBS =   $(LIB_PATH) $(MALLOC) -lm /home/d/daniel/fftw332/lib/libfftw3.a

