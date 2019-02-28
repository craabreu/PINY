# user specific directories
#--------------------------
ROOT   = /home/d/daniel/RESPA2_DANIEL
CODE   = $(ROOT)
DCODE  = $(ROOT)
ECODE  = $(ROOT)/runable

# general directories
#--------------------------
INCLUDES = $(DCODE)/include/pentium_nopar
EXE      = $(ECODE)/piny_md_pentium_fftw3_opt
CMALLOC = 

# HP compiler
#--------------------------
FC = gfortran -DLINUX -DFFTW3
CC = gcc -DLINUX -DFFTW3
OPT = -O2 
OPT_CARE = -O2 
OPT_GRP = -O2
CFLAGS =  -I /home/d/daniel/fftw332/include 
FFLAGS =  
LIBS =   $(LIB_PATH) $(MALLOC) /opt/intel/fc/10.1.018/lib/libsvml.so -lm /home/d/daniel/fftw332/lib64/libfftw3.a

