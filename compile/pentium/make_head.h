# user specific directories
#--------------------------
ROOT   = /home/tuck/PINY/work_DANIEL_RESPA2
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
FC = ifort -DLINUX -DFFTW3
CC = icc -DLINUX -DFFTW3
OPT = -O2 
OPT_CARE = -O2 
OPT_GRP = -O2
CFLAGS =  -I /home/tuck/PROGS/FFTW/fftw-3.3.3/include
FFLAGS =  -nofor_main
LIBS =   $(LIB_PATH) $(MALLOC)  -lm /home/tuck/PROGS/FFTW/fftw-3.3.3/lib64/libfftw3.a

