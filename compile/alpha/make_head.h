# user specific directories
#--------------------------
ROOT   = /home/dcode
CODE   = $(ROOT)/work_PARA_VER_3_JULY_03
DCODE  = $(ROOT)/work_PARA_VER_3_JULY_03
ECODE  = $(DCODE)/runable

# general directories
#--------------------------
INCLUDES = $(DCODE)/include/dec
EXE      = $(ECODE)/piny_md_alpha.x
CMALLOC = 

#  alpha compiler
#--------------------------
FC = f77
CC = cc
OPT = -O4
OPT_CARE = -O4
OPT_GRP = -O4
CFLAGS =
FFLAGS = -nofor_main
LIBS =  -lm
                                                 


