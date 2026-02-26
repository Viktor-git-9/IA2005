EXEC_PATH = ../
EXEC = $(EXEC_PATH)/executable
COMPILER ?= gfortran
MACHINE ?= local

ifeq ($(MACHINE),local)
    	ifeq ($(COMPILER),gfortran)
        CC = gfortran
        FFLAGS = -ffree-form -ffree-line-length-1024 -fopenmp -fno-stack-arrays -O3 -march=native
    	endif

else ifeq ($(MACHINE),cluster)
    	ifeq ($(COMPILER),gfortran)
        CC = gfortran
        FFLAGS = -ffree-form -ffree-line-length-1024 -fopenmp -O3 -march=znver2 -mtune=znver2 \
                 -ftree-parallelize-loops=4 -funroll-loops -mfma -mavx2
    	else ifeq ($(COMPILER),nvfortran)
        CC = nvfortran
        FFLAGS = -acc -ta:tesla:cc80 -fPIC \
                 -I/cm/shared/apps/nvidia-hpc-sdk-11.1/Linux_x86_64/20.11/cuda/11.1/include/ \
                 -Minfo=accel
    	endif
endif

F90_OBJS := $(patsubst %.f90,%.o,$(wildcard *.f90))
F_OBJS   := $(patsubst %.f,%.o,$(wildcard *.f))
OBJS     := $(F90_OBJS) $(F_OBJS)

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(FFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(EXEC) $(OBJS) *.mod

.SUFFIXES:
.SUFFIXES: .f90 .o .f .o

%.o: %.f90
	$(CC) $(FFLAGS) -c $<

%.o: %.f
	$(CC) $(FFLAGS) -c $<

include makefile.depend
