# Makefile
SHELL := /bin/bash

default:
	cd ./Numerics_with_Fortran; make -j 16 COMPILER=gfortran
gpu: 
	cd ./Numerics_with_Fortran; make -j 16 COMPILER=nvfortran LAPACK=mkl MACHINE=cluster
clean:
	cd ./Numerics_with_Fortran; make clean
