## Makefile for mcmcrun project

PROJECT=mcmcrun

SRC=mcmcrun.F90

#FC = gfortran
#FFLAGS = -O3 -mtune=native -frecord-marker=4 -fexternal-blas
#FFLAGS = -g -Wall -fbounds-check -frecord-marker=4 -fexternal-blas
# In Linux, we need lapack and blas
# In Mac blas and lapack are provided in Accelerate framework
#LIBS = -lmcmcrun -framework Accelerate

# cray
#FC = pgf90
#FFLAGS = -g
#LIBS = -lmcmcrun -lacml

# cray with intel
#FC = ftn
# Intel
#FFLAGS = -g -mkl -lmkl_core -lmcmcrun
# tai pgi
#FFLAGS = -g
#LIBS = -lmcmcrun -lacml

FC=ifort
FFLAGS = -g -mkl
LIBS = -lmcmcrun

# Need to know where libmcmcrun.a and the *.mod files are
FFLAGS += -I . -L .


all: $(PROJECT)

run: mcmcrun
	./mcmcrun
 
mcmcrun: driver.o acdc.o system.o vode.o vodea.o mcmcprec.o matutils.o mcmcrun.o
	ifort $(FFLAGS) $^ -o $@ $(LIBS)

mcmcprec.o: mcmcprec.F90
	ifort -O3 -c $< -o $@

matutils.o: matutils.F90
	ifort -O3 -c $< -o $@

mcmcrun.o: mcmcrun.F90 system.o driver.o
	ifort $(FFLAGS) -c $< -o $@

driver.o: driver_acdc.f90 system.o 
	ifort -O3 -c $< -o $@

acdc.o: acdc_equations292K_2.f90 
	ifort -O3 -c $< -o $@

system.o: acdc_system292K.f90
	ifort -O3 -c $< -o $@

vode.o: vode.f
	ifort -O3 -c $< -o $@

vodea.o: vodea.f
	ifort -O3 -c $< -o $@

.PHONY: clean

clean:
	rm -f *.o *~ core *seed.dat chain.dat chain.mat \
        sschain.dat sschain.mat s2chain.dat s2chain.mat $(PROJECT) \
	mcmccovf.dat mcmcparf.dat mcmcsigma2f.dat mcmcinit_final.nml \
        mcmcmean.dat
	rm -rf *.dSYM


