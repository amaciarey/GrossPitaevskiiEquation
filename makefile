#f90=gfortran -Wall
#f90=gfortran -O2 -Wall -fopenmp -fcheck=bounds #-heap-arrays #-qopenmp
f90=ifort -qopenmp -heap-arrays -mmic
FFT_LIB=/home/amacia/local/FFTW/3.3.5/lib
FFT_INC=/home/amacia/local/FFTW/3.3.5/include

gpe3d: gpe3d.o simpson.o tridiagonal.o
	$(f90) -o gpe3d *.o 

gpe3d.o: gpe3d.f90
	$(f90) -c gpe3d.f90 

simpson.o: simpson.f90
	$(f90) -c simpson.f90

tridiagonal.o: tridiagonal.f90
	$(f90) -c tridiagonal.f90

clean:
	rm -f *.o 
