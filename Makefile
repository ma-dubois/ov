all:
	gfortran -fpic *.c *.f -c -O3
	ar crs ov.a *.o
	gfortran -o test *.o -lm -lgfortran -O3

clean:
	rm -f *.o *.a *.so* *.exe *~ test
