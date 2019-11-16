FC = ifort
FCFLAGS = -traceback -debug all -O3
LIBPATH = -Wl,  /usr/local/lib/liblapack.a /usr/local/lib/libblas.a 

overlap:	overlap.f90
	$(FC) -o overlap $(FCFLAGS) overlap.f90 $(LIBPATH)

clean:	
	rm -rf overlap overlap.dSYM 
