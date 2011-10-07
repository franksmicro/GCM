gcm:	gcm.o netcdf.o
	cc -o gcm gcm.o netcdf.o -lm -lnetcdf -lfftw3

gcm.o:	gcm.c griddefs.h
	cc -c gcm.c

netcdf.o: netcdf.c griddefs.h
	cc -c netcdf.c

