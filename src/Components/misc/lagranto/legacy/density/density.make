F77     = ${FORTRAN}
FFLAGS  = -O
OBJS    = density.o ${LAGRANTO}/lib/times.a ${LAGRANTO}/lib/iotra.a ${LAGRANTO}/lib/libcdfio.a ${LAGRANTO}/lib/libcdfplus.a
INCS    = ${NETCDF_INC}
LIBS    = ${NETCDF_LIB}

.f.o: $*.f
	${F77} -c ${FFLAGS} ${INCS} $*.f

density:	$(OBJS)
	${F77} -o density $(OBJS) ${INCS} $(LIBS)
