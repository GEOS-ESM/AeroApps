F77	= ${FORTRAN}
FFLAGS	= -O
OBJS	= getvars.o ${LAGRANTO}/lib/ioinp.a ${LAGRANTO}/lib/libcdfio.a ${LAGRANTO}/lib/libcdfplus.a
INCS	= ${NETCDF_INC}
LIBS	= ${NETCDF_LIB}

.f.o: $*.f
	${F77} -c ${FFLAGS} ${INCS} $*.f

getvars:	$(OBJS)
	${F77} -o getvars $(OBJS) ${INCS} $(LIBS)
