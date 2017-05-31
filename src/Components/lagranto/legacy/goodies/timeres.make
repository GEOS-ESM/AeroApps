F77	= ${FORTRAN}
FFLAGS	= -O
OBJS	= timeres.o ${LAGRANTO}/lib/iotra.a ${LAGRANTO}/lib/times.a ${LAGRANTO}/lib/libcdfio.a ${LAGRANTO}/lib/libcdfplus.a
INCS	= ${NETCDF_INC}
LIBS	= ${NETCDF_LIB} 

.f.o: $*.f
	${F77} -c ${FFLAGS} ${INCS} $*.f

timeres:	$(OBJS)
	${F77} -o timeres $(OBJS) ${INCS} $(LIBS)
