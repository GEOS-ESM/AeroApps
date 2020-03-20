F77	= ${FORTRAN}
FFLAGS	= -O
OBJS	= datelist.o ${LAGRANTO}/lib/iotra.a ${LAGRANTO}/lib/times.a ${LAGRANTO}/lib/libcdfio.a ${LAGRANTO}/lib/libcdfplus.a
INCS	= ${NETCDF_INC}
LIBS	= ${NETCDF_LIB} 

.f.o: $*.f
	${F77} -c ${FFLAGS} ${INCS} $*.f

datelist:	$(OBJS)
	${F77} -o datelist $(OBJS) ${INCS} $(LIBS)
