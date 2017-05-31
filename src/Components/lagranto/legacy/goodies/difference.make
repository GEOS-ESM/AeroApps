F77	= ${FORTRAN}
FFLAGS	= -O
OBJS	= difference.o ${LAGRANTO}/lib//iotra.a  ${LAGRANTO}/lib/libcdfio.a ${LAGRANTO}/lib/libcdfplus.a
INCS	= ${NETCDF_INC}
LIBS	= ${NETCDF_LIB} 

.f.o: $*.f
	${F77} -c ${FFLAGS} ${INCS} $*.f

difference:	$(OBJS)
	${F77} -o difference $(OBJS) ${INCS} $(LIBS)
