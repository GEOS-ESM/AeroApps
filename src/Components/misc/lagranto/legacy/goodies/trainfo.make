F77	= ${FORTRAN}
FFLAGS	= -O
OBJS	= trainfo.o ${LAGRANTO}/lib/iotra.a  ${LAGRANTO}/lib/times.a ${LAGRANTO}/lib/libcdfio.a ${LAGRANTO}/lib/libcdfplus.a
INCS	= ${NETCDF_INC}
LIBS	= ${NETCDF_LIB} 

.f.o: $*.f
	${F77} -c ${FFLAGS} ${INCS} $*.f

trainfo:	$(OBJS)
	${F77} -o trainfo $(OBJS) ${INCS} $(LIBS)
