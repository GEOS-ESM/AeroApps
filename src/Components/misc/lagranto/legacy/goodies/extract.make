F77	= ${FORTRAN}
FFLAGS	= -O
OBJS	= extract.o ${LAGRANTO}/lib/iotra.a ${LAGRANTO}/lib/libcdfio.a ${LAGRANTO}/lib/libcdfplus.a
INCS	= ${NETCDF_INC}
LIBS	= ${NETCDF_LIB} 

.f.o: $*.f
	${F77} -c ${FFLAGS} ${INCS} $*.f

extract:	$(OBJS)
	${F77} -o extract $(OBJS) ${INCS}  $(LIBS)
