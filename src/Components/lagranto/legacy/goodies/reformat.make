F77	= ${FORTRAN}
FFLAGS	= -O
OBJS	= reformat.o ${LAGRANTO}/lib//iotra.a  ${LAGRANTO}/lib/libcdfio.a ${LAGRANTO}/lib/libcdfplus.a
INCS	= ${NETCDF_INC}
LIBS	= ${NETCDF_LIB} 

.f.o: $*.f 
	${F77} -c ${FFLAGS} ${INCS} $*.f

reformat:	$(OBJS)
	${F77} -o reformat $(OBJS) ${INCS} $(LIBS) 
