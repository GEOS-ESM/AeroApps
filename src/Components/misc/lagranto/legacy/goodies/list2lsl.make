F77	= ${FORTRAN}
FFLAGS	= -O
OBJS	= list2lsl.o ${LAGRANTO}/lib/iotra.a  ${LAGRANTO}/lib/libcdfio.a ${LAGRANTO}/lib/libcdfplus.a
INCS	= ${NETCDF_INC}
LIBS	= ${NETCDF_LIB} 

.f.o: $*.f
	${F77} -c ${FFLAGS} ${INCS} $*.f

list2lsl:	$(OBJS)
	${F77} -o list2lsl $(OBJS) ${INCS} $(LIBS)
