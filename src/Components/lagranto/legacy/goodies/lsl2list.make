F77	= ${FORTRAN}
FFLAGS	= -O
OBJS	= lsl2list.o ${LAGRANTO}/lib/iotra.a  ${LAGRANTO}/lib/libcdfio.a ${LAGRANTO}/lib/libcdfplus.a
INCS	= ${NETCDF_INC}
LIBS	= ${NETCDF_LIB} 

.f.o: $*.f
	${F77} -c ${FFLAGS} ${INCS} $*.f

lsl2list:	$(OBJS)
	${F77} -o lsl2list $(OBJS) ${INCS} $(LIBS)
