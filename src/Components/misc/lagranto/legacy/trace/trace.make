F77	= ${FORTRAN}
FFLAGS	= -O
OBJS	= ${LAGRANTO}/lib/times.a ${LAGRANTO}/lib/iotra.a ${LAGRANTO}/lib/ioinp.a ${LAGRANTO}/lib/inter.a ${LAGRANTO}/lib/libcdfio.a ${LAGRANTO}/lib/libcdfplus.a
INCS	= ${NETCDF_INC}
LIBS	= ${NETCDF_LIB} 

trace:	$(OBJS)
	${F77} -o trace trace.f90 calvar.f $(OBJS) ${INCS} $(LIBS)
