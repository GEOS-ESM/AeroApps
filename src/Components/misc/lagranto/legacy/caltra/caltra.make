F77	= ${FORTRAN}
FFLAGS	= -O
OBJS	= caltra.o ${LAGRANTO}/lib/times.a ${LAGRANTO}/lib/iotra.a ${LAGRANTO}/lib/ioinp.a ${LAGRANTO}/lib/inter.a ${LAGRANTO}/lib/libcdfio.a ${LAGRANTO}/lib/libcdfplus.a
INCS	= ${NETCDF_INC}
LIBS	= ${NETCDF_LIB} 

.f.o: $*.f 
	${F77} -c ${FFLAGS} ${INCS} $*.f

caltra:	$(OBJS)
	${F77} -o caltra $(OBJS) ${INCS} $(LIBS)
