F77	= ${FORTRAN}
FFLAGS	= -O
OBJS	= lidar.o ${LAGRANTO}/lib/times.a ${LAGRANTO}/lib/iotra.a ${LAGRANTO}/lib/ioinp.a ${LAGRANTO}/lib/inter.a ${LAGRANTO}/lib/libcdfio.a ${LAGRANTO}/lib/libcdfplus.a
INCS	= ${NETCDF_INC}
LIBS	= ${NETCDF_LIB} 

.f.o: $*.f 
	${F77} -c ${FFLAGS} ${INCS} $*.f

lidar:	$(OBJS)
	${F77} -o lidar $(OBJS) ${INCS} $(LIBS)
