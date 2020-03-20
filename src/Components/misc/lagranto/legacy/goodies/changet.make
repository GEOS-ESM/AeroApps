F77     = ${FORTRAN}
FFLAGS  = -O
OBJS    = changet.o ${LAGRANTO}/lib/libcdfio.a ${LAGRANTO}/lib/libcdfplus.a
INCS    = ${NETCDF_INC}
LIBS    = ${NETCDF_LIB} 

.f.o: $*.f
	${F77} -c ${FFLAGS} ${INCS} $*.f

changet:	$(OBJS)
	${F77} -o changet $(OBJS) ${INCS} $(LIBS)
