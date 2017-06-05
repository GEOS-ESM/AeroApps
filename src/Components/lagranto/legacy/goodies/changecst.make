F77     = ${FORTRAN}
FFLAGS  = -O
OBJS    = changecst.o ${LAGRANTO}/lib/libcdfio.a ${LAGRANTO}/lib/libcdfplus.a
INCS    = ${NETCDF_INC}
LIBS    = ${NETCDF_LIB} 

.f.o: $*.f 
	${F77} -c ${FFLAGS} ${INCS} $*.f

changecst:	$(OBJS)
	${F77} -o changecst $(OBJS) ${INCS} $(LIBS)
