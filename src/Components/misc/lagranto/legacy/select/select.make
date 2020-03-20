F77	= ${FORTRAN}
FFLAGS	= -O
OBJS	= special.o select.o ${LAGRANTO}/lib/iotra.a ${LAGRANTO}/lib/libcdfio.a ${LAGRANTO}/lib/libcdfplus.a
INCS	= ${NETCDF_INC}
LIBS	= ${NETCDF_LIB} 

.f.o: $*.f
	${F77} -c ${FFLAGS} ${INCS} $*.f

select:	$(OBJS)
	${F77} -o select $(OBJS) ${INCS} $(LIBS)
