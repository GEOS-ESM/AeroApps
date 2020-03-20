F77	= ${FORTRAN}
FFLAGS	= -O
OBJS	= gettidiff.o ${LAGRANTO}/lib/times.a
LIBS	= 

.f.o: $*.f $(INCS)
	${F77} -c ${FFLAGS} $*.f

gettidiff:	$(OBJS)
	${F77} -o gettidiff $(OBJS) $(LIBS)
