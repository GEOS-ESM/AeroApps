#ifdef _DLAPACK_
# define _SPPTRF	DPPTRF
# define _SPPTRS	DPPTRS
# define _SCOPY		DCOPY
# define _SAXPY		DAXPY
# define _SDOT		DDOT
# define _SNRM2		DNRM2
# define _SSPMV		DSPMV
# define _SGEMV		DGEMV
# define _SGETRF	DGETRF
# define _SGETRS	DGETRS
#else
# define _SPPTRF	SPPTRF
# define _SPPTRS	SPPTRS
# define _SCOPY		SCOPY
# define _SAXPY		SAXPY
# define _SDOT		SDOT
# define _SNRM2		SNRM2
# define _SSPMV		SSPMV
# define _SGEMV		SGEMV
# define _SGETRF	SGETRF
# define _SGETRS	SGETRS
#endif
