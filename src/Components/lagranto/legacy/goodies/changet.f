      program changetime

C     Changes the time value of a NetCDF file.


C-----declarations------------------------------------------------

      integer   ierr,ntimes
      real      tstart,time
      integer   cdfid,varid
      character*30  filnam

C-----start of program--------------------------------------------

      include 'netcdf.inc'

      call ncpopt(NCVERBOS)

      read(9,10)filnam
   10 format(a30)
      read(9,*)time

C     Open the data file

      cdfid=ncopn(filnam,NCWRITE,ierr)

C     Get time value from file

      call gettimes(cdfid,tstart,ntimes,ierr)

C     Get index for time variable

      varid=ncvid(cdfid,'time',ierr)

C     Overwrite time-value

      call ncvpt1(cdfid,varid,1,time,ierr)

      write(*,20)'file ',trim(filnam),' time value changed from ',
     >           tstart,' to ',time
   20 format(a,a,a,f8.2,a,f8.2)

C     Close open NetCDF files

      call clscdf(cdfid,ierr)

      end
