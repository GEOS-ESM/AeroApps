      program changecstname

C     Changes the constants file name of a NetCDF file.

      implicit none

C-----declarations------------------------------------------------

      integer   ierr,lenstr
      integer   cdfid,nctype
      character*80 filnam,cstnam,ocstnam

      integer   iargc
      character*(80) arg

      integer   strbeg,strend

C-----start of program--------------------------------------------

      include 'netcdf.inc'

      call ncpopt(NCVERBOS)

c     check for sufficient requested arguments
      if (iargc().ne.2) then
        print*,'USAGE: changecst filename new_cst_filename'
        call exit(1)
      endif

c     read and transform input
      call getarg(1,arg)
      filnam=trim(arg)

      call getarg(2,arg)
      cstnam=trim(arg)

C     Open the data file

      cdfid=ncopn(trim(filnam),NCWRITE,ierr)

C     Put file into define mode

      call ncredf(cdfid,ierr)

C     Get current value of the constants_file_name attribute

      call ncainq(cdfid,NCGLOBAL,'constants_file_name',nctype,
     >            lenstr,ierr)
      call ncagtc(cdfid,NCGLOBAL,'constants_file_name',ocstnam,
     >            lenstr,ierr)

C     Overwrite constants_file_name
      lenstr = len_trim(cstnam)
      call ncaptc(cdfid,NCGLOBAL,'constants_file_name',NCCHAR,
     >            lenstr,cstnam(1:lenstr),ierr)

      print*,'file ',trim(filnam),
     >       ' cst_file_name changed from ',trim(ocstnam),
     >       ' to ',trim(cstnam)

C     Close open NetCDF files

      call clscdf(cdfid,ierr)

      end
