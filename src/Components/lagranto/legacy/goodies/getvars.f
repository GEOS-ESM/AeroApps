      program getvarnames
      
c     ***********************************************************************
c     * Get list of variables on netCDF file                                *
c     * Michael Sprenger / Spring, summer 2016                              *
c     ***********************************************************************

      use netcdf
      
      implicit none
 
      integer      nvars,ierr
      integer      cdfid,i
      character*80 cdfname
      character*80 vnam(200)
 
      integer   iargc
      character*(80) arg

c     check for sufficient requested arguments
      if (iargc().ne.1) then
         print*,'USAGE: getvars NetCDF-filename'
         call exit(1)
      endif
 
c     read and transform input
      call getarg(1,arg)
      cdfname=trim(arg)
      
c     Get list and write to standard out      
      call input_open    (cdfid,cdfname)
      call input_getvars (cdfid,vnam,nvars)
      call input_close   (cdfid)
      do i=1,nvars
        write(*,*)  vnam(i)
      enddo

      end 
