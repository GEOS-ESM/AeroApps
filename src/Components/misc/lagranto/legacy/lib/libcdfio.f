      subroutine clscdf (cdfid, error)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine closes an open netCDF file.
c     Aguments :
c        cdfid  int  input   the id of the file to be closed.
c        error  int  output  indicates possible errors found in this
c                            routine.
c                            error = 0   no errors detected.
c                            error = 1   error detected.
c     History:
c        Nov. 91  PPM  UW  Created.
c-----------------------------------------------------------------------

      include "netcdf.inc"

c     Argument declarations.
      integer      cdfid, error

c     Local variable declarations.
      integer      ncopts

c     Get current value of error options.
      call ncgopt (ncopts)

c     Make sure netCDF errors do not abort execution.
      call ncpopt (NCVERBOS)

c     Close requested file.
      call ncclos (cdfid, error)

c     Reset error options.
      call ncpopt (ncopts)

      end


      subroutine crecdf (filnam, cdfid, phymin, phymax, ndim, cfn, 
     &                   error) 
c-----------------------------------------------------------------------
c     Purpose:
c        This routine is called to create a netCDF file for use with 
c        the UWGAP plotting package.
c           Any netCDF file written to must be closed with the call
c        'call clscdf(cdfid,error)', where cdfid and error are
c        as in the argumentlist below. 
c     Arguments:
c        filnam  char  input   the user-supplied netCDF file name.
c        cdfid   int   output  the file-identifier
c        phymin  real  input   the minimum physical dimension of the
c                              entire physical domain along each axis.
c                              phymin is dimensioned (ndim)
c        phymax  real  input   the maximum physical dimension of the
c                              entire physical domain along each axis.
c                              phymax is dimensioned (ndim)
c        ndim    int   input   the number of dimensions in the file
c                              (i.e. number of elements in phymin,
c                              phymax)
c        cfn     char  input   constants file name 
c                              ('0' = no constants file).
c        error   int   output  indicates possible errors found in this
c                              routine.
c                              error = 0   no errors detected.
c                              error = 1   error detected.
c     History:
c        Nov. 91  PPM  UW  Created cr3df.
c        Jan. 92  CS   UW  Created crecdf.
c-----------------------------------------------------------------------

      include "netcdf.inc"

c     Argument declarations.
      integer        MAXDIM
      parameter      (MAXDIM=4)
      integer        ndim, error
      character *(*) filnam,cfn
      real           phymin(*), phymax(*)

c     Local variable declarations.
      character *(20) attnam
      character *(1)  chrid(MAXDIM)
      integer         cdfid, k, ibeg, iend, lenfil, ncopts
      data            chrid/'x','y','z','a'/

c     Get current value of error options, and make sure netCDF-errors do
c     not abort execution
      call ncgopt (ncopts)
      call ncpopt(NCVERBOS)

c     Initially set error to indicate no errors.
      error = 0

c     create the netCDF file
      cdfid = nccre (trim(filnam), NCCLOB, error)
      if (error.ne.0) go to 920

c     define global attributes
      do k=1,ndim
        attnam(1:3)='dom'
        attnam(4:4)=chrid(k)
        attnam(5:7)='min'
        attnam=attnam(1:7)
        call ncapt(cdfid,NCGLOBAL,attnam,NCFLOAT,1,phymin(k),error)
        if (error.gt.0) goto 920

        attnam(1:3)='dom'
        attnam(4:4)=chrid(k)
        attnam(5:7)='max'
        attnam=attnam(1:7)
        call ncapt(cdfid,NCGLOBAL,attnam,NCFLOAT,1,phymax(k),error)
        if (error.gt.0) goto 920
      enddo

c     define constants file name
      if (cfn.ne.'0') then
        call ncaptc (cdfid, NCGLOBAL, 'constants_file_name',
c    &             NCCHAR, len_trim(cfn)+1, cfn // char(0) , error)
     &             NCCHAR, len_trim(cfn), cfn , error)
        if (error.gt.0) goto 920
      endif

c     End variable definitions.
      call ncendf (cdfid, error)
      if (error.gt.0) goto 920

c     normal exit
      call ncpopt (ncopts)
      return

c     error exit
 920  write (6, *) 'ERROR: An error occurred while attempting to ',
     &             'create the data file in subroutine crecdf.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      end

      subroutine opncdf(filnam, cdfid,
     &                 phymin, phymax, ndim, varnam, nvar, cfn, error) 
c-----------------------------------------------------------------------
c     Purpose:
c        This routine is called to open a netCDF file for read and write
c        with the UWGAP plotting package. 
c     Arguments:
c        filnam  char  input   the user-supplied netCDF file name.
c        cdfid   int   output  the file-identifier
c        phymin  real  output  the minimum physical dimension of the
c                              entire physical domain along each axis.
c                              phymin is dimensioned (ndim)
c        phymax  real  output  the maximum physical dimension of the
c                              entire physical domain along each axis.
c                              phymax is dimensioned (ndim)
c        ndim    int   output  the number of dimensions in the file
c                              (i.e. number of elements in phymin,
c                              phymax)
c        varnam  char  output  an array containing the variable names.
c                              varnam is dimensioned (nvar).
c        nvar    int   output  the number of variables in the file
c        cfn     char  output  constants file name 
c                              ('0'=no constants file).
c        error   int   output  indicates possible errors found in this
c                              routine.
c                              error = 0   no errors detected.
c                              error = 1   error detected.
c     History:
c        Nov. 91  PPM  UW  Created cr3df.
c        Jan. 92  CS   UW  Created opncdf.
c-----------------------------------------------------------------------

      include "netcdf.inc"

c     Argument declarations.
      integer        MAXDIM
      parameter      (MAXDIM=4)
      integer        ndim, nvar, error
      character *(*) filnam, varnam(*),cfn
      real           phymin(*), phymax(*)

c     Local variable declarations.
      character *(20) attnam,vnam
      character *(1)  chrid(MAXDIM)
      integer         cdfid, i,k
      integer         ncopts, ndims,ngatts,recdim
      integer         nvdims,vartyp,nvatts,vardim(MAXDIM)
      real            attval
      integer         lenstr
      data            chrid/'x','y','z','a'/
      data            lenstr/80/

c     Get current value of error options and make sure netCDF-errors do 
c     not abort execution
      call ncgopt (ncopts)
      call ncpopt(NCVERBOS)

c     Initially set error to indicate no errors.
      error = 0

c     open the netCDF file for write
      cdfid = ncopn (trim(filnam), NCWRITE, error)
      if (error.ne.0) then
c       try to open the netCDF file for read
        cdfid = ncopn (trim(filnam), NCNOWRIT, error)
        if (error.ne.0) go to 920
      endif

c     inquire for number of variables
      call ncinq(cdfid,ndims,nvar,ngatts,recdim,error)
      if (error.eq.1) goto 920

c     read the variables
      do i=1,nvar
        call ncvinq(cdfid,i,varnam(i),vartyp,nvdims,vardim,
     &                         nvatts,error)
        if (vartyp.ne.NCFLOAT) error=1
        if (error.gt.0) goto 920
      enddo

c     get global attributes
      k=0
  100 continue
        k=k+1
        attnam(1:3)='dom'
        attnam(4:4)=chrid(k)
        attnam(5:7)='min'
        attnam=attnam(1:7)

c       switch off error message
        call ncpopt(0)

c       check whether dimension k is present
        call ncagt(cdfid,NCGLOBAL,attnam,attval,error)
        if (error.gt.0) goto 110
        phymin(k)=attval

        attnam(1:3)='dom'
        attnam(4:4)=chrid(k)
        attnam(5:7)='max'
        attnam=attnam(1:7)
        call ncagt(cdfid,NCGLOBAL,attnam,attval,error)
        if (error.gt.0) goto 920
        phymax(k)=attval
      if (k.lt.3) goto 100
      k=k+1

c     define ndim-parameter
 110  continue
      ndim=k-1
      error=0

c     switch on error messages
      call ncpopt(NCVERBOS)

c     get constants file name
c      call ncagt(cdfid,NCGLOBAL,'constants_file_name',cfn,error)
c     ! chrigel
      call ncagtc(cdfid,NCGLOBAL,'constants_file_name',cfn,lenstr,error)
      if (error.gt.0) cfn='0'

c     normal exit
      call ncpopt (ncopts)
      return

c     error exit
 920  write (6, *) 'ERROR: An error occurred while attempting to ',
     &             'read the data file in subroutine opncdf.'
      call ncclos (cdfid, error)
      call ncpopt (ncopts)
      end


      subroutine readcdf(filnam, cdfid,
     &                 phymin, phymax, ndim, varnam, nvar, cfn, error) 
c-----------------------------------------------------------------------
c     Purpose:
c        This routine is called to open a netCDF file for read 
c        with the UWGAP plotting package. 
c     Arguments:
c        filnam  char  input   the user-supplied netCDF file name.
c        cdfid   int   output  the file-identifier
c        phymin  real  output  the minimum physical dimension of the
c                              entire physical domain along each axis.
c                              phymin is dimensioned (ndim)
c        phymax  real  output  the maximum physical dimension of the
c                              entire physical domain along each axis.
c                              phymax is dimensioned (ndim)
c        ndim    int   output  the number of dimensions in the file
c                              (i.e. number of elements in phymin,
c                              phymax)
c        varnam  char  output  an array containing the variable names.
c                              varnam is dimensioned (nvar).
c        nvar    int   output  the number of variables in the file
c        cfn     char  output  constants file name 
c                              ('0'=no constants file).
c        error   int   output  indicates possible errors found in this
c                              routine.
c                              error = 0   no errors detected.
c                              error = 1   error detected.
c     History:
c        Nov. 91  PPM  UW  Created cr3df.
c        Jan. 92  CS   UW  Created opncdf.
c-----------------------------------------------------------------------

      include "netcdf.inc"

c     Argument declarations.
      integer        MAXDIM
      parameter      (MAXDIM=4)
      integer        ndim, nvar, error
      character *(*) filnam, varnam(*),cfn
      real           phymin(*), phymax(*)


c     Local variable declarations.
      character *(20) attnam
      character *(1)  chrid(MAXDIM)
      integer         cdfid, i,k
      integer         ncopts, ndims,ngatts,recdim
      integer         nvdims,vartyp,nvatts,vardim(MAXDIM)
      real            attval
      integer         lenstr
      data            chrid/'x','y','z','a'/
      data            lenstr/80/

c     Get current value of error options.
      call ncgopt (ncopts)

c     make sure netCDF-errors do not abort execution
      call ncpopt(NCVERBOS)

c     Initially set error to indicate no errors.
      error = 0

c     open the netCDF file for read
      cdfid = ncopn (trim(filnam), NCNOWRIT, error)
      if (error.ne.0) go to 920

c     inquire for number of variables
      call ncinq(cdfid,ndims,nvar,ngatts,recdim,error)
      if (error.eq.1) goto 920

c     read the variables
      do i=1,nvar
        call ncvinq(cdfid,i,varnam(i),vartyp,nvdims,vardim,
     &                         nvatts,error)
        if (vartyp.ne.NCFLOAT) error=1
c       print *,varnam(i),nvdims,nvatts
        if (error.gt.0) goto 920
      enddo

c     get global attributes
      k=0
  100 continue
        k=k+1
        attnam(1:3)='dom'
        attnam(4:4)=chrid(k)
        attnam(5:7)='min'
        attnam=attnam(1:7)

c       switch off error message
        call ncpopt(0)

c       check whether dimension k is present
        call ncagt(cdfid,NCGLOBAL,attnam,attval,error)
        if (error.gt.0) goto 110
        phymin(k)=attval

        attnam(1:3)='dom'
        attnam(4:4)=chrid(k)
        attnam(5:7)='max'
        attnam=attnam(1:7)
        call ncagt(cdfid,NCGLOBAL,attnam,attval,error)
        if (error.gt.0) goto 920
        phymax(k)=attval
      if (k.lt.4) goto 100
      k=k+1

c     define ndim-parameter
 110  continue
      ndim=k-1
      error=0

c     switch on error messages
      call ncpopt(NCVERBOS)

c     get constants file name
c      call ncagt(cdfid,NCGLOBAL,'constants_file_name',cfn,error)
c     ! chrigel
      call ncagtc(cdfid,NCGLOBAL,'constants_file_name',cfn,lenstr,error)
      if (error.gt.0) cfn='0'
c     print *,cfn

c     normal exit
      call ncpopt (ncopts)
      return

c     error exit
 920  write (6, *) 'ERROR: An error occurred while attempting to ',
     &             'read the data file in subroutine opncdf.'
      call ncclos (cdfid, error)
      call ncpopt (ncopts)
      end



      subroutine getcdf (cdfid, varnam, ndim, misdat, 
     &                       vardim, varmin, varmax, stag, dat, error)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine is called to get a variable and its attributes
c        from a netCDF file for use with the UWGAP plotting package. 
c        It is assumed that the data is floating-point data. Prior to
c        calling this routine, the file must be opened with a call to
c        opncdf.
c     Arguments:
c        cdfid   int   input   file-identifier
c                              (can be obtained by calling routine 
c                              opncdf)
c        varnam  char  input   the user-supplied variable name.
c                              (can be obtained by calling routine 
c                              opncdf)
c        ndim    int   output  the number of dimensions (ndim<=4)
c        misdat  real  output  missing data value for the variable. 
c        vardim  int   output  the dimensions of the variable.
c                              is dimensioned at least (ndim). 
c        varmin  real  output  the location in physical space of the
c                              origin of each variable.
c                              is dimensioned at least Min(ndim,3). 
c        varmax  real  output  the extent of each variable in physical
c                              space.
c                              is dimensioned at least Min(ndim,3). 
c        stag    real  output  the grid staggering for each variable.
c                              is dimensioned at least Min(ndim,3). 
c        dat     real  output  data-array dimensioned suffiecently 
c                              large, at least 
c                              vardim(1)* ... vardim(ndim)
c        error   int   output  indicates possible errors found in this
c                              routine.
c                              error = 0   no errors detected.
c                              error = 1   error detected.
c     History:
c        Nov. 91  PPM  UW  Created cr3df.
c        Jan. 92  CS   UW  Created getcdf.
c-----------------------------------------------------------------------

      include "netcdf.inc"

c     Argument declarations.
      integer        MAXDIM
      parameter      (MAXDIM=4)
      character *(*) varnam
      integer        vardim(*), ndim, error, cdfid
      real           misdat,  stag(*), varmin(*), varmax(*), dat(*)

c     Local variable declarations.
      character *(20) dimnam(100),attnam
      character *(1)  chrid(MAXDIM)
      integer         id,i,k,corner(MAXDIM)
      integer         ndims,nvars,ngatts,recdim,dimsiz(100)
      integer         vartyp,nvatts, ncopts
      data            chrid/'x','y','z','a'/
      data            corner/1,1,1,1/

c     Get current value of error options, and make sure netCDF-errors do 
c     not abort execution 
      call ncgopt (ncopts)
      call ncpopt(NCVERBOS)

c     Initially set error to indicate no errors.
      error = 0

c     inquire for number of dimensions
      call ncinq(cdfid,ndims,nvars,ngatts,recdim,error)
      if (error.eq.1) goto 920

c     read dimension-table
      do i=1,ndims 
        call ncdinq(cdfid,i,dimnam(i),dimsiz(i),error)
        if (error.gt.0) goto 920
      enddo

c     get id of the variable
      id=ncvid(cdfid,varnam,error)
      if (error.eq.1) goto 910

c     inquire about variable
      call ncvinq(cdfid,id,varnam,vartyp,ndim,vardim,nvatts,error)
      if (vartyp.ne.NCFLOAT) error=1
      if (error.gt.0) goto 920

c     Make sure ndim <= MAXDIM.
      if (ndim.gt.MAXDIM) then
         error = 1
         go to 900
      endif

c     get dimensions from dimension-table
      do k=1,ndim 
        vardim(k)=dimsiz(vardim(k))
      enddo

c     get attributes
      do k=1,min0(ndim,3)
c       get staggering
        attnam(1:1)=chrid(k)
        attnam(2:5)='stag'
        attnam=attnam(1:5)
        call ncagt(cdfid,id,attnam,stag(k),error)
        if (error.gt.0) goto 920
c       get min postion
        attnam(1:1)=chrid(k)
        attnam(2:4)='min'
        attnam=attnam(1:4)
        call ncagt(cdfid,id,attnam,varmin(k),error)
        if (error.gt.0) goto 920
c       get max position     
        attnam(1:1)=chrid(k)
        attnam(2:4)='max'
        attnam=attnam(1:4)
        call ncagt(cdfid,id,attnam,varmax(k),error)
        if (error.gt.0) goto 920     
      enddo

c     get missing data value
      call ncagt(cdfid,id,'missing_data',misdat,error)
      if (error.gt.0) goto 920     

c     get data
      call ncvgt(cdfid,id,corner,vardim,dat,error)
      if (error.gt.0) goto 920

c     normal exit
      call ncpopt (ncopts)
      return


c     Error exits.
 900  write (6, *) 'ERROR: When calling getcdf, the number of ',
     &             'variable dimensions must be less or equal 4.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      return

 910  write (6, *) 'ERROR: The selected variable could not be found ',       
     &             'in the file by getcdf.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      return

 920  write (6, *) 'ERROR: An error occurred while attempting to ',
     &             'read the data file in subroutine getcdf.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      return

      end


      subroutine putcdf (cdfid, varnam, ndim, misdat, 
     &                       vardim, varmin, varmax, stag, dat, error)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine is called to put a variable and its attributes
c        onto a netCDF file for use with the UWGAP plotting package. 
c        It is assumed that the data is floating-point data. Prior to 
c        calling this routine, the file must be created (crecdf) or
c        opened (opncdf). 
c           Any netCDF file written to must be closed with the call
c        call ncclos(cdfid,error), where cdfid and error are
c        as in the argumentlist below.
c     Arguments:
c        cdfid   int   input   file-identifier
c                              (can be obtained by calling routine 
c                              opncdf)
c        varnam  char  input   the user-supplied variable name.
c                              (can be obtained by calling routine 
c                              opncdf)
c        ndim    int   input   the number of dimensions (ndim<=4)
c        misdat  real  input   missing data value for the variable. 
c        vardim  int   input   the dimensions of the variable.
c                              is dimensioned at least (ndim). 
c        varmin  real  input   the location in physical space of the
c                              origin of each variable.
c                              is dimensioned at least Min(ndim,3). 
c        varmax  real  input   the extent of each variable in physical
c                              space.
c                              is dimensioned at least Min(ndim,3).  
c        stag    real  input   the grid staggering for each variable.
c                              is dimensioned at least Min(ndim,3).  
c        dat     real  input   data-array dimensioned suffiecently 
c                              large, at least 
c                              vardim(1)* ... vardim(ndim)
c        error   int   output  indicates possible errors found in this
c                              routine.
c                              error = 0   no errors detected.
c                              error = 1   error detected.
c     History:
c        Nov. 91  PPM  UW  Created cr3df, wr3df.
c        Jan. 92  CS   UW  Created putcdf.
c-----------------------------------------------------------------------

      include "netcdf.inc"

c     Argument declarations.
      integer        MAXDIM
      parameter      (MAXDIM=4)
      character *(*) varnam
      integer        vardim(*), ndim, error, cdfid
      real           misdat,  stag(*), varmin(*), varmax(*), dat(*)

c     Local variable declarations.
      character *(20) dimnam,attnam,dimchk
      character *(1)  chrid(MAXDIM)
      character *(20) dimnams(MAXNCDIM)
      integer         dimvals(MAXNCDIM)
      integer         numdims,numvars,numgats,dimulim
      integer         id,did(MAXDIM),i,k,corner(MAXDIM)
      integer         ncopts
      integer         ibeg,iend
      data            chrid/'x','y','z','t'/
      data            corner/1,1,1,1/

c     Get current value of error options.
      call ncgopt (ncopts)

c     make sure netCDF-errors do not abort execution
      call ncpopt(NCVERBOS)

c     Initially set error to indicate no errors.
      error = 0

c     Make sure ndim <= MAXDIM.
      if (ndim.gt.MAXDIM) then
         error = 1
         go to 900
      endif

c     Read existing dimensions-declarations from the file
      call ncinq(cdfid,numdims,numvars,numgats,dimulim,error)
      if (error.ne.0) numdims=0
      if (numdims.gt.0) then
        do i=1,numdims
          call ncdinq(cdfid,i,dimnams(i),dimvals(i),error)
c          print *,dimnams(i),dimvals(i)
        enddo
      endif

c     put file into define mode
      call ncredf(cdfid,error)
      if (error.ne.0) goto 920

c     define the dimension
      do k=1,ndim
c       define the dimension-name
        dimnam(1:3)='dim'
        dimnam(4:4)=chrid(k)
        dimnam(5:5)='_'
        dimnam(6:5+len_trim(varnam))=trim(varnam)
        dimnam=dimnam(1:5+len_trim(varnam))
        did(k)=-1
        if (numdims.gt.0) then
c         check if an existing dimension-declaration can be used
c         instead of defining a nuw dimension
          do i=1,numdims
            dimchk=dimnams(i)
            if ((vardim(k).eq.dimvals(i)).and.
     &        (dimnam(1:4).eq.dimchk(1:4))) then 
              did(k)=i
              goto 100
            endif
          enddo
 100      continue
        endif
        if (did(k).lt.0) then
c         define the dimension
          did(k)=ncddef(cdfid,dimnam,vardim(k),error)
          if (error.ne.0) goto 920
        endif
      enddo

c     define variable
      id=ncvdef(cdfid,varnam,NCFLOAT,ndim,did,error)
      if (error.ne.0) goto 920

c     define attributes
      do k=1,min0(ndim,3)
c       staggering
        attnam(1:1)=chrid(k)
        attnam(2:5)='stag'
        attnam=attnam(1:5)
        call ncapt(cdfid,id,attnam,NCFLOAT,1,stag(k),error)
        if (error.gt.0) goto 920
c       min postion
        attnam(1:1)=chrid(k)
        attnam(2:4)='min'
        attnam=attnam(1:4)
        call ncapt(cdfid,id,attnam,NCFLOAT,1,varmin(k),error)
        if (error.gt.0) goto 920
c       max position     
        attnam(1:1)=chrid(k)
        attnam(2:4)='max'
        attnam=attnam(1:4)
        call ncapt(cdfid,id,attnam,NCFLOAT,1,varmax(k),error)
        if (error.gt.0) goto 920     
      enddo

c     define missing data value
      call ncapt(cdfid,id,'missing_data',NCFLOAT,1,misdat,error)
      if (error.gt.0) goto 920     

c     leave define mode
      call ncendf(cdfid,error)
      if (error.gt.0) goto 920     

c     define data
      call ncvpt(cdfid,id,corner,vardim,dat,error)
      if (error.gt.0) goto 920

c     synchronyse output to disk and exit
      call ncsnc (cdfid,error)
      call ncpopt (ncopts)
      return

c     Error exits.
 900  write (6, *) 'ERROR: When calling putcdf, the number of ',
     &             'variable dimensions must be less or equal 4.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      return

 920  write (6, *) 'ERROR: An error occurred while attempting to ',
     &             'write the data file in subroutine putcdf.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      return
      end
c
c
      subroutine getdef (cdfid, varnam, ndim, misdat, 
     &                              vardim, varmin, varmax, stag, error)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine is called to get the dimensions and attributes of 
c        a variable from an IVE-NetCDF file for use with the IVE plotting
c        package. Prior to calling this routine, the file must be opened
c        with a call to opncdf.
c     Arguments:
c        cdfid   int   input   file-identifier
c                              (can be obtained by calling routine 
c                              opncdf)
c        varnam  char  input   the user-supplied variable name.
c                              (can be obtained by calling routine 
c                              opncdf)
c        ndim    int   output  the number of dimensions (ndim<=4)
c        misdat  real  output  missing data value for the variable. 
c        vardim  int   output  the dimensions of the variable.
c                              Is dimensioned at least (ndim). 
c        varmin  real  output  the location in physical space of the
c                              origin of each variable. 
c                              Is dimensioned at least Min(3,ndim). 
c        varmax  real  output  the extend of each variable in physical
c                              space.
c                              Is dimensioned at least Min(3,ndim). 
c        stag    real  output  the grid staggering for each variable.
c                              Is dimensioned at least Min(3,ndim). 
c        error   int   output  indicates possible errors found in this
c                              routine.
c                              error = 0   no errors detected.
c                              error = 1   the variable is not on the file.
c                              error =10   other errors.
c     History:
c       Apr. 93    Christoph Schaer (ETHZ)     Created.
c-----------------------------------------------------------------------

      include "netcdf.inc"

c     Argument declarations.
      integer        MAXDIM
      parameter      (MAXDIM=4)
      character *(*) varnam
      integer        vardim(*), ndim, error, cdfid
      real           misdat,  stag(*), varmin(*), varmax(*)

c     Local variable declarations.
      character *(20) dimnam(MAXNCDIM),attnam,vnam
      character *(1)  chrid(MAXDIM)
      integer         id,i,k
      integer         ndims,nvars,ngatts,recdim,dimsiz(MAXNCDIM)
      integer         vartyp,nvatts, ncopts
      data            chrid/'x','y','z','t'/

c     Get current value of error options.
      call ncgopt (ncopts)

c     make sure NetCDF-errors do not abort execution
      call ncpopt(NCVERBOS)

c     Initially set error to indicate no errors.
      error = 0

c     inquire for number of dimensions
      call ncinq(cdfid,ndims,nvars,ngatts,recdim,error)
      if (error.eq.1) goto 920

c     read dimension-table
      do i=1,ndims 
        call ncdinq(cdfid,i,dimnam(i),dimsiz(i),error)
        if (error.gt.0) goto 920
      enddo

c     get id of the variable
      id=ncvid(cdfid,varnam,error)
      if (error.eq.1) goto 910

c     inquire about variable
      call ncvinq(cdfid,id,vnam,vartyp,ndim,vardim,nvatts,error)
      if (vartyp.ne.NCFLOAT) error=1
      if (error.gt.0) goto 920

c     Make sure ndim <= MAXDIM.
      if (ndim.gt.MAXDIM) then
         error = 1
         go to 900
      endif

c     get dimensions from dimension-table
      do k=1,ndim 
        vardim(k)=dimsiz(vardim(k))
      enddo

c     get attributes
      do k=1,min0(3,ndim)
c       get min postion
        attnam(1:1)=chrid(k)
        attnam(2:4)='min'
        attnam=attnam(1:4)
        call ncagt(cdfid,id,attnam,varmin(k),error)
        if (error.gt.0) goto 920
c       get max position     
        attnam(1:1)=chrid(k)
        attnam(2:4)='max'
        attnam=attnam(1:4)
        call ncagt(cdfid,id,attnam,varmax(k),error)
        if (error.gt.0) goto 920     
c       get staggering
        attnam(1:1)=chrid(k)
        attnam(2:5)='stag'
        attnam=attnam(1:5)
        call ncagt(cdfid,id,attnam,stag(k),error)
        if (error.gt.0) goto 920
      enddo

c     get missing data value
      call ncagt(cdfid,id,'missing_data',misdat,error)
      if (error.gt.0) goto 920     

c     normal exit
      call ncpopt (ncopts)
      return

c     Error exits.
 900  write (6, *) '*ERROR*: When calling getcdf, the number of ',
     &             'variable dimensions must be less or equal 4.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      return

 910  write (6, *) '*ERROR*: The selected variable could not be found ',       
     &             'in the file by getcdf.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      return

 920  write (6, *) '*ERROR*: An error occurred while attempting to ',
     &             'read the data file in subroutine getcdf.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      end



      subroutine getdat(cdfid, varnam, time, level, dat, error)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine is called to read the data of a variable
c        from an IVE-NetCDF file for use with the IVE plotting package. 
c        Prior to calling this routine, the file must be opened with 
c        a call to opncdf (for extension) or crecdf (for creation) or
c        readcdf (for readonly).
c     Arguments:
c        cdfid   int   input   file-identifier
c                              (must be obtained by calling routine 
c                              opncdf,readcdf  or crecdf)
c        varnam  char  input   the user-supplied variable name (must 
c                              previously be defined with a call to
c                              putdef)
c        time    real  input   the user-supplied time-level of the
c                              data to be read from the file (the time-
c                              levels stored in the file can be obtained
c                              with a call to gettimes). 
c        level   int   input   the horizontal level(s) to be read 
c                              to the NetCDF file. Suppose that the
c                              variable is defined as (nx,ny,nz,nt).
c                              level>0: the call reads the subdomain
c                                       (1:nx,1:ny,level,itimes)
c                              level=0: the call reads the subdomain
c                                       (1:nx,1:ny,1:nz,itimes)
c                              Here itimes is the time-index corresponding
c                              to the value of 'time'. 
c        dat     real  output  data-array dimensioned sufficiently 
c                              large. The dimensions (nx,ny,nz)
c                              of the variable must previously be defined
c                              with a call to putdef. No previous 
c                              definition of the time-dimension is
c                              required.
c        error   int   output  indicates possible errors found in this
c                              routine.
c                              error = 0   no errors detected.
c                              error = 1   the variable is not present on
c                                          the file.
c                              error = 2   the value of 'time' is not
c                                          known.to the file.
c                              error = 3   inconsistent value of level
c                              error =10   another error.
c     History:
c       March 93    Heini Wernli (ETHZ)      Created wr2cdf.
c       April 93    Bettina Messmer (ETHZ)   Created putdat.
c       June  93    Christoph Schaer (ETHZ)  Created getdat
c       Sept. 07    Johannes Jenkner (ETHZ)  Integer and double times
c-----------------------------------------------------------------------

      include "netcdf.inc"

C     Declaration of local variables
      character*(*) varnam
      character*(20) chars,dimnam
      integer cdfid

      real        dat(*)
      real        misdat,varmin(3),varmax(3),stag(3)
      real        time, timeval

      integer     corner(4),edgeln(4),didtim,vardim(4),ndims
      integer     error, ierr
      integer     level,ntime
      integer     idtime,idvar,iflag
      integer     i

      integer     vtyp,dn,nat
      integer     dims(4)
      integer     inttime
      double precision doubletime

      call ncpopt(NCVERBOS)

c     access the variable
      call getdef (cdfid, trim(varnam), ndims, misdat, 
     &                           vardim, varmin, varmax, stag, ierr)
      if (ierr.ne.0) then
        print *,'*ERROR* in getdef in getdat'
        error=1
        return
      endif
      idvar=ncvid(cdfid,trim(varnam),ierr)
      if (ierr.ne.0) then
        print *,'*ERROR* in ncvid in getdat'
        error=1
        return
      endif

C     Get times-array
      didtim=ncdid(cdfid,'time',ierr)
      if (ierr.ne.0) then
        print *,'*ERROR* didtim in getdat'
        error=10
        return
      endif
      call ncdinq(cdfid,didtim,chars,ntime,ierr)
      if (ierr.ne.0) then
        print *,'*ERROR* in ncdinq in getdat'
        error=10
        return
      endif
      idtime=ncvid(cdfid,'time',ierr)
      if (ierr.ne.0) then
        print *,'*ERROR* in ncvid for time in getdat'
        error=10
        return
      endif
c     find appropriate time-index
      iflag=0
      call ncvinq(cdfid,idtime,dimnam,vtyp,dn,dims,natt,ierr)
      do i=1,ntime
        if (vtyp.eq.5) then
          call ncvgt1(cdfid,idtime,i,timeval,ierr)
        elseif (vtyp.eq.4) then ! integer version 
          call ncvgt1(cdfid,idtime,i,inttime,ierr) 
          timeval=real(inttime)
        elseif (vtyp.eq.6) then ! double precision version
           call ncvgt1(cdfid,idtime,i,doubletime,ierr) 
           timeval=real(doubletime)   
        endif      
        if (ierr.ne.0) print *,'*ERROR* in ncvgt1 in getdat'
        if (time.eq.timeval) iflag=i
      enddo
      if (iflag.eq.0) then
        error=2
        print *,'Error: Unknown time in getdat'
        return
      endif

C     Define data volume to be written (index space)
      corner(1)=1
      corner(2)=1
      edgeln(1)=vardim(1)
      edgeln(2)=vardim(2)
      if (level.eq.0) then
        corner(3)=1
        edgeln(3)=vardim(3)
      else if ((level.le.vardim(3)).and.(level.ge.1)) then
        corner(3)=level
        edgeln(3)=1
      else
        error=3
        return
      endif
      corner(4)=iflag
      edgeln(4)=1

C     Read data from NetCDF file
c      print *,'getdat vor Aufruf ncvgt'
c      print *,'cdfid ',cdfid
c      print *,'idvar ',idvar
c      print *,'corner ',corner
c      print *,'edgeln ',edgeln

      call ncvgt(cdfid,idvar,corner,edgeln,dat,error)
      if (error.ne.0) then
        print *, '*ERROR* in ncvgt in getdat'
        error=10
      endif
      end


      subroutine putdat(cdfid, varnam, time, level, dat, error)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine is called to write the data of a variable
c        to an IVE-NetCDF file for use with the IVE plotting package. 
c        Prior to calling this routine, the file must be opened with 
c        a call to opncdf (for extension) or crecdf (for creation), the 
c        variable must be defined with a call to putdef.
c     Arguments:
c        cdfid   int   input   file-identifier
c                              (must be obtained by calling routine 
c                              opncdf or crecdf)
c        varnam  char  input   the user-supplied variable name (must 
c                              previously be defined with a call to
c                              putdef)
c        time    real  input   the user-supplied time-level of the
c                              data to be written to the file (the time-
c                              levels stored in the file can be obtained
c                              with a call to gettimes). If 'time' is not
c                              yet known to the file, a knew time-level is
c                              allocated and appended to the times-array.
c        level   int input     the horizontal level(s) to be written 
c                              to the NetCDF file. Suppose that the
c                              variable is defined as (nx,ny,nz,nt).
c                              level>0: the call writes the subdomain
c                                       (1:nx,1:ny,level,itimes)
c                              level=0: the call writes the subdomain
c                                       (1:nx,1:ny,1:nz,itimes)
c                              Here itimes is the time-index corresponding
c                              to the value of 'time'. 
c        dat     real  output  data-array dimensioned sufficiently 
c                              large. The dimensions (nx,ny,nz)
c                              of the variable must previously be defined
c                              with a call to putdef. No previous 
c                              definition of the time-dimension is
c                              required.
c        error   int output    indicates possible errors found in this
c                              routine.
c                              error = 0   no errors detected.
c                              error = 1   the variable is not present on
c                                          the file.
c                              error = 2   the value of 'time' is new, but
c                                          appending it would yield a non
c                                          ascending times-array.
c                              error = 3   inconsistent value of level
c                              error =10   another error.
c     History:
c       March 93    Heini Wernli (ETHZ)      Created wr2cdf.
c       April 93    Bettina Messmer (ETHZ)   Created putdat.
c-----------------------------------------------------------------------

      include "netcdf.inc"

C     Declaration of local variables

      character*(*) varnam
      character*(20) chars
      integer cdfid


      real 	dat(*)
      real	misdat,varmin(3),varmax(3),stag(3)
      real	time, timeval
      data	stag/0.,0.,0./

      integer	corner(4),edgeln(4),did(4),vardim(4),ndims
      integer	error, ierr
      integer	level,ntime
      integer	idtime,idvar,iflag
      integer	i

      call ncpopt(NCVERBOS)

c     get definitions of data
      call getdef (cdfid, trim(varnam), ndims, misdat, 
     &                           vardim, varmin, varmax, stag, ierr)
      if (ierr.ne.0)  print *,'*ERROR* in getdef in putdat'

c     get id of variable
      idvar=ncvid(cdfid,trim(varnam),ierr)
      if (ierr.ne.0) print *,'*ERROR* in ncvid in putdat'

c     get times-array 
      did(4)=ncdid(cdfid,'time',ierr)
      if (ierr.ne.0) print *,'*ERROR* did(4) in putdat'
      call ncdinq(cdfid,did(4),chars,ntime,ierr)
      if (ierr.ne.0) print *,'*ERROR* in ncdinq in putdat'
      idtime=ncvid(cdfid,'time',ierr)
      if (ierr.ne.0) print *,'*ERROR* in ncvid in putdat'
C     Check if a new time step is starting
      iflag=0
      do i=1,ntime
        call ncvgt1(cdfid,idtime,i,timeval,ierr)
        if (ierr.ne.0) print *,'*ERROR* in ncvgt1 in putdat'
        if (time.eq.timeval) iflag=i
      enddo
      if (iflag.eq.0) then		! new time step
        ntime=ntime+1
        iflag=ntime
        idtime=ncvid(cdfid,'time',ierr)
        if (ierr.ne.0) print *, '*ERROR* in ncvid in putdat'
        call ncvpt1(cdfid,idtime,ntime,time,ierr)
        if (ierr.ne.0) print *, '*ERROR* in ncvpt1 in putdat'
      endif

C     Define data volume to write on the NetCDF file in index space
      corner(1)=1               ! starting corner of data volume
      corner(2)=1
      edgeln(1)=vardim(1)       ! edge lengthes of data volume
      edgeln(2)=vardim(2)
      if (level.eq.0) then
        corner(3)=1
        edgeln(3)=vardim(3)
      else
        corner(3)=level
        edgeln(3)=1
      endif
      corner(4)=iflag
      edgeln(4)=1
      
C     Put data on NetCDF file

c      print *,'vor Aufruf ncvpt d.h. Daten schreiben in putdat '
c      print *,'cdfid ',cdfid
c      print *,'idvar ',idvar
c      print *,'corner ',corner
c      print *,'edgeln ',edgeln

      call ncvpt(cdfid,idvar,corner,edgeln,dat,error)
      if (error.ne.0) then
        print *, '*ERROR* in ncvpt in putdat - Put data on NetCDF file'
      endif

C     Synchronize output to disk and close the files

      call ncsnc(cdfid,ierr)
      if (ierr.ne.0) print *, '*ERROR* in ncsnc in putdat'
      end






      subroutine putdef (cdfid, varnam, ndim, misdat, 
     &                            vardim, varmin, varmax, stag, error)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine is called to define the dimensions and the
c        attributes of a variable on an IVE-NetCDF file for use with the
c        IVE plotting package. Prior to calling this routine, the file must
c        be opened with a call to opncdf (extend an existing file) or
c        crecdf (create a new file).
c     Arguments:
c        cdfid   int   input   file-identifier
c                              (can be obtained by calling routine 
c                              opncdf)
c        varnam  char  input   the user-supplied variable name.
c        ndim    int   input   the number of dimensions (ndim<=4). 
c                              Upon ndim=4, the fourth dimension of the
c                              variable is specified as 'unlimited'
c                              on the file (time-dimension). It can 
c                              later be extended to arbitrary length.
c        misdat  real  input   missing data value for the variable. 
c        vardim  int   input   the dimensions of the variable.
c                              Is dimensioned at least Min(3,ndim). 
c        varmin  real  input   the location in physical space of the
c                              origin of each variable.
c                              Is dimensioned at least Min(3,ndim). 
c        varmax  real  input   the extent of each variable in physical
c                              space.
c                              Is dimensioned at least Min(ndim). 
c        stag    real  input   the grid staggering for each variable.
c                              Is dimensioned at least Min(3,ndim). 
c        error   int   output  indicates possible errors found in this
c                              routine.
c                              error = 0   no errors detected.
c                              error =10   other errors detected.
c     History:
c       Apr. 93    Christoph Schaer (ETHZ)     Created.
c-----------------------------------------------------------------------

      include "netcdf.inc"

c     Argument declarations.
      integer        MAXDIM
      parameter      (MAXDIM=4)
      character *(*) varnam
      integer        vardim(*), ndim, error, cdfid
      real           misdat,  stag(*), varmin(*), varmax(*)

c     Local variable declarations.
      character *(20) dimnam,attnam,dimchk
      character *(1)  chrid(MAXDIM)
      character *(20) dimnams(MAXNCDIM)
      integer         dimvals(MAXNCDIM)
      integer         numdims,numvars,numgats,dimulim
      integer         id,did(MAXDIM),idtime,i,k,ierr
      integer         ncopts
      integer         ibeg,iend
      data            chrid/'x','y','z','t'/

c     Get current value of error options.
      call ncgopt (ncopts)

c     make sure NetCDF-errors do not abort execution
      call ncpopt(NCVERBOS)

c     Initially set error to indicate no errors.
      error = 0

c     Make sure ndim <= MAXDIM.
      if (ndim.gt.MAXDIM) then
         error = 10
         go to 900
      endif

c     Read existing dimensions-declarations from the file
      call ncinq(cdfid,numdims,numvars,numgats,dimulim,error)
      if (numdims.gt.0) then
        do i=1,numdims
          call ncdinq(cdfid,i,dimnams(i),dimvals(i),error)
c         print *,dimnams(i),dimvals(i)
        enddo
      endif

c     put file into define mode
      call ncredf(cdfid,error)
      if (error.ne.0) goto 920

c     define spatial dimensions
      do k=1,min0(3,ndim)
c       define the default dimension-name
        dimnam(1:3)='dim'
        dimnam(4:4)=chrid(k)
        dimnam(5:5)='_'
        dimnam(6:5+len_trim(varnam))=trim(varnam)
        dimnam=dimnam(1:5+len_trim(varnam))
        did(k)=-1
        if (numdims.gt.0) then
c         check if an existing dimension-declaration can be used
c         instead of defining a new dimension
          do i=1,numdims
            dimchk=dimnams(i)
            if ((vardim(k).eq.dimvals(i)).and.
     &        (dimnam(1:4).eq.dimchk(1:4))) then 
              did(k)=i
              goto 100
            endif
          enddo
 100      continue
        endif
        if (did(k).lt.0) then
c         define the dimension
          did(k)=ncddef(cdfid,dimnam,vardim(k),error)
          if (error.ne.0) goto 920
        endif
      enddo

c     define the times-array
      if (ndim.eq.4) then
c       define dimension and variable 'time'
        if (numdims.ge.4) then
          did(4)=ncdid(cdfid,'time',ierr)
          idtime=ncvid(cdfid,'time',ierr)
        else
c         this dimension must first be defined
          did(4) = ncddef (cdfid,'time',NCUNLIM,ierr)
          idtime = ncvdef (cdfid,'time',NCFLOAT,1,did(4),ierr)
        endif
      endif

c     define variable
      id=ncvdef(cdfid,varnam,NCFLOAT,ndim,did,error)
      if (error.ne.0) goto 920

c     define attributes
      do k=1,min0(ndim,3)
c       min postion
        attnam(1:1)=chrid(k)
        attnam(2:4)='min'
        attnam=attnam(1:4)
        call ncapt(cdfid,id,attnam,NCFLOAT,1,varmin(k),error)
        if (error.gt.0) goto 920
c       max position     
        attnam(1:1)=chrid(k)
        attnam(2:4)='max'
        attnam=attnam(1:4)
        call ncapt(cdfid,id,attnam,NCFLOAT,1,varmax(k),error)
        if (error.gt.0) goto 920     
c       staggering
        attnam(1:1)=chrid(k)
        attnam(2:5)='stag'
        attnam=attnam(1:5)
        call ncapt(cdfid,id,attnam,NCFLOAT,1,stag(k),error)
        if (error.gt.0) goto 920
      enddo

c     define missing data value
      call ncapt(cdfid,id,'missing_data',NCFLOAT,1,misdat,error)
      if (error.gt.0) goto 920     

c     leave define mode
      call ncendf(cdfid,error)
      if (error.gt.0) goto 920     

c     synchronyse output to disk and exit
      call ncsnc (cdfid,error)
      call ncpopt (ncopts)
      return

c     Error exits.
 900  write (6, *) '*ERROR*: When calling putcdf, the number of ',
     &             'variable dimensions must be less or equal 4.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      return

 920  write (6, *) '*ERROR*: An error occurred while attempting to ',
     &             'write the data file in subroutine putcdf.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      return
      end


      subroutine puttimes(cdfid,times,ntimes,ierr)
C------------------------------------------------------------------------
C     Purpose:
C        Redefine all times on the specified NetCDF file
C     Arguments: 
C        cdfid  int  input   identifier for NetCDF file
C        times	real input   array contains all time values on the file,
C                            dimensioned at least times(ntimes)
C        ntimes int  input   number of times on the file
C        error  int  output  errorflag 
C     History:
C        Heini Wernli, ETHZ  
C        Christoph Schaer, ETHZ
C        Johannes Jenkner, ETHZ (adjustment for integer and double times)
C     Note:
C        This preliminary version does not define the times-array, but only
C        overwrites or extends an existing times-array.
C------------------------------------------------------------------------

      integer	ierr,i
      real times(*)
      integer didtim,ntimes

      integer	cdfid,idtime,nfiltim
      integer	ncdid,ncvid

      integer   vtyp,dn,nat
      integer dims(4)

      idtime=ncvid(cdfid,'time',ierr)   ! inquire id for time array
      if (ierr.ne.0) return
      didtim=ncdid(cdfid,'time',ierr)	! inquire id for time dimension
      if (ierr.ne.0) return

      call ncdinq(cdfid,didtim,'time',nfiltim,ierr)   ! inquire # of times
      if (ierr.ne.0) return

      call ncvinq(cdfid,idtime,dimnam,vtyp,dn,dims,natt,ierr)

      if (nfiltim.lt.ntimes) then
        print *,'Warning: puttimes is extending times-array'
      else if (nfiltim.gt.ntimes) then
        print *,'Warning: puttimes does not cover range of times-array'
      endif

      if (vtyp.eq.5) then
        do i=1,ntimes
          call ncvpt1(cdfid,idtime,i,times(i),ierr)
          if (ierr.ne.0) return
        enddo
      elseif (vtyp.eq.4) then ! integer version  
        do i=1,ntimes
          call ncvpt1(cdfid,idtime,i,int(times(i)),ierr) 
          if (ierr.ne.0) return
        enddo
      elseif (vtyp.eq.6) then ! double precision version
        do i=1,ntimes
          call ncvgt1(cdfid,idtime,i,dble(times(i)),ierr) 
          if (ierr.ne.0) return
        enddo        
      else
        return
      endif

      end



      subroutine gettimes(cdfid,times,ntimes,ierr)
C------------------------------------------------------------------------
C     Purpose:
C        Get all times on the specified NetCDF file
C     Arguments: 
C        cdfid  int  input   identifier for NetCDF file
C        times	real output  array contains all time values on the file,
C                            dimensioned at least times(ntimes)
C        ntimes int  output  number of times on the file
C        error  int  output  errorflag 
C     History:
C        Heini Wernli, ETHZ  
C        Johannes Jenkner, ETHZ (adjustment for integer and double times)
C------------------------------------------------------------------------

      include "netcdf.inc"

      integer	ierr,i
      real times(*)
      integer didtim,ntimes

      integer	cdfid,idtime
      integer	ncopts
      character*(20) dimnam

      integer   vtyp,dn,nat
      integer dims(4)

      integer,dimension(:),allocatable :: inttimes
      double precision,dimension(:),allocatable :: doubletimes

c     Get current value of error options, and make sure netCDF-errors do 
c     not abort execution
      call ncgopt (ncopts)
      call ncpopt(NCVERBOS)

      didtim=ncdid(cdfid,'time',ierr)	! inquire id for time dimension
      if (ierr.ne.0) goto 900
      idtime=ncvid(cdfid,'time',ierr)   ! inquire id for time array
      if (ierr.ne.0) goto 900
      call ncdinq(cdfid,didtim,dimnam,ntimes,ierr)      ! inquire # of times
      if (ierr.ne.0) goto 900
  
      call ncvinq(cdfid,idtime,dimnam,vtyp,dn,dims,natt,ierr)

      if (vtyp.eq.5) then
        do i=1,ntimes
          call ncvgt1(cdfid,idtime,i,times(i),ierr) ! get times
          if (ierr.ne.0) goto 900
        enddo   
      elseif (vtyp.eq.4) then ! integer version  
        allocate(inttimes(ntimes))
        do i=1,ntimes
          call ncvgt1(cdfid,idtime,i,inttimes(i),ierr) ! get times
          if (ierr.ne.0) goto 900
        enddo
        times(1:ntimes)=real(inttimes(1:ntimes))
      elseif (vtyp.eq.6) then ! double precision version
        allocate(doubletimes(ntimes))
        do i=1,ntimes
          call ncvgt1(cdfid,idtime,i,doubletimes(i),ierr) ! get times
          if (ierr.ne.0) goto 900
        enddo
        times(1:ntimes)=real(doubletimes(1:ntimes))        
      else
        goto 900
      endif

c     normal exit
      call ncpopt (ncopts)
      return

c     error exit
 900  ntimes=1
      times(1)=0.
      call ncpopt (ncopts)
      end




      subroutine cpp_crecdf(filnam,filnam_len,cdfid,phymin,phymax,ndim,
     &     cfn,cfn_len,error)
C------------------------------------------------------------------------
C     Purpose:
C        allows to call crecdf from c++
C     Arguments: 
C        see crecdf
C        additionally: fname_len and cfn_len, the length of the 
C           strings
C        
C        
C     History:
C        981221  Mark A. Liniger ETHZ
C        
C     Note:
C        
C        
C------------------------------------------------------------------------
      integer        filnam_len,ndim,cfn_len,error,cdfid
      character *(*) filnam,cfn
      real           phymin(*),phymax(*)

      call crecdf (filnam(1:filnam_len),cdfid,phymin,phymax,ndim,
     &     cfn(1:cfn_len),error)

      end


      subroutine cpp_putdef(cdfid,varnam,varnam_len,ndim,misdat,
     &     vardim,varmin,varmax,stag,error)
C------------------------------------------------------------------------
C     Purpose:
C        allows to call putdef from c++
C     Arguments: 
C        see crecdf
C        additionally: varnam_len, the length of the 
C           strings
C        
C        
C     History:
C        981221  Mark A. Liniger ETHZ
C        
C     Note:
C        
C        
C------------------------------------------------------------------------
      integer        varnam_len,ndim,error,vardim(*),cdfid
      character *(*) varnam
      real           misdat,stag(*),varmin(*),varmax(*)

      call putdef (cdfid, varnam(1:varnam_len), ndim, misdat, 
     &     vardim, varmin, varmax, stag, error)

      end


      subroutine cpp_putdat(cdfid, varnam,varnam_len, 
     &     time, level, dat, error)
C------------------------------------------------------------------------
C     Purpose:
C        allows to call putdef from c++
C     Arguments: 
C        see crecdf
C        additionally: varnam_len, the length of the 
C           strings
C        
C        
C     History:
C        981221  Mark A. Liniger ETHZ
C        
C     Note:
C        
C        
C------------------------------------------------------------------------
      integer        varnam_len,cdfid,error,level
      character *(*) varnam
      real           dat(*)
      real           time

      call putdat(cdfid, varnam(1:varnam_len), time, level, dat, error)



      end
