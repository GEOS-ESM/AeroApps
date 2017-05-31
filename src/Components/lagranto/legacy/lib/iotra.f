c     ****************************************************************
c     * This package provides IO routines for trajectories. A file   *
c     * is characterised by the filename <filename> and the file     *
c     * identifier <fid>. Different modes <mode> are supported:      *
c     *     mode=1: ascii, sorted by trajectory;                     *
c     *     mode=2: ascii, sorted by time;                           *
c     *     mode=3: fortran (unformatted)                            * 
c     *     mode=4: IVE netcdf (for compatibiltzy reasons)           *
c     * A trajectory set is given as 3d array <tra(ntra,ntim,ncol)>  *
c     * where <ntra> is the number of trajectories, <ntim> the       *
c     * number of times of each trajectory and <ncol> the number of  *
c     * columns of the trajectory. The first 4 columns are: time,    *
c     * longitude, latitude, pressure. The other columns are traced  *
c     * fields. The complete list of all columns is given in the     *
c     * array <vars(ncol)>. Finally, the reference date is given in  *
c     * the array <time(6)=year,month,day,hour,time length of the    *
c     * trajectory (hour,min)>.                                      *
c     *                                                              *
c     * Author: Michael Sprenger, September 2008                     *
c     ****************************************************************

c     ----------------------------------------------------------------
c     Open a trajectory file for reading
c     ----------------------------------------------------------------
      
      subroutine ropen_tra(fid,filename,ntra,ntim,ncol,time,vars,mode)

      implicit none
      
c     Declaration of subroutine parameters
      integer      fid
      character*80 filename
      integer      mode
      integer      ntra,ntim,ncol
      integer      time(6)
      character*80 vars(ncol)

c     Auxiliary variables
      integer      vardim(4)
      real         varmin(4),varmax(4),stag(4)
      real         mdv
      character*80 cfn
      integer      ierr
      integer      i
      integer      nvars

c     Open file
      if (mode.eq.1) then
         fid = 10
         open(fid,file=filename)
      elseif (mode.eq.2) then
         fid = 10
         open(fid,file=filename)
      elseif (mode.eq.3) then
         open(fid,file=filename,form='unformatted')
      elseif (mode.eq.4) then
         call cdfopn(filename,fid,ierr)
      elseif (mode.eq.5) then
         print*,' ERROR: Reading KML not supported'
         stop
      endif

c     Read header information
      call read_hea(fid,time,vars,ntra,ntim,ncol,mode)

      end
      

c     ----------------------------------------------------------------
c     Open a trajectory file for wrinting
c     ----------------------------------------------------------------
      
      subroutine wopen_tra(fid,filename,ntra,ntim,ncol,time,vars,mode)
      
      implicit none
      
c     Declaration of subroutine parameters
      integer      fid
      character*80 filename
      integer      mode
      integer      ntra,ntim,ncol
      integer      time(6) 
      character*80 vars(ncol)

c     Auxiliary variables
      integer      vardim(4)
      real         varmin(4),varmax(4),stag(4)
      real         mdv
      character*80 cfn
      integer      ierr
      integer      i
      character*80 varname
      real         rtime(6)

c     Open file
      if (mode.eq.1) then
         fid = 10
         open(fid,file=filename)
      elseif (mode.eq.2) then
         fid = 10
         open(fid,file=filename)
      elseif (mode.eq.3) then
         open(fid,file=filename,form='unformatted')
      elseif (mode.eq.4) then
         vardim(1)=ntra
         vardim(2)=1
         vardim(3)=1
         vardim(4)=1
         cfn      =trim(filename)//'_cst'
         mdv      =-999.98999
         call crecdf(filename,fid,varmin,varmax,3,cfn,ierr)
      elseif (mode.eq.5) then
         fid = 10
         open(fid,file=filename)
      endif

c     Write header information
      call write_hea(fid,time,vars,ntra,ntim,ncol,mode)

      end


c     ----------------------------------------------------------------
c     Read a trajectory
c     ----------------------------------------------------------------

      subroutine read_tra(fid,tra,ntra,ntim,ncol,mode)

      implicit none

c     Declaration of subroutine parameters
      integer   fid
      integer   ntim
      integer   ncol
      integer   ntra
      real      tra(ntra,ntim,ncol)
      integer   mode

c     Auxiliary variables
      integer      i,j,n
      real         arr(ntra)
      integer      ntimes
      real         times(10000)
      integer      ierr
      character*80 vars(ncol+2)
      integer      nvars

c     Read ascii mode, sorted by trajectory (mode=1)
      if (mode.eq.1) then
         read(fid,*,end=100)
         do n=1,ntra
            do i=1,ntim
               read(fid,*,end=110) (tra(n,i,j),j=1,ncol)
            enddo
         enddo

c     Read ascii mode, sorted by time (mode=2)
      elseif (mode.eq.2) then
         read(fid,*,end=100)
         do i=1,ntim
            do n=1,ntra
               read(fid,*,end=100) (tra(n,i,j),j=1,ncol)
            enddo
         enddo

c     Read fortran mode (mode=3)
      elseif (mode.eq.3) then
         read(fid) tra

c     Read IVE netcdf mode (mode=4)
      elseif (mode.eq.4) then
          call gettimes(fid,times,ntimes,ierr)
          call getvars(fid,nvars,vars,ierr)
          do i=1,ntim
             do j=1,ncol
                if (j.eq.1) then
                   do n=1,ntra
                      tra(n,i,1)=times(i)
                   enddo
                else
                   call getdat(fid,vars(j),times(i),0,arr,ierr)
                   do n=1,ntra
                      tra(n,i,j)=arr(n)
                   enddo
                endif
            enddo
         enddo

      endif

      return

c     End of file has been reached: set negative <fid>
 100  fid=-fid
      return

c     Error: incomplete trajectory
 110  print*,'<read_tra>: Incomplete trajectory... Stop'
      stop
      
      end


c     ----------------------------------------------------------------
c     Write a trajectory
c     ----------------------------------------------------------------

      subroutine write_tra(fid,tra,ntra,ntim,ncol,mode)

      implicit none

c     Declaration of subroutine parameters
      integer   fid
      integer   ntim
      integer   ncol
      integer   ntra
      real      tra(ntra,ntim,ncol)
      integer   mode

c     Auxiliary variables
      integer      i,j,n
      real         arr(ntra)
      integer      ierr
      real         time
      character*80 vars(ncol+2)
      integer      nvars
      character*20 lonstr,latstr,levstr
      character*80 outstr
      real         ref_z(3000),ref_p(3000),ref_t(3000)
      real         lev
      character*80 path

c     Write ascii mode, sorted by trajectory (mode=1)
      if (mode.eq.1) then
         do n=1,ntra
            write(fid,*)
            do i=1,ntim

c              Avoid ugly *s or missing space in output
               do j=5,ncol
                  if ( abs(tra(n,i,j)).gt.9999.) then                     
                     print*,'Format problem : ',tra(n,i,j),' -> -999.99'
                     tra(n,i,j) = -999.99
                  endif
               enddo

               write(fid,'(1f7.2,f9.2,f8.2,i6,100f10.3)') 
     >               (tra(n,i,j),j=1,3),             ! time, lon, lat
     >               nint(tra(n,i,4)),               ! p
     >               (tra(n,i,j),j=5,ncol)           ! fields
            enddo
         enddo

c     Write ascii mode, sorted by time (mode=2)
      elseif (mode.eq.2) then
         do i=1,ntim
            write(fid,*)
            do n=1,ntra

c              Avoid ugly *s or missing space in output
               do j=5,ncol
                  if ( abs(tra(n,i,j)).gt.9999.) then                     
                     print*,'Format problem : ',tra(n,i,j),' -> -999.99'
                     tra(n,i,j) = -999.99
                  endif
               enddo

               write(fid,'(1f7.2,f9.2,f8.2,i6,100f10.3)') 
     >               (tra(n,i,j),j=1,3),             ! time, lon, lat
     >               nint(tra(n,i,4)),               ! p
     >               (tra(n,i,j),j=5,ncol)           ! fields
            enddo
         enddo

c     Write fortran mode (mode=3)
      elseif (mode.eq.3) then
         write(fid) tra                              

c     Write netcdf mode (mode=4)
      elseif (mode.eq.4) then
         call getvars(fid,nvars,vars,ierr)
         do i=1,ntim
            time=tra(1,i,1)
            do j=2,ncol
               do n=1,ntra
                  arr(n)=tra(n,i,j)
               enddo
               call putdat(fid,vars(j),time,0,arr,ierr)
            enddo
         enddo

c     Write KML mode (mode=5)
      elseif (mode.eq.5) then

         call getenv('DYN_TOOLS',path)
         path = trim(path)//'/lagranto.ecmwf/goodies/'

         open(fid+1,file=trim(path)//'reformat.refprof')
           
           do n=1,6
              read(fid+1,*)
           enddo
           do n=1,3000
              read(fid+1,*)  ref_z(n),ref_t(n),ref_p(n)
              ref_p(n) = 0.01 * ref_p(n)
           enddo

         close(fid+1)

         do n=1,ntra
           write(fid,"(A)") '<Placemark>'
           write(fid,"(A)") '<name>Absolute Extruded</name>'
           write(fid,"(A)") '<styleUrl>#yellowkLineGreenPoly</styleUrl>'
           write(fid,"(A)") '<LineString>'
           write(fid,"(A)") '<extrude>1</extrude>'
           write(fid,"(A)") '<tessellate>1</tessellate>'
           write(fid,"(A)") '<altitudeMode>absolute</altitudeMode>'
           write(fid,"(A)") '<coordinates>'

           do i=1,ntim
             write(lonstr,*) tra(n,i,2)
             write(latstr,*) tra(n,i,3)
             
             call binary(lev,tra(n,i,4),ref_z,ref_p)
             write(levstr,*) lev

             outstr = trim(adjustl(lonstr))//','//
     >                trim(adjustl(latstr))//','//
     >                trim(adjustl(levstr))

             write(fid,"(A)") outstr

           enddo

           write(fid,*) '</coordinates>'
           write(fid,*) '</LineString>'
           write(fid,*) '</Placemark>'
         enddo

      endif

      end


c     ----------------------------------------------------------------
c     Read header from trajectory file
c     ----------------------------------------------------------------

      subroutine read_hea(fid,time,vars,ntra,ntim,ncol,mode)

      implicit none
      
c     Declaration of subroutine parameters
      integer       fid
      integer       time(6)
      integer       ntra,ntim,ncol
      character*80  vars(ncol)
      integer       mode

c     Auxiliary variables
      integer       i
      character     ch(500)
      character*500 str
      integer       ich(500)
      integer       isstr,ileft,iright
      character*80  varname
      real          rtime(6)
      integer       ierr
      integer       nvars
      character*15  str1
      character     str2
      character*13  str3
      character*4   str4
      character*80  linestr
      integer       itmp1,itmp2

c     Read ascii format (mode=1,2)
      if ( (mode.eq.1).or.(mode.eq.2) ) then

c        Read the time specification (old and new format)
         read(fid,'(a80)') linestr
         
         if ( linestr(1:14).eq.'Reference date' ) then
            read(linestr,'(a15,i4,i2,i2,a1,i2,i2,a13,i8,a4)') 
     >           str1,
     >           time(1),time(2),time(3),str2,time(4),time(5),
     >           str3,time(6),str4
                       
         elseif ( linestr(1:11).eq.'time period' ) then
            read(linestr,'(a12,i4,i2,i2,a1,i2,a4,i6,a1,i2,a2)') 
     >           str1,
     >           time(1),time(2),time(3),str2,time(4),
     >           str3,itmp1,str3,itmp2,str4
            time(5) = 0
            time(6) = itmp1 * 60 + itmp2

         endif

c        Skip the empty line and read field names
         read(fid,*)
         read(fid,'(a500)',end=100) str
         do i=1,500
            ch(i)=str(i:i)
         enddo

c        Split the input string
         isstr=0
         nvars=0
         do i=1,500
            if ( (isstr.eq.0).and.(ch(i).ne.' ') ) then
               isstr=1
               ileft=i
            elseif ( (isstr.eq.1).and.(ch(i).eq.' ') ) then
               isstr=0
               iright=i-1
               nvars=nvars+1
               vars(nvars)=str(ileft:iright)
            endif
         enddo

c        Skip the empty line
         read(fid,*,end=100)

c     Read fortran mode (mode=3)
      elseif (mode.eq.3) then
         read(fid) ntra,ntim,ncol
         read(fid) time
         read(fid) vars

c     Read IVE netcdf mode (mode=4)
      elseif (mode.eq.4) then
         call getvars(fid,nvars,vars,ierr)  
         varname='BASEDATE'
         call getdat(fid,varname,0.,0,rtime,ierr)
         do i=1,6
            time(i)=nint(rtime(i))
         enddo
      endif

      return

c     End of file has been reached
 100  fid=-fid
      return

c     Excetion handling
 110  print*,'<read_hea>: Unexspected time format.... Stop'
      stop
   
      end
      

c     ----------------------------------------------------------------
c     Write header to trajectory file (in ascii mode)
c     ----------------------------------------------------------------

      subroutine write_hea(fid,time,vars,ntra,ntim,ncol,mode)

      implicit none
      
c     Declaration of subroutine parameters
      integer       fid
      integer       time(6)
      integer       ntra,ntim,ncol
      character*80  vars(ncol)
      integer       mode

c     Auxiliary variables
      integer       i
      character*500 str
      character*4   str1
      character*2   str2,str3,str4,str5,str6
      integer       vardim(4)
      real          varmin(4),varmax(4),stag(4)
      real          mdv
      integer       ierr
      character*80  varname
      real          rtime(6)
      integer       nvars

c     Write ascii format (mode=1,2)
      if ( (mode.eq.1).or.(mode.eq.2) ) then

c        Get the strings for output
         write(str1,'(i4)') time(1)
         write(str2,'(i2)') time(2)
         write(str3,'(i2)') time(3)
         write(str4,'(i2)') time(4)
         write(str5,'(i2)') time(5)
         if (time(2).eq. 0) str2(1:1)='0'
         if (time(3).eq. 0) str3(1:1)='0'
         if (time(4).eq. 0) str4(1:1)='0'
         if (time(5).eq. 0) str5(1:1)='0'
         if (time(2).lt.10) str2(1:1)='0'
         if (time(3).lt.10) str3(1:1)='0'
         if (time(4).lt.10) str4(1:1)='0'
         if (time(5).lt.10) str5(1:1)='0'

c        Write the time specification
         write(fid,'(a15,a4,a2,a2,a1,a2,a2,a13,i8,a4)') 
     >          'Reference date ',
     >           str1,str2,str3,'_',str4,str5,
     >          ' / Time range',time(6), ' min'
         write(fid,*)

c        Write variable names
         str=''
         do i=1,ncol
            varname = vars(i)
            vars(i) = varname(1:9)
            str     = trim(str)//trim(vars(i))
         enddo
         write(fid,'(a6,a9,a8,a6,100a10)') (trim(vars(i)),i=1,ncol)
         write(fid,'(a6,a9,a8,a6,100a10)') 
     >              '------','---------','--------','------',
     >              ('----------',i=5,ncol)

c     Write fortran mode (mode=3)
      elseif (mode.eq.3) then
         write(fid) ntra,ntim,ncol
         write(fid) time
         write(fid) vars

c     Write IVE netcdf format (mode=4)
      elseif (mode.eq.4) then
         vardim(1)=ntra
         vardim(2)=1
         vardim(3)=1
         vardim(4)=1
         mdv      =-999.98999
         do i=2,ncol
            call putdef(fid,vars(i),4,mdv,vardim,
     >                  varmin,varmax,stag,ierr)
         enddo
         varname='BASEDATE'
         vardim(1)=6
         call putdef(fid,varname,4,mdv,vardim,
     >               varmin,varmax,stag,ierr)
         do i=1,6
            rtime(i)=real(time(i))
         enddo
         call putdat(fid,varname,0.,0,rtime,ierr)

c     Write KML format (mode=5)
      elseif (mode.eq.5) then

      write(fid,"(A)") '<?xml version="1.0" encoding="UTF-8"?>'
      write(fid,"(A)") '<kml xmlns="http://www.opengis.net/kml/2.2">'
      write(fid,"(A)") '<Document>'
      write(fid,"(A)") '<name>Paths</name>'
      write(fid,"(A)") '<Style id="yellowLineGreenPoly">'
      write(fid,"(A)") '<LineStyle>'
c      write(fid,*) '<color>7f00ffff</color>'    ! Yellow
      write(fid,"(A)") '<color>500A0A0A</color>'     ! Black
      write(fid,"(A)") '<width>4</width>'
      write(fid,"(A)") '</LineStyle>'
      write(fid,"(A)") '<PolyStyle>'
      write(fid,"(A)") '<color>7f00ff00</color>'
      write(fid,"(A)") '</PolyStyle>'
      write(fid,"(A)") '</Style>'

      endif

      end
      
      
c     ----------------------------------------------------------------
c     Close a trajectory file
c     ----------------------------------------------------------------
      
      subroutine close_tra(fid,mode)

      implicit none
      
c     Declaration of subroutine parameters
      integer      fid
      integer      mode
      
c     Auxiliary variables
      integer      ierr

c     Close file
      if (mode.eq.1) then
         close(abs(fid))

      elseif (mode.eq.2) then
         close(abs(fid))

      elseif (mode.eq.3) then
         close(fid)

      elseif (mode.eq.4) then
         call clscdf(fid,ierr)

      elseif (mode.eq.5) then
          write(fid,"(A)") '</Document>'
          write(fid,"(A)") '</kml>'
         close(abs(fid))
      endif

      end
      

c     ----------------------------------------------------------------
c     Determine the mode of a trajectory file
c     ----------------------------------------------------------------

      subroutine mode_tra(mode,filename)
      
      implicit none

c     Declaration of subroutine parameters
      integer        mode
      character*80   filename

c     Auxiliary variables
      integer        len
      character      char0,char1,char2,char3,char4

c     Get mode
      mode=-1
      
      len  = len_trim(filename)

c     Mode specified by number
      char0 = filename((len-1):(len-1))
      char1 = filename(len:len)

      if ( (char0.eq.'.').and.(char1.eq.'1') ) mode=1
      if ( (char0.eq.'.').and.(char1.eq.'2') ) mode=2
      if ( (char0.eq.'.').and.(char1.eq.'3') ) mode=3
      if ( (char0.eq.'.').and.(char1.eq.'4') ) mode=4
      if ( (char0.eq.'.').and.(char1.eq.'5') ) mode=5

      if ( mode.gt.0 ) return

c     Mode specified by appendix
      char0 = filename((len-3):(len-3))
      char1 = filename((len-2):(len-2))
      char2 = filename((len-1):(len-1))
      char3 = filename(len:len)
      if ( (char1.eq.'.').and.(char2.eq.'l').and.(char3.eq.'s') ) mode=1
      if ( (char1.eq.'.').and.(char2.eq.'t').and.(char3.eq.'i') ) mode=2
      if ( (char1.eq.'.').and.(char2.eq.'d').and.(char3.eq.'u') ) mode=3
      if ( (char1.eq.'.').and.(char2.eq.'n').and.(char3.eq.'c') ) mode=4

      if ( (char0.eq.'.').and.(char1.eq.'k').and.
     >                        (char2.eq.'m').and.
     >                        (char3.eq.'l') ) mode = 5

      end


c     ----------------------------------------------------------------
c     Get dimension of a trajectory file
c     ----------------------------------------------------------------
    
      subroutine info_tra(filename,ntra,ntim,ncol,mode)

      implicit none
      
c     Declaration of subroutine parameters
      integer      fid
      character*80 filename
      integer      mode
      integer      ntra,ntim,ncol

c     Auxiliary variables
      integer       vardim(4)
      real          varmin(4),varmax(4),stag(4)
      real          mdv
      character*80  cfn
      integer       ierr
      integer       i,ndim
      integer       nvars
      integer       ntimes
      real          times(100)
      character*500 str
      integer       nline0,nline1,nline2
      integer       isstr,isok
      character     ch
      character*80  vars(200)
      integer       ios

c     Open file
      if (mode.eq.1) then
         fid=10
         open(fid,file=filename)
      elseif (mode.eq.2) then
         fid=10
         open(fid,file=filename)
      elseif (mode.eq.3) then
         fid=10
         open(fid,file=filename,form='unformatted')
      elseif (mode.eq.4) then
         call cdfopn(filename,fid,ierr)
      endif

c     Get dimension information
      if ( (mode.eq.1).or.(mode.eq.2) ) then
         read(fid,*)
         read(fid,*)
         read(fid,'(a500)') str
         read(fid,*)

c        Get the number of columns
         isstr=0
         ncol =0
         do i=1,500
            ch = str(i:i)
            if ( (isstr.eq.0).and.(ch.ne.' ') ) then
               isstr=1
            elseif ( (isstr.eq.1).and.(ch.eq.' ') ) then
               isstr=0
               ncol=ncol+1
            endif
         enddo

c        Init the dimensions for empty trajectory files
         ntra = 0
         ntim = 0

c        Get the first data block
         nline0  = 5
         nline1  = 5
         read(fid,*,end=140,iostat=ios)
 100     read(fid,'(a500)',end=110,iostat=ios) str
         if (str.ne.'') then
            nline1 = nline1 + 1
            goto 100
         endif
 110     continue
         
c        Get the total numbers of lines in the data block
         nline2 = nline1
         if ( ios.eq.0 ) then
 120        read(fid,*,end=130)
            nline2 = nline2 + 1
            goto 120
 130        nline2 = nline2 + 1
         endif

c        Set the dimensions
         if (mode.eq.1) then
            ntim = nline1 - nline0
            ntra = (nline2-nline0+1)/(ntim+1)
         else
            ntra = nline1 - nline0
            ntim = (nline2-nline0+1)/(ntra+1)
         endif

      elseif (mode.eq.3) then
         read(fid) ntra,ntim,ncol

      elseif (mode.eq.4) then
         call gettimes(fid,times,ntimes,ierr)
         call getvars(fid,nvars,vars,ierr)
         call getdef(fid,trim(vars(2)),ndim,mdv,vardim,
     >               varmin,varmax,stag,ierr)
         ntra = vardim(1)
         ntim = ntimes
         ncol = nvars-1

      endif

c     Close file
 140  continue
      if (mode.eq.1) then
         close(fid)
      elseif (mode.eq.2) then
         close(fid)
      elseif (mode.eq.3) then
         close(fid)
      elseif (mode.eq.4) then
         call clscdf(fid,ierr)
      endif

      end
     

c     ----------------------------------------------------------------
c     Binary search algorithm
c     ----------------------------------------------------------------
 
      subroutine binary (z,p,ref_z,ref_p)

      implicit none

c     Declaration of subroutine parameters
      real   z
      real   p
      real   ref_z(3000)
      real   ref_p(3000)

c     Auxiliary variables
      integer i0,i1,im

c     Binary search
      i0 = 1
      i1 = 3000
 100  continue
        im = (i0 + i1) / 2
        if ( p.lt.ref_p(im) ) then
           i0 = im
        else 
           i1 = im
        endif
      if ( (i1-i0).gt.1 ) goto 100

c     Linear interpolation in between
      z = ref_z(i0) + ( p - ref_p(i0) ) / ( ref_p(i1) - ref_p(i0) ) *
     >                ( ref_z(i1) - ref_z(i0) ) 

      end





