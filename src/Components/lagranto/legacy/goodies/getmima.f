      PROGRAM getmima

c     ***********************************************************************
c     * Get minimum and maximum value of a varaiable on pressure            *
c     * or isentropic surface                                               *
c     * Michael Sprenger / Spring, summer 2016                              *
c     ***********************************************************************

      use netcdf

      implicit none

c     ----------------------------------------------------------------------
c     Declaration of parameters and variables
c     ----------------------------------------------------------------------
      
c     Interpolation method
      integer   ipom
      parameter (ipom=0)    
      
c     Flag for timecheck
      character*80 timecheck
      parameter    (timecheck = 'no' )        
      
c     netCDF fields      
      real,allocatable, dimension (:,:)   :: sp,varf
      real,allocatable, dimension (:,:,:) :: field,varl,tt,p3
      real,allocatable, dimension (:)     :: ak,bk
      integer                                stat
      character*80                           cdfnam,varnam
      character*80                           vnam(200)
      integer                                nvars
      integer	                             cdfid,ierr
      real	                                 mdv
      real                                   stagz
      
c     Grid description      
      real      xmin,xmax       ! Zonal grid extension
      real      ymin,ymax       ! Meridional grid extension
      integer   nx,ny,nz        ! Grid dimensions
      real      dx,dy           ! Horizontal grid resolution
      real	    pollon,pollat   ! Pole position

c     Output variables
      real      min,max                      ! Minimum & maximum value
      real      lonmin,lonmax,latmin,latmax  ! Position of min & max

c     Auxiliary variables      
      integer        iargc
      character*(80) arg
      real           rd
      integer        i,j
      integer        minindx,minindy,maxindx,maxindy
      character*80   tvar
      character*1    mode
      character*80   clev
      real	         xphys,yphys
      real	         level
      integer        flag
       
c     Externals      
      real      lmstolm,phstoph
      external  lmstolm,phstoph

c     ----------------------------------------------------------------------
c     Preparation - argument handling, grid paraemters, memor allocation
c     ----------------------------------------------------------------------

c     Check for sufficient requested arguments
      if (iargc().lt.2) then
        print*,'USAGE: getmima filename var ',
     >         '[lev in the form Pval or Tval]'
        call exit(1)
      endif
 
c     Read and transform input
      call getarg(1,arg)
      cdfnam=trim(arg)
      call getarg(2,arg)
      varnam=trim(arg)

      if (iargc().eq.3) then
        call getarg(3,arg)
        mode=arg(1:1)
        clev=arg(2:len(arg)) 
        call checkchar(clev,".",flag)
        if (flag.eq.0) clev=trim(clev)//"."
        read(clev,'(f10.2)') level
      else
        mode='X'
        level=0.
      endif

c     Init level type and minimum,maximum
      min= 1.e19
      max=-1.e19

c     Open netCDF file      
      call input_open (cdfid,cdfnam)

c     Get list of variables on file
      call input_getvars (cdfid,vnam,nvars)

C     Get infos about data domain - in particular: dimensions
      nx       = 1
      ny       = 1
      nz       = 1
      call input_grid (-cdfid,varnam,xmin,xmax,ymin,ymax,dx,dy,nx,ny,
     >               0.,pollon,pollat,rd,rd,nz,rd,rd,rd,stagz,timecheck)

C     Get memory for dynamic arrays
      allocate(sp(nx,ny),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array sp ***'
      allocate(varf(nx,ny),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array varf ***'
      allocate(varl(nx,ny,nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array varl ***'
      allocate(tt(nx,ny,nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array tt ***'
      allocate(field(nx,ny,nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array field ***'
      allocate(ak(nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array ak ***'
      allocate(bk(nz),stat=stat)
      if (stat.ne.0) print*,'*** error allocating array bk ***'
   
c     ----------------------------------------------------------------------
c     Load data
c     ----------------------------------------------------------------------

c     Get grid info for variable - non-constant part  (P, PS, AK, BK)   
      call input_grid                               
     >       (cdfid,varnam,xmin,xmax,ymin,ymax,dx,dy,nx,ny,
     >        0.,pollon,pollat,varl,sp,nz,ak,bk,stagz,timecheck)
    
c     Load Variable    
      call input_wind (cdfid,varnam,field,0.,stagz,mdv,
     >                 xmin,xmax,ymin,ymax,dx,dy,nx,ny,nz,timecheck)

c     In case of isentropic level, read temperature
      if (mode.eq.'T') then
            tvar="T"
            do i=1,nvars
              if (vnam(i).eq."THETA") tvar="THETA"
              if (vnam(i).eq."TH")    tvar="TH"
            enddo
            if (tvar.eq.'T') then
                 call input_wind (cdfid,tvar,tt,0.,stagz,mdv,xmin,xmax,
     >                            ymin,ymax,dx,dy,nx,ny,nz,timecheck)
                 call pottemp(varl,tt,sp,nx,ny,nz,ak,bk)
            else
                 call input_wind (cdfid,tvar,tt,0.,stagz,mdv,xmin,xmax,
     >                            ymin,ymax,dx,dy,nx,ny,nz,timecheck)
            endif
      endif

c     ----------------------------------------------------------------------
c     Do interpolation and get minimum, maximum
c     ----------------------------------------------------------------------

C     If required do interpolation on pressure or theta level
      if (mode.eq.'P') then
        call vipo(field,varl,level,varf,nx,ny,nz,mdv,ipom)
      elseif (mode.eq.'T') then
        call thipo(field,varl,level,varf,nx,ny,nz,mdv,ipom)
      elseif (mode.eq.'L') then
        do i=1,nx
           do j=1,ny
               varf(i,j) = field(i,j,nint(level))
           enddo
        enddo
      elseif (mode.eq.'X') then
         do i=1,nx
            do j=1,ny
                varf(i,j) = field(i,j,1)
            enddo
         enddo
      endif
 
C     Determine minimum and maximum
      do i=1,nx
        do j=1,ny

            if ((varf(i,j).ne.mdv).and.(varf(i,j).lt.min)) then
               min     = varf(i,j)
               minindx = i
               minindy = j
            endif
            if ((varf(i,j).ne.mdv).and.(varf(i,j).gt.max)) then
               max     = varf(i,j)
               maxindx = i
               maxindy = j
          endif
        enddo
      enddo
      lonmin = xmin + real(minindx-1) * dx
      latmin = ymin + real(minindy-1) * dy
      lonmax = xmin + real(maxindx-1) * dx
      latmax = ymin + real(maxindy-1) * dy

c     Rotate position to true longitude/latitude 
      if ((pollon.ne.0.).or.(pollat.ne.90.)) then
        xphys=lmstolm(latmin,lonmin,pollat,pollon)
        yphys=phstoph(latmin,lonmin,pollat,pollon)
        lonmin=xphys
        latmin=yphys
        xphys=lmstolm(latmax,lonmax,pollat,pollon)
        yphys=phstoph(latmax,lonmax,pollat,pollon)
        lonmax=xphys
        latmax=yphys
      endif
 
c     Write output
      write(*,101) min,max,lonmin,latmin,nint(level),
     >                     lonmax,latmax,nint(level)
  101 format(2f10.3,2f8.2,i6,2f8.2,i6)
 
c     Close netCDF file
      call input_close(cdfid)
 
      end
      
c     ----------------------------------------------------------------------
C     Interpolates the 3d variable var3d on the pressure surface defined
C     by lev of the variable varl. The interpolated field is
C     returned as var.
C     ipom determines the way of vertical interpolation:
C       ipom=0 is for linear interpolation
C       ipom=1 is for logarithmic interpolation
c     ----------------------------------------------------------------------

      subroutine vipo(var3d,varl,lev,var,nx,ny,nz,mdv,ipom)

      integer   nx,ny,nz,ipom
      real      lev,mdv
      real      var3d(nx,ny,nz),varl(nx,ny,nz),var(nx,ny)
 
      integer   i,j,k
      real      kind
      real      int3dm
 
      do i=1,nx
      do j=1,ny
 
        do k=1,nz-1
          if ((varl(i,j,k).ge.lev).and.(varl(i,j,k+1).le.lev)) then
            kind=float(k)+(varl(i,j,k)-lev)/
     >                   (varl(i,j,k)-varl(i,j,k+1))
            goto 100
          endif
        enddo
 100    continue
 
        var(i,j)=int3dm(var3d,nx,ny,nz,float(i),float(j),kind,mdv)
 
      enddo
      enddo
 
      return
      end
      
c     ----------------------------------------------------------------------
C     Interpolates the 3d variable var3d on the isentropic surface defined
C     by lev of the variable varl. The interpolated field is
C     returned as var.
C     ipom determines the way of vertical interpolation:
C       ipom=0 is for linear interpolation
C       ipom=1 is for logarithmic interpolation
c     ----------------------------------------------------------------------

      subroutine thipo(var3d,th3d,lev,var,nx,ny,nz,mdv,mode)
 
      integer   nx,ny,nz,mode
      real      lev,mdv
      real      var3d(nx,ny,nz),th3d(nx,ny,nz),var(nx,ny)
 
      integer   i,j,k
      real      kind
      real      int3dm
 
      do i=1,nx
      do j=1,ny
 
        do k=1,nz-1
          if ((th3d(i,j,k).le.lev).and.(th3d(i,j,k+1).ge.lev)) then
            kind=float(k)+(th3d(i,j,k)-lev)/
     >                   (th3d(i,j,k)-th3d(i,j,k+1))
            goto 100
          endif
        enddo
 100    continue
 
        var(i,j)=int3dm(var3d,nx,ny,nz,float(i),float(j),kind,mdv)
 
      enddo
      enddo
 
      return
      end
      
c     ----------------------------------------------------------------------
c     Check whether character appears within string
c     ----------------------------------------------------------------------

      subroutine checkchar(string,char,flag)
 
      character*(*)     string
      character*(1)     char
      integer           n,flag
 
      flag=0
      do n=1,len(string)
        if (string(n:n).eq.char) then
          flag=n
          return
        endif
      enddo
      end
      
c     ----------------------------------------------------------------------
c     3D interpolation
c     ----------------------------------------------------------------------
      
      real function int3d(ar,n1,n2,n3,rid,rjd,rkd)

c     Purpose:
c        This subroutine interpolates a 3d-array to an arbitrary
c        location within the grid.
c     Arguments:
c        ar        real  input   surface pressure, define as ar(n1,n2,n3)
c        n1,n2,n3  int   input   dimensions of ar
c        ri,rj,rk  real  input   grid location to be interpolated to
c     History:

c     argument declarations
      integer   n1,n2,n3
      real      ar(n1,n2,n3), rid,rjd,rkd

c     local declarations
      integer   i,j,k,ip1,jp1,kp1,ih,jh,kh
      real      frac0i,frac0j,frac0k,frac1i,frac1j,frac1k,ri,rj,rk

c     do linear interpolation
      ri=amax1(1.,amin1(float(n1),rid))
      rj=amax1(1.,amin1(float(n2),rjd))
      rk=amax1(1.,amin1(float(n3),rkd))
      ih=nint(ri)
      jh=nint(rj)
      kh=nint(rk)

c     Check for interpolation in i
*     if (abs(float(ih)-ri).lt.1.e-3) then
*       i  =ih
*       ip1=ih
*     else
        i =min0(int(ri),n1-1)
        ip1=i+1
*     endif

c     Check for interpolation in j
      if (abs(float(jh)-rj).lt.1.e-3) then
        j  =jh
        jp1=jh
      else
        j =min0(int(rj),n2-1)
        jp1=j+1
      endif

c     Check for interpolation in k
*     if (abs(float(kh)-rk).lt.1.e-3) then
*       k  =kh
*       kp1=kh
*     else
        k =min0(int(rk),n3-1)
        kp1=k+1
*     endif

      if (k.eq.kp1) then
c       no interpolation in k
        if ((i.eq.ip1).and.(j.eq.jp1)) then
c          no interpolation at all
           int3d=ar(i,j,k)
c          print *,'int3d 00: ',rid,rjd,rkd,int3d
        else
c          horizontal interpolation only
           frac0i=ri-float(i)
           frac0j=rj-float(j)
           frac1i=1.-frac0i
           frac1j=1.-frac0j
           int3d = ar(i  ,j  ,k  ) * frac1i * frac1j
     &           + ar(i  ,jp1,k  ) * frac1i * frac0j
     &           + ar(ip1,j  ,k  ) * frac0i * frac1j
     &           + ar(ip1,jp1,k  ) * frac0i * frac0j
c          print *,'int3d 10: ',rid,rjd,rkd,int3d
        endif
      else 
        frac0k=rk-float(k)
        frac1k=1.-frac0k
        if ((i.eq.ip1).and.(j.eq.jp1)) then
c          vertical interpolation only
           int3d = ar(i  ,j  ,k  ) * frac1k
     &           + ar(i  ,j  ,kp1) * frac0k
c          print *,'int3d 01: ',rid,rjd,rkd,int3d
        else
c          full 3d interpolation
           frac0i=ri-float(i)
           frac0j=rj-float(j)
           frac1i=1.-frac0i
           frac1j=1.-frac0j
           int3d = ar(i  ,j  ,k  ) * frac1i * frac1j * frac1k
     &           + ar(i  ,jp1,k  ) * frac1i * frac0j * frac1k 
     &           + ar(ip1,j  ,k  ) * frac0i * frac1j * frac1k
     &           + ar(ip1,jp1,k  ) * frac0i * frac0j * frac1k
     &           + ar(i  ,j  ,kp1) * frac1i * frac1j * frac0k
     &           + ar(i  ,jp1,kp1) * frac1i * frac0j * frac0k 
     &           + ar(ip1,j  ,kp1) * frac0i * frac1j * frac0k
     &           + ar(ip1,jp1,kp1) * frac0i * frac0j * frac0k
c          print *,'int3d 11: ',rid,rjd,rkd,int3d
        endif
      endif
      end
      
c     ----------------------------------------------------------------------
c     3D interpolation including missing data check
c     ----------------------------------------------------------------------
  
      real function int3dm(ar,n1,n2,n3,rid,rjd,rkd,misdat)

c     Purpose:
c        This subroutine interpolates a 3d-array to an arbitrary
c        location within the grid. The interpolation includes the 
c        testing of the missing data flag 'misdat'. 
c     Arguments:
c        ar        real  input   surface pressure, define as ar(n1,n2,n3)
c        n1,n2,n3  int   input   dimensions of ar
c        ri,rj,rk  real  input   grid location to be interpolated to
c        misdat    real  input   missing data flag (on if misdat<>0)
c     Warning:
c        This routine has not yet been seriously tested
c     History:

c     argument declarations
      integer   n1,n2,n3
      real      ar(n1,n2,n3), rid,rjd,rkd, misdat

c     local declarations
      integer   i,j,k,ip1,jp1,kp1,ih,jh,kh
      real      frac0i,frac0j,frac0k,frac1i,frac1j,frac1k,ri,rj,rk,int3d

c     check if routine without missing data checking can be called instead
      if (misdat.eq.0.) then
        int3dm=int3d(ar,n1,n2,n3,rid,rjd,rkd)
        return
      endif

c     do linear interpolation
      ri=amax1(1.,amin1(float(n1),rid))
      rj=amax1(1.,amin1(float(n2),rjd))
      rk=amax1(1.,amin1(float(n3),rkd))
      ih=nint(ri)
      jh=nint(rj)
      kh=nint(rk)

c     Check for interpolation in i
*     if (abs(float(ih)-ri).lt.1.e-3) then
*       i  =ih
*       ip1=ih
*     else
        i =min0(int(ri),n1-1)
        ip1=i+1
*     endif

c     Check for interpolation in j
*     if (abs(float(jh)-rj).lt.1.e-3) then
*       j  =jh
*       jp1=jh
*     else
        j =min0(int(rj),n2-1)
        jp1=j+1
*     endif

c     Check for interpolation in k
*     if (abs(float(kh)-rk).lt.1.e-3) then
*       k  =kh
*       kp1=kh
*     else
        k =min0(int(rk),n3-1)
        kp1=k+1
*     endif

      if (k.eq.kp1) then
c       no interpolation in k
        if ((i.eq.ip1).and.(j.eq.jp1)) then
c          no interpolation at all
           if (misdat.eq.ar(i,j,k)) then
             int3dm=misdat
           else
             int3dm=ar(i,j,k)
           endif
c          print *,'int3dm 00: ',rid,rjd,rkd,int3dm
        else
c          horizontal interpolation only
           if ((misdat.eq.ar(i  ,j  ,k  )).or.
     &         (misdat.eq.ar(i  ,jp1,k  )).or.
     &         (misdat.eq.ar(ip1,j  ,k  )).or.
     &         (misdat.eq.ar(ip1,jp1,k  ))) then
             int3dm=misdat
           else
             frac0i=ri-float(i)
             frac0j=rj-float(j)
             frac1i=1.-frac0i
             frac1j=1.-frac0j
             int3dm = ar(i  ,j  ,k  ) * frac1i * frac1j
     &              + ar(i  ,jp1,k  ) * frac1i * frac0j
     &              + ar(ip1,j  ,k  ) * frac0i * frac1j
     &              + ar(ip1,jp1,k  ) * frac0i * frac0j
c            print *,'int3dm 10: ',rid,rjd,rkd,int3dm
           endif
        endif
      else 
        frac0k=rk-float(k)
        frac1k=1.-frac0k
        if ((i.eq.ip1).and.(j.eq.jp1)) then
c          vertical interpolation only
           if ((misdat.eq.ar(i  ,j  ,k  )).or.
     &         (misdat.eq.ar(i  ,j  ,kp1))) then
             int3dm=misdat
           else
             int3dm = ar(i  ,j  ,k  ) * frac1k
     &              + ar(i  ,j  ,kp1) * frac0k
c            print *,'int3dm 01: ',rid,rjd,rkd,int3dm
           endif
        else
c          full 3d interpolation
           if ((misdat.eq.ar(i  ,j  ,k  )).or.
     &         (misdat.eq.ar(i  ,jp1,k  )).or.
     &         (misdat.eq.ar(ip1,j  ,k  )).or.
     &         (misdat.eq.ar(ip1,jp1,k  )).or.
     &         (misdat.eq.ar(i  ,j  ,kp1)).or.
     &         (misdat.eq.ar(i  ,jp1,kp1)).or.
     &         (misdat.eq.ar(ip1,j  ,kp1)).or.
     &         (misdat.eq.ar(ip1,jp1,kp1))) then
             int3dm=misdat
           else
             frac0i=ri-float(i)
             frac0j=rj-float(j)
             frac1i=1.-frac0i
             frac1j=1.-frac0j
             int3dm = ar(i  ,j  ,k  ) * frac1i * frac1j * frac1k
     &              + ar(i  ,jp1,k  ) * frac1i * frac0j * frac1k
     &              + ar(ip1,j  ,k  ) * frac0i * frac1j * frac1k
     &              + ar(ip1,jp1,k  ) * frac0i * frac0j * frac1k
     &              + ar(i  ,j  ,kp1) * frac1i * frac1j * frac0k
     &              + ar(i  ,jp1,kp1) * frac1i * frac0j * frac0k
     &              + ar(ip1,j  ,kp1) * frac0i * frac1j * frac0k
     &              + ar(ip1,jp1,kp1) * frac0i * frac0j * frac0k
c            print *,'int3dm 11: ',rid,rjd,rkd,int3dm
           endif
        endif
      endif
      end

c     ----------------------------------------------------------------------
c     Calculate potential temperature
c     ----------------------------------------------------------------------

      subroutine pottemp(pt,t,sp,ie,je,ke,ak,bk)
 
c     argument declaration
      integer   ie,je,ke
      real      pt(ie,je,ke),t(ie,je,ke),sp(ie,je),
     >     	ak(ke),bk(ke)
 
c     variable declaration
      integer   i,j,k
      real      rdcp,tzero
      data      rdcp,tzero /0.286,273.15/
 
c     statement-functions for the computation of pressure
      real      pp,psrf
      integer   is
      pp(is)=ak(is)+bk(is)*psrf
 
c     computation of potential temperature
      do i=1,ie
      do j=1,je
        psrf=sp(i,j)
        do k=1,ke
c     distinction of temperature in K and deg C
          if (t(i,j,k).lt.100.) then
            pt(i,j,k)=(t(i,j,k)+tzero)*( (1000./pp(k))**rdcp )
          else
            pt(i,j,k)=t(i,j,k)*( (1000./pp(k))**rdcp )
          endif
        enddo
      enddo
      enddo
      end
      
c     ----------------------------------------------------------------------
c     Calculate 3D pressure
c     ----------------------------------------------------------------------

      subroutine pres(pr,sp,ie,je,ke,ak,bk)

c     argument declaration
      integer  ie,je,ke
      real,intent(OUT) :: pr(ie,je,ke)
      real,intent(IN)  :: sp(ie,je)
      real,intent(IN)  :: ak(ke),bk(ke)

c     variable declaration
      integer  i,j,k

c     computation pressure
      do i=1,ie
      do j=1,je
      do k=1,ke
        pr(i,j,k)=ak(k)+bk(k)*sp(i,j)
      enddo
      enddo
      enddo
      end
      
c     ----------------------------------------------------------------------
c     Coordinate rotations from COSMO
c     ----------------------------------------------------------------------      
      
      REAL FUNCTION PHTOPHS (PHI, LAM, POLPHI, POLLAM)
C
C%Z% Modul %M%, V%I% vom %G%, extrahiert am %H%
C
C**** PHTOPHS  -   FC:UMRECHNUNG DER WAHREN GEOGRAPHISCHEN BREITE PHI
C****                 AUF EINEM PUNKT MIT DEN KOORDINATEN (PHIS, LAMS)
C****                 IM ROTIERTEN SYSTEM. DER NORDPOL DES SYSTEMS HAT
C****                 DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
C**   AUFRUF   :   PHI = PHTOPHS (PHI, LAM, POLPHI, POLLAM)
C**   ENTRIES  :   KEINE
C**   ZWECK    :   UMRECHNUNG DER WAHREN GEOGRAPHISCHEN BREITE PHI AUF
C**                EINEM PUNKT MIT DEN KOORDINATEN (PHIS, LAMS) IM
C**                ROTIERTEN SYSTEM. DER NORDPOL DIESES SYSTEMS HAT
C**                DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
C**   VERSIONS-
C**   DATUM    :   03.05.90
C**
C**   EXTERNALS:   KEINE
C**   EINGABE-
C**   PARAMETER:   PHI    REAL BREITE DES PUNKTES IM GEOGR. SYSTEM
C**                LAM    REAL LAENGE DES PUNKTES IM GEOGR. SYSTEM
C**                POLPHI REAL GEOGR.BREITE DES N-POLS DES ROT. SYSTEMS
C**                POLLAM REAL GEOGR.LAENGE DES N-POLS DES ROT. SYSTEMS
C**   AUSGABE-
C**   PARAMETER:   ROTIERTE BREITE PHIS ALS WERT DER FUNKTION
C**                ALLE WINKEL IN GRAD (NORDEN>0, OSTEN>0)
C**
C**   COMMON-
C**   BLOECKE  :   KEINE
C**
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**   VERFASSER:   G. DE MORSIER
 
      REAL        LAM,PHI,POLPHI,POLLAM
 
      DATA        ZRPI18 , ZPIR18  / 57.2957795 , 0.0174532925 /
 
      ZSINPOL = SIN(ZPIR18*POLPHI)
      ZCOSPOL = COS(ZPIR18*POLPHI)
      ZLAMPOL = ZPIR18*POLLAM
      ZPHI    = ZPIR18*PHI
      ZLAM    = LAM
      IF(ZLAM.GT.180.0) ZLAM = ZLAM - 360.0
      ZLAM    = ZPIR18*ZLAM
      ZARG    = ZCOSPOL*COS(ZPHI)*COS(ZLAM-ZLAMPOL) + ZSINPOL*SIN(ZPHI)
 
      PHTOPHS = ZRPI18*ASIN(ZARG)
 
      RETURN
      END
      
      
      REAL FUNCTION PHSTOPH (PHIS, LAMS, POLPHI, POLLAM)
C
C%Z% Modul %M%, V%I% vom %G%, extrahiert am %H%
C
C**** PHSTOPH  -   FC:BERECHNUNG DER WAHREN GEOGRAPHISCHEN BREITE FUER
C****                 EINEN PUNKT MIT DEN KOORDINATEN (PHIS, LAMS) IM
C****                 ROTIERTEN SYSTEM. DER NORDPOL DIESES SYSTEMS HAT
C****                 DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
C**   AUFRUF   :   PHI = PHSTOPH (PHIS, LAMS, POLPHI, POLLAM)
C**   ENTRIES  :   KEINE
C**   ZWECK    :   BERECHNUNG DER WAHREN GEOGRAPHISCHEN BREITE FUER
C**                EINEN PUNKT MIT DEN KOORDINATEN (PHIS, LAMS) IM
C**                ROTIERTEN SYSTEM. DER NORDPOL DIESES SYSTEMS HAT
C**                DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
C**   VERSIONS-
C**   DATUM    :   03.05.90
C**
C**   EXTERNALS:   KEINE
C**   EINGABE-
C**   PARAMETER:   PHIS     REAL   GEOGR. BREITE DES PUNKTES IM ROT.SYS.
C**                LAMS     REAL   GEOGR. LAENGE DES PUNKTES IM ROT.SYS.
C**                POLPHI   REAL   WAHRE GEOGR. BREITE DES NORDPOLS
C**                POLLAM   REAL   WAHRE GEOGR. LAENGE DES NORDPOLS
C**   AUSGABE-
C**   PARAMETER:   WAHRE GEOGRAPHISCHE BREITE ALS WERT DER FUNKTION
C**                ALLE WINKEL IN GRAD (NORDEN>0, OSTEN>0)
C**
C**   COMMON-
C**   BLOECKE  :   KEINE
C**
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**   VERFASSER:   D.MAJEWSKI
 
      REAL        LAMS,PHIS,POLPHI,POLLAM
 
      DATA        ZRPI18 , ZPIR18  / 57.2957795 , 0.0174532925 /
 
      SINPOL = SIN(ZPIR18*POLPHI)
      COSPOL = COS(ZPIR18*POLPHI)
      ZPHIS  = ZPIR18*PHIS
      ZLAMS  = LAMS
      IF(ZLAMS.GT.180.0) ZLAMS = ZLAMS - 360.0
      ZLAMS  = ZPIR18*ZLAMS
      ARG     = COSPOL*COS(ZPHIS)*COS(ZLAMS) + SINPOL*SIN(ZPHIS)
 
      PHSTOPH = ZRPI18*ASIN(ARG)
 
      RETURN
      END
      REAL FUNCTION LMTOLMS (PHI, LAM, POLPHI, POLLAM)
C
C%Z% Modul %M%, V%I% vom %G%, extrahiert am %H%
C
C**** LMTOLMS  -   FC:UMRECHNUNG DER WAHREN GEOGRAPHISCHEN LAENGE LAM
C****                 AUF EINEM PUNKT MIT DEN KOORDINATEN (PHIS, LAMS)
C****                 IM ROTIERTEN SYSTEM. DER NORDPOL DES SYSTEMS HAT
C****                 DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
C**   AUFRUF   :   LAM = LMTOLMS (PHI, LAM, POLPHI, POLLAM)
C**   ENTRIES  :   KEINE
C**   ZWECK    :   UMRECHNUNG DER WAHREN GEOGRAPHISCHEN LAENGE LAM AUF
C**                EINEM PUNKT MIT DEN KOORDINATEN (PHIS, LAMS) IM
C**                ROTIERTEN SYSTEM. DER NORDPOL DIESES SYSTEMS HAT
C**                DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
C**   VERSIONS-
C**   DATUM    :   03.05.90
C**
C**   EXTERNALS:   KEINE
C**   EINGABE-
C**   PARAMETER:   PHI    REAL BREITE DES PUNKTES IM GEOGR. SYSTEM
C**                LAM    REAL LAENGE DES PUNKTES IM GEOGR. SYSTEM
C**                POLPHI REAL GEOGR.BREITE DES N-POLS DES ROT. SYSTEMS
C**                POLLAM REAL GEOGR.LAENGE DES N-POLS DES ROT. SYSTEMS
C**   AUSGABE-
C**   PARAMETER:   WAHRE GEOGRAPHISCHE LAENGE ALS WERT DER FUNKTION
C**                ALLE WINKEL IN GRAD (NORDEN>0, OSTEN>0)
C**
C**   COMMON-
C**   BLOECKE  :   KEINE
C**
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**   VERFASSER:   G. DE MORSIER
 
      REAL        LAM,PHI,POLPHI,POLLAM
 
      DATA        ZRPI18 , ZPIR18  / 57.2957795 , 0.0174532925 /
 
      ZSINPOL = SIN(ZPIR18*POLPHI)
      ZCOSPOL = COS(ZPIR18*POLPHI)
      ZLAMPOL =     ZPIR18*POLLAM
      ZPHI    =     ZPIR18*PHI
      ZLAM    = LAM
      IF(ZLAM.GT.180.0) ZLAM = ZLAM - 360.0
      ZLAM    = ZPIR18*ZLAM
 
      ZARG1   = - SIN(ZLAM-ZLAMPOL)*COS(ZPHI)
      ZARG2   = - ZSINPOL*COS(ZPHI)*COS(ZLAM-ZLAMPOL)+ZCOSPOL*SIN(ZPHI)
      IF (ABS(ZARG2).LT.1.E-30) THEN
        IF (ABS(ZARG1).LT.1.E-30) THEN
          LMTOLMS =   0.0
        ELSEIF (ZARG1.GT.0.) THEN
              LMTOLMS =  90.0
            ELSE
              LMTOLMS = -90.0
            ENDIF
      ELSE
        LMTOLMS = ZRPI18*ATAN2(ZARG1,ZARG2)
      ENDIF
 
      RETURN
      END
      
      
      REAL FUNCTION LMSTOLM (PHIS, LAMS, POLPHI, POLLAM)
C
C%Z% Modul %M%, V%I% vom %G%, extrahiert am %H%
C
C**** LMSTOLM  -   FC:BERECHNUNG DER WAHREN GEOGRAPHISCHEN LAENGE FUER
C****                 EINEN PUNKT MIT DEN KOORDINATEN (PHIS, LAMS)
C****                 IM ROTIERTEN SYSTEM. DER NORDPOL DES SYSTEMS HAT
C****                 DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
C**   AUFRUF   :   LAM = LMSTOLM (PHIS, LAMS, POLPHI, POLLAM)
C**   ENTRIES  :   KEINE
C**   ZWECK    :   BERECHNUNG DER WAHREN GEOGRAPHISCHEN LAENGE FUER
C**                EINEN PUNKT MIT DEN KOORDINATEN (PHIS, LAMS)
C**                IM ROTIERTEN SYSTEM. DER NORDPOL DIESES SYSTEMS HAT
C**                DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
C**   VERSIONS-
C**   DATUM    :   03.05.90
C**
C**   EXTERNALS:   KEINE
C**   EINGABE-
C**   PARAMETER:   PHIS     REAL   GEOGR. BREITE DES PUNKTES IM ROT.SYS.
C**                LAMS     REAL   GEOGR. LAENGE DES PUNKTES IM ROT.SYS.
C**                POLPHI   REAL   WAHRE GEOGR. BREITE DES NORDPOLS
C**                POLLAM   REAL   WAHRE GEOGR. LAENGE DES NORDPOLS
C**   AUSGABE-
C**   PARAMETER:   WAHRE GEOGRAPHISCHE LAENGE ALS WERT DER FUNKTION
C**                ALLE WINKEL IN GRAD (NORDEN>0, OSTEN>0)
C**
C**   COMMON-
C**   BLOECKE  :   KEINE
C**
C**   FEHLERBE-
C**   HANDLUNG :   KEINE
C**   VERFASSER:   D.MAJEWSKI
 
      REAL        LAMS,PHIS,POLPHI,POLLAM
 
      DATA        ZRPI18 , ZPIR18  / 57.2957795 , 0.0174532925 /
 
      ZSINPOL = SIN(ZPIR18*POLPHI)
      ZCOSPOL = COS(ZPIR18*POLPHI)
      ZLAMPOL = ZPIR18*POLLAM
      ZPHIS   = ZPIR18*PHIS
      ZLAMS   = LAMS
      IF(ZLAMS.GT.180.0) ZLAMS = ZLAMS - 360.0
      ZLAMS   = ZPIR18*ZLAMS
 
      ZARG1   = SIN(ZLAMPOL)*(- ZSINPOL*COS(ZLAMS)*COS(ZPHIS)  +
     1                          ZCOSPOL*           SIN(ZPHIS)) -
     2          COS(ZLAMPOL)*           SIN(ZLAMS)*COS(ZPHIS)
      ZARG2   = COS(ZLAMPOL)*(- ZSINPOL*COS(ZLAMS)*COS(ZPHIS)  +
     1                          ZCOSPOL*           SIN(ZPHIS)) +
     2          SIN(ZLAMPOL)*           SIN(ZLAMS)*COS(ZPHIS)
      IF (ABS(ZARG2).LT.1.E-30) THEN
        IF (ABS(ZARG1).LT.1.E-30) THEN
          LMSTOLM =   0.0
        ELSEIF (ZARG1.GT.0.) THEN
              LMSTOLAM =  90.0
            ELSE
              LMSTOLAM = -90.0
            ENDIF
      ELSE
        LMSTOLM = ZRPI18*ATAN2(ZARG1,ZARG2)
      ENDIF
 
      RETURN
      END
