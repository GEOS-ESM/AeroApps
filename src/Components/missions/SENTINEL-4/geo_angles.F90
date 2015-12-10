!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    geo_angles
! PURPOSE
!     Creates files of satellite and solar position angles for a geostationary satellite
! INPUT
!     date  : string of variable name
!     time  : file to be read
!     inst  : instrument name
!     indir : main directory for input data
!     outdir: directory for output data
! OUTPUT
!     None
!  HISTORY
!     29 May 2015 P. Castellanos 
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
#  include "MAPL_Generic.h"
#  include "MAPL_ErrLogMain.h"
program geo_angles

  use ESMF                         ! ESMF modules
  use MAPL_Mod
  use MAPL_ShmemMod                ! The SHMEM infrastructure
  use netcdf                       ! for reading the NR files
  !use mpi
  use mp_netcdf_Mod                ! Module with netcdf routines
  use netcdf_Mod
  use GeoAngles                    ! Module with geostationary satellite algorithms for scene geometry

  implicit none
  include "mpif.h"

! ESMF Objects
!----------------
  type(ESMF_Config)       :: cf
  type(ESMF_VM)           :: vm   

! RC Inputs 
! --------
  character(len=256)                    :: arg
  character(len=8)                      :: date 
  character(len=2)                      :: time 
  character(len=256)                    :: instname, indir, outdir


! File names
! ----------
  character(len=256)                    :: MET_file, INV_file, OUT_file 

! Global, 3D arrays to be allocated using SHMEM
! ---------------------------------------------
  real*8, pointer                       :: CLON(:,:) => null()
  real*8, pointer                       :: CLAT(:,:) => null()
  real*8, pointer                       :: SCANTIME(:) => null()
  real, pointer                         :: SZA(:,:) => null()
  real, pointer                         :: VZA(:,:) => null()
  real, pointer                         :: SAA(:,:) => null()
  real, pointer                         :: VAA(:,:) => null()

! Working variables
!------------------------------
  real                                  :: sat_lon                   ! satellite longitude
  real                                  :: sat_lat                   ! satellite latitude
  real                                  :: sat_alt                   ! satellite altitude


! Satellite domain variables
!------------------------------
  integer                               :: im, jm                    ! size of TEMPO domain

! netcdf variables
!----------------------  
  integer                               :: ncid, szaVarID, vzaVarID, saaVarID, vaaVarID

! Miscellaneous
! -------------
  integer                               :: ierr                      ! MPI error message
  integer                               :: myid, npet, CoresPerNode  ! MPI dimensions and processor id
  integer                               :: p                         ! i-processor
  real                                  :: progress 
  real                                  :: MISSING

! ESMF variables
  integer                               :: rc, status          ! MPI error message

! System tracking variables
! -----------------------------
  integer*8                             :: t1, t2, clock_max
  real*8                                :: clock_rate
  character(len=*), parameter           :: Iam = 'geo_angles'

! Start Timing
! -----------------
  call system_clock ( t1, clock_rate, clock_max )
  progress = -1


! Initialize MPI with ESMF
! --------------
  call ESMF_Initialize (logkindflag=ESMF_LOGKIND_NONE, vm=vm, __RC__)

  call ESMF_VMGet(vm, localPET=myid, PETcount=npet) 
  if ( MAPL_am_I_root() ) write(*,'(A,I4,A)')'Starting MPI on ',npet, ' processors'


! Initialize SHMEM
! ----------------
  CoresPerNode = MAPL_CoresPerNodeGet(MPI_COMM_WORLD,rc=ierr) ! a must
  call MAPL_InitializeShmem(rc=ierr)

! Parse Resource file for input info 
! ----------------------------------------------------------------------------
  CALL getarg(1, arg)

  cf = ESMF_ConfigCreate()
  call ESMF_ConfigLoadFile(cf, fileName=trim(arg), __RC__)

  ! Read in variables
  ! -----------------------
  call ESMF_ConfigGetAttribute(cf, date, label = 'DATE:',__RC__)
  call ESMF_ConfigGetAttribute(cf, time, label = 'TIME:',__RC__)
  call ESMF_ConfigGetAttribute(cf, instname, label = 'INSTNAME:',__RC__)
  call ESMF_ConfigGetAttribute(cf, indir, label = 'INDIR:',__RC__)
  call ESMF_ConfigGetAttribute(cf, outdir, label = 'OUTDIR:',default=indir)

! Write out settings to use
! -------------------------------
  if (myid == 0) then
    write(*,*) 'Processing ', lower_to_upper(trim(instname)),' domain on ',date,' ', time, 'Z'
    write(*,*) 'Input directory: ',trim(indir)
    write(*,*) 'Output directory: ',trim(outdir)
    write(*,*) ' '
  end if 

! INFILES & OUTFILE NAMES
!-------------------------
  call filenames()

! Query for domain dimensions, satellite constants, and missing value
!--------------------------------------------------------------
  call mp_readDim("ew", MET_file, im)
  call mp_readDim("ns", MET_file, jm)
  call mp_readGattr("sat_lat", INV_FILE, sat_lat)
  call mp_readGattr("sat_lon", INV_FILE, sat_lon)
  call mp_readGattr("sat_alt", INV_FILE, sat_alt)
  call mp_readVattr("missing_value", INV_FILE, "clon", MISSING)

! Create OUTFILE
!-------------------------
  if ( MAPL_am_I_root() ) call create_outfile(date,time,ncid,szaVarID, vzaVarID, saaVarID, vaaVarID)

! Allocate the Global arrays using SHMEM
! It will be available on all processors
! ---------------------------------------------------------
  call allocate_shared(im, jm)

! Read in position data to shared memeory
! --------------------------------------------
  call mp_readvar1Dchunk('scanTime', MET_file, (/im/), 1, npet, myid, SCANTIME)
  call mp_readvar2Dchunk('clon', INV_file, (/im, jm/), 1, npet, myid, CLON)
  call mp_readvar2Dchunk('clat', INV_file, (/im, jm/), 1, npet, myid, CLAT)

! Wait for everyone to finish reading 
! ------------------------------------------------------------------  
  call MAPL_SyncSharedMemory(rc=ierr)  

! Calculate the viewing geometry - these are shared memory arrays
! Distribute calculations over processors
!-----------------------------------------------
  call get_geometry()


  if (MAPL_am_I_root()) then
    write(*,*) '<> Finished calculating sensor and solar angles for '//trim(lower_to_upper(instname))//' domain'
    call system_clock ( t2, clock_rate, clock_max )
    write ( *, * ) 'Elapsed real time = ', real ( t2 - t1 ) / real ( clock_rate )
    write(*,*) ' '
  end if

! All done
! --------
  call MAPL_SyncSharedMemory(rc=ierr)
  call shutdown()

! -----------------------------------------------------------------------------------

  contains


!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     get_geometry
! PURPOSE
!     calculate angles
! INPUT
!     none
! OUTPUT
!     none
!  HISTORY
!     29 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
subroutine get_geometry()
  integer                              :: yr, mo, day, hr, min 
  real                                 :: sec
  integer                              :: i, j, p, startj, countj, endj
  integer, dimension(npet)             :: nchunk                             ! how many chunks each processor works on                
  real,dimension(4)                    :: sat_angles


  ! Parse out date and time
  ! -------------------------
  read(time,*) hr
  read(date(1:4),*)  yr
  read(date(5:6),*)  mo
  read(date(7:8),*) day

  ! Everyone Figure out how many indeces each PE has to read
  ! -----------------------------
  nchunk = 0
  if (npet >= jm) then
    nchunk(1:npet) = 1
  else if (npet < jm) then
    nchunk(1:npet) = jm/npet
    nchunk(npet)   = nchunk(npet) + mod(jm,npet)
  end if 


  do i = 1, im
    min = SCANTIME(i)/60
    sec = SCANTIME(i) - min*60.0
    if (myid == 0) then
      startj = 1
    else
      startj = sum(nchunk(1:myid))+1
    end if
    countj = nchunk(myid+1)
    endj   = startj + countj - 1

    do j = startj, endj  
      if (CLAT(i,j) .ne. MISSING) then      
        sat_angles = satellite_angles(yr,mo,day,hr,min,dble(sec),dble(0.0),dble(CLAT(i,j)),dble(CLON(i,j)),dble(sat_lat),dble(sat_lon),dble(sat_alt))

        ! Store in shared memeory
        SZA(i,j) = sat_angles(4)
        VZA(i,j) = sat_angles(2)
        SAA(i,j) = sat_angles(3)
        VAA(i,j) = sat_angles(1)

      end if
    end do
  end do
  call MAPL_SyncSharedMemory(rc=ierr)
  if (MAPL_am_I_root()) then
    call check( nf90_open(OUT_file, nf90_write, ncid), "opening file " // OUT_file )
    call check(nf90_put_var(ncid,szaVarID,SZA), "writing out sza")
    call check(nf90_put_var(ncid,vzaVarID,VZA), "writing out vza")
    call check(nf90_put_var(ncid,saaVarID,SAA), "writing out saa")
    call check(nf90_put_var(ncid,vaaVarID,VAA), "writing out vaa")
    call check( nf90_close(ncid), "close outfile" )
  end if

end subroutine get_geometry

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     filenames
! PURPOSE
!     populate file name varaibles
! INPUT
!     none
! OUTPUT
!     none
!  HISTORY
!     26 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
subroutine filenames()

  ! INFILES
  write(MET_file,'(14A)') trim(indir),'/LevelB/Y',date(1:4),'/M',date(5:6),'/D',date(7:8),'/',trim(instname),'-g5nr.lb2.met_Nv.',date,'_',time,'z.nc4'
  write(INV_file,'(4A)')  trim(indir),'/LevelG/invariant/',trim(instname),'.lg1.invariant.nc4'

! OUTFILES
  write(OUT_file,'(14A)') trim(outdir),'/LevelB/Y',date(1:4),'/M',date(5:6),'/D',date(7:8),'/',trim(instname),'.lb2.angles.',date,'_',time,'z.nc4'
end subroutine filenames

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     lower_to_upper
! PURPOSE
!     converts lower case characters to uppercase
! INPUT
!     word
! OUTPUT
!     lower_to_upper
!  HISTORY
!     18 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;   

  function lower_to_upper ( word )
    implicit none

    character(len=*), intent(in)   :: word
    character(len=:),allocatable   :: lower_to_upper
    integer                        :: i,n

    character(*), parameter :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
    character(*), parameter :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    allocate(character(len=len(word)) :: lower_to_upper)

    lower_to_upper = word

    ! Loop over string elements
    do i = 1, LEN(word)
    ! Find location of letter in lower case constant string
      n = INDEX(LOWER_CASE, word( i:i ))
    ! If current substring is a lower case letter, make it upper case
      if ( n /= 0 ) lower_to_upper( i:i ) = UPPER_CASE( n:n )
    end do

  end function lower_to_upper

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     allocate_shared
! PURPOSE
!     allocates all the shared memory arrays that I need
! INPUT
!     im, jm   :: domain size
! OUTPUT
!     none
!  HISTORY
!     29 May 2015 P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  subroutine allocate_shared(im, jm)
    integer, intent(in)    :: im, jm

    call MAPL_AllocNodeArray(SCANTIME,(/im/),rc=ierr)
    call MAPL_AllocNodeArray(CLON,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(CLAT,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(SZA,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(VZA,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(SAA,(/im,jm/),rc=ierr)
    call MAPL_AllocNodeArray(VAA,(/im,jm/),rc=ierr)

  end subroutine allocate_shared


!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!     create_outfile
! PURPOSE
!     creates a netcdf4 file that can be written to my multiple processors independently
! INPUT
!     none
! OUTPUT
!     none
!  HISTORY
!     6 May 2015 P. Castellanos
!    29 May 2015 P. Castellanos - modified for satellite angles only
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  subroutine create_outfile(date, time, ncid, szaVarID, vzaVarID, saaVarID, vaaVarID)
    character(len=*)                   :: date, time
    integer,intent(out)                :: ncid    
    integer,intent(out)                :: szaVarID, vzaVarID, saaVarID, vaaVarID    
    
    integer, dimension(4)              :: chunk_size
    integer                            :: ewDimID, nsDimID
    integer                            :: clonVarID, clatVarID, timeVarID, ewVarID, nsVarID
    real*8,allocatable,dimension(:,:)  :: clon, clat
    real*8,allocatable,dimension(:)    :: scantime, ew, ns


    call check(nf90_create(OUT_file, IOR(nf90_netcdf4, nf90_clobber), ncid), "creating file " // OUT_file)
    call check(nf90_def_dim(ncid, "ew", im, ewDimID), "creating ew dimension") !im
    call check(nf90_def_dim(ncid, "ns", jm, nsDimID), "creating ns dimension") !jm

    call check(nf90_put_att(ncid,NF90_GLOBAL,'title','Satellite and solar position angles for '//lower_to_upper(instname)),"title attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,'institution','NASA/Goddard Space Flight Center'),"institution attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,'source','Global Model and Assimilation Office'),"source attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,'history','Created from geo_angles.x'),"history attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,'references','n/a'),"references attr")    
    call check(nf90_put_att(ncid,NF90_GLOBAL,'comment','This file contains SZA, VZA, SAA and VAA for the TEMPO geostationary grid '),"comment attr")   
    call check(nf90_put_att(ncid,NF90_GLOBAL,"contact","Arlindo da Silva <arlindo.dasilva@nasa.gov>"),"contact attr")
    call check(nf90_put_att(ncid,NF90_GLOBAL,"Conventions","cf"),"conventions attr")


    call check(nf90_def_var(ncid,'scanTime',nf90_float,(/ewDimID/),timeVarID),"create scanTime var")
    call check(nf90_def_var(ncid,'ew',nf90_float,(/ewDimID/),ewVarID),"create ew var")
    call check(nf90_def_var(ncid,'ns',nf90_float,(/nsDimID/),nsVarID),"create ns var")
    call check(nf90_def_var(ncid,'clon',nf90_float,(/ewDimID,nsDimID/),clonVarID),"create clon var")
    call check(nf90_def_var(ncid,'clat',nf90_float,(/ewDimID,nsDimID/),clatVarID),"create clat var")
    call check(nf90_def_var(ncid,'solar_zenith',nf90_float,(/ewDimID,nsDimID/),szaVarID),"create solar_zenith var")
    call check(nf90_def_var(ncid,'sensor_zenith',nf90_float,(/ewDimID,nsDimID/),vzaVarID),"create sensor_zenith var")
    call check(nf90_def_var(ncid,'solar_azimuth',nf90_float,(/ewDimID,nsDimID/),saaVarID),"create solar_azimuth var")
    call check(nf90_def_var(ncid,'sensor_azimuth',nf90_float,(/ewDimID,nsDimID/),vaaVarID),"create sensor_azimuth var")    

    call check(nf90_put_att(ncid,timeVarID,'long_name','Initial Time of Scan'),"long_name attr")
    call check(nf90_put_att(ncid,timeVarID,'units','seconds since '//date(1:4)//'-'//date(5:6)//'-'//date(7:8)//' '// &
                                                  time//':00:00'),"units attr")

    call check(nf90_put_att(ncid,ewVarID,'long_name','pseudo longitude'),"long_name attr")
    call check(nf90_put_att(ncid,ewVarID,'units','degrees_east'),"units attr")
    call check(nf90_put_att(ncid,nsVarID,'long_name','pseudo latitude'),"long_name attr")
    call check(nf90_put_att(ncid,nsVarID,'units','degrees_north'),"units attr")      

    call check(nf90_put_att(ncid,clonVarID,'long_name','pixel center longitude'),"long_name attr")
    call check(nf90_put_att(ncid,clonVarID,'missing_value',real(MISSING)),"missing_value attr")
    call check(nf90_put_att(ncid,clatVarID,'long_name','pixel center latitude'),"long_name attr")
    call check(nf90_put_att(ncid,clatVarID,'missing_value',real(MISSING)),"missing_value attr")  

    call check(nf90_put_att(ncid,szaVarID,'long_name','solar zenith angle (SZA)'),"long_name attr")
    call check(nf90_put_att(ncid,szaVarID,'missing_value',real(MISSING)),"missing_value attr")
    call check(nf90_put_att(ncid,szaVarID,'units','degrees'),"units attr")

    call check(nf90_put_att(ncid,vzaVarID,'long_name','sensor zenith angle (VZA)'),"long_name attr")
    call check(nf90_put_att(ncid,vzaVarID,'missing_value',real(MISSING)),"missing_value attr")
    call check(nf90_put_att(ncid,vzaVarID,'units','degrees'),"units attr") 

    call check(nf90_put_att(ncid,saaVarID,'long_name','solar azimuth angle (SAA)'),"long_name attr")
    call check(nf90_put_att(ncid,saaVarID,'missing_value',real(MISSING)),"missing_value attr")
    call check(nf90_put_att(ncid,saaVarID,'units','degrees clockwise from North'),"units attr")  

    call check(nf90_put_att(ncid,vaaVarID,'long_name','sensor viewing azimuth angle (VAA)'),"long_name attr")
    call check(nf90_put_att(ncid,vaaVarID,'missing_value',real(MISSING)),"missing_value attr")
    call check(nf90_put_att(ncid,vaaVarID,'units','degrees clockwise from North'),"units attr")  

    !Leave define mode
    call check(nf90_enddef(ncid),"leaving define mode")

    ! write out clon, clat, scantime
    allocate (scantime(im))
    allocate (clon(im, jm))
    allocate (clat(im, jm))
    allocate (ew(im))
    allocate (ns(jm))

    call readvar1D("scanTime", MET_file, scantime)
    call check(nf90_put_var(ncid,timeVarID,scantime), "writing out scantime")

    call readvar1D("ew", MET_file, ew)
    call check(nf90_put_var(ncid,ewVarID,ew), "writing out ew")

    call readvar1D("ns", MET_file, ns)
    call check(nf90_put_var(ncid,nsVarID,ns), "writing out ns")    

    call readvar2D("clon", INV_file, clon)
    call check(nf90_put_var(ncid,clonVarID,clon), "writing out clon")

    call readvar2D("clat", INV_file, clat)
    call check(nf90_put_var(ncid,clatVarID,clat), "writing out clat")

    deallocate (clon)
    deallocate (clat)  
    deallocate (scantime)
    deallocate (ns)
    deallocate (ew)

    call check( nf90_close(ncid), "close outfile" )
  end subroutine create_outfile  

!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
! NAME
!    shutdown
! PURPOSE
!     Clean up allocated arrays and finalize MPI
! INPUT
!     None
! OUTPUT
!     None
!  HISTORY
!     27 April P. Castellanos
!;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  subroutine shutdown()         

    ! shmem must deallocate shared memory arrays
    call MAPL_DeallocNodeArray(CLON,rc=ierr)
    call MAPL_DeallocNodeArray(CLAT,rc=ierr)
    call MAPL_DeallocNodeArray(SCANTIME,rc=ierr)

    call MAPL_FinalizeShmem (rc=ierr)

    call ESMF_ConfigDestroy(cf)

    call ESMF_Finalize(__RC__)

  end subroutine shutdown


end program geo_angles
