!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !PROGRAM:  saber2ods.x --- Reads a SABER data files and writes ODS files
! 
! !DESCRIPTION: 
!  This module ingests SABER data files and writes ODS files.  It uses 
!  the m_saber.f module.
! 
! !USAGE: saber2ods.x [-o template] saber_ncfile(s)
!
!         where:  template            template for ODS output file(s)
!                 saber_ncfile(s)     Input netcdf SABER file(s)"
!
!
! !REVISION HISTORY: 
!
! 18Mar2003  T. King      First crack.
!
!EOP
!-------------------------------------------------------------------------

      program saber2ods

      use m_saber

      use m_ods

      implicit none

      integer :: i,j,k,l
      integer :: num_options,nhms(5),nymd,iargc,pos
      integer :: rc_saber_get,rc_ods_clean,rc_ods_init,
     $     rc_ods_put,rc_ods_get,rc_ods_open,rc_ods_close
      integer :: ierr_write,ierr_app
      integer :: ncid
      integer :: fyear,fdoy,fjul_day,fcal_day,fjul_day_p1,fcal_day_p1,
     $     jul_day
      integer :: ODS_Julian,ODS_Caldat
      integer :: len_template
      character*180 fnames(200),tempfname,arg,template,fn,odsname
      character*8 cnymd
      character*4 cfyear
      character*3 cfdoy
      character*9 ftype
      character*10 s_string
      logical :: present,append,hit
      type(ods_vect)  ::  ods_struct,ods_old ! ODS vector
      integer :: tnobs,fc

! Initialize variables
      tnobs=450000
      ncid=11
      ftype='pre_anal'
      nhms(1)=000000
      nhms(2)=060000
      nhms(3)=120000
      nhms(4)=180000
      nhms(5)=000000
      i=1
      fc=1
      s_string='SABER_L2A_'
! Process argument list
      num_options = iargc()
      if (num_options .eq. 0 ) then
         print *, "usage: saber2ods.x [-o template] saber_ncfile(s)"
         print *, "template         template for ODS output file(s)"
         print *, "saber_ncfile(s)  Input netcdf SABER file(s)"
         stop
      endif
      template=''
      do while (i .le. num_options)
         call getarg (i,arg)
         pos=index(arg,'-o')
         if (i.eq.1.and.pos.gt.0)then
            call getarg (2,arg)
            template=arg
            i=3
         else
            if(len(trim(template)).eq.0)then
               template='saber.l2a.obs.%y4%m2%d2.ods'
            endif
            fnames(fc)=trim(arg)
            fc=fc+1
            i=i+1
         endif
      enddo

! Begin reading files.
      do j=1,fc-1
         
! Search for date string in the filename.
         pos=index(fnames(j),s_string)
         tempfname=fnames(j)
         cfyear=tempfname(pos+10:pos+13)
         cfdoy=tempfname(pos+14:pos+16)
         read(cfyear,204)fyear
 204     format(i4)
         read(cfdoy,204)fdoy
 205     format(i3)
         fyear=(10000*fyear)+101

         fjul_day=ODS_Julian(fyear)+fdoy-1
         fjul_day_p1=ODS_Julian(fyear)+fdoy

         fcal_day=ODS_Caldat(fjul_day)
         fcal_day_p1=ODS_Caldat(fjul_day_p1)
         hit=.false.
         do k=1,5 ! Cycle through the possible synoptic times in the file
            present=.false.
            append=.false.
            nymd=fcal_day
            if(k.eq.5)then
               nymd=fcal_day_p1
            endif

! Initialize ods vector
            call ODS_Init(ods_struct,tnobs,rc_ods_init)
            call SABER_Get(fnames(j),nymd,nhms(k),ods_struct,
     $           rc_saber_get)

! If data matching the selected synoptic time were found, continue.
            if(rc_saber_get.eq.0)then ! If something was read from the file, continue
               hit=.true.

! Create ods file name from template
               pos=index(template,'%y4%m2%d2')
               len_template=len(trim(template))
               write(cnymd,'(i8)')nymd
               odsname=template(1:pos-1)//cnymd//template(pos+9:len_template)
               print*,'Writing data at ',nymd,nhms(k),'Z to ',odsname 

! Check to see if ods file already exists.
               inquire(file=odsname,exist=present)
               if(present)then

! If the file is present read it and get the number of obs
! and then add those nobs to the ks values so ks is unique
! throughout the synoptic time.  We also need to know if 
! there are already obs written for the current synoptic 
! time.  If so, we have to set append=.true.

                  call ODS_Get(odsname,nymd,nhms(k),ftype,ods_old,rc_ods_get)
                  print*, 'rc_ods_get = ',rc_ods_get
                  print*,'number of obs found for ',nhms(k),'Z is ',ods_old%data%nobs
                  if(ods_old%data%nobs.eq.0)then
                     print*,'Writing data to a new synopitic time'
                     call ODS_Put(odsname,ftype,nymd,nhms(k),ods_struct,
     $                    rc_ods_put)
                     print*,'rc_ods_put = ',rc_ods_put
                  else
                     print*, 'Appending data to ',nhms(k),' Z'
                     do l=1,ods_struct%data%nobs
                        ods_struct%data%ks(l)=ods_old%data%ks(ods_old%data%nobs)
     $                       +ods_struct%data%ks(l)
                     enddo
! Append data to the current synoptic time.
                     call ODS_Open(ncid,odsname,'w',rc_ods_open)
                     print*,'rc_ods_open=',rc_ods_open
                     call ODS_Append(ncid,ods_struct%data%nobs,ierr_app)
                     print*,'ierr_app=',ierr_app
                     jul_day=ODS_Julian(nymd)
                     call ODS_PutR(ncid,'lat',jul_day,nhms(k)/10000,ods_struct%data%nobs,
     $              ods_struct%data%lat,ierr_write)
                     print*,'ierr_write=',ierr_write
                     call ODS_PutR(ncid,'lon',jul_day,nhms(k)/10000,ods_struct%data%nobs,
     $              ods_struct%data%lon,ierr_write)
                     print*,'ierr_write=',ierr_write
                     call ODS_PutR(ncid,'lev',jul_day,nhms(k)/10000,ods_struct%data%nobs,
     $              ods_struct%data%lev,ierr_write)
                     print*,'ierr_write=',ierr_write
                     call ODS_PutI(ncid,'time',jul_day,nhms(k)/10000,ods_struct%data%nobs,
     $              ods_struct%data%time,ierr_write)
                     print*,'ierr_write=',ierr_write
                     call ODS_PutI(ncid,'kt',jul_day,nhms(k)/10000,ods_struct%data%nobs,
     $              ods_struct%data%kt,ierr_write)
                     print*,'ierr_write=',ierr_write
                     call ODS_PutI(ncid,'kx',jul_day,nhms(k)/10000,ods_struct%data%nobs,
     $              ods_struct%data%kx,ierr_write)
                     print*,'ierr_write=',ierr_write
                     call ODS_PutI(ncid,'ks',jul_day,nhms(k)/10000,ods_struct%data%nobs,
     $              ods_struct%data%ks,ierr_write)
                     print*,'ierr_write=',ierr_write
                     call ODS_PutR(ncid,'xm',jul_day,nhms(k)/10000,ods_struct%data%nobs,
     $              ods_struct%data%xm,ierr_write)
                     print*,'ierr_write=',ierr_write
                     call ODS_PutI(ncid,'qcexcl',jul_day,nhms(k)/10000,ods_struct%data%nobs,
     $              ods_struct%data%qcexcl,ierr_write)
                     print*,'ierr_write=',ierr_write
                     call ODS_PutI(ncid,'qchist',jul_day,nhms(k)/10000,ods_struct%data%nobs,
     $              ods_struct%data%qchist,ierr_write)
                     print*,'ierr_write=',ierr_write
                     call ODS_PutR(ncid,'obs',jul_day,nhms(k)/10000,ods_struct%data%nobs,
     $              ods_struct%data%obs,ierr_write)
                     print*,'ierr_write=',ierr_write

                     call ODS_Close(ncid,'SABER',rc_ods_close)
                     print*,'rc_ods_close=',rc_ods_close
                  endif
               else
! Write data to an entirely NEW ods file
                  call ODS_Put(odsname,ftype,nymd,nhms(k),ods_struct,rc_ods_put)
                  print*,'rc_ods_put = ',rc_ods_put
               endif
            else
               print*,'No data at ',nymd,nhms(k) 
            endif
! Deallocate the ods vector
            call ODS_Clean(ods_struct,rc_ods_clean)
         enddo 
         if(.not.hit)then
            print*, 'Warning! No data matching the files date and times were found.'
         endif
      enddo
      stop
      end

