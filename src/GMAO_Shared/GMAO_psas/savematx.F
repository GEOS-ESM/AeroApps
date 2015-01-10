	subroutine savematx(prefix,kri,kti,lni,krj,ktj,lnj,cor)

	implicit none
	character(len=*),intent(in) :: prefix
	integer, intent(in) :: kri,krj
	integer, intent(in) :: kti,ktj
	integer, intent(in) :: lni,lnj
	real   , intent(in) :: cor(lni,lnj)

! Locals.

	integer  luavail,lnblnk
	external luavail,lnblnk

	integer lumat,i,j,ios

	character*80 fnmat
	
		! Set the filename

	  j=max(lnblnk(prefix),1)
	  write(fnmat,'(a,4(a,i2.2),2(a,i5.5),a)') prefix(1:j),
     &	    '.',kri,'-',krj,'.v',kti,'-',ktj,'.',lni,'x',lnj,'.mat'
	  do i=1,j
	    if(fnmat(i:i).eq.' ') fnmat(i:i)='_'
	  end do

		! Check the existence

	  lumat=luavail()
	  open(lumat,file=fnmat,status='new',iostat=ios,
     &	    form='unformatted',access='direct',recl=lni)

		! Save the file
	  if(ios.eq.0) then

	    do j=1,lnj
	      write(lumat,rec=j,iostat=ios) (cor(i,j),i=1,lni)
	    end do

	    close(lumat)
	  endif

	end
!.
