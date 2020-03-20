	subroutine read_prop
     +	(nlay,ndat,ntypes,nspars,nleg,depol,scatfile,
     +	 wavnum,gasdat,taudp,omega,coefsr,coefsa,ncoefsa,
     +	 spars,heights,theta,theta0,phi,fr,fa)

	implicit none

c  Inputs

	integer nlay,ndat,ntypes,nspars,nleg
	double precision depol
	character*15 scatfile(ntypes)

c  Outputs

	integer ncoefsa(ntypes,2)
	double precision wavnum(ndat)
	double precision gasdat(nlay,ndat)
	double precision taudp(nlay,ndat)
	double precision omega(nlay,ndat)
	double precision coefsr(0:2,6)
	double precision coefsa(ntypes,0:nleg-1,6,2)
	double precision spars(nspars)
	double precision heights(0:nlay)
	double precision theta,theta0,phi
	double precision fr(nlay,ndat)
	double precision fa(ntypes,nlay,ndat)

c  Local variables

	integer i,j,n
	character*2 string
	double precision raydat(nlay,ndat)
	double precision aerdat(ntypes,nlay,ndat)
	double precision omegaa(ntypes,nlay,ndat)
	double precision ray,aer(ntypes)

c  Read scattering matrix data

	call read_scat
     +	(ntypes,nleg,scatfile,depol,
     +   ncoefsa,coefsr,coefsa)

c  Read extinction and ssa profiles

	open(1,file='../data/od_ray_01.dat',status='OLD')
	do i = 1, ndat
	  read(1,'(61(e18.10,1x))') 
     +		wavnum(i),(raydat(n,i),n=nlay,1,-1)
	enddo
	close(1)

	do j= 1, ntypes
	  write(string,fmt='(i2)') j
	  open(1,file='../data/od_aero_0'//trim(adjustl(string))//
     +	      '_01.dat',status='OLD')
	  do i = 1, ndat
	    read(1,'(61(e18.10,1x))') wavnum(i),
     +		(aerdat(j,n,i),n=nlay,1,-1)
	  enddo
	  close(1)
	enddo

	open(1,file='../data/od_gas_01.dat',status='OLD')
	do i = 1, ndat
	  read(1,'(61(e18.10,1x))') wavnum(i),(gasdat(n,i),n=nlay,1,-1)
	enddo
	close(1)

	open(1,file='../data/ssa_01.dat',status='OLD')
	do i = 1, ndat
	  read(1,'(61(e18.10,1x))') wavnum(i),(omega(n,i),n=nlay,1,-1)
	enddo
	close(1)

	do j= 1, ntypes
	  write(string,fmt='(i2)') j
	  open(1,file='../data/ssa_aero_0'//trim(adjustl(string))//
     +	      '_01.dat',status='OLD')
	  do i = 1, ndat
	    read(1,'(61(e18.10,1x))') wavnum(i),
     +		(omegaa(j,n,i),n=nlay,1,-1)
	  enddo
	  close(1)
	enddo

	open(4,file='../data/geometry.dat',status='OLD')
	read(4,*) theta0
	read(4,*) theta
	read(4,*) phi
	close(4)

	open(4,file='../data/alb_01.dat',status='OLD')
	read(4,*) spars(1)
	close(4)

	do n = 0, nlay
	  heights(n) = dble(nlay)-dble(n)
	enddo

	do i = 1, ndat
	  do n = 1, nlay
	    taudp(n,i) = sum(aerdat(:,n,i))+gasdat(n,i)+raydat(n,i)
	  enddo
	enddo

	do i = 1, ndat
	  do n = 1, nlay
            ray = raydat(n,i)
	    do j = 1, ntypes
	      aer(j) = aerdat(j,n,i)*omegaa(j,n,i)
	    enddo
	    fr(n,i) = ray/(ray+sum(aer))
	    do j = 1, ntypes
	      fa(j,n,i) = aer(j)/(ray+sum(aer))
	    enddo
	  enddo
	enddo

c  Finish

	return
	end subroutine read_prop

	subroutine read_scat
     +	(ntypes,nleg,scatfile,depol,
     +   ncoefsa,coefsr,coefsa)

	implicit none

c  Inputs

	integer ntypes,nleg
	character*15 scatfile(ntypes)
	double precision depol

c  Outputs

	integer ncoefsa(ntypes,2)
	double precision coefsr(0:2,6)
	double precision coefsa(ntypes,0:nleg-1,6,2)

c  Local variables

	integer i,j,l,o1
	character*15 fname
	double precision temp

c  Rayleigh Greek matrix

	coefsr(0:,:) = 0.d0
	coefsr(0,1) = 1.d0
	coefsr(1,4) = 3.d0*(1.d0-2.d0*depol)/(2.d0+depol)
	coefsr(2,1) = (1.d0-depol)/(2.d0+depol)
	coefsr(2,5) = sqrt(6.d0)*(1.d0-depol)/(2.d0+depol)
	coefsr(2,2) = 6.d0*(1.d0-depol)/(2.d0+depol)

c  Read aerosol Greek matrices

	do i = 1, ntypes
	  fname = scatfile(i)
	  open(15,file='../data/mom_mie/'//trim(fname),status='OLD')
	  do j = 1,1
	    read(15,*)
	  enddo
	  read(15,*) temp,ncoefsa(i,1)
	  do l = 0, ncoefsa(i,1)
	    read(15,*) coefsa(i,l,2,1),coefsa(i,l,1,1),coefsa(i,l,5,1),
     +	       coefsa(i,l,4,1),coefsa(i,l,6,1),coefsa(i,l,3,1)
	    coefsa(i,l,5,1) = -coefsa(i,l,5,1)
	  enddo
	  do l = ncoefsa(i,1)+1, nleg-1
	    do o1 = 1, 6
	      coefsa(i,l,o1,1) = 0.d0
	    enddo
	  enddo
	  read(15,*) temp,ncoefsa(i,2)
	  do l = 0, ncoefsa(i,2)
	    read(15,*) coefsa(i,l,2,2),coefsa(i,l,1,2),coefsa(i,l,5,2),
     +	       coefsa(i,l,4,2),coefsa(i,l,6,2),coefsa(i,l,3,2)
	    coefsa(i,l,5,2) = -coefsa(i,l,5,2)
	  enddo
	  do l = ncoefsa(i,2)+1, nleg-1
	    do o1 = 1, 6
	      coefsa(i,l,o1,2) = 0.d0
	    enddo
	  enddo
	  close(15)
	enddo

c  Finish

	return
	end subroutine read_scat
