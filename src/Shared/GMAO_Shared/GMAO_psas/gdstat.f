	subroutine GDSTAT (lu,mx,my,mz,a,h,atype,htype,amiss,header,inc)

c	..Print statistics of one 3-d variable
	implicit none

	integer lu		! Output unit
	integer mx,my,mz	! Array sizes
	real a(mx,my,mz)	! The array
	real h(mz)		! The argument(levels)
	character*4 atype	! Type of the variable
	character*4 htype	! Typf of the levels
	real amiss		! missing value flag of a
	character*(*) header	! A header message
	integer inc		! order of the listing

	integer i,j,k
	integer kfr,kto,kinc
	integer imx,imn,jmx,jmn
	integer knt
	real amx,amn
	real avg,dev,d
	logical first

c	..A practical value for the magnitude of the fraction of a real
c	number.

	real rfrcval
	parameter(rfrcval=1.e-5)

c	..function

	logical spv
	real aspv
	spv(aspv)=abs((aspv-amiss)/amiss).le.rfrcval

	write(lu,'(/a)') header
	if(htype.eq.'PRES') then
	  write(lu,'(a,3x,a,2x,a,5x,a,6x,a,9x,a,15x,a)') 'lvl','mbar',
     >	    'count','mean','stdv','maxi','mini'
	elseif(htype.eq.'HGHT') then
	  write(lu,'(a,2x,a,2x,a,5x,a,6x,a,9x,a,15x,a)') 'lvl','meter',
     >	    'count','mean','stdv','maxi','mini'
	elseif(htype.eq.'TEMP') then
	  write(lu,'(a,4x,a,4x,a,5x,a,6x,a,9x,a,15x,a)') 'lvl','K',
     >	    'count','mean','stdv','maxi','mini'
	else
	  write(lu,'(a,4x,a,4x,a,5x,a,6x,a,9x,a,15x,a)') 'lvl',htype,
     >	    'count','mean','stdv','maxi','mini'
	endif

c	..Check the order of the listing, increase or decrease
	if(inc.ge.0) then
	  kfr=1
	  kto=mz
	  kinc=1
	else
	  kfr=mz
	  kto=1
	  kinc=-1
	endif

	do k=kfr,kto,kinc
	  knt=0
	  avg=0.
	  do j=1,my
	    do i=1,mx
	      if(.not.spv(a(i,j,k))) then
		knt=knt+1
	      	avg=avg+a(i,j,k)
	      endif
	    end do
	  end do
	  avg=avg/max(1,knt)

	  dev=0.
	  do j=1,my
	    do i=1,mx
	      if(.not.spv(a(i,j,k))) then
		d=a(i,j,k)-avg
		dev=dev+d*d
	      endif
	    end do
	  end do
	  dev=sqrt(dev/max(1,knt-1))

	  amx=a(1,1,k)
	  amn=a(1,1,k)
	  first=.true.
	  do j=1,my
	    do i=1,mx
	      if(.not.spv(a(i,j,k))) then
		if(first) then
		  imx=i
		  imn=i
		  jmx=j
		  jmn=j
		  amx=a(imx,jmx,k)
		  amn=a(imn,jmn,k)
		  first=.false.
		else
		  if(a(i,j,k).gt.amx) then
		    amx=a(i,j,k)
		    imx=i
		    jmx=j
		  endif
	      	  if(a(i,j,k).lt.amn) then
		    amn=a(i,j,k)
		    imn=i
		    jmn=j
		  endif
		endif
	      endif
	    end do
	  end do

	  if(atype.eq.'RELH') then
	    avg=avg*100.
	    dev=dev*100.
	    amx=amx*100.
	    amn=amn*100.
	  endif

	  if(htype.eq.'PRES'.or.htype.eq.'HGHT') then
	    if(atype.eq.'HGHT'.or.atype.eq.'STRM') then
	      write(lu,'(i3,2i7,2i10,'//
     &		'2(i10,a,i3,a,i3,a))')
     &		k,nint(h(k)),knt,nint(avg),nint(dev),
     &		nint(amx),'(',imx,',',jmx,')',
     &		nint(amn),'(',imn,',',jmn,')'
	    elseif(atype.eq.'TEMP'.or.atype.eq.'PRES'.or.
     &	      atype.eq.'WIND'.or.atype.eq.'%REH'.or.
     &	      atype.eq.'RELH'.or.atype.eq.'MIXR') then
	      write(lu,'(i3,2i7,2f10.2,'//
     &		'2(f10.2,a,i3,a,i3,a))')
     &		k,nint(h(k)),knt,avg,dev,
     &		amx,'(',imx,',',jmx,')',
     &		amn,'(',imn,',',jmn,')'
	    elseif(atype.eq.'NORM') then
	      write(lu,'(i3,2i7,2f10.4,'//
     &		'2(f10.4,a,i3,a,i3,a))')
     &		k,nint(h(k)),knt,avg,dev,
     &		amx,'(',imx,',',jmx,')',
     &		amn,'(',imn,',',jmn,')'
	    else
	      write(lu,'(i3,2i7,1p,2e10.3e1,0p,'//
     &		'2(1p,e10.3e1,0p,a,i3,a,i3,a))')
     &		k,nint(h(k)),knt,avg,dev,
     &		amx,'(',imx,',',jmx,')',
     &		amn,'(',imn,',',jmn,')'
	    endif

	  elseif(htype.eq.'TEMP') then
	    if(atype.eq.'HGHT'.or.atype.eq.'STRM') then
	      write(lu,'(i3,f7.2,i7,2i10,'//
     &		'2(i10,a,i3,a,i3,a))')
     &		k,h(k),knt,nint(avg),nint(dev),
     &		nint(amx),'(',imx,',',jmx,')',
     &		nint(amn),'(',imn,',',jmn,')'
	    elseif(atype.eq.'TEMP'.or.atype.eq.'PRES'.or.
     &	      atype.eq.'WIND'.or.atype.eq.'%REH'.or.
     &	      atype.eq.'RELH'.or.atype.eq.'MIXR') then
	      write(lu,'(i3,f7.2,i7,2f10.2,'//
     &		'2(f10.2,a,i3,a,i3,a))')
     &		k,h(k),knt,avg,dev,
     &		amx,'(',imx,',',jmx,')',
     &		amn,'(',imn,',',jmn,')'
	    elseif(atype.eq.'NORM') then
	      write(lu,'(i3,f7.2,i7,2f10.4,'//
     &		'2(f10.4,a,i3,a,i3,a))')
     &		k,h(k),knt,avg,dev,
     &		amx,'(',imx,',',jmx,')',
     &		amn,'(',imn,',',jmn,')'
	    else
	      write(lu,'(i3,f7.2,i7,1p,2e10.3e1,0p,'//
     &		'2(1p,e10.3e1,0p,a,i3,a,i3,a))')
     &		k,h(k),knt,avg,dev,
     &		amx,'(',imx,',',jmx,')',
     &		amn,'(',imn,',',jmn,')'
	    endif

	  else
	    if(atype.eq.'HGHT'.or.atype.eq.'STRM') then
	      write(lu,'(i3,1p,e10.3e1,0p,i7,2i10,'//
     &		'2(i10,a,i3,a,i3,a))')
     &		k,h(k),knt,nint(avg),nint(dev),
     &		nint(amx),'(',imx,',',jmx,')',
     &		nint(amn),'(',imn,',',jmn,')'
	    elseif(atype.eq.'TEMP'.or.atype.eq.'PRES'.or.
     &	      atype.eq.'WIND'.or.atype.eq.'%REH'.or.
     &	      atype.eq.'RELH'.or.atype.eq.'MIXR') then
	      write(lu,'(i3,1p,e10.3e1,0p,i7,2f10.2,'//
     &		'2(f10.2,a,i3,a,i3,a))')
     &		k,h(k),knt,avg,dev,
     &		amx,'(',imx,',',jmx,')',
     &		amn,'(',imn,',',jmn,')'
	    elseif(atype.eq.'NORM') then
	      write(lu,'(i3,1p,e10.3e1,0p,i7,2f10.4,'//
     &		'2(f10.4,a,i3,a,i3,a))')
     &		k,h(k),knt,avg,dev,
     &		amx,'(',imx,',',jmx,')',
     &		amn,'(',imn,',',jmn,')'
	    else
	      write(lu,'(i3,1p,e10.3e1,0p,i7,1p,2e10.3e1,0p,'//
     &		'2(1p,e10.3e1,0p,a,i3,a,i3,a))')
     &		k,h(k),knt,avg,dev,
     &		amx,'(',imx,',',jmx,')',
     &		amn,'(',imn,',',jmn,')'
	    endif

	  endif
	end do		! k=kfr,kto,kinc
	
	end	! gdstat()
