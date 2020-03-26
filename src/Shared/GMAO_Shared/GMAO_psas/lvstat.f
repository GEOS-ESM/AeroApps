	subroutine LVSTAT (lu,mx,my,a,h,atype,htype,amiss,flag)

c	..Print statistics of one 2-d variable
	implicit none

	integer lu		! Output unit
	integer mx,my		! Array sizes
	real a(mx,my)		! The array
	real h			! The argument
	character*4 atype	! Type of the variable(array)
	character*4 htype	! Type of the level(argument)
	real amiss		! missing value flag of a
	character*(*) flag

	integer i,j
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

	knt=0
	avg=0.
	do j=1,my
	  do i=1,mx
	    if(.not.spv(a(i,j))) then
	      avg=avg+a(i,j)
	      knt=knt+1
	    endif
	  end do
	end do
	avg=avg/max(1,knt)

	dev=0.
	do j=1,my
	  do i=1,mx
	    if(.not.spv(a(i,j))) then
	      d=a(i,j)-avg
	      dev=dev+d*d
	    endif
	  end do
	end do
	dev=sqrt(dev/max(1,knt-1))

	amx=a(1,1)
	amn=a(1,1)
	first=.true.
	do j=1,my
	  do i=1,mx
	    if(.not.spv(a(i,j))) then
	      if(first) then
		imx=i
		imn=i
		jmx=j
		jmn=j
		amx=a(imx,jmx)
		amn=a(imn,jmn)
		first=.false.
	      else
		if(a(i,j).gt.amx) then
		  amx=a(i,j)
		  imx=i
		  jmx=j
		endif
		if(a(i,j).lt.amn) then
		  amn=a(i,j)
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

	if(htype.eq.'NULL'.or.htype.eq.'SRFC') then
	  if(atype.eq.'HGHT'.or.atype.eq.'STRM') then
	    write(lu,'(2a,i6,i7,a,i6,2(a,i7,a,i3,a,i3,a))')
     &	      flag,':',knt,nint(avg),'+',nint(dev),
     &	      ' mx=',nint(amx),'(',imx,',',jmx,')',
     &	      ' mi=',nint(amn),'(',imn,',',jmn,')'
	  elseif(atype.eq.'PRES'.or.atype.eq.'TEMP'.or.
     &	    atype.eq.'WIND'.or.atype.eq.'%REH'.or.
     &	    atype.eq.'RELH'.or.atype.eq.'MIXR') then
	    write(lu,'(2a,i6,f7.1,a,f6.1,2(a,f7.1,a,i3,a,i3,a))')
     &	      flag,':',knt,avg,'+',dev,
     &	      ' mx=',amx,'(',imx,',',jmx,')',
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  elseif(atype.eq.'NORM') then
	    write(lu,'(2a,i6,f7.3,a,f6.3,2(a,f7.3,a,i3,a,i3,a))')
     &	      flag,':',knt,avg,'+',dev,
     &	      ' mx=',amx,'(',imx,',',jmx,')',
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  else
	    write(lu,'(2a,i6,1p,e10.3e1,a,e9.3e1,0p,'//
     &	      '2(a,1p,e10.3e1,0p,a,i3,a,i3,a))')
     &	      flag,':',knt,avg,'+',dev,
     &	      ' mx=',amx,'(',imx,',',jmx,')',
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  endif

	elseif(htype.eq.'PRES') then
	  if(atype.eq.'HGHT'.or.atype.eq.'STRM') then
	    write(lu,'(2a,i4,a,i6,i7,a,i6,2(a,i7,a,i3,a,i3,a))')
     &	      flag,'(',nint(h),')',knt,nint(avg),'+',nint(dev),
     &	      ' mx=',nint(amx),'(',imx,',',jmx,')',
     &	      ' mi=',nint(amn),'(',imn,',',jmn,')'
	  elseif(atype.eq.'PRES'.or.atype.eq.'TEMP'.or.
     &	    atype.eq.'WIND'.or.atype.eq.'%REH'.or.
     &	    atype.eq.'RELH'.or.atype.eq.'MIXR') then
	    write(lu,'(2a,i4,a,i6,f7.1,a,f6.1,'//
     &	      '2(a,f7.1,a,i3,a,i3,a))')
     &	      flag,'(',nint(h),')',knt,avg,'+',dev,
     &	      ' mx=',amx,'(',imx,',',jmx,')',
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  elseif(atype.eq.'NORM') then
	    write(lu,'(2a,i4,a,i6,f7.3,a,f6.3,'//
     &	      '2(a,f7.3,a,i3,a,i3,a))')
     &	      flag,'(',nint(h),')',knt,avg,'+',dev,
     &	      ' mx=',amx,'(',imx,',',jmx,')',
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  else
	    write(lu,'(2a,i4,a,i6,1p,e10.3e1,a,e9.3e1,0p,'//
     &	      '2(a,e10.3e1,a,i3,a,i3,a))')
     &	      flag,'(',nint(h),')',knt,avg,'+',dev,
     &	      ' mx=',amx,'(',imx,',',jmx,')',
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  endif

	elseif(htype.eq.'HGHT') then
	  if(atype.eq.'HGHT'.or.atype.eq.'STRM') then
	    write(lu,'(2a,i5,a,i6,i7,a,i6,2(a,i7,a,i3,a,i3,a))')
     &	      flag,'(',nint(h),')',knt,nint(avg),'+',nint(dev),
     &	      ' mx=',nint(amx),'(',imx,',',jmx,')',
     &	      ' mi=',nint(amn),'(',imn,',',jmn,')'
	  elseif(atype.eq.'PRES'.or.atype.eq.'TEMP'.or.
     &	    atype.eq.'WIND'.or.atype.eq.'%REH'.or.
     &	    atype.eq.'RELH'.or.atype.eq.'MIXR') then
	    write(lu,'(2a,i5,a,i6,f7.1,a,f6.1,'//
     &	      '2(a,f7.1,a,i3,a,i3,a))')
     &	      flag,'(',nint(h),')',knt,avg,'+',dev,
     &	      ' mx=',amx,'(',imx,',',jmx,')',
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  elseif(atype.eq.'NORM') then
	    write(lu,'(2a,i5,a,i6,f7.3,a,f6.3,'//
     &	      '2(a,f7.3,a,i3,a,i3,a))')
     &	      flag,'(',nint(h),')',knt,avg,'+',dev,
     &	      ' mx=',amx,'(',imx,',',jmx,')',
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  else
	    write(lu,'(2a,i5,a,i6,1p,e10.3e1,a,e9.3e1,0p,'//
     &	      '2(a,e10.3e1,a,i3,a,i3,a))')
     &	      flag,'(',nint(h),')',knt,avg,'+',dev,
     &	      ' mx=',amx,'(',imx,',',jmx,')',
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  endif

	elseif(htype.eq.'TEMP') then
	  if(atype.eq.'HGHT'.or.atype.eq.'STRM') then
	    write(lu,'(2a,f7.2,a,i6,i7,a,i6,2(a,i7,a,i3,a,i3,a))')
     &	      flag,'(',h,')',knt,nint(avg),'+',nint(dev),
     &	      ' mx=',nint(amx),'(',imx,',',jmx,')',
     &	      ' mi=',nint(amn),'(',imn,',',jmn,')'
	  elseif(atype.eq.'PRES'.or.atype.eq.'TEMP'.or.
     &	    atype.eq.'WIND'.or.atype.eq.'%REH'.or.
     &	    atype.eq.'RELH'.or.atype.eq.'MIXR') then
	    write(lu,'(2a,f7.2,a,i6,f7.1,a,f6.1,'//
     &	      '2(a,f7.1,a,i3,a,i3,a))')
     &	      flag,'(',h,')',knt,avg,'+',dev,
     &	      ' mx=',amx,'(',imx,',',jmx,')',
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  elseif(atype.eq.'NORM') then
	    write(lu,'(2a,f7.2,a,i6,f7.3,a,f6.3,'//
     &	      '2(a,f7.3,a,i3,a,i3,a))')
     &	      flag,'(',h,')',knt,avg,'+',dev,
     &	      ' mx=',amx,'(',imx,',',jmx,')',
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  else
	    write(lu,'(2a,f7.2,a,i6,1p,e10.3e1,a,e9.3e1,0p,'//
     &	      '2(a,e10.3e1,a,i3,a,i3,a))')
     &	      flag,'(',h,')',knt,avg,'+',dev,
     &	      ' mx=',amx,'(',imx,',',jmx,')',
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  endif

	else
	  if(atype.eq.'HGHT'.or.atype.eq.'STRM') then
	    write(lu,'(2a,1p,e10.3e1,0p,a,i6,i7,a,i6,'//
     &	      '2(a,i7,a,i3,a,i3,a))')
     &	      flag,'(',h,')',knt,nint(avg),'+',nint(dev),
     &	      ' mx=',nint(amx),'(',imx,',',jmx,')',
     &	      ' mi=',nint(amn),'(',imn,',',jmn,')'
	  elseif(atype.eq.'PRES'.or.atype.eq.'TEMP'.or.
     &	    atype.eq.'WIND'.or.atype.eq.'%REH'.or.
     &	    atype.eq.'RELH'.or.atype.eq.'MIXR') then
	    write(lu,'(2a,1p,e10.3e1,0p,a,i6,f7.1,a,f6.1,'//
     &	      '2(a,f7.1,a,i3,a,i3,a))')
     &	      flag,'(',h,')',knt,avg,'+',dev,
     &	      ' mx=',amx,'(',imx,',',jmx,')',
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  elseif(atype.eq.'NORM') then
	    write(lu,'(2a,1p,e10.3e1,0p,a,i6,f7.3,a,f6.3,'//
     &	      '2(a,f7.3,a,i3,a,i3,a))')
     &	      flag,'(',h,')',knt,avg,'+',dev,
     &	      ' mx=',amx,'(',imx,',',jmx,')',
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  else
	    write(lu,'(2a,1p,e10.3e1,0p,a,i6,1p,e10.3e1,a,e9.3e1,0p,'//
     &	      '2(a,e10.3e1,a,i3,a,i3,a))')
     &	      flag,'(',h,')',knt,avg,'+',dev,
     &	      ' mx=',amx,'(',imx,',',jmx,')',
     &	      ' mi=',amn,'(',imn,',',jmn,')'
	  endif

	endif

	end
