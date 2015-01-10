!	(hfecQQ.h)
!  18Sep95  - Jing G.	- Created multiple copies of statistical tables
!  25Jul94  - Meta Sienkiewicz	- get moist tables based on cosine of
!				  separation angle, cos(dist/rade) 

!  The resolution is uniquely determined by a number, nQQtab.  Others
!  are derived form it, based on an assumption of an exp() function.

	integer   nQQtab	! ntab_fQQ
	parameter(nQQtab = 1800 +1)
	
	integer   nQQtb1
	parameter(nQQtb1=nQQtab/2)
	integer   nQQtb2
	parameter(nQQtb2=nQQtab-1-nQQtb1)

	real hfecQQ
	common /hfecQQ0/	hfecQQ(MXveclev,nQQtab)

	real QQbeg2,Qcoslim,qxQtb1,qxQtb2
	common /hfecQQ1/	QQbeg2,Qcoslim,qxQtb1,qxQtb2

!. end-of-include
