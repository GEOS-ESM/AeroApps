!                                                    02/10/93 - bands.h
!.......................................................................
      integer nbands
      parameter  ( nbands = 5 )

      real		seplim
      real              cosseplim
      integer		maxpass
      real		criter
      character*8	cgname
      integer		msmall
      integer		minpass
      integer		method
      integer		nbandcg
      real		diagscal
      logical		cgverb

      common/bandscom0/seplim(nbands),criter(nbands,2),maxpass(nbands)
      common/bandscom1/cgname(nbands),msmall,         minpass(nbands)
      common/bandscom2/method,nbandcg,diagscal,       cgverb(nbands)
      common/bandscom3/cosseplim(nbands)


!.......................................................................
!.... end-of-include
