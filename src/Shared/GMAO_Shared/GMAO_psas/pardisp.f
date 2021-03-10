C***********************************************************************
C              subroutine pardisp
C***********************************************************************
C
C WRITTEN BY:  Jim Pfaendtner
C
C!DESCRIPTION: Display the parameter settings for current run
C
C CALLED FROM: analdrv
C
C SYSTEM ROUTINES USED: none
C
C SUBROUTINES CALLED: none
C
C!INPUT PARAMETERS:
C  c8time   - character*8 ASCII current time
C  c9date   - character*9 ASCII current date
C  cgname   - character*8 vector of length nbands with subroutine 
C             identifiers for the conjugate gradients stratification
C  criter   - real vector of length nbands with convergence criterion 
C             for the levels of the conjugate gradients stratification
C  delfile  - character*40 name of 'del' file to read
C  expid    - character*8 experiment identifier
C  idelprb  - integer index of begining innovation component to print
C  idelpre  - integer index of maximum innovation component to print
C  idelpri  - integer increment for innovation components to print
C  ktmax    - integer maximum number of data types
C  kxmax    - integer maximum number of data sources
C  lu       - integer unit number for parameter display
C  ludelf   - integer unit number for 'del' file
C  luverb   - integer unit number for more comprehensive output
C  maxpass  - integer vector of length nbands with maximum number of 
C             iterations at each level of the conjugate gradients 
C             stratification
C  minpass   - integer vector of length nbands with number of iterations
C             at each level of the conjugate gradients stratification
C             in which the convergence measure can increase
C  msmall   - integer size of diagonal blocks to be treated via direct
C             methods
C  nbands   - integer number of levels in the matrix stratification
C  nobsmax  - integer maximum number of observation
C  nhms     - integer GMT time of center of data window (hhmmss format)
C  ntwidth  - integer width of data window (hhmmss format)
C  nymd     - integer date (yymmdd format) of center of data window
C  seplim   - real vector of length nbands with the region separations
C             (in degrees) defining the conjugate gradients 
C             stratification
C  verbose  - logical switch to turn on more comprehensive output
C
C!OUTPUT PARAMETERS: none
C
C!INPUT/OUTPUT PARAMETERS: none
C
C FILES USED: none
C
C!REVISION HISTORY:
!  25jul96 - Jing Guo - modified some variables, minmax to minpass and
!			criter(:) to criter(:,2)
C  22jun94  - Jim Pf.  - Added header
C
C***********************************************************************

      subroutine pardisp( lu,
     1                    expid,   c9date, c8time,
     1                    nobsmax, kxmax,  ktmax,
     1                    verbose, luverb, idelprb, idelpre, idelpri,
     1                    delfile, ludelf, nymd,    nhms,    ntwidth,
     1                    nbands,  msmall,
     1                    cgname,  seplim, criter,  minpass, maxpass   )

c.......................................................................
c.... Argument declarations.

      integer       lu
      character*(*)   expid
      character*(*)   c9date
      character*(*)   c8time
      integer       nobsmax, kxmax,  ktmax
      logical       verbose
      integer       luverb
      integer       idelprb, idelpre, idelpri
      character*(*)  delfile
      integer       ludelf, nymd, nhms, ntwidth
      integer       nbands,  msmall
      character*(*)   cgname(nbands)
      real          seplim(nbands), criter(nbands,2)
      integer       minpass(nbands), maxpass(nbands)

C***********************************************************************

      write(lu,900)
      write(lu,902)
      write(lu,904)

      write(lu,900)
      write(lu,906) 'expid   = ', expid
      write(lu,907) 'c9date  = ', c9date
      write(lu,906) 'c8time  = ', c8time

      write(lu,900)
      write(lu,908) 'nobsmax = ', nobsmax
      write(lu,908) 'kxmax   = ', kxmax
      write(lu,908) 'ktmax   = ', ktmax

      write(lu,900)
      write(lu,910) 'verbose = ', verbose
      write(lu,908) 'luverb  = ', luverb
      write(lu,908) 'idelprb = ', idelprb
      write(lu,908) 'idelpre = ', idelpre
      write(lu,908) 'idelpri = ', idelpri

      write(lu,900)
      write(lu,912) 'delfile = ', delfile
      write(lu,908) 'ludelf  = ', ludelf
      write(lu,908) 'nymd    = ', nymd
      write(lu,908) 'nhms    = ', nhms
      write(lu,908) 'ntwidth = ', ntwidth

      write(lu,900)
      write(lu,908) 'nbands  = ', nbands
      write(lu,908) 'msmall  = ', msmall
      do 100 kband = 1, nbands
         write(lu,900)
         write(lu,906) 'cgname     = ', cgname(kband)
         write(lu,914) 'seplim     = ', seplim(kband)
         write(lu,914) 'max_criter = ', criter(kband,1)
         write(lu,914) 'min_criter = ', criter(kband,2)
         write(lu,908) 'minpass    = ', minpass(kband)
         write(lu,908) 'maxpass    = ', maxpass(kband)
  100 continue

      write(lu,900)
      write(lu,904)
      write(lu,900)

  900 format('  pardisp:')
  902 format('  pardisp: Run time parameters were set as follows:')
  904 format('  pardisp: ----------------------------------------')
  906 format('  pardisp: ',A10,2X,A8)
  907 format('  pardisp: ',A10,1X,A9)
  908 format('  pardisp: ',A10,I10)
  910 format('  pardisp: ',A10,L10)
  912 format('  pardisp: ',A10,A32)
  914 format('  pardisp: ',A10,1X,E13.5)

      return

C***********************************************************************
      end
