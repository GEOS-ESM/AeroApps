        module m_sigAdef

        implicit none


!       Define eigen-problem configuration
!       ----------------------------------
	character*2, save :: which
!	parameter ( which = 'SM' )	! smallest eigenvalues
!	parameter ( which = 'LM' )	! largest eigenvalues
!	parameter ( which = 'BE' )	! half largest/half smallest evals

        real   , parameter :: r_inc0   = 0.3	! percentage of nobs used as
						!   increments for no. of
						!   modes in eigendecomp.
        real   , parameter :: r_tol    = 1.e-2	! relative tolerance

!       Define eigen-problem case
!       -------------------------
	character, save ::   bmat	 ! define eigen-problem type in ARPACK sense
!	parameter ( bmat = 'I' )	 ! standard eigen-problem
!	parameter ( bmat = 'G' )	 ! generalized eigen-problem
 
	integer, save :: iparam(11)      ! array of specific options for ARPACK
	integer, save :: n_chk(3)        ! if not zero, check on eigen-decomposition
                                         ! (calculate residuals)
	real   , save ::   tol           ! relative accuracy of eigenvalues
        integer, save ::   perc_var      ! percentage variance captured in innov cov
	logical, save ::   rvec          ! calculate eigenvectors
	logical, save ::   eigchk        ! check on accuracy of eigen-decomposition
	logical, save ::   eig_corr      ! when .t. solve eigen-problem for corr matrix
        integer, save ::   njump         ! number of modes to skip between two 
                                         !  consecutive estimates of sigA

        character*4, save :: eigtype     ! type of eigen-problem to solve:
                                         !   lpr  - for Lar Peter Riishojgaard's proposal
                                         !   cohn - for Steve Cohn's generalized eigen-problem
                                         !   scrt - for R. Todling's simplification of S. Cohn's scheme

        real, parameter :: pma = 0.2     ! percentage of modes added by PfH'

        end module m_sigAdef
