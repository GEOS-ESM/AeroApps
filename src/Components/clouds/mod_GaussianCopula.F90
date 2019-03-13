! file: mod_GaussianCopula.f90
! Gaussian Copula class

module mod_GaussianCopula

  use mod_mkl
  use mod_utils, only : myError
  implicit none

  private
  public GaussianCopula
  public GaussianCopula_create, GaussianCopula_destroy
  public GaussianCopula_generate
  public GaussianCopula_copy

  type GaussianCopula
    private
    integer :: ndim
    real*4, allocatable :: zeroMean (:)  !(ndim)
    real*4, allocatable :: uprCholFacOfCorrMat (:,:)  !(ndim,ndim)
  end type GaussianCopula

contains

  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine GaussianCopula_create (self, ndim, &
    upperCholeskyFactorOfValidCorrMatrixWithLowerTriZeroed)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    implicit none
    type (GaussianCopula), intent(out) :: self
    integer, intent(in) :: ndim
    real*4,  intent(in), dimension (ndim,ndim) :: &
      upperCholeskyFactorOfValidCorrMatrixWithLowerTriZeroed

    call GaussianCopula_destroy(self)
    self%ndim = ndim
    allocate(self%zeroMean(ndim)); self%zeroMean = 0.
    allocate(self%uprCholFacOfCorrMat(ndim,ndim))
    self%uprCholFacOfCorrMat = &
      upperCholeskyFactorOfValidCorrMatrixWithLowerTriZeroed

  end subroutine GaussianCopula_create

  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine GaussianCopula_destroy (self)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    implicit none
    type (GaussianCopula), intent(out) :: self
    self%ndim = -1
    if (allocated(self%zeroMean)) &
      deallocate(self%zeroMean)
    if (allocated(self%uprCholFacOfCorrMat)) &
      deallocate(self%uprCholFacOfCorrMat)
  end subroutine GaussianCopula_destroy


  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine GaussianCopula_copy (self, new)
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    implicit none
    type (GaussianCopula), intent(in)  :: self
    type (GaussianCopula), intent(out) :: new
    call GaussianCopula_create(new, self%ndim, self%uprCholFacOfCorrMat)
  end subroutine GaussianCopula_copy


  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  ! generate nsamples rank vectors from the GCOP
  ! the result has dimensions (nsamples, self%ndim)
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine GaussianCopula_generate ( &
    self, nsamples, brngStream, rank)

    use MKL_VSL_TYPE
    use MKL_VSL
    implicit none

    type (GaussianCopula), intent(in) :: self
    integer, intent(in) :: nsamples
    type (VSL_STREAM_STATE), intent(inout) :: brngStream
    real*4, intent(out) :: rank (nsamples, self%ndim)
     
    integer :: rc

    ! generate (nsamples, ndim) multi-variate Gaussian samples

    ! this method generates univariate Gaussian random numbers
    ! and does the Cholesky multiply explicitly. It seems to
    ! be a little faster than the vsrnggaussianmv() approach
    ! below. Note that for getting a (nsamples,ndim) matrix, 
    !   (LZ)^T = Z^T L^T = "randn"(nsamples,ndim) U.
    rc = vsrnggaussian( &
      VSL_RNG_METHOD_GAUSSIAN_ICDF, brngStream, &
      size(rank), rank, 0., 1.)
    if (rc /= VSL_STATUS_OK) &
      call myError('using vsrnggaussian')
    call strmm('R','U','N','N', &
      nsamples, self%ndim, 1., &
      self%uprCholFacOfCorrMat, self%ndim, &
      rank, nsamples)

!   ! method using vsRngGaussianMV
!   ! (see the VSL notes. It needs an uprChol in fortran,
!   ! due to the vsl code being in C, but it still produces
!   ! an (ndim,nsamples) matrix needing a transpose.
!   real*4 :: rnmvn (self%ndim, nsamples)
!   rc = vsrnggaussianmv( &
!     VSL_METHOD_SGAUSSIANMV_BOXMULLER2, brngStream, &
!     nsamples, rnmvn, self%ndim, &
!     VSL_MATRIX_STORAGE_FULL, &
!     self%zeroMean, self%uprCholFacOfCorrMat)
!   if (rc /= VSL_STATUS_OK) &
!     call myError('using vsrnggaussianmv')
!   rank = transpose(rnmvn)

    ! transform back to ranks with CDF of Normal
    call vscdfnorm(size(rank),rank,rank)

  end subroutine GaussianCopula_generate

end module mod_GaussianCopula
