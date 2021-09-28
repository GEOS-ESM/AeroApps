module ErrorHandleMod

  integer, parameter :: MAXSTRING = 255

  integer, parameter :: SUCCESS =  0
  integer, parameter :: FAILURE = -1

contains

  logical function MY_RTRN(A,iam,line,rc)
    integer,           intent(IN ) :: A
    character*(*),     intent(IN ) :: iam
    integer,           intent(IN ) :: line
    integer, optional, intent(OUT) :: RC
    MY_RTRN = .true.
    if(A/=SUCCESS)print'(A40,I10)',Iam,line
    if(present(RC)) RC=A
  end function MY_RTRN

  logical function MY_VRFY(A,iam,line,rc)
    integer,           intent(IN ) :: A
    character*(*),     intent(IN ) :: iam
    integer,           intent(IN ) :: line
    integer, optional, intent(OUT) :: RC
    MY_VRFY = A/=SUCCESS
    if(MY_VRFY)then
      if(present(RC)) then
        print'(A40,I10)',Iam,line
        RC=A
      endif
    endif
  end function MY_VRFY

  logical function MY_ASRT(A,iam,line,rc)
    logical,           intent(IN ) :: A
    character*(*),     intent(IN ) :: iam
    integer,           intent(IN ) :: line
    integer, optional, intent(OUT) :: RC
    MY_ASRT = .not.A
    if(MY_ASRT)then
      if(present(RC))then
        print'(A40,I10)',Iam,LINE
        RC=FAILURE
      endif
    endif
  end function MY_ASRT

  logical function NC_VRFY(A,iam,line,rc)
    use netcdf
    integer,           intent(IN ) :: A
    character*(*),     intent(IN ) :: iam
    integer,           intent(IN ) :: line
    integer, optional, intent(OUT) :: RC
    NC_VRFY = (A /= nf90_noerr)
    if(NC_VRFY)then
      if(present(RC))then
        print *, trim(nf90_strerror(A))
        print'(A40,I10)',Iam,LINE
        RC=FAILURE
      endif
    endif
  end function NC_VRFY

end module ErrorHandleMod
