subroutine tidy (a)
  !
  !.. Implicit Declarations ..
    implicit none
  !
  !.. Formal Arguments ..
    character (len=120), intent (inout) :: a
  !
  !.. Local Scalars ..
    logical :: first = .true.
    integer :: i, icapa, icapz, idiff, ismala, j
  !
  !.. Intrinsic Functions ..
    intrinsic Char, Ichar
  !
  ! ... Executable Statements ...
  !
    if (first) then
      icapa = Ichar ("A")
      ismala = Ichar ("a")
      icapz = Ichar ("Z")
      idiff = icapa - ismala
    end if
    do i = 1, 120
      j = Ichar (a(i:i))
      if (j >= icapa .and. j <= icapz) then
        a(i:i) = Char(j-idiff)
      end if
    end do
end subroutine tidy
