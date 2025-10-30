subroutine sd (hr, info, i, k)
  !
  !.. Implicit Declarations ..
    implicit none
  !
  !.. Formal Arguments ..
    character, intent (inout) :: hr
    character (len=80), intent (in) :: info
    integer, intent (in) :: i
    integer, intent (inout) :: k
  !
  !.. Local Scalars ..
    integer :: j
  !
  !.. Intrinsic Functions ..
    intrinsic Index
  !
  ! ... Executable Statements ...
  !
    j = Index (info(i+3:), " ") + i + 3
    if (Index (info(i:j), ",") == 0) then
      if (hr == "*") then
        hr = "@"
      else
        hr = "+"
      end if
    else
      k = k + 1
    end if
end subroutine sd
