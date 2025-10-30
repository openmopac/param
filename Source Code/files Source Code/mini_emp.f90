subroutine mini_emp (nosort, base2, first)
    use common_elements, only : maxele
!
!  On exit, nosort = maximum number of each type of element
!           xscale = unique number indicating priority of compound
! 
!
!
  !
  !.. Implicit Declarations ..
    implicit none
  !
  !.. Formal Arguments ..
    integer, dimension (107), intent (in) :: nosort
    double precision, dimension (maxele), intent (inout) :: base2
    logical, intent (inout) :: first 
  !
  !.. Local Scalars ..
    integer :: i
  !
  !
  !
  ! ... Executable Statements ...
  !
      if (first) then
        first = .false.
        do i = 1, maxele
          base2(i) = 0
        end do
        go to 1400
      end if !**********************************************************************
    !
    !
1400    do i = 1, maxele
          if (base2(i) < nosort(i)) then
          base2(i) = nosort(i)
        end if
      end do
end subroutine mini_emp
