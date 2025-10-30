subroutine mini_emp2 (xscale, nosort, base2, order)
    use common_elements, only : maxele
!
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
    double precision, intent (inout) :: xscale
    integer, dimension (107), intent (in) :: nosort
    double precision, dimension (90), intent (in) :: base2
  !
  !.. Local Scalars ..
    integer :: i, j
    double precision :: order
    double precision :: sum
  !
  !.. Local Arrays ..
    double precision, dimension (90) :: base
  !
    !**********************************************************************
    !
    !  REDUCE IT TO AN EMPIRICAL FORMULA
    !
      xscale = nosort(maxele+1)*1.d-20
    !
    !   WORK OUT A SEQUENCE-NUMBER
    !
    !   EACH MOLECULE IS GIVEN A SMALL BIAS SO AS TO ALLOW MOLECULES
    !   WITH THE SAME EMPIRICAL FORMULA TO BE DISTINGUISHED
    !
      if (order < 0.5d0) then
        sum = 4.d-20
        base = 0.d0
        base(1) = sum

          j = 1
          do i = 2, maxele 
            if (base2(i) < 0.1d0) cycle
            sum = sum * (base2(j)+1.d0)
            base(i) = sum
            j = i
          end do
          base(maxele+1) = 1.d-20
      end if
      xscale = xscale + order / 20000
      order = order + 1.d-20
      sum = 0
      do i = 1, maxele + 1
        if (nosort(i) /= 0) then
          xscale = xscale + nosort(i) * base(i) + sum
        end if
      end do
end subroutine mini_emp2
