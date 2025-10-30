  subroutine getsym (l_sym, ir, locpar, idepfn, locdep, depmul) 
    use molkst_C, only: natoms, ndep, id, line
    implicit none
!
!  Read in all symmetry data, if any.
! 
    logical, intent (in) :: l_sym
    integer, intent (in) :: ir
    integer, dimension (3*natoms), intent (out) :: idepfn, locpar, locdep
    double precision, dimension (natoms), intent (out) :: depmul
    integer :: n, nvalue, i, ll, j, l, nerror, i_sym(4), ivalue(100)
    double precision :: sum, value(100)
    save i_sym
    data i_sym/3, 6, 8, 12/  
!
    ndep = 0
    if (.not. l_sym) return
    nerror = 0
    n = 0 
    depmul(1) = 0.D0 
  20 continue 
    read (ir, '(A)', end=90) line 
    call upcase(line, len_trim(line)) 
    call nuchar (line, len_trim(line), value, nvalue) 
    do i = 1, nvalue 
      ivalue(i) = nint(value(i)) 
    end do 
!   FILL THE LOCDEP ARRAY
    if (nvalue==0 .or. Abs(value(3)) < 1.d-20) go to 90 
    if (ivalue(2) == 19) then
      do i = 4, nvalue
        if (ivalue(i) == 0) exit
        ndep = ndep + 1
        locdep(ndep) = ivalue(i)
        locpar(ndep) = ivalue(1)
        idepfn(ndep) = 19
        n = n + 1
!
!  Check:  Is multiplier a square root of a rational ratio?
!
        sum = value(3) ** 2
        do ll = 1, 4
          l = i_sym(ll)
          j = Nint (sum*l)
          if (Abs (j-sum*l) < l*1.d-4) then
            value(3) = Sqrt ((1.d0*j)/l)
            exit
          end if
        end do
        depmul(n) = value(3)
      end do
    else
      do i = 3, nvalue
        if (ivalue(i) == 0) exit
        ndep = ndep + 1
        locdep(ndep) = ivalue(i)
        locpar(ndep) = ivalue(1)
        idepfn(ndep) = ivalue(2)
        if (ivalue(i) > natoms + id) then
          nerror = 1
        end if
      end do
    end if
    ll = i - 1
    if (nerror == 1) then
      call mopend("A SYMMETRY FUNCTION IS USED TO DEFINE A NON-EXISTENT ATOM")
      return
    end if
  goto 20
!
! CLEAN UP
  90 continue 
    return  
end subroutine getsym 
  