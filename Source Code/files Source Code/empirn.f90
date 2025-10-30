subroutine empirn (allmol_i, xscale, mode, nosort, TV, name)
    use common_molkst, only : numat, nat, ioffset
    use common_keywrd, only : keywrd
    use common_elements, only : numbrs, elemnt, inorg, lowel, maxele

!
!  On exit, nosort = maximum number of each type of element
!           xscale = unique number indicating priority of compound
! 
!
!  empirn puts the empirical formula into allmol_i(1:12).
!  The sequence used is either Cox & Pilcher or JANAF, depending on
!  the value of "inorg" (here, hard-wired to Cox & Pilcher)
!
  !
  !.. Implicit Declarations ..
    implicit none
  !
  !.. Formal Arguments ..
    character (len=110), intent (inout) :: allmol_i
    double precision, intent (inout) :: xscale
    character :: name*300
    integer, intent (in) :: mode
    integer, dimension (107), intent (inout) :: nosort
    logical :: tv
  !
  !.. Local Scalars ..
    logical :: first = .true., first1 = .true.
    character (len=10) :: swap
    integer :: i, j, l, Z, nel(107), k, m, mers(3)
    integer :: nc
    double precision :: order, sum
    save :: order, numz
  !
  !.. Local Arrays ..
    character (len=3), dimension (107) :: presnt
    integer, dimension (2000) :: numn
    integer, dimension (107) :: nos, ntype
    integer, dimension (107, 3) :: nmap, numz
    double precision, dimension (90) :: base, base2
    intrinsic Abs, Index, Nint
    double precision, external :: reada
  !
  !
  ! ... Executable Statements ...
  !
    if (mode == 1) then
      if (first) then
        first = .false.
        order = 1
        do i = 1, maxele
          base2(i) = 0
        end do
      !
      !  CONSTRUCT THE MAP TO DETERMINE THE ORDER IN WHICH THE ELEMENTS WILL
      !  BE PRINTED IN THE EMPIRIC FORMULA
      !
        do i = 1, 107
          do j = 1, 107
            if (elemnt(j) == lowel(i, 1)) go to 1000
          end do
          go to 1300
1000      nmap(i, 1) = j
          do j = 1, 107
            if (elemnt(j) == lowel(i, 2)) go to 1100
          end do
          go to 1200
1100      nmap(i, 2) = j
          do j = 1, 107
            if (elemnt(j) == lowel(i, 3)) go to 1101
          end do
          go to 1200
1101      nmap(i, 3) = j
        end do
        do j = 1, 3
          do i = 1, maxele
            if (lowel(i, j) (2:2) == " ") then
              numz(i, j) = 1
            else
              numz(i, j) = 2
            end if
          end do
        end do
        numn(1) = 0
        do i = 1, 199
          if (i < 10) then
            write(numbrs(i),"(i1)")i
          else if (i < 100) then
            write(numbrs(i),"(i2)")i
          else
            write(numbrs(i),"(i3)")i
            end if
        end do
        do i = 2, 199          
            numn(i) = len_trim(numbrs(i))
        end do
        go to 1400
1200    write (7,*) "   FAULT:", lowel (i, 2), i
        stop
1300    write (7,*) "   FAULT:", lowel (i, 1), i
        stop
      end if 
    !**********************************************************************
    !
    !  INITIALIZE ELEMENT ARRAYS
    !
1400  do i = 1, 107
        nos(i) = 0
      end do
    !
    !  READ FORMULA AND WORK OUT ELEMENTAL COMPOSITION
    !
      do i = 1, numat
        j = nat(i)
        nos(j) = nos(j) + 1
      end do
      nos(99) = 0
    !
    !
    !   SEQUENCE ELEMENTS TO LOOK NICE
    !
   
      do i = 1, maxele
        nosort(i) = nos(nmap(i, inorg))
      end do
      if (TV) then
!
!  Reduce to simplest empirical formula
!
        call upcase(keywrd, len_trim(keywrd))
        keywrd = " "//trim(keywrd)//" "
        if (index(keywrd, " MERS") /= 0)then
           j = index(keywrd," MERS")
           mers = 0
           k = 0
           i = Index (keywrd(j + 1:), " ") + j
           do l = 1, 3
             j = j + k
             if (l > 1 .and. k == 0) exit
             mers(l) = Nint (reada (keywrd(j:), 1))
             k = Index (keywrd(j:i), ",")
           end do
           i = mers(1)
           do l = 2, 3
             if (mers(l) == 0) exit
             i = i*mers(l)
           end do
        else
           mers = 0
        end if
        i = Index(keywrd, " Z=")
        if (i /= 0) then
          i = Nint(reada(keywrd,i))
          z = mers(1)*mers(2)*mers(3)*i
          if (Index (keywrd, " BCC") /= 0) z = z/2
        else
          nel = 0
          nel(:100) = nosort(:100)
          j = 0
          do i = 1, 100
            if (nel(i) > 0) then
              j = j + 1
              nel(j) = nel(i)
             end if
          end do
          k = 1000
          do i = 1, j
            if(nel(i) < k) k = nel(i)
          end do
!
!  k is the smallest number of atoms of any element in the formula
!
          do i = 1, 100
            m = 0
            do l = 1, j
              if (Abs((i*nel(l))/k - (i*1.d0*nel(l))/k) > 1.d-5) m = 1 
            end do
            if (m == 0) exit
          end do
!
!  Number of empirical units  = k/i
!
          z = k/i
        end if
        do i = 1, 106
          if (nosort(i) > 0) then
            nosort(i) = nosort(i)/z
            if (nosort(i) == 0) then
              write(*,*) "Z incorrect in file: '"//trim(name)//"'"
            end if
          end if
        end do
      end if
      do i = 1, maxele
        if (base2(i) < nosort(i)) then
          base2(i) = nosort(i)
        end if
      end do
      l = 0
      i = Index (keywrd, "CHARGE") + Index (keywrd, "charge")
      if (i /= 0) then
        sum = Nint (reada (keywrd, i))
        nosort(maxele+1) = -Nint (sum)
        if (base2(maxele+1) < Abs (nosort(maxele+1))) then
          base2(maxele+1) = Abs (nosort(maxele+1))
        end if
      else
        nosort(maxele+1) = 0
      end if
    else
    !**********************************************************************
    !
    !  REDUCE IT TO AN EMPIRICAL FORMULA
    !
      xscale = nosort(maxele+1)*1.d-300
    !
    !   WORK OUT A SEQUENCE-NUMBER
    !
    !   EACH MOLECULE IS GIVEN A SMALL BIAS SO AS TO ALLOW MOLECULES
    !   WITH THE SAME EMPIRICAL FORMULA TO BE DISTINGUISHED
    !
      if (first1) then
        write (7,*) " MAXIMUM NUMBERS OF EACH TYPE OF ELEMENT IN ANY REFERENCE DATA SET"
        write (7, "(4(10X,A2,I4))") (lowel(j, inorg), Nint(base2(j)), j=1, &
       & maxele)
        first1 = .false.
        !
        !  ASSUME RANGE OF CHARGE IS +/- 3
        !
          sum = 4.d-300
          base(maxele+1) = 1.d-300
          base(1) = sum + 1.d-300
          j = 1
          do i = 2, maxele
            if(base2(i) < 0.5d0) cycle
            sum = sum*(base2(j) + 1.d0)
            base(i) = sum
            j = i
          end do
          select case (inorg)
      case (1)
      write (7, "(/,'   SEQUENCE OF COMPOUNDS AS IN JANAF',&
               &' THERMOCHEMICAL TABLES',/)")      
      case (2)
      write (7, "(/,'   SEQUENCE OF COMPOUNDS AS IN COX',&
               &' & PILCHER THERMOCHEMICAL TABLES')")      
      case (3)
      write (7, "(/,'   SEQUENCE OF COMPOUNDS IN PERIODIC TABLE ORDER',/)")
      end select
       l = ioffset + 8
          write (7, "(' Emp. Formula               File name  ', &
          & "//numbrs(l)//"x,'    H.o.F. Dipole I.P. Elements')")
      end if
      xscale = xscale + order*1.d-300 / 20000
      if (Index(allmol_i,",") /= 0) xscale = xscale + 2.5*1.d-300 / 20000
      order = order + 1.d0
      
      if (nosort(maxele+1) /= 0) l = l - 1
    !#         XSCALE=XSCALE+XSTORE
    !#      WRITE(7,'(F23.5)')XSCALE
    !
    !   CONVERT IT TO A CHARACTER FORMULA
    !
      allmol_i (:9) = " "
      
      l = 0
      sum = 0
      do i = 1, maxele + 1
        if (nosort(i) /= 0) then
          xscale = xscale + nosort(i) * base(i) !+ sum
          sum = 0
          l = l + 1
          ntype(l) = i
          presnt(l) = lowel(i, inorg)
        else if (base2(i) > 1.d-1) then
          sum = sum + base(i) * (base2(i) + 1.d0)
        end if
      end do
!
!  ntype(i): atomic number of elements in formula, in order of printing
!
      nc = 2
      do i = 1, l
        if (nc > 28) exit
        if (ntype(i) > maxele) exit
        allmol_i(nc:28) = presnt(i) (:numz(ntype(i), inorg))
        nc = nc + numz(ntype(i), inorg)
        if (nosort(ntype(i)) /= 1) then
          allmol_i(nc:28) = numbrs(nosort(ntype(i)))
          nc = nc + numn(nosort(ntype(i)))
        end if
      end do
!
!  If first two elements are H and C, then reverse the order.
!
      if (allmol_i(2:2) == "H") then
        do i = 3,10
          if (allmol_i(i:i) < "0" .or. allmol_i(i:i) > "9") then
            if (allmol_i(i:i) == "C" .and. (allmol_i(i+1:i+1) < "a" .or. allmol_i(i+1:i+1) > "z")) then
              if (allmol_i(i-1:i-1) /= "H" .and. allmol_i(i-1:i-1) >= "0" .and. allmol_i(i-1:i-1) <= "9") then
                do j = i+1,10
                  if (allmol_i(j:j) < "0" .or. allmol_i(j:j) > "9") exit
                end do
  !
  ! Characters 2 - i-1 need to be swapped with characters i - j-1
  !
              swap = allmol_i(i:j-1)
              swap(j-i+1:) = allmol_i(2:i-1)
              allmol_i(2:j-1) = swap
              end if
            end if
          end if
        end do
      end if
      return
    end if
end subroutine empirn
