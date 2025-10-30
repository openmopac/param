subroutine sort (a, ia, n, asort, mark, mend, iasort)
  !
  !.. Implicit Declarations ..
    implicit none
  !
  !.. Formal Arguments ..
    integer, intent (in) :: n
    integer, dimension (n), intent (out) :: ia, iasort, mark, mend
    double precision, dimension (n), intent (out) :: a, asort
  !
  !.. Local Scalars ..
    integer :: i, ii, iloop1, iloop2, ipull, j, m, nbins, nlim, nn
  !
  !.. Intrinsic Functions ..
    intrinsic Log, Min, Nint
  !
  ! ... Executable Statements ...
  !
  !
  !  GENERAL SORTING PROCEDURE.  SORT WILL 'SORT' THE CONTENTS OF THE
  !  ARRAY 'A' INTO INCREASING ORDER.  THE MAP WILL BE RETURNED IN
  !  ARRAY 'IA'.  ASORT, MARK, MEND, AND IASORT ARE WORKSPACE.
  !  THE TIME REQUIRED INCREASES AS N.LOG(N).
  !
  !  SET UP INITIAL MAP ARRAY (IA(I) = I)
  !
    m = 1
    do i = 1, n
      ia(i) = i
    end do
  !
  !  NLIM IS THE NUMBER OF PAIRWISE SORTS NECESSARY TO FULLY SEQUENCE 'A'
  !
    nlim = Nint (Log(n*1.d0)/Log(2.d0)+1.d0)
    do ii = 1, nlim
    !#      WRITE(6,*)' SORT',II
    !#      WRITE(6,*)(IA(I),A(I),I=1,N)
    !
    !  SORT BINS TOGETHER, TWO BINS AT A TIME.  THE NUMBER OF BINS IS EQUAL
    !  TO THE NUMBER OF BINS THAT ARE FULLY FILLED, PLUS ONE IF EVEN ONE
    !  NUMBER IS STILL LEFT OVER (THIS IS THE REASON FOR THE (M-1))
    !
      nbins = (n+m-1) / m
    !
    !   MARK HOLDS THE START ADDRESS OF EACH BIN
    !   MEND HOLDS THE END ADDRESS OF EACH BIN
    !
      do i = 1, nbins
        mark(i) = (i-1) * m + 1
        mend(i) = i * m
      end do
    !
    !  THE LAST BIN MIGHT BE PARTLY EMPTY
    !
      mend(nbins) = Min (n, mend(nbins))
    !
    !   GO THROUGH THE BINS IN A PAIRWISE MANNER
    !
      nn = 0
      do iloop1 = 1, nbins, 2
        iloop2 = iloop1 + 1
        if (iloop2 <= nbins) then
          do i = 1, m * 2
            if (mark(iloop1) <= mend(iloop1) .and. mark(iloop2) <= &
           & mend(iloop2)) then
            !
            !   OPTION 'A': BOTH BINS HAVE NUMBERS LEFT IN THEM
            !
            !
            !   DECIDE WHICH BIN TO PULL A NUMBER OUT OF
            !
              if (a(mark(iloop1)) > a(mark(iloop2))) then
                ipull = iloop2
              else
                ipull = iloop1
              end if
            !
            !   STORE THE PULLED NUMBER IN ASORT AND IASORT
            !
              nn = nn + 1
              asort(nn) = a(mark(ipull))
              iasort(nn) = ia(mark(ipull))
              mark(ipull) = mark(ipull) + 1
            else
              go to 1000
            end if
          end do
          cycle
1000      if (mark(iloop2) <= mend(iloop2)) then
            ipull = iloop2
          else
            ipull = iloop1
          end if
        !
        !   OPTION 'B': ONE OR OTHER BIN IS EMPTY
        !
        !
        !   THE REST OF THE REMAINING NUMBERS HAVE A 'FREE PASS'
        !  - NO SORT NECESSARY AS THE SECOND BIN IS EMPTY
        !   COPY ALL OF THE REMAINING NUMBERS INTO ASORT AND IASORT
        !
          do j = mark(ipull), mend(ipull)
            nn = nn + 1
            asort(nn) = a(j)
            iasort(nn) = ia(j)
          end do
        else
        !
        !   OPTION 'C': THE SECOND BIN NEVER HAD ANYTHING IN IT.
        !
        !
        !   A 'FREE PASS' - NO SORT NECESSARY AS THE SECOND BIN IS MISSING
        !
          do i = nn + 1, n
            asort(i) = a(i)
            iasort(i) = ia(i)
          end do
          nn = n
        end if
      end do
    !
    !  COPY SORTED NUMBERS BACK INTO 'A' AND 'IA'
    !
      do i = 1, nn
        a(i) = asort(i)
        ia(i) = iasort(i)
      end do
    !
    !  FINALLY, SET BIN SIZE FOR NEXT SORT
    !
      m = m + m
    end do
    return
end subroutine sort
