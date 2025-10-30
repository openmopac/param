      subroutine geout(mode1) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      use common_arrays_C, only : labels, na, nb, nc, geo, loc 
      USE molkst_C, ONLY: natoms, ndep, line

      USE symmetry_C, ONLY: depmul, locpar, idepfn, locdep 
      USE elemts_C, only : elemnt
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: mode1 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: mode, iprt, i, j, n, ia, ii, k, igui
      double precision :: degree, w, x
      logical :: cart
      character , dimension(3) :: q*2 
      character :: flag1*2, flag0*2, blank*80, fmt1*4, fmt2*4, fmt3*4
      mode = mode1 
      igui = -10 ! Set to impossible value
      flag1 = '+1' 
      flag0 = '+0' 
      iprt = abs(mode) 
      degree = 57.29577951308232D0       
      cart = .true.
      do i = 1, natoms
        if (na(i) > 0) cart = .false.
      end do
      fmt1  = "13.8"
      fmt2  = "13.7"
      fmt3  = "13.7"
      n = 1 
      ia = loc(1,1) 
      ii = 0 
      blank = " "
      do i = 1, natoms 
        do j = 1, 3 
          q(j) = flag0 
          if (ia /= i) cycle  
          if (j /= loc(2,n)) cycle  
          q(j) = flag1 
          n = n + 1 
          ia = loc(1,n) 
        end do         
        if (na(i) > 0) then
          w = geo(2,i)*degree 
          x = geo(3,i)*degree  
          if (i > 3) then
!
!  CONSTRAIN ANGLE TO DOMAIN 0 - 180 DEGREES
!
            w = w - aint(w/360.D0)*360.D0 
            if (w < -1.d-6) w = w + 360.D0 
            if (w > 180.000001D0) then 
              x = x + 180.D0 
              w = 360.D0 - w 
            endif 
!
!  CONSTRAIN DIHEDRAL TO DOMAIN -180 - 180 DEGREES
!
            x = x - aint(x/360.D0 + sign(0.5D0 - 1.D-9,x) - 1.D-9)*360.D0 
          end if
        else
          w = geo(2,i) 
          x = geo(3,i)  
        endif 
        blank = elemnt(labels(i))
        if (labels(i) /= 99 .and. labels(i) /= 107) ii = ii + 1 
        if (labels(i) == 0) cycle  
        if (na(i) == igui .or. cart) then 
          if (labels(i)/=99 .and. labels(i)/=107) then    
              write(line,'(a,f8.4)')blank(41:59 + k)
            write (iprt, '(1X,A,F'//fmt1//',1X,A2,F'//fmt2//',1X,A2,F'//fmt3//',1X, A2, A)') &
                blank(:j), geo(1,i), q(1), w, q(2), x, q(3),trim(line)
          else 
            write (iprt, '(1X,A,F'//fmt1//',1X,A2,F'//fmt2//',1X,A2,F'//fmt3//',1X,A2,a)') &
            blank(:j), geo(1,i), q(1), w, q(2), x, q(3), " "
          endif  
        else 
          if (labels(i) /= 107) then 
            write(line,'(3i6)')na(i), nb(i), nc(i)                
            write(line(len_trim(line) + 1:),'(a,f8.4)') blank(41:41 + k)
            write (iprt, &
              '(1X,A,F'//fmt1//',1X,A2,F'//fmt2//',1X,A2,F'//fmt3//',1X,A2,A)') &
              & blank(:j), geo(1,i), q(1), w, q(2), x, q(3), trim(line) 
          else 
            write (iprt, '(1X,A,F'//fmt1//',1X,A2,F'//fmt2//',1X,A2,F'//fmt3//',1X,A2,3I6)') &
              blank(:j), geo(1,i), q(1), w, q(2), x, q(3), na(i), nb(i), nc(i) 
          endif 
        endif 
      end do 
      write (iprt, *) 
      if (ndep /= 0) then  
!
!   OUTPUT SYMMETRY DATA.
!
        n = 1
        i = 1
        outer_loop: do
          j = i
          do
            if (j == ndep) exit outer_loop
             !
             !  Group together symmetry functions of the same type
             !  (same reference atom, same reference function, same multiplier,
             !   if function 18 or 19)
             !  (Maximum number of dependent atoms on a line: 9)
             !
            if (locpar(j) /= locpar(j+1) .or. idepfn(j) /= idepfn(j+1) .or. &
                 & j-i >= 9) exit
            if (idepfn(i) == 18 .or. idepfn(i) == 19) then
              if (Abs(depmul(n) - depmul(n+1)) > 1.d-10) exit
              n = n + 1
            end if
            j = j + 1
          end do
          if (idepfn(i) == 18 .or. idepfn(i) == 19) then
            write (iprt, "(I4,I3,F13.9,10I5)") locpar (i), idepfn (i), &
                 & depmul(n), (locdep(k), k=i, j)
            n = n + 1
          else
            write (iprt, "(I4,I3,10I5)") locpar (i), idepfn (i), (locdep(k), k=i, j)
          end if
          i = j + 1
        end do outer_loop
        if (idepfn(i) == 19 .or. idepfn(i) == 18) then
          write (iprt, "(I4,I3,F13.9,10I5)") locpar (i), idepfn (i), &
               & depmul(n), (locdep(k), k=i, j)
        else
          write (iprt, "(I4,I3,10I5)") locpar (i), idepfn (i), (locdep(k), k=i, j)
        end if
        write (iprt,*)
      end if
      end subroutine geout 
