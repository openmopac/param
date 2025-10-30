      subroutine upcase(keywrd, n) 
!...Translated by Pacific-Sierra Research 77to90  4.4G  08:36:11  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: n 
      character , intent(inout) :: keywrd*(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: icapa, ilowa, ilowz, i, iline, j 
      character :: keybuf*2000 
!-----------------------------------------------
!
!  UPCASE WILL TAKE A CHARACTER STRING, IN KEYWRD, AND PUT IT INTO
!  UPPER CASE.  KEYWRD IS LIMITED TO 80 CHARACTERS
!
      icapa = ichar('A') 
      ilowa = ichar('a') 
      ilowz = ichar('z') 
      keybuf = keywrd 
      do i = 1, n 
        iline = ichar(keywrd(i:i)) 
        if (iline>=ilowa .and. iline<=ilowz) &
        keywrd(i:i) = char(iline + icapa - ilowa) 
        if (iline /= 9) cycle  
!
!  Change tabs to spaces.  A tab is ASCII character 9.
!
        keywrd(i:i) = ' ' 
      end do 
!
!   If the word EXTERNAL is present, do NOT change case of the following
!   character string.
!
      i = index(keywrd,'EXTERNAL=') 
      if (i /= 0) then 
        j = index(keywrd(i+1:),' ') + i 
        keywrd(i+9:j) = keybuf(i+9:j) 
      endif 
      return  
  end subroutine upcase 
  subroutine mopend(txt)  
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE molkst_C, only : errtxt 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character , intent(in) :: txt*(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!----------------------------------------------- 
      errtxt = txt 
      return  
  end subroutine mopend 
        subroutine wrttxt(iprt)  
      use molkst_C, only : koment, title, refkey, keywrd, line
      implicit none
      integer , intent(in) :: iprt 
      integer :: i, j, k
      logical :: l_chains = .false., l_start = .false.
!-----------------------------------------------
!
! Is CHAINS or START_RES present in the data set?
!
      do i = 1, 6
        if (index(refkey(i), " NULL") /= 0) exit
        line = " "//trim(refkey(i))
        call upcase(line, len_trim(line))
        if (.not. l_chains) l_chains = (index(line, " CHAINS") /= 0)
        if (.not. l_start)  l_start  = (index(line, " START_RES") /= 0)
      end do
!
!  Is CHAINS present in the keyword?
!
      i = index(keywrd," CHAINS")
      if(i /= 0 .and. .not. l_chains) then
!
!  CHAINS is present in keyword, but was not present in the data-set,
!  so add CHAINS keyword to refkey(1)
!
        j = index(keywrd(i + 7:), ")") + i + 7
        refkey(1) = keywrd(i:j)//trim(refkey(1))
      end if
!
!  Is START_RES present in the keyword?
!
      i = index(keywrd," START_RES")
      if(i /= 0 .and. .not.l_start) then
!
!  START_RES is present in keyword, but was not present in the data-set,
!  so add START_RES keyword to refkey(1)
!
        j = index(keywrd(i + 10:), ")") + i + 10
        refkey(1) = keywrd(i:j)//trim(refkey(1))
      end if   
      if (refkey(2) == " ") then
        i = index(refkey(1), " +")
        if (i /= 0) then
          refkey(1)(i:i + 1) = " "
          refkey(2) = " NULL"
        end if
      end if        
      k = 0
      do i = 1, 6
        if (index(refkey(i), " NULL") /= 0) exit
        if (index(refkey(i), " +") == 0) k = k + 1
        write(iprt,'(a)', iostat = j)trim(refkey(i))
        if (j /= 0) then
          call mopend("ERROR DETECTED WHILE TRYING TO WRITE KEYWORDS TO A FILE")
          return
        end if
      end do
      if (index(koment, " NULL") == 0 .and. k < 3) write (iprt, '(A)') trim(koment) 
      if (index(koment, " NULL") == 0 .and. k < 4) write (iprt, '(A)') trim(title) 
      return  
  end subroutine wrttxt 
  

