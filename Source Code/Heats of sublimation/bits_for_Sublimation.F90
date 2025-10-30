 module Common_C
    character*100 :: mol_names(3000), sol_names(3000), all_mol_names(10000), all_sol_names(3000), &
      sorted_mol_names(3000), mol_temp_names(3000), sol_temp_names(3000)
    character*63 :: location_of_mols, location_of_solids, list_of_compounds, &
      header, one_line, footer
    data location_of_mols /"M:\PARAM\HTML Files\table_of_heats.html                        "/, &
      location_of_solids  /"M:\PARAM\HTML Files\Survey_of_Solids.html                      "/, &
      list_of_compounds   /"M:\PARAM\Source Code\Heats of sublimation\List of compounds.txt"/, &
      header              /"M:\PARAM\Source Code\Heats of sublimation\header.txt           "/, &
      footer              /"M:\PARAM\Source Code\Heats of sublimation\footer.txt           "/
    integer :: ir = 5, iw = 6, ir_mol = 7, ir_sol = 8
    double precision :: heat_of_sub(2000), all_PM6_D3H4_g(10000), all_PM7_g(10000), all_PM6_D3H4_s(3000), &
      all_PM7_s(3000), ref_heat_sub(3000), calc_sub_PM7(3000), calc_sub_PM6_D3H4(3000), sorted_ref_heat(3000)
    character :: refs(3000)*10    
  end module Common_C
  
  
  real(kind(0.0d0)) function reada (string, istart) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.4G  09:33:41  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: istart 
      character  :: string*(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i0, i9, idot, ineg, ipos, icapd, icape, ismld, ismle, l, i, &
        iadd, n, j 
      logical :: expnnt 
      double precision, external :: digit
!-----------------------------------------------
!     FORTRAN FUNCTION TO EXTRACT NUMBER FROM STRING
!
!
!     DEFINE ASCII VALUES OF NUMERIC FIELD CHARACTERS
      i0 = ichar('0') 
      i9 = ichar('9') 
      idot = ichar('.') 
      ineg = ichar('-') 
      ipos = ichar('+') 
      icapd = ichar('D') 
      icape = ichar('E') 
      ismld = ichar('d') 
      ismle = ichar('e') 
!
      l = len(string) 
!
!     FIND THE START OF THE NUMERIC FIELD
      do i = istart, l 
        iadd = 0 
        n = ichar(string(i:i)) 
!
!       SIGNAL START OF NUMERIC FIELD IF DIGIT FOUND
        if (n>=i0 .and. n<=i9) go to 20 
!
!       ACCOUNT FOR CONSECUTIVE SIGNS [- AND(OR) +]
        if (n==ineg .or. n==ipos) then 
          iadd = iadd + 1 
          if (i + iadd > l) go to 50 
          n = ichar(string(i+iadd:i+iadd)) 
          if (n>=i0 .and. n<=i9) go to 20 
        endif 
!
!       ACCOUNT FOR CONSECUTIVE DECIMAL POINTS (.)
        if (n /= idot) cycle  
        iadd = iadd + 1 
        if (i + iadd > l) go to 50 
        n = ichar(string(i+iadd:i+iadd)) 
        if (n>=i0 .and. n<=i9) go to 20 
      end do 
      go to 50 
!
!     FIND THE END OF THE NUMERIC FIELD
   20 continue 
      expnnt = .FALSE. 
      do j = i + 1, l 
        iadd = 0 
        n = ichar(string(j:j)) 
!
!       CONTINUE SEARCH FOR END IF DIGIT FOUND
        if (n>=i0 .and. n<=i9) cycle  
!
!       CONTINUE SEARCH FOR END IF SIGN FOUND AND EXPNNT TRUE
        if (n==ineg .or. n==ipos) then 
          if (.not.expnnt) go to 40 
          iadd = iadd + 1 
          if (j + iadd > l) go to 40 
          n = ichar(string(j+iadd:j+iadd)) 
          if (n>=i0 .and. n<=i9) cycle  
        endif 
        if (n == idot) then 
          iadd = iadd + 1 
          if (j + iadd > l) go to 40 
          n = ichar(string(j+iadd:j+iadd)) 
          if (n>=i0 .and. n<=i9) cycle  
          if (n==icape .or. n==ismle .or. n==icapd .or. n==ismld) cycle  
        endif 
        if (n==icape .or. n==ismle .or. n==icapd .or. n==ismld) then 
          if (expnnt) go to 40 
          expnnt = .TRUE. 
          cycle  
        endif 
        go to 40 
      end do 
      j = l + 1 
   40 continue 
      n = ichar(string(j-1:j-1)) 
      if (n==icape .or. n==ismle .or. n==icapd .or. n==ismld) j = j - 1 
!
!     FOUND THE END OF THE NUMERIC FIELD (IT RUNS 'I' THRU 'J-1')
      n = 0 
      n = n + index(string(i:j-1),'e') 
      n = n + index(string(i:j-1),'E') 
      n = n + index(string(i:j-1),'d') 
      n = n + index(string(i:j-1),'D') 
      if (n == 0) then 
        reada = digit(string(i:j-1),1) 
      else 
        reada = digit(string(:i+n-2),i)*1.D1**digit(string(:j-1),i+n) 
      endif 
      return  
!
!     DEFAULT VALUE RETURNED BECAUSE NO NUMERIC FIELD FOUND
   50 continue 
      reada = 0.D0 
      return  
  end function reada 
  
  !     ******************************************************************
      real(kind(0.0d0)) function digit (string, istart) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.4G  09:34:53  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: istart 
      character , intent(in) :: string*(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i0, i9, ineg, ipos, idot, ispc, l, idig, i, n, j 
      double precision :: c1, c2, deciml 
      logical :: sign 
!-----------------------------------------------
!     FORTRAN FUNCTION TO CONVERT NUMERIC FIELD TO DOUBLE PRECISION
!     NUMBER.  THE STRING IS ASSUMED TO BE CLEAN (NO INVALID DIGIT
!     OR CHARACTER COMBINATIONS FROM ISTART TO THE FIRST NONSPACE,
!     NONDIGIT, NONSIGN, AND NONDECIMAL POINT CHARACTER).
!
!
!     DEFINE ASCII VALUES OF NUMERIC FIELD CHARACTERS
      i0 = ichar('0') 
      i9 = ichar('9') 
      ineg = ichar('-') 
      ipos = ichar('+') 
      idot = ichar('.') 
      ispc = ichar(' ') 
!
      c1 = 0.D0 
      c2 = 0.D0 
      sign = .TRUE. 
      l = len(string) 
!
!     DETERMINE THE CONTRIBUTION TO THE NUMBER GREATER THAN ONE
      idig = 0 
      do i = istart, l 
        n = ichar(string(i:i)) 
        if (n>=i0 .and. n<=i9) then 
          idig = idig + 1 
          c1 = c1*1.D1 + n - i0 
        else if (n==ineg .or. n==ipos .or. n==ispc) then 
          if (n == ineg) sign = .FALSE. 
        else if (n == idot) then 
          exit  
        else 
          go to 40 
        endif 
      end do 
!
!     DETERMINE THE CONTRIBUTION TO THE NUMBER LESS THAN THAN ONE
      deciml = 1.D0 
      do j = i + 1, l 
        n = ichar(string(j:j)) 
        if (n>=i0 .and. n<=i9) then 
          deciml = deciml/1.D1 
          c2 = c2 + (n - i0)*deciml 
        else if (n /= ispc) then 
          exit  
        endif 
      end do 
!
!     PUT THE PIECES TOGETHER
   40 continue 
      digit = c1 + c2 
      if (.not.sign) digit = -digit 
      return  
      end function digit 

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
 
  integer function numb(value, n_decimals)
    double precision, intent (in) :: value
    integer, intent (in) :: n_decimals
    character :: num, text*(20)
    integer :: i
!
!  Work out how many characters are needed in order
!  to print the number correctly
!
    num = char(ichar("0") + n_decimals)
    write(text, '(f13.'//num//')')value
    do i = 1, 20
      if (text(i:i) /= " ") exit
    end do
    numb = 11 - i    
    return
  end function numb




  

