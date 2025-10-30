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
      use molkst_C, only : moperr
      implicit none
      character , intent(in) :: txt*(*) 
      moperr = .true.
      if (txt /= "dummy") return
      return  
  end subroutine mopend 
  subroutine tidy1(line, len)
  implicit none
  integer :: len
  character (len = len + 1) :: line
  integer :: i, j
!
!  Clean keyword line to remove unwanted key-words
!
    do i = 1, len + 1
      if (line(i:i) >= "a" .and. line(i:i) <= "z") line(i:i) = char(ichar(line(i:i)) + ichar("A") - ichar("a")) 
    end do
    i = index(line,"EXTERNAL")
    if (i > 0) line(i:i + index(line(i:)," ") - 1) = " "
    i = index(line,"OLDENS")
    if (i > 0) line(i:i + 6) = " "
    i = index(line,"DENOUT")
    if (i > 0) line(i:i + 6) = " "
    i = index(line,"RESTART")
    if (i > 0) line(i:i + 6) = " "
    i = index(line,"DEBUG")
    if (i > 0) line(i:i + 4) = " "
    i = index(line," PL")
    if (i > 0) line(i:i + 2) = " "
    i = index(line," XYZ")
    if (i > 0) line(i:i + 3) = " "
    i = index(line,"H= ")
    if (i > 0) line(i:i + 11) = " " 
    j = 0 
    do i = 1, len 
      j = j + 1
      if (line(i:i + 1) == "  ")then
        j = j - 1
      else
        line(j:j) = line(i:i)
      end if
    end do  
    line(j + 1:) = " "
end subroutine tidy1
  subroutine sub_rab(iels, niels) 
    use molkst_C, only : natoms
    use common_arrays_C, only : nat, labels
    implicit none
    integer :: iels(2,107), niels
    integer :: i, j, k
!
!  Work out the empirical formula
!
      niels= 0
      iels = 0
      j = 0
      do i = 1, natoms - 3
        if (labels(i) > 90) cycle
        j = j + 1
        nat(j) = labels(i)
        do k = 1, niels
          if (iels(1,k) == labels(i)) exit
        end do  
        if (k > niels) then
          niels = niels + 1
          iels(1,k) = labels(i)
        end if
        iels(2,k) = iels(2,k) + 1
      end do
    return
  end subroutine sub_rab
  
 
   subroutine nuchar(line, l_line, value, nvalue)
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: l_line
      integer , intent(out) :: nvalue 
      character  :: line*(*)
      double precision , intent(out) :: value(40) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(40) :: istart 
      integer :: i 
      logical :: leadsp 
      character :: tab, comma, space 
      double precision, external :: reada

      save comma, space 
!-----------------------------------------------
!***********************************************************************
!
!   NUCHAR  DETERMINS AND RETURNS THE REAL VALUES OF ALL NUMBERS
!           FOUND IN 'LINE'. ALL CONNECTED SUBSTRINGS ARE ASSUMED
!           TO CONTAIN NUMBERS
!   ON ENTRY LINE    = CHARACTER STRING
!   ON EXIT  VALUE   = ARRAY OF NVALUE REAL VALUES
!
!***********************************************************************
      data comma, space/ ',', ' '/  
      tab = char(9) 
!
! CLEAN OUT TABS AND COMMAS
!
      do i = 1, l_line 
        if (line(i:i)/=tab .and. line(i:i)/=comma) cycle  
        line(i:i) = space 
      end do 
!
! FIND INITIAL DIGIT OF ALL NUMBERS, CHECK FOR LEADING SPACES FOLLOWED
!     BY A CHARACTER
!
      leadsp = .TRUE. 
      nvalue = 0 
      do i = 1, l_line 
        if (leadsp .and. line(i:i)/=space) then 
          nvalue = nvalue + 1 
          istart(nvalue) = i 
        endif 
        leadsp = line(i:i) == space 
      end do 
!
! FILL NUMBER ARRAY
!
      do i = 1, nvalue 
        value(i) = reada(line,istart(i)) 
      end do 
      return  
  end subroutine nuchar 
  
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




  

