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

      implicit none
      character , intent(in) :: txt*(*) 
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
  


  subroutine getdata(natoms, gnorm, hof, rms, idata, cutoff, charge)
    implicit none
    integer :: natoms, idata
    double precision :: gnorm, hof, rms, cutoff, charge
    integer :: i, j, io_stat
    character :: line*400
    double precision, external :: reada
    gnorm = 0.d0
    hof = 0.d0
    rms = 0.d0
    rewind (idata)
    do i = 1, 40
      read(idata,'(a)', iostat = io_stat)line
      if (io_stat /= 0) return
      j = index(line, "HEAT OF")
      if (j /= 0) then
        hof = reada(line, j + 20)
      end if
    end  do
    rewind (idata)
    do i = 1, 40
      read(idata,'(a)')line
      j = index(line, "NT NO")
      if (j /= 0) then
        gnorm = reada(line, j + 15)
        exit
      end if
    end  do
    rewind (idata)
    do i = 1, 40
      read(idata,'(a)')line
      j = index(line, "RMS DIS")
      if (j /= 0) then
        rms = reada(line, j + 20)
        exit
      end if
    end  do
    rewind (idata)
    do i = 1, 40
      read(idata,'(a)')line
      j = index(line, "Empir")
      if (j /= 0) then
        j = index(line, " = ")
        natoms  = nint(reada(line, j + 3))
        exit
      end if
    end  do
    do i = 1, 40
      read(idata,'(a)')line
      j = index(line, "CUTOFF=")
      if (j /= 0) then
        j = index(line, "CUTOFF=")
        cutoff  = reada(line, j + 3)
      end if
      j = index(line, "CHARGE=")
      if (j /= 0) then
        j = index(line, "CHARGE=")
        charge  = reada(line, j + 3)
      end if
    end  do
  end subroutine getdata
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
    double precision function meci()
    meci = 0.d0
  return
  end function meci
  subroutine fock2
  end subroutine fock2
  subroutine to_screen
  end subroutine to_screen
  subroutine chrge_for_MOZYME
  end subroutine chrge_for_MOZYME

  

