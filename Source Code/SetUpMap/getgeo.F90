      subroutine getgeo(iread, labels, geo, xyz, lopt, na, nb, nc, int)  
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
!
      use parameters_C, only : ams
!
      use molkst_C, only : natoms, keywrd, numat, maxtxt, line, moperr
!
      use chanel_C, only : iw, ir, input_fn
!
      use common_arrays_C, only :atmass, simbol, txtatm, nat
!
!
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  10:47:17  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: iread 
      integer , intent(out) :: labels(*) 
      integer , intent(out) :: lopt(3,*)
      integer, intent(out)  :: na(*) 
      integer, intent(out)  :: nb(*) 
      integer, intent(out)  :: nc(*) 
      double precision, intent(out)  :: geo(3,*), xyz(3,*)
      logical :: int, lmop, solid, exists

!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(40) :: istart 
      integer :: i, icapa, icapz, iserr, k, icomma, khar, nvalue, label, j, &
      jj, ltl, max_atoms
      double precision :: weight, real, sum 
      logical :: lxyz, leadsp, ircdrc, saddle, lturn, mini
      character , dimension(107) :: elemnt*2 
      character :: space, nine, zero, comma, string*120, ele*2, turn, no
      double precision, external :: reada
      save elemnt, space, nine, zero, comma 
!-----------------------------------------------
!***********************************************************************
!
!   GETGEO READS IN THE GEOMETRY. THE ELEMENT IS SPECIFIED BY IT'S
!          CHEMICAL SYMBOL, OR, OPTIONALLY, BY IT'S ATOMIC NUMBER.
!
!  ON INPUT   IREAD  = CHANNEL NUMBER FOR READ, NORMALLY 5
!             AMS    = DEFAULT ATOMIC MASSES.
!
! ON OUTPUT LABELS = ATOMIC NUMBERS OF ALL ATOMS, INCLUDING DUMMIES.
!           GEO    = INTERNAL COORDINATES, IN ANGSTROMS, AND DEGREES.
!                    OR CARTESIAN COORDINATES, DEPENDING ON WHETHER
!                    KEYWORD ' XYZ' IS PRESENT
!           LOPT   = INTEGER ARRAY, A '1' MEANS OPTIMIZE THIS PARAMETER,
!                    '0' MEANS DO NOT OPTIMIZE, AND A '-1' LABELS THE
!                    REACTION COORDINATE.
!           NA     = INTEGER ARRAY OF ATOMS (SEE DATA INPUT)
!           NB     = INTEGER ARRAY OF ATOMS (SEE DATA INPUT)
!           NC     = INTEGER ARRAY OF ATOMS (SEE DATA INPUT)
!           ATMASS = ATOMIC MASSES OF ATOMS.
!***********************************************************************
      data (elemnt(i),i=1,107)/ 'H', 'HE', 'LI', 'BE', 'B', 'C', 'N', 'O', 'F'&
        , 'NE', 'NA', 'MG', 'AL', 'SI', 'P', 'S', 'CL', 'AR', 'K', 'CA', 'SC', &
        'TI', 'V', 'CR', 'MN', 'FE', 'CO', 'NI', 'CU', 'ZN', 'GA', 'GE', 'AS', &
        'SE', 'BR', 'KR', 'RB', 'SR', 'Y', 'ZR', 'NB', 'MO', 'TC', 'RU', 'RH', &
        'PD', 'AG', 'CD', 'IN', 'SN', 'SB', 'TE', 'I', 'XE', 'CS', 'BA', 'LA', &
        'CE', 'PR', 'ND', 'PM', 'SM', 'EU', 'GD', 'TB', 'DY', 'HO', 'ER', 'TM'&
        , 'YB', 'LU', 'HF', 'TA', 'W', 'RE', 'OS', 'IR', 'PT', 'AU', 'HG', 'TL'&
        , 'PB', 'BI', 'PO', 'AT', 'RN', 'FR', 'RA', 'AC', 'TH', 'PA', 'U', 'NP'&
        , 'PU', 'AM', 'CM', 'BK', 'MI', 'XX', '+T', '-T', 'CB', '++', '+', '--'&
        , '-', 'TV'/  
      data comma, space, nine, zero/ ',', ' ', '9', '0'/  
      ircdrc = index(keywrd,' IRC') + index(keywrd,' DRC') /= 0 
      if (index(keywrd," FORCETS") /= 0) ircdrc = .false.
      icapa = ichar('A') 
      icapz = ichar('Z') 
      lturn = .true.
      saddle = (index(keywrd, "SADDLE") > 0)
      lxyz = (index(keywrd, " XYZ") > 0 .or. saddle .or. (index(keywrd, " LOCATE-TS") /= 0))
      int = (index(keywrd, " INT ") > 0)
      lmop = (Index (keywrd, " MOPAC") /= 0)
      max_atoms = size(txtatm)
      maxtxt = 0 
      istart(1:1) = shape(simbol)
      if (istart(1) < natoms*3) then
        natoms = 0 
        numat = 0 
        return
      end if
      simbol(:natoms*3) = '---------' 
      natoms = 0 
      numat = 0 
     ! call dbreak()
      iserr = 0      
      mini = (index(keywrd, " MINI") /= 0)
   20 continue 
      read (iread, '(A241)', end=120, err=210) line 
      if (line == " ") goto 120
      ltl = len_trim(line)
      icomma = ichar(comma) 
      do i = 1, ltl 
        khar = ichar(line(i:i)) 
        if (khar == icomma .or. khar == 9) line(i:i) = space 
      end do 
      if (natoms == max_atoms) then
        write(iw,"(//10x,a)")" Maximum number of atoms exceeded"
        write(iw,"(10x,a,i6)")" Maximum allowed:", max_atoms
        call mopend("Maximum number of atoms exceeded")
        return
      end if
!
!   SEE IF TEXT IS ASSOCIATED WITH THIS ELEMENT
!
      i = index(line,'(') 
      if (i /= 0) then 
!
!  YES, ELEMENT IS LABELLED.
!
        k = index(line,')')
        if (k == 0) then
          write(iw,"(a,i5,a)")" Atom",natoms," has an opening parenthesis but no closing parenthesis"
          write(iw,"(a)")" Line :'"//line(:len_trim(line))//"'"
          write (iw, '(/,A)') ' GEOMETRY IS FAULTY.  GEOMETRY READ IN IS' 
          atmass(1) = -1.D0 
          nat(1) = 0
          call geout (iw) 
          call mopend ('GEOMETRY IS FAULTY') 
          return  
        end if
        if (line(k:k) == ")") then
          jj = k - 1
        else
          jj = k
        end if
        txtatm(natoms + 1) = line(i + 1:jj)         
        if (maxtxt == 27 .and. k - i - 1 /= 27) then
          if (index(keywrd, " GEO-OK") == 0) then
            i = natoms + 1
            no = char(Nint(log10(i + 1.1)) + ichar("1"))
            write(iw,'(//10x,a,i'//no//',a,i3,a)') "Atom label length of 27 detected, but atom ", natoms + 1, " has a label of", &
            k - i - 1, " characters."
            write(iw,'(/,a)')"Faulty line = '"//trim(line)//"'"
            write (iw,'(/10x,a)') "If this is intended, add keyword ""GEO-OK"""
            write(line,'(a,i'//no//')')"Fault detected in atom number ", natoms + 1
            call mopend (trim(line))
            stop
          end if
        end if
        maxtxt = max(maxtxt,k - i - 1) 
        if (maxtxt > 38) then
          if (index(line,"ATOM") + index(line,"HETATM") + index(line,"TITLE") + &
            index(line,"HEADER") + index(line,"ANISOU") + index(line,"COMPND") + &
            index(line,"SOURCE") + index(line,"KEYWDS") + index(line,"USER ")  + &
            index(line,"HELIX") + index(line,"SHEET") /= 0)  goto 70
          write(iw,"(a)")" Atom labels must not exceed 38 characters"
          string = " "
          write(iw,"(a)")string(:i)//"                1         2         3         4"
          write(iw,"(a)")string(:i)//"       1234567890123456789012345678901234567890"
          write(iw,"(a)")"Line :"""//line(:len_trim(line))//""""
          write (iw, '(/,A)') ' GEOMETRY IS FAULTY.  GEOMETRY READ IN IS' 
          atmass(1) = -1.D0 
          nat(1) = 0
          call geout (iw) 
          call mopend ('GEOMETRY IS FAULTY') 
          return  
        end if
        string = line(1:i-1)//line(k+1:) 
        line = string 
      else 
        txtatm(natoms + 1) = ' ' 
      endif 
!   CLEAN THE INPUT DATA
      call upcase (line, ltl) 
!
!   INITIALIZE ISTART TO INTERPRET BLANKS AS ZERO'S
      istart(:10) = len_trim(line) + 1
!
! FIND INITIAL DIGIT OF ALL NUMBERS, CHECK FOR LEADING SPACES FOLLOWED
!     BY A CHARACTER AND STORE IN ISTART
      leadsp = .TRUE. 
      nvalue = 0 
      i = 0
      do jj = 1, ltl 
        i = i + 1
        if (leadsp .and. line(i:i) /= space) then 
          nvalue = nvalue + 1 
          istart(nvalue) = i  
          if (i > 2) then
            if (line(i - 1:i - 1) ==  "-" .and. line(i - 2:i - 2) /=  " ") istart(nvalue) = i - 1
          endif 
          if (line(i:i) == '"') then
!
!  Run to other end of quoted text
!
            do j = 1, 40
              i = i + 1
              if (line(i:i) == '"') exit
            end do
          end if
        end if
!
!  set leadsp true if a space is detected, or if a "-" sign is found and it's part of a new number
!
        leadsp = (line(i:i) == space .or. (i /= istart(max(nvalue,1)) .and. line(i:i) ==  "-"))
        if (i > 1) then
          if ((line(i-1:i) == "D-" .or. line(i-1:i) == "E-")) leadsp = .false.
        end if
      end do 
      if (nvalue == 4) then !  Cartesian coordinates without optimization flags
!
! Check: is the fourth datum only 1 or 2 characters long?
!
        i = index(line(istart(4):), " ") + istart(4) - 2
        k = 0
        if (i > istart(4) + 1) k = 1
!
! Is is alphabetic?
!
        do j = istart(4), i
          if (line(j:j) < "A" .or. line(j:j) > "Z") k = 1        
        end do
!
! Is the first datum definitely not pure alphabetic?
!
        i = index(line(istart(1):), " ") + istart(1) - 2
        do j = istart(1), i
          if (line(j:j) >= "A" .and. line(j:j) <= "Z") exit        
        end do
        if (j > i .and. k == 0) then
          natoms = natoms + 1
          geo(1,natoms) = reada(line,istart(1)) 
          geo(2,natoms) = reada(line,istart(2)) 
          geo(3,natoms) = reada(line,istart(3)) 
          lopt(1,natoms) = 1
          lopt(2,natoms) = 1
          lopt(3,natoms) = 1 
          na(natoms) = 0
          nb(natoms) = 0
          nc(natoms) = 0       
          do i = 1, 107 
            if (line(istart(4): istart(4) + 1) /= elemnt(i)) cycle  
            labels(natoms) = i 
            atmass(natoms) = ams(i) 
            exit
          end do 
          goto 20
        end if
      end if
!
! ESTABLISH THE ELEMENT'S NAME AND ISOTOPE, CHECK FOR ERRORS OR E.O.DATA
!
      weight = 0.D0 
      string = line(istart(1):istart(2)-1) 
      if (string == "+3") string = "+T"
      if (string == "-3") string = "-T"
      if (string(1:1) >= zero .and. string(1:1) <= nine) then 
!  ATOMIC NUMBER USED: NO ISOTOPE ALLOWED
        label = nint(reada(string,1)) 
        if (label == 0) go to 120 
        if (label < 0 .or. label > 107) then 
          write (iw, '(''  ILLEGAL ATOMIC NUMBER'')') 
          go to 210 
        endif 
        go to 70 
      endif 
!  ATOMIC SYMBOL USED
      real = abs(reada(string,1)) 
      if (real < 1.D-15) then 
!   NO ISOTOPE
        ele = string(1:2) 
      else 
        weight = real 
        if (string(2:2) >= zero .and. string(2:2) <= nine) then 
          ele = string(1:1) 
        else 
          ele = string(1:2) 
        endif 
      endif 
!   CHECK FOR ERROR IN ATOMIC SYMBOL
      if (ele(1:1)=='-' .and. ele(2:2)/='-') ele(2:2) = ' ' 
      do i = 1, 107 
        if (ele /= elemnt(i)) cycle  
        label = i 
        go to 70 
      end do 
      if (ele(1:1) == 'X') then 
        label = 99 
        go to 70 
      endif 
      if (ele == "D ") then
        label = 1 
        weight = 2.014d0
        go to 70 
      else if (ele == "T ") then
        label = 1 
        weight = 3.016d0
        go to 70 
      else
         if (index(line,"ATOM") + index(line,"HETATM") + index(line,"TITLE") + index(line,"HEADER") + &
        index(line,"ANISOU") + index(line,"COMPND") + index(line,"SOURCE") + index(line,"KEYWDS") + &
        index(line,"HELIX") + index(line,"SHEET") + index(line,"REMARK") + index(line,"USER ")  + &
        index(line, "SEQRES") /= 0) goto 70
        if (trim(line) == " *                    *") return          
        write (iw, '(''  UNRECOGNIZED ELEMENT NAME: ('',A,'')'')') ele 
        write(iw,'(/,"  Faulty line: """,a,"""",/)')trim(line)
        go to 210 
      end if
!
! ALL O.K.
!
 70   continue 
      natoms = natoms + 1 
      nb(natoms) = 0
      nc(natoms) = 0
      if (label /= 99) numat = numat + 1 
      if (weight /= 0.D0) then 
        atmass(numat) = weight 
      else 
        if (label /= 99) atmass(numat) = ams(label) 
      endif 
      labels(natoms) = label 
      if (nvalue == 4) then !  Cartesian coordinates without optimization flags
        geo(1,natoms) = reada(line,istart(2)) 
        geo(2,natoms) = reada(line,istart(3)) 
        geo(3,natoms) = reada(line,istart(4)) 
        lopt(1,natoms) = 1
        lopt(2,natoms) = 1
        lopt(3,natoms) = 1 
      else
        geo(1,natoms) = reada(line,istart(2)) 
        geo(2,natoms) = reada(line,istart(4)) 
        geo(3,natoms) = reada(line,istart(6)) 
        if (.not. mini .and. ircdrc) then 
          turn = line(istart(3):istart(3)) 
          if (turn == 'T') then 
            lopt(1,natoms) = 1 
            if (lturn) then
              write (iw, '(A)') &
              ' IN DRC MONITOR POTENTIAL ENERGY TURNING POINTS' 
              lturn = .false.
            end if
          else 
            lopt(1,natoms) = 0 
          endif 
          turn = line(istart(5):istart(5)) 
          if (turn == 'T') then 
            lopt(2,natoms) = 1 
          else 
            lopt(2,natoms) = 0 
          endif 
          turn = line(istart(7):istart(7)) 
          if (turn == 'T') then 
            lopt(3,natoms) = 1 
          else 
            lopt(3,natoms) = 0 
          endif 
        else 
          lopt(1,natoms) = nint(reada(line,istart(3))) 
          lopt(2,natoms) = nint(reada(line,istart(5))) 
          lopt(3,natoms) = nint(reada(line,istart(7))) 
          do i = 3, 7, 2 
            if (.not.(ichar(line(istart(i):istart(i)))>=icapa .and. ichar(line(&
              istart(i):istart(i)))<=icapz .and. natoms>1)) cycle  
            iserr = 1 
          end do 
        endif 
      end if
      string = trim(line)
      if (moperr) then
        write(iw,'(/10x,a, i5)')" Error detected in definition of atom no.:", natoms
        write(iw,'(/,a)')" Text of faulty atom: '"//trim(elemnt(labels(natoms)))//&
        & "("//txtatm(natoms)//") "//trim(string(istart(2):))//"'"
        return
      end if
      j = len_trim(line)
      if (ltl /= j) then
        i = index(line(istart(8):), " ") + istart(8)
        nvalue = 8
        leadsp = .true.
        do i = i, j
          if (leadsp .and. line(i:i) /= space) then 
            nvalue = nvalue + 1 
            istart(nvalue) = i  
          end if
          leadsp = (line(i:i) == " ")
        end do 
      end if
      sum = reada(line,istart(8)) 
      i = index(line(istart(8):), " ") + istart(8) ! Find end of 8'th datum
      if (index(line(istart(8):i), ".") /= 0) sum = 0.D0 ! if 8'th datum contains a decimal, then it's not a connectivity
      na(natoms) = nint(sum) 
      if (natoms == 1) na(1) = 0
      if (lmop .and. natoms == 2) then
        na(2) = 1  
      else if (lmop .and. natoms == 3 .and. na(3) == 0) then
        na(3) = 2 
        nb(3) = 1 
        geo(2,3) = geo(2,3)*1.7453292519943D-02 
        geo(3,3) = 0.D0 
        lopt(3,3) = 0 
      else if (na(natoms) > 0) then
        nb(natoms) = nint(reada(line,istart(9))) 
        nc(natoms) = nint(reada(line,istart(10))) 
        geo(2:3,natoms) = geo(2:3,natoms)*1.7453292519943D-02 
      end if
      
!
!  SPECIAL CASE OF USERS FORGETTING TO ADD DIHEDRAL DATA FOR ATOM 3
!
      if (natoms == 3) then 
        if (lopt(1,3) /= 2 .and. lopt(2,3) /= 2 .and. lopt(3,3) == 2) then 
          na(3) = 1 
          nb(3) = 2
          geo(3,3) = 0.D0 
          lopt(3,3) = 0 
        else if (lopt(3,3)==1 .and. line(istart(6):istart(6) + 1) == "2 ") then 
          na(3) = 2 
          nb(3) = 1 
          geo(3,3) = 0.D0 
          geo(2,3) = geo(2,3)*1.7453292519943D-02 
          lopt(3,3) = 0 
        endif 
      endif 
      if (iserr == 1) then 
      if (index(line,"ATOM") + index(line,"HETATM") + index(line,"TITLE") + index(line,"HEADER") + &
        index(line,"ANISOU") + index(line,"COMPND") + index(line,"SOURCE") + index(line,"KEYWDS") + &
        index(line,"HELIX") + index(line,"SHEET") + index(line,"REMARK") + index(line,"USER ")  + &
        index(line, "SEQRES") /= 0) then
!
!  Geometry is definitely PDB
!
        natoms = -2
        return
      end if
!
!  MUST BE GAUSSIAN GEOMETRY INPUT
!
        do i = 2, natoms 
          do k = 1, 3 
            j = nint(geo(k,i)) 
            if (abs(geo(k,i)-j) <= 1.D-5) cycle  
!
!   GEOMETRY CANNOT BE GAUSSIAN
!
           write (iw, '(A)') ' GEOMETRY IS FAULTY.  GEOMETRY READ IN IS' 
            atmass(1) = -1.D0 
            nat(1) = 0
            call geout (iw) 
            call mopend ('GEOMETRY IS FAULTY') 
            return  
          end do 
        end do 
        natoms = -1 
        return  
      endif 
      go to 20 
!***********************************************************************
! ALL DATA READ IN, CLEAN UP AND RETURN
!***********************************************************************
  120 continue 
      if (maxtxt > 0) then
        do i = 1, natoms
          j = max(1, len_trim(txtatm(i)))
          if (j /= maxtxt) then           !  Pad out text with blanks
            txtatm(i)(j + 1:) = " "       !  (Removes rubbish like "null" characters)
          end if
        end do
      end if
      na(1) = 0
      nb(1) = 0
      nc(1) = 0 
      nb(2) = 0
      nc(2) = 0
      nc(3) = 0    
!
!     READ IN VELOCITY VECTOR, IF PRESENT
!
      solid = .false.
      do i = 1, natoms 
        if (labels(i) <= 0) then 
          write (iw, '('' ATOMIC NUMBER OF '',I3,'' ?'')') labels(i) 
          if (i == 1) then 
            write (iw, '(A)') ' THIS WAS THE FIRST ATOM' 
          else 
            write (iw, '(A)') &
              '    GEOMETRY UP TO, BUT NOT INCLUDING, THE FAULTY ATOM' 
            natoms = i - 1  
            call geout (iw) 
          endif  
          call mopend ('Error in READMO') 
          return  
        endif 
        if (labels(i) == 107) solid = .true.
        if (solid .and. labels(i) < 99) then
          call mopend ('TRANSLATION VECTORS MUST BE AT THE END OF THE DATA SET') 
          return
        end if
        if (i == 1 .or. na(i) == 0) cycle
        j = 0
        if (na(i) == nb(i).and. i > 1 .or. &
          (na(i) == nc(i) .or. nb(i) == nc(i)) .and. i > 2 .or. &
          nb(i)*nc(i) == 0 .and. i > 3) j = 1  !  Error condition
        if (na(i) >= i .or. nb(i) >= i .or. nc(i) >= i) then 
!
! An atom is being defined using a connectivity involving one or more atoms
! whose positions have not yet been defined.  This is likely to be an error.
!
          j = j + 1
!
!  Test: The atom is defined in internal coordinates, so check that
!  the dependent atoms are defined in Cartesian coordinates
! 
          if ((na(i) < i .or. na(max(1,na(i))) == 0) .and. & 
              (nb(i) < i .or. na(max(1,nb(i))) == 0) .and. &
              (nc(i) < i .or. na(max(1,nc(i))) == 0) ) j = j - 1 
        end if
        if (j == 0) cycle
        j = max(i, na(i), nb(i), nc(i))
        no = char(Nint(log10(j + 1.1)) + ichar("1"))
        write (line, '(" ATOM NUMBER ",I'//no//'," IS ILL-DEFINED")') i 
        write (iw, '(//10x, a)') trim(line) 
        write(iw,'(/10x," Connectivity of atom ",i'//no//',": NA=",i'//no//',", NB=",i'//no//',", NC=",i'//no//')') &
        i, na(i), nb(i), nc(i)
        if (na(i) > i .and. na(na(i)) /= 0) &
          write(iw,'(/10x,a)')" NA is defined using internal coordinates" 
        if (nb(i) > i .and. na(nb(i)) /= 0) &
          write(iw,'(/10x,a)')" NB is defined using internal coordinates" 
        if (nc(i) > i .and. na(nc(i)) /= 0) &
          write(iw,'(/10x,a)')" NC is defined using internal coordinates" 
          call web_message(iw,"geometry_specification.html")
        if (i == 1) then 
          return  
        endif 
        write (0, '(//10x, a)') trim(line)
        write (iw, '(/,''  GEOMETRY READ IN'',/)')
        nat = 0 
        call geout (iw) 
        call mopend (trim(line) )  
        inquire (file = input_fn, exist = exists)
        if (exists) close(ir, status = 'delete', err = 99) 
99      stop
      end do 
      if (natoms > 0) then
        call gmetry (geo, xyz) 
      else
        call mopend ("No atoms!")
      end if
      if (moperr) return
!
!  Switch for converting between coordinate systems.
!
      if (lxyz .or. int) then 
!
!    Coordinates should all be Cartesian or all be internal. 
!    First, unconditionally convert to Cartesian.
!
        k = 0
        do i = 1, natoms
          do j = 1, 3
            k = k + lopt(j,i)
          end do
        end do
!
!  Get rid of dummy atoms
!
        numat = 0 
        do i = 1, natoms 
          if (labels(i) /= 99) then 
            numat = numat + 1 
            labels(numat) = labels(i)
            txtatm(numat) = txtatm(i)
            lopt(:,numat) = lopt(:,i)
            na(numat) = 0
          endif 
          geo(:,i) = xyz(:,i) 
        end do 
!
!   If everything is marked for optimization then unconditionally mark the first
!   three atoms for optimization
!
        if (k >= 3*numat - 6) lopt(:,:min(3, numat)) = 1
        natoms = numat
        if (saddle .or. (index(keywrd, " LOCATE-TS") /= 0)) then
          lopt(:,:numat) = 1  ! In a saddle or locate-ts calculation, all parameters must be optimizable.
          na(:natoms) = 0
        end if
      end if 
      return  
! ERROR CONDITIONS
  210 continue 
      j = natoms - 1 
      write (iw, '('' DATA CURRENTLY READ IN ARE: '',/)') 
      do k = 1, j 
        if (na(k) > 0) geo(2:3,k) = geo(2:3,k)/1.7453292519943D-02 
        write (iw, '(3x,a,2x,3(f10.5,2x,i2,2x),3(i2,1x))') &
        elemnt(labels(k)), (geo(jj,k),lopt(jj,k),jj=1,3), &
        na(k), nb(k), nc(k) 
      end do 
      natoms = 0
      return
      end subroutine getgeo 
      subroutine web_message(channel, txt)
        implicit none
        character (*) ::txt
        integer :: channel
        write(channel,'(/10x,a,/)')"For more information, see: HTTP://OpenMOPAC.net/Manual/"//trim(txt)
        return
    end subroutine web_message
    subroutine getsym (locpar, idepfn, locdep, depmul)
    use molkst_C, only: natoms, ndep, id
    use chanel_C, only : iw, ir
    use common_arrays_C, only : na
    use molkst_C, only : line
   !
   !.. Implicit Declarations ..
    implicit none
   !
   !.. Formal Arguments ..
    integer, dimension (3*natoms), intent (inout) :: idepfn, locpar
    integer, dimension (3*natoms), intent (inout) :: locdep
    double precision, dimension (natoms), intent (out) :: depmul
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  10:47:17  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------


!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer , dimension(100) :: ivalue 
      integer :: n, nvalue, i, ll, j, l, nerror, n_used, i_sym(4)
      double precision :: sum
      double precision, dimension(100) :: value 
      character , dimension(19) :: texti*60, textx*60 
      character, dimension(19,2) :: text*60 
      character, dimension(38) :: used*60
      logical :: ok

      save text, i_sym
!-----------------------------------------------
      
      equivalence (text(1,1), texti), (text(1,2), textx) 
      data i_sym/3, 6, 8, 12/
      data texti/ &
        ' BOND LENGTH    IS SET EQUAL TO THE REFERENCE BOND LENGTH   ', &
        ' BOND ANGLE     IS SET EQUAL TO THE REFERENCE BOND ANGLE    ', &
        ' DIHEDRAL ANGLE IS SET EQUAL TO THE REFERENCE DIHEDRAL ANGLE', &
        ' DIHEDRAL ANGLE VARIES AS  90 DEGREES - REFERENCE DIHEDRAL  ', &
        ' DIHEDRAL ANGLE VARIES AS  90 DEGREES + REFERENCE DIHEDRAL  ', &
        ' DIHEDRAL ANGLE VARIES AS 120 DEGREES - REFERENCE DIHEDRAL  ', &
        ' DIHEDRAL ANGLE VARIES AS 120 DEGREES + REFERENCE DIHEDRAL  ', &
        ' DIHEDRAL ANGLE VARIES AS 180 DEGREES - REFERENCE DIHEDRAL  ', &
        ' DIHEDRAL ANGLE VARIES AS 180 DEGREES + REFERENCE DIHEDRAL  ', &
        ' DIHEDRAL ANGLE VARIES AS 240 DEGREES - REFERENCE DIHEDRAL  ', &
        ' DIHEDRAL ANGLE VARIES AS 240 DEGREES + REFERENCE DIHEDRAL  ', &
        ' DIHEDRAL ANGLE VARIES AS 270 DEGREES - REFERENCE DIHEDRAL  ', &
        ' DIHEDRAL ANGLE VARIES AS 270 DEGREES + REFERENCE DIHEDRAL  ', &
        ' DIHEDRAL ANGLE VARIES AS - REFERENCE DIHEDRAL              ', &
        ' BOND LENGTH VARIES AS HALF THE REFERENCE BOND LENGTH       ', &
        ' BOND ANGLE VARIES AS HALF THE REFERENCE BOND ANGLE         ', &
        ' BOND ANGLE VARIES AS 180 DEGREES - REFERENCE BOND ANGLE    ', &
        ' DO NOT USE - USE SYMMETY FUNCTION 19 INSTEAD               ', &
        ' BOND LENGTH IS A MULTIPLE OF THE REFERENCE BOND LENGTH     '/  
      data textx/ &
        ' X COORDINATE IS SET EQUAL TO   THE REFERENCE X COORDINATE  ', &
        ' Y COORDINATE IS SET EQUAL TO   THE REFERENCE Y COORDINATE  ', &
        ' Z COORDINATE IS SET EQUAL TO   THE REFERENCE Z COORDINATE  ', &
        ' X COORDINATE IS SET EQUAL TO - THE REFERENCE X COORDINATE  ', &
        ' Y COORDINATE IS SET EQUAL TO - THE REFERENCE Y COORDINATE  ', &
        ' Z COORDINATE IS SET EQUAL TO - THE REFERENCE Z COORDINATE  ', &
        ' X COORDINATE IS SET EQUAL TO   THE REFERENCE Y COORDINATE  ', &
        ' Y COORDINATE IS SET EQUAL TO   THE REFERENCE Z COORDINATE  ', &
        ' Z COORDINATE IS SET EQUAL TO   THE REFERENCE X COORDINATE  ', &
        ' X COORDINATE IS SET EQUAL TO - THE REFERENCE Y COORDINATE  ', &
        ' Y COORDINATE IS SET EQUAL TO - THE REFERENCE Z COORDINATE  ', &
        ' Z COORDINATE IS SET EQUAL TO - THE REFERENCE X COORDINATE  ', &
        ' X COORDINATE IS SET EQUAL TO   THE REFERENCE Z COORDINATE  ', &
        ' Y COORDINATE IS SET EQUAL TO   THE REFERENCE X COORDINATE  ', &
        ' Z COORDINATE IS SET EQUAL TO   THE REFERENCE Y COORDINATE  ', &
        ' X COORDINATE IS SET EQUAL TO - THE REFERENCE Z COORDINATE  ', &
        ' Y COORDINATE IS SET EQUAL TO - THE REFERENCE X COORDINATE  ', &
        ' Z COORDINATE IS SET EQUAL TO - THE REFERENCE Y COORDINATE  ', &
        ' NOT USED                                                   '/  
      n_used = 0
      nerror = 0
!
! TITLE OUTPUT
!
!
! INPUT SYMMETRY : FUNCTION, REFERANCE PARAMETER, AND DEPENDENT ATOMS
!
      n = 0 
      depmul(1) = 0.D0 
      ndep = 0
   20 continue 
      read (ir, '(A)', end=90) line 
      call upcase(line, len_trim(line))
      call nuchar (line, len_trim(line), value, nvalue) 
!   INTEGER VALUES
      do i = 1, nvalue 
        ivalue(i) = nint(value(i)) 
      end do 
!   FILL THE LOCDEP ARRAY
      if (nvalue==0 .or. Abs(value(3)) < 1.d-20) go to 90 
      if (ivalue(2) == 19) then
        if (na(ivalue(1)) == 0) then
            !
            !  Not allowed: a Cartesian coordinate cannot use function 19
            !
          write (iw,*) "Atom ", ivalue (1), " is Cartesian.  " // &
                     & "Function 19 cannot be used here."
          call mopend ("Error in Symmetry Data")
          return
        end if
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
        if (na(ivalue(1)) /= 0 .and. ivalue(2) == 18) then
            !
            !  Not allowed: an internal coordinate cannot use function 18
            !
          write (iw,*) "Atom ", ivalue (1), " is internal.  " // &
                     & "Function 18 cannot be used here."
          call mopend ("Error in Symmetry Data")
          return
        end if
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
      if (na(ivalue(1)) == 0) then
        i = 2
      else
        i = 1
      end if
      line = text(ivalue(2), i)
      do i = 1, n_used
        if (used(i) == line) exit
      end do
      if (i > n_used) then
        n_used = n_used + 1
        used(n_used) = line(:len_trim(line))
      end if
      if (nerror == 1) then
        call mopend("A SYMMETRY FUNCTION IS USED TO DEFINE A NON-EXISTENT ATOM")
        return
      end if
    goto 20
!
! CLEAN UP
   90 continue 
      do j = 1, 18 
        do i = 1, n_used 
          if (used(i) == texti(j)) go to 120 
        end do 
        cycle  
  120   continue 
      end do 
      do j = 1, 18 
        do i = 1, n_used 
          if (used(i) == textx(j)) go to 121 
        end do 
        cycle  
  121   continue 
      end do 
      write(iw,*)
      return  
  end subroutine getsym 
        subroutine nuchar(line, l_line, value, nvalue) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.4G  10:47:31  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: l_line
      integer , intent(out) :: nvalue 
      character  :: line*(*)
      double precision, intent(out) :: value(40) 
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


