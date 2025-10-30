      subroutine getgeo(iread, labels, geo, xyz, lopt, na, nb, nc)  
!
!  Read in a geometry
!
      use molkst_C, only : natoms, keywrd, numat, line, moperr, loop
!
      use chanel_C, only : iw
!
      use common_arrays_C, only : nat 
!
      use sets_of_names_C, only : data_set_name
!
      integer, intent(in) :: iread 
      integer, intent(out) :: labels(*), lopt(3,*), na(*), nb(*), nc(*) 
      double precision, intent(out)  :: geo(3,*), xyz(3,*)
!
!  Local
!
      logical :: solid, leadsp  
      integer :: i, k, icomma, nvalue, label, j, &
      jj, ltl, ii, istart(40), maxtxt
      double precision :: real, sum 
      character :: space, nine, zero, comma, string*120, ele*2, no, elemnt(107)*2
      double precision, external :: reada
      save elemnt, space, nine, zero, comma 
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
      natoms = 0 
      numat = 0      
!
!  Read in atoms, one at a time.
!
      ii = 0
   20 continue 
      read (iread, '(A241)', end=120, err=210) line 
      if (line(1:1) == '*') go to 20
      if (line == ' ') then
        if(natoms == 0) then
!
!  Check:  Is this an ARC file?
!
          rewind (iread)
          sum = 0.d0
          do i = 1, 10000
            read (iread, '(A)', end=120, err=210) line 
            if (index(line, "HEAT OF FORMATION") > 0) sum = reada(line,20)
            if (index(line, "FINAL GEOMETRY OBTAINED") > 0) exit   
            if (index(line, "GEOMETRY IN CARTESIAN COORDINATE") > 0) exit 
            if (index(line, "GEOMETRY IN MOPAC Z-MATRIX") > 0) exit         
          end do
!
!  Yes, it's an arc file
! 
          ii = 0
          if (index(line, "FINAL GEOMETRY OBTAINED") > 0) then
            do
              read (iread, '(A)', end=120, err=210) line 
              if (index(line, " &") + index(line, " +") == 0 .and. line(1:1) /="*") ii = ii + 1
              if (ii == 3) then
                read (iread, '(A)', end=120, err=210) line 
                exit
              end if
            end do
          else
            natoms = -3
            return
          end if
        end if
        if (ii == 0) goto 120
        ii = 0
      end if
      ltl = len_trim(line)
      icomma = ichar(comma) 
      do i = 1, ltl 
        k = ichar(line(i:i)) 
        if (k == icomma .or. k == 9) line(i:i) = space 
      end do 
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
        string = line(1:i-1)//line(k+1:) 
        line = string 
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
            exit
          end do           
          goto 20
        end if
      end if
!
! ESTABLISH THE ELEMENT'S NAME, CHECK FOR ERRORS OR E.O.DATA
!
      string = line(istart(1):istart(2)-1) 
      real = abs(reada(string,1)) 
      if (real < 1.D-15) then 
!   NO ISOTOPE
        ele = string(1:2) 
      else 
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
      if (ele == "D ") label = 1 
!
! ALL O.K.
!
 70   continue 
      natoms = natoms + 1 
      nb(natoms) = 0
      nc(natoms) = 0
      if (label /= 99) numat = numat + 1 
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
        lopt(1,natoms) = nint(reada(line,istart(3))) 
        lopt(2,natoms) = nint(reada(line,istart(5))) 
        lopt(3,natoms) = nint(reada(line,istart(7)))   
      end if
      string = trim(line)
      if (moperr) then
        write(iw,'(/10x,a, i5)')" Error detected in definition of atom no.:", natoms
        write(iw,'(/,a)')" Text of faulty atom: '"//trim(elemnt(labels(natoms)))//trim(string(istart(2):))//"'"
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
      if (na(natoms) > 0) then
        nb(natoms) = nint(reada(line,istart(9))) 
        nc(natoms) = nint(reada(line,istart(10))) 
        geo(2:3,natoms) = geo(2:3,natoms)*1.7453292519943D-02 
      end if
      go to 20 
!***********************************************************************
! ALL DATA READ IN, CLEAN UP AND RETURN
!***********************************************************************
  120 continue 
      if (natoms == 0) return  
      na(1) = 0
      nb(1) = 0
      nc(1) = 0 
      nb(2) = 0
      nc(2) = 0
      nc(3) = 0    
      solid = .false.
      do i = 1, natoms 
        if (labels(i) <= 0) then 
          write(iw,'(a)')"In system """//trim(data_set_name(loop))//""""
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
        if (i == 1) then 
          return  
        endif 
        write (0, '(//10x, a)') trim(line)
        write (iw, '(/,''  GEOMETRY READ IN'',/)')
        nat = 0 
        call geout (iw) 
        call mopend (trim(line) )  
        stop
      end do 
      if (natoms > 0) then
        call gmetry (geo, xyz)         
      else
        call mopend ("No atoms!")
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
    