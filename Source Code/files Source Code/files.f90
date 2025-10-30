program files
!
!  Given a list of filenames, sort the list into a chemically-sensible order,
!  H - He -Li, etc.
!
!  The list can be solids, in which case the ARC file is used,
!  or molecules, in which case the MOP file is used
!  
!  
    use common_molkst, only : numat, nat, ioffset
    use common_keywrd, only : keywrd, all_keywords
    use common_elements, only : lowel, maxele, numbrs
  !
  !.. Implicit Declarations ..
    implicit none
  !
  !.. Local Scalars ..
    character :: ch, hr
    character (len=120) :: title
    character (len=100) :: atom
    character (len=300) :: name, name_store, info, line
    logical :: exists, fileok, first = .true., is_ref, is_root, xxx = .true., TV, safe, l_sort
    integer :: i, idipo, igeom, iheat, ii, iip, iloop, j, jdipo, jgeom, jheat, &
   & jip, k, kdipo, kheat, kip, l, ll, m, n, nrefs, io_stat, icapa, icapz, idiff, ismala
    double precision :: order, sum
  !
  !.. Local Arrays ..
    character (len=2), dimension (107) :: elemnt
    character (len=3), dimension (20) :: presnt
    character, allocatable :: refs(:)*400, mnemonic(:)*8, allmol(:,:)*300, allb(:,:)*300
    logical, dimension (107) :: lelem
    integer, dimension (50000) :: ia, iw1, iw2, iw3
    integer, dimension (107, 6) :: ndatas
    integer, dimension (107, 50000) :: nosort, nob
    double precision, dimension (50000) :: work1, xmols, xmolb
    integer, external :: iargc
    double precision, external :: reada
  !
  !.. External Functions ..
  !
  !.. Intrinsic Functions ..
    intrinsic Char, Ichar, Index, Max, Real
  !
  !.. Data Declarations ..
    data (elemnt(i), i=1, 107) / "H", "HE", "LI", "BE", "B", "C", "N", "O", &
         &"F", "NE", "NA", "MG", "AL", "SI", "P", "S", "CL", "AR", "K", &
         &"CA", "SC", "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", &
         &"GA", "GE", "AS", "SE", "BR", "KR", "RB", "SR", "Y", "ZR", "NB", &
         &"MO", "TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN", "SB", "TE", &
         &"I", "XE", "CS", "BA", "LA", "CE", "PR", "ND", "PM", "SM", "EU", &
         &"GD", "TB", "DY", "HO", "ER", "TM", "YB", "LU", "HF", "TA", "W", &
         &"RE", "OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO", "AT", &
         &"RN", "FR", "RA", "AC", "TH", "PA", "U", "NP", "PU", "AM", "CM", &
         &"BK", "CF", "XX", "FM", "MD", "NO", "++", "+", "--", "-", "TV" /
  !
  ! ... Executable Statements ...
  allocate (refs(10000), mnemonic(10000), allmol(2,50000), allb(2,50000))
  allmol = " "
  allb = " "
  !
  !       Read in journal references
  !
    ioffset=148 
    refs = " "
    i = iargc()
    name = " "
    if (i == 2) call getarg(2, name)
    if (i >  0) call getarg(1, name(len_trim(name) + 1:))
    call upcase(name, len_trim(name))
      
    TV = (index(name, "SOLID") /= 0)
    l_sort = (index(name, "NOSORT") == 0)
    safe =  (index(name, "SAFE") /= 0)
    if (.not. TV) then
    open(unit=13, file="../REFERENCES.txt")
      do i = 1, 10000
        read(13,"(a8)", iostat=io_stat)mnemonic(i)
        if (io_stat /= 0) exit
        if (first) then
          icapa = Ichar ("A")
          ismala = Ichar ("a")
          icapz = Ichar ("Z")
          idiff = icapa - ismala
        end if
        do k = 1, 8
          j = Ichar (mnemonic(i)(k:k))
          if (j >= icapa .and. j <= icapz) then
            mnemonic(i)(k:k) = Char(j-idiff)
          end if
        end do
        read(13,"(a120)", iostat=io_stat)refs(i)(:120)
        if (io_stat /= 0) then
          if (mnemonic(i) == " ") exit
          write(*,*)" An error in the REFERENCES file was detected in reference",i
          write(*,*)" The mnemonic was "//trim(mnemonic(i))
          exit
        end if
        read(13,"(a120)", iostat=io_stat)refs(i)(121:240)
        read(13,"(a120)", iostat=io_stat)refs(i)(241:360)
        do
          j = index(refs(i)(:len_trim(refs(i))),"  ") 
          if ( j /= 0) then
            refs(i)(j:) = refs(i)(j + 1:)
          else
            exit
          end if
        end do      
      end do
    endif
    nrefs = i - 1
  5 inquire (file="all.txt", exist=exists)
    if ( .not. exists) then
      write (*,*)
      write (*,*) " The input data file " // "all" // ".txt does not exist"
      write (*,*)
      stop
    end if
    open (unit=5, file="all"//".txt", status="OLD")
    rewind (5)
    open (unit=7, file="files.out")
    open (unit=11, file="..\files.txt")
    open (unit=8, file="bits.txt")
    open (unit=9, form="FORMATTED", status="SCRATCH")
    iheat = 0
    jheat = 0
    kheat = 0
    idipo = 0
    jdipo = 0
    kdipo = 0
    iip = 0
    jip = 0
    kip = 0
    igeom = 0
    iloop = 0
    jgeom = 0
!
!   Read in all the data sets, one at a time
!
    do
      !
      !   EXTRACT NAME OF FILE
      !
      read (5, "(A)", END=1500, err=10) name
10    continue
      if (name(1:1) == "*") cycle
      iloop = iloop + 1
      do i = 1, 15
        if (name(1:1) == " ") then
          name = name(2:)
        end if
      end do
      name_store = name
      do i = 1,300
        if (name(i:i) >= "a" .and. name(i:i) <= "z") &
        & name (i:i) = Char (Ichar (name(i:i))+Ichar ("A")-Ichar ("a"))
      end do
      i = Index (name, ".MOP") + Index (name, ".ARC")
      if (i == 0) then
        iloop = iloop - 1
      else
        name_store = name_store(1: i - 1)
!
!  Open a specific reference data set file
!
        open (14, status="UNKNOWN", err=1300, iostat=j, file=name)
        rewind (14)
        title = name(:i)
        if (Index (name, ".ARC") /= 0) then
          do 
            read (14, "(A)", END=1200, err=20) line
            if (index(line, " FINAL GEOMETRY OBTAINED") /= 0) exit
          end do
        end if
  19    read (14, "(A)", END=1200, err=20) keywrd
        if (keywrd(1:1) == "*") goto 19
        all_keywords(iloop) = keywrd
        read (14, "(A)", END=1200, err=20) name, info

20      continue
!
!   Information MUST be written in reverse order for the next step to work correctly
!
!   Order: 
!       (last)Elements present, 
!             I.P., 
!             dipole, 
!             H.o.F, 
!             system-name, 
!             data-set name
!             No. of elements present
!
!
!   Work out the types of elements in the molecule
!
        do i = 1, 107
          lelem(i) = .false.
        end do
        do i = 1, 1000
            read (14, "(A)", END=1100, err=30) atom
30        continue
          do j = 1, 8
            if (atom(j:j) /= " ") exit
          end do
          line = trim(atom(j:))
          atom = trim(line)
          call upcase(atom,2)
          if (atom(2:2) == "(") atom(2:2) = " "
          do j = 1, 107
            if (elemnt(j) (1:2) == atom(1:2)) go to 1000
          end do
          exit
1000      lelem(j) = .true.
          nat(i) = j
        end do        
        lelem(99) =.false.
        lelem(107) =.false.
1100    numat = i - 1
        call empirn (allmol(1,iloop), xmols(iloop), 1, nosort(1, iloop), TV, name_store)
        if (.not. TV .and. nat(numat) == 107) then
          TV = .true.
          goto 5        
        end if
        l = 0
        do i = 1, 106
          ndatas(i, 1) = 0
          if (lelem(i)) then
            l = l + 1
            ndatas(i, 1) = 1
            presnt(l) = lowel(i,3)
          end if
        end do
        i = ioffset + 41
        write (9, "('+',"//numbrs(i)//"X,1X,21A)") (presnt(i), i=1, l)
!
!  Convert all reference information into lower case
!
        call tidy (info)
        if (.not. TV) then
!
!            Test for Geometric Data
!
          fileok = (Index (info, " georef")+Index (info, " gr=") /= 0)
          i = Index (info, " georef") + Index (info, " gr=")
          if (i /= 0) then
          !
          !  COUNT THE NUMBER OF GEOMETRIC REFERENCES
          !
            m = 1
            ll = 0
            do l = 1, 10
              ii = Index (info(m:), "<")
              if (ii == 0) exit
              m = m + ii
              ll = ll + 1
            !
            !  DON'T COUNT NULL GEOMETRIC PARAMETERS
            !
              if (info(m:m) == " ") then
                ll = ll - 1
              end if
              if (info(m:m) == ">") then
                ll = ll - 1
              end if
            end do
            do k = 1, 107
              ndatas(k, 5) = ndatas(k, 5) + ndatas(k, 1) * ll
            end do
            igeom = igeom + 1
            hr = " "
            if (Index (info, "gr=") == 0) then
              hr = "*"
              jgeom = jgeom + 1
            end if
            i = ioffset + 36
            write (9, "('+',"//numbrs(i)//"X,'YES',A)") hr
          end if
  !
  !            Test for Ionization Potential
  !
          i = Index (info, " i")
          ch = info(i+2:i+2)
          if (i /= 0 .and. (ch == "a" .or. ch == "s" .or. ch == "e" &
               &.or. ch > "0" .and. ch < "9" .or. ch == "=")) then
            fileok = .true.
            iip = iip + 1
            do k = 1, 107
              ndatas(k, 3) = ndatas(k, 3) + ndatas(k, 1)
            end do
            hr = " "   
            if (Index (info, "ir=") == 0) then
              hr = "*"
              jip = jip + 1
            end if
            call sd (hr, info, i, kip)
            j = ioffset + 35
            write (9, "('+',"//numbrs(j)//"X,F5.1,A)") reada (info, i+3), hr
  !
  !  Used 100 characters so-far
  !  
          end if
  
  !
  !            Test for Dipole Moment
  !
          i = Index (info, "d=")
          if (i /= 0) then
            do k = 1, 107
              ndatas(k, 4) = ndatas(k, 4) + ndatas(k, 1)
            end do
            idipo = idipo + 1
            fileok = .true.
            hr = " "
            if (Index (info, "dr=") == 0) then
              hr = "*"
              jdipo = jdipo + 1
            end if
            call sd (hr, info, i, kdipo)
             j = ioffset + 29
            write (9, "('+',"//numbrs(j)//"X,F5.1,A)") reada (info, i), hr
  !
  !  Used 84 characters so-far
  !  
          end if
  
  !
  !            Test for Heat of Formation
  !
      
          i = Index (info, "h=") + index(info, "pka")
          if (i /= 0) then
            if (safe .and. &
            index(info, '+"') + index(info, 'hr=ref') /= 0) then
              iloop = iloop - 1
              cycle
            end if
            do k = 1, 107
              ndatas(k, 2) = ndatas(k, 2) + ndatas(k, 1)
            end do
            iheat = iheat + 1
            fileok = .true.
            hr = " "
            if (Index (info, "hr=") == 0) then
              hr = "*"
              jheat = jheat + 1
            end if
            call sd (hr, info, i, kheat)
            j = ioffset + 21
            write (9, "('+',"//numbrs(j)//"X,F7.1,A)") reada (info, i), hr
  !
  !  Used 88 characters so-far
  !  
          end if
        end if  ! for Tv
        i = Index (info, "pm6")
        if (i /= 0) then
          fileok = .true.
        end if
!
! Name of system in data set can be  long
!
        write (9, "('+',A,1X,A)")trim(name_store)
        do i = 1, 15
          if (name(1:1) == " ") then
            name = name(2:)
          end if
        end do
        allmol(2,iloop) = name
!
!  Used 81 characters so-far
!
        if ((iloop/10)*10 == iloop) then
          write (9, "(1X)")
        end if
        write (9, "(1X)")
   

        rewind (9)
        i = Index (info, " type")
        ch = " "
        if (i /= 0) then
          ch = info(i+6:i+6)
        end if
        info = " "
        name = " "
        if(iloop == 700) then
        iloop = 700
      end if
      do 40 j = 1, 10
        read (9, "(A)", END=40, err=50) info
50      continue
        if (info(1:1) == "+") then
!
!   Write reference data, starting from the RIGHT hand side, writing to the LEFT side.
!
          do i = 300, 3, -1
            if (info(i:i) /= " ") exit
          end do
          name(2:i) = info(2:i)
        else if (fileok .or. TV) then
          if (name /= " " .and. Index (name, "title") == 0 .and. &
               & name(1:1) /= "1") then
            call tidy(keywrd)
            allmol(1,iloop) (28:) = trim(name)
            allmol(1,iloop) (27:27) = ch
          end if
          name = info
        else
          write (7, "(A)") " NO REFERENCE DATA:" // trim(title)
          name = info
        end if
40    continue
      rewind (9)
      cycle
1200  write (7,*) " KEYWRD, NAME, OR INFO MISSING:" // title
      go to 1400
1300  write (7,*) " FILE DOES NOT EXIST:" // trim(name (:i))
1400  rewind (9)
      iloop = iloop - 1
    end if     
  end do
!
!  At this point, all files have been read in.
!
1500 iloop = iloop - 1
    iheat = Max (1, iheat)
    idipo = Max (1, idipo)
    iip = Max (1, iip)
    igeom = Max (1, igeom)
    jheat = iheat - jheat
    jdipo = idipo - jdipo
    jip = iip - jip
    jgeom = igeom - jgeom
  !
  !  PRINT NUMBER OF REFERENCE DATA PER ELEMENT
  !
    write (7,*) "     ELEMENT  TOT  H.O.F.    I.P.   DIPOLE   GEOMETRY"
    do i = 1, 97
      ndatas(i, 6) = ndatas(i, 2) + ndatas(i, 3) + ndatas(i, 4) + ndatas(i, 5)
    end do
    do i = 1, 97
      k = ndatas(i, 6)
      l = i
      do j = 1, 97
        if (ndatas(j, 6) > k) then
          k = ndatas(j, 6)
          l = j
        end if
      end do
      if (k == 0) exit
      write (7, "(I4,4X,A2,I7, I7,4I9)") i, elemnt (l), k, (ndatas(l, j), &
           &j=2, 5)
      do j = 1, 6
        ndatas(l, j) = 0
      end do
    end do
    if (iloop == 0) then
      write (*,*) " NO REFERENCE DATA FOUND!"
      stop
    end if
!
! At this point, allmol contains data on individual molecules, unsorted, in the format:
!
!  (1:18) - unused
!  (20:)  - Name of the data-set
!
    do i = 1, iloop
      keywrd = all_keywords(i)
      call empirn (allmol(1,i), xmols(i), 2, nosort(1, i), TV, name_store)
    end do
    if (.not. l_sort) then
      allb(:,:iloop) = allmol(:,:iloop)
      xmolb(:iloop) = xmols(:iloop)
      goto 99
    end if
!
! At this point, allmol contains data on individual molecules, unsorted, in the format:
!
!  (1:18) - empirical formula
!  (20:)  - Name of the data-set
!
!
!  SORT BY FORMULA
!
!
    call sort (xmols, ia, iloop, work1, iw1, iw2, iw3)
!
!   Put compounds into the correct sequence
!
    k = 0
    do j = 1, iloop
      i = ia(j)
      if (Index (allmol(1,i) (19:), " ") /= 0) then
        do n = 120, 2, -1
          if (allmol(1,i) (n:n) /= " ") exit
        end do
        if (allmol(1,i) (18:) /= " ") then
          k = k + 1
          allb(:,k) = allmol(:,i) 
          xmolb(k) = xmols(j)
          do l = 1,90
          nob(l,k) = nosort(l,i)
          end do
        end if
      end if
    end do
!
! At this point, allb(1,i)  = text, in sorted sequence
!                nob(l,i) = number of elements of atomic number l in system allb(i)
!                xmolb(i) = number defined by empirical formula of nob(i)
!
    allmol(:,1:k) = allb(:,1:k)
    xmols(1:k) = xmolb(1:k)
    nosort(1:90,1:k) = nob(1:90,1:k)                
!
! Check: do 2 or more compounds have the same priority
!
    iloop = k
    i = 1
    do j = 1, iloop
      i = Max(i,j)
      if (j /= i) cycle
      if (j > 6950) then 
        ii = j
      end if
      do i = j + 1, iloop
        sum = Abs(xmolb(i) - xmolb(j))/xmolb(j)
        if( sum > 1.d-22) exit
      end do
      i = i - 1
      if(i /= j) then
!
!  Two or more molecules have the same flags
!
!
!  Work out "local" empirical formula
!
        first = .true.
        do k = j,i
          call mini_emp (nosort(1, k), work1, first)
        end do
        work1(maxele+1) = 4
!
!  Re-prioritize using local formula
!
        order=0
        if (j >= 8142 .and. xxx) then
          xxx = .false.
        end if
        do k = j,i
            call mini_emp2 (xmols(k), nosort(1, k), work1, order)
        end do
!
!  work1 contains the maximum numbers of each element in compounds j-i
!
!
!  SORT BY FORMULA
!
!
        call sort (xmols(j), ia(j), i-j+1, work1, iw1, iw2, iw3)
!

        do k = j,i
          allb(:,k) = allmol(:,ia(k)+j-1)
          xmolb(k) = xmols(ia(k)+j-1)
        end do
      end if
    end do
    if (TV) then
      do j = 1, iloop
        write(8, '(a)')trim(allmol(1,j)(29:))
      end do
      stop
    end if  
99  continue 
    do j = 1, iloop
      do i = 300,2,-1
        if (allb(1,j)(i:i) .ne. " ") exit
      end do
      k = ioffset + 70
      if (allb(1,j)(k:k + 10) /= " ") write(7,"(a)")allb(1,j)(1:i)
      do i = 198,2,-1
        if (allb(1,j)(i:i) .ne. " ") exit
      end do
      if (ichar(allb(1,j)(29:29)) > 32 .and. ichar(allb(1,j)(29:29)) < 123) then
        write(8,"(a)")allb(1,j)(29:i)
      else
        continue
      end if
      write(11,"(a)")allb(1,j)(2:19)
      write(11,"(A)")allb(1,j)(21:20 + len_trim(allb(1,j)(21:198)))
!
!  Read reference data
!
      info = allb(1,j)(21:20+len_trim(allb(1,j)(21:190)))//".mop"     
 98     open (14, status="OLD", file=info, iostat = l)
        if ( l /= 0) then
  !
  !  The file exists, but is not currently accessible
  !
            write(*,*) "Problem with", info(:len_trim(info))
            call sleep(3)
            goto 98
        end if
        rewind (14)
     if (.not. TV) then
       do
        read (14, "(A)", iostat = i) keywrd
        if (i /= 0 .or. keywrd(1:1) /= "*") exit
       end do
       do
        read (14, "(A)", iostat = i) name
        if (i /= 0 .or. name(1:1) /= "*") exit
       end do
       do
        read (14, "(A)", iostat = i) info
        if (i /= 0 .or. info(1:1) /= "*") exit
       end do
        keywrd = trim(info)
        call tidy (info)
        if (index(allb(1,j), "YES") /= 0) then
          allb(2,j)(len_trim(allb(2,j)) + 1:) = " (Geo)"
          i = index(info, " gr=")
          if (i /= 0) then
            name = info(i + 4:index(info(i + 5:)," ") + i + 4)
            do i = 1, nrefs
              if (index(name(:8),mnemonic(i)(2:)) /= 0) exit
            end do
            if (i <= nrefs) then
                write(11,"(2a, f7.1,5a)")"REF: <B> REF: </B>",refs(i)(:len_trim(refs(i)))
            end if
          end if
        else
  !
  !  Check for Heat of Formation
  !
          is_ref = .false.
          is_root = .false.
          i = index(info," h=")
          if (i /= 0) then
            sum = reada(info, i+3)
            i = index(info, " hr=")
            if (i /= 0) then
              name = info(i + 4:index(info(i + 5:)," ") + i + 4)
              do i = 1, nrefs
                if (index(name(:8),mnemonic(i)(2:)) /= 0) exit
              end do
              is_ref = (i <= nrefs) 
            end if
            l = index(info, " root=")
            if (l /= 0) then
              is_root = .true.
              k = index(info(l + 6:)," ") + l + 5
            end if 
            if (is_root) then
              write(11,"(5a)")"REF: For electronic state ", keywrd(l+6:k)
            end if
            if (is_ref) then
              write(11,"(2a, f7.1,5a)")"REF: ", '<B><span class="b">&Delta;</span>H<sub>f</sub>:</B>',&
              sum," kcal/mol, <b>&nbsp;&nbsp;&nbsp; REF: </b>",refs(i)(:len_trim(refs(i)))
            else
              write(11,"(2a, f7.1,5a)")"REF: ", '<B><span class="b">&Delta;</span>H<sub>f</sub>:</B>',&
              sum," kcal/mol"
            end if
          end if
  !
  !  Check for dipole 
  !
          i = index(info," d=")
          if (i /= 0) then
            sum = reada(info, i+3)
            i = index(info, " dr=")
            if (i /= 0) then
              name = info(i + 4:index(info(i + 5:)," ") + i + 4)
              do i = 1, nrefs
                if (index(name(:8),mnemonic(i)(2:)) /= 0) exit
              end do
              write(11,"(2a, f7.1,5a)")"REF: ", '<B>Dipole:</B>',&
              sum," Debye, <b>&nbsp;&nbsp;&nbsp; REF: </b>",refs(i)(:len_trim(refs(i)))
            else
              write(11,"(2a, f7.1,5a)")"REF: ", '<B><span class="b">&Delta;</span>H<sub>f</sub>:</B>',&
              sum," kcal/mol"
            end if 
          end if
  !
  !  Check for I.P. 
  !
          i = index(info," i=")
          if (i /= 0) then
            sum = reada(info, i+3)
            i = index(info, " ir=")
            if (i /= 0) then
              name = info(i + 4:index(info(i + 5:)," ") + i + 4)
              do i = 1, nrefs
                if (index(name(:8),mnemonic(i)(2:)) /= 0) exit
              end do
              write(11,"(2a, f7.1,5a)")"REF: ", '<B>I.P.:</B>',&
              sum," eV, <b>&nbsp;&nbsp;&nbsp; REF: </b>",refs(i)(:len_trim(refs(i)))
            else
              write(11,"(2a, f7.1,5a)")"REF: ", '<B>I.P.:</B>',&
              sum," eV"
            end if 
          end if
        end if
      end if
!
!  Check for PM6 geometry 
!
      i = index(info," pm6")
      if (i /= 0) then
        write(11,"(2a, f7.1,5a)")"REF:  <B>Geometry calculated using PM6 - for information only. </B>"
      end if
      write(11,"(a)")allb(2,j)(:len_trim(allb(2,j)))      
    end do
    write (7, "(' TOTAL NUMBER ',50X)")
    write (7, "(' OF FILES:             ',I4)") iloop
    write (7, "(' OF HEATS OF FORMATION:',i4,',  WITH',' REFS:',i4,' =',&
         &f6.2,'%   AND SD''S:',i4,' =',f6.2)") iheat, jheat, &
         &real (jheat) / iheat * 100, kheat, (100.d0*kheat) / iheat
    write (7, "(' OF DIPOLE MOMENTS:    ',i4,',  WITH',' REFS:',i4,' ='&
         &,f6.2,'%   AND SD''S:',i4,' =',f6.2)") idipo, jdipo, &
         &real (jdipo) / idipo * 100, kdipo, (100.d0*kdipo) / idipo
    write (7, "(' OF I.P.s:             ',i4,',  WITH',' REFS:',i4,' ='&
         &,f6.2,'%   AND SD''S:',i4,' =',f6.2)") iip, jip, &
         &real (jip) / iip * 100, kip, (100.d0*kip) / iip * 100.d0
    write (7, "(' OF GEOMETRY FILES:    ',i4,',  WITH',' REFS:',i4,' =',&
         &f6.2,'%')") igeom, jgeom, real (jgeom) / igeom * 100
end program files
