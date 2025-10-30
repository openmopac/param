 Program SetUpMap
      use chanel_C, only : ir, iw

!
      USE symmetry_C, ONLY: idepfn, locdep, depmul, locpar 
!
      use molkst_C, only : ndep, numat, koment, title, natoms, keywrd, nvar, refkey
!
      use elemts_C, only : elemnt
!
      use common_arrays_C, only : loc, labels, na, nb, nc, &
      & geo, coord, nat
!
      use dfport, only : system
!
!
!  Build a simple HTML document to list all the molecules
!
!  Start in M:/data   Use "molecules.txt" as the argument
!  "Molecules.txt" is made by adding "header_for_molecules.txt" to the start of "files.txt" 
!  "files.txt" is made by running "files" in M:/data/data_normal
!
  implicit none
  character (len=340), dimension (-20:10000) :: system_name, data_set_name, formula, clean_data_set_name
  character (len=340), dimension (-20:10000,4) :: ref_data = " "
  character (len=240) :: line, line1, set, header(200), working, spaces, location, folder
  integer ::  file_name = 10, input = 11, PM7_file = 14, &
  io_stat, counter, i, all_mols = 13, nmols, j, keys, notes = 16, nheader, &
  jmol = 17, lopt(3,1000), k
  logical :: exists, prt, pre, int, ts, l_intermolecular
  integer, external :: iargc
  location = "M:\PARAM\HTML Files\data_molecules\"
  inquire (DIRECTORY = trim(location), exist = exists)
  if ( .not. exists) then
    line = "mkdir """//trim(location)//""""
    i = SYSTEM(trim(line))
    if (i /= 0) then
    write(*,'(///10x,a,///)') "   Folder """//trim(location)//""" does not exist"
    call sleep(60)
    stop
    end if
  end if
  refkey = " NULL"
  koment = " NULL"
  title  = " NULL"
  ir = file_name
  write(spaces,"(30a)")('&nbsp;', i = 1,30)
  natoms = 1000
  call setup_mopac_arrays(natoms)
  i = iargc()
  if (i >= 1) then
    call getarg(1, set)
    call getarg(2, folder)   
  else
    write(*,*) "Needs argument = name of list of molecules"
    call sleep(1000)
  end if
  l_intermolecular = (index(set, "Interm") /= 0)
  i = len_trim(set)
  if (set(i - 3:i - 3) ==".") set(i - 3:) = " "
  inquire (file="files.txt", exist = exists)
  if (.not. exists) then
    write(*,*)" File 'files.txt' does not exist"
    stop
  end if
  line = trim(set)
  call upcase(line, len_trim(line))
  ts = (index(line, "TRANSITION") /= 0)
  open(unit=input, file="files.txt", status='UNKNOWN', form='FORMATTED') 
  open(unit=all_mols, file="M:/PARAM/HTML Files/"//trim(set)//".html")
!
!  Text at the top of all_molecules.html
!
  write(all_mols,*)" <HTML>"
  write(all_mols,*)" <BODY>"
  call fdate(line) 
      write(all_mols,'(a)')" Date:"//trim(line)
!  write(all_mols,*)'Jmol requires Java and ActiveX control to be installed.<BR>' 
!
!  Header for set of molecules file
!
  do i = 1,1000
    read(input,"(a240)", iostat=io_stat)line
    if (io_stat /= 0) exit
    if (index(line,"++++") /= 0) exit
    write(all_mols,"(a)")trim(line)
  end do
  if (i > 1000) then
    write(*,*)" No lines containing '++++' found"
    stop
  end if
  if (i == 1) then
    write(*,*)" File contains no lines!"
    stop
  end if
!
!  Header for each molecule
!
  do nheader = 1, 100
    read(input,"(a180)", iostat=io_stat)header(nheader)
    if (index(header(nheader),"++++") /= 0) exit
  end do
  nheader = nheader - 1
  counter = 1
  do j = 1, 100000
    read(input,"(a120)", iostat=io_stat)formula(counter)
    read(input,"(a120)", iostat=io_stat)data_set_name(counter)
    i = len_trim(data_set_name(counter))
    if (i > 1) then
      do
        if (data_set_name(counter)(1:1) /= " ") exit
        data_set_name(counter) = data_set_name(counter)(2:)
      end do
    end if
    if (.not. ts) then
      line = trim(data_set_name(counter))
      call upcase(line, len_trim(line))
      if (index(line," TS") /= 0) cycle
      if (index(line," TRANSI") /= 0) cycle
    end if
    do i = 1, 4
      read(input,"(a340)", iostat=io_stat)system_name(counter)
      if (system_name(counter)(1:4) ==  "REF:") then
        ref_data(counter,i) = system_name(counter)(5:)
      else
        exit
      end if
    end do
    if (io_stat /= 0) exit
    i = Index(system_name(counter),">") 
    if (i /= 0) system_name(counter)(i:i) = "-"
    i = Index(system_name(counter),"""") 
    if (i /= 0) system_name(counter)(i:i) = " "
    counter = counter + 1
  end do 
  nmols = counter - 1
  do i = 1,11
    system_name(nmols + i) = system_name(i)
    data_set_name(nmols + i) = data_set_name(i)
    formula(nmols + i) = formula(i)
    ref_data(nmols + i,:) = ref_data(i,:)
    system_name(-11 + i) = system_name(nmols - 11 + i)
    data_set_name(-11 + i) = data_set_name(nmols - 11 + i)
    formula(-11 + i) = formula(nmols - 11 + i)   
    ref_data(-11 + i,:) = ref_data(nmols - 11 + i,:) 
  end do

  nmols = counter - 1
!
! Create "Round Robin"
  do i = 1,11
    data_set_name(nmols + i) = data_set_name(i)
    data_set_name(-11 + i) = data_set_name(nmols - 11 + i)
  end do
  do counter = -10, nmols + 11
!
! make "clean" names for use by HTML
!
    line = trim(data_set_name(counter))
    line1 = " "
    j = 0
    do i = 1, len_trim(line)
      if (line(i:i) == "+") then
        j = j + 4
        line1(j - 3:j) = "plus"
      else if (line(i:i) == "[") then
        j = j + 1
        line1(j:j) = "("
       else if (line(i:i) == "]") then
        j = j + 1
        line1(j:j) = ")"
      else
        j = j + 1
        line1(j:j) = line(i:i)
      end if
    end do
    clean_data_set_name(counter) = line1
  end do
  iw = jmol
  do counter = 1, nmols
  !
  !  Write master list of all molecules
  !
   if (counter < 10) then
      i = 4
    else if (counter < 100) then
      i = 3
    else if (counter < 1000) then
      i = 2
    else
      i = 1
    end if
    write(all_mols,"(a,i5,1x,20(a,/))")'&nbsp; &nbsp;'//spaces(:i*6), counter, &
    & ' &nbsp; &nbsp;', &
    & '<a href="./data_molecules/'//trim(clean_data_set_name(counter))//'_jmol.html"> (JSmol)</a> &nbsp; &nbsp;', &
    & trim(system_name(counter)), &
    & '<a name="'//trim(clean_data_set_name(counter))//'"></a><BR>'
  !
  !  Write map for each molecule
  !
    working = trim(clean_data_set_name(counter))
    do i = 1, len_trim(working)
      if (working(i:i) == " ")working(i:i) = "_"
      if (working(i:i) == ",")working(i:i) = "_"
      if (working(i:i) == "(")working(i:i) = "_"
      if (working(i:i) == ")")working(i:i) = "_"
      if (working(i:i) == "'")working(i:i) = "_"
    end do 

    open(unit=jmol, file=trim(location)//trim(clean_data_set_name(counter))//"_jmol.html")
!
!  Write common header
!
    write(jmol,"(a)") '<!DOCTYPE html>', &
    ' <HTML>',(trim(header(i)),i = 1,nheader),'<H3>', &
    ' <html><title>'//trim(clean_data_set_name(counter))//'</title> <head>', &
    '<meta charset="utf-8"> <script type="text/javascript" src="../../jsmol/JSmol.min.js"></script>', &
    '<script type="text/javascript"> ', &
    '$(document).ready(function() {Info = {', &
    'width: 500, ', &
    'height: 500, ', &
    'color: "0x000000", ', &
    'disableInitialConsole: true, ', &
    'addSelectionOptions: false, ', & 
    'j2sPath: "../../jsmol/j2s",', &
    'jarPath: "../../jsmol/java",', &
    'use: "HTML5", script: ', &
    '"load 	\'
    write(jmol,"(a)") trim(working)//'.xyz; \'
    write(jmol,"(a)") &
    '	set measurementUnits ANGSTROMS; \', &
    '	hBonds calculate; \', &
    '	connect 0.8  1.5 (hydrogen) (phosphorus) create; \', &
    '	set zoomLarge false;"', &
    '} ', &
    '$("#mydiv").html(Jmol.getAppletHtml("jmolApplet0",Info))}); ', &
    '</script>', &
    '</head>', &
    '<body>'
      if (counter < 10) then
        i = 4
      else if (counter < 100) then
        i = 3
      else if (counter < 1000) then
        i = 2
      else
        i = 1
      end if
      write(jmol,"(a,i5,1x,30(a,/))") &
      spaces(:i*6), counter, trim(system_name(counter)), &
   '  </H3>', &
   '  <a href="'//trim(clean_data_set_name(counter -  1))//'_Jmol.html"> ', &
   ' (Previous)</a> &nbsp;&nbsp;&nbsp;', &
   ' <a href="../'//set(:len_trim(set))//'.html#',trim(clean_data_set_name(counter))//'">', &
   ' (Back)</a>  &nbsp;&nbsp;&nbsp;', &
   '  <a href="'//trim(clean_data_set_name(counter +  1))//'_Jmol.html"> ', &
   ' (Next)</a> &nbsp; &nbsp; &nbsp; &nbsp; '
      if (l_intermolecular) then
        if (index(data_set_name(counter),  "water dimer") /= 0) then
          write(jmol,'(a)') '<B>Geometry from G. S. Tschumper, M. L. Leininger, B. C. Hoffman, '// &
            'E. F. Valeev, H. F. Schaefer III,  M. Quack, J. Chem. Phys. 116, 690 (2002)</B><BR>'
        else
          write(jmol,'(a)') '<B>Geometry from the BEGDB - the Benchmark Energy and Geometry Database</B><BR>'
        end if
      else
        if (index(system_name(counter), "(Geo)") == 0) write(jmol,'(a)') '<B>Geometry predicted using PM7</B><BR>'
      end if
!
!  See if notes exist
!
     line1 = data_set_name(counter)(1:len_trim(data_set_name(counter)))
     do i = 1, len_trim(line1)
        if (line1(i:i) ==  " ") line1(i:i) = "_"
      end do
     inquire (file = trim(line1)//".txt", exist = exists)
      if (exists) then
        close(notes)
  !
  !  Write out notes on this system
  !
        pre = .false.
        open(unit=notes, file=trim(line1)//".txt")
        do
          read(notes,"(a120)", iostat=io_stat) line       
          if (io_stat /= 0) exit
          if (index(line,"<PRE>") /= 0) pre = .true.
          if (index(line,"</PRE>") /= 0) pre = .false.
          if (pre) then
             write(jmol,"(3a)")trim(line)
          else
             write(jmol,"(3a)")trim(line),'<BR>'
          end if
         
        end do
      end if
      write(jmol,"(a)")'<table border="0" cellpadding="0" cellspacing="0" &
    &style="border-collapse: collapse" bordercolor="#111111" width=800><tr>'
!
! Start of big table (1by2)
!   
!
!
!  Jmol script
!
    write(jmol,"(30(a,/))") &
' <TD><h3 align="left">', &
' <span id=mydiv></span><a href="javascript:Jmol.script(jmolApplet0)"></a></h3>'
    open(file_name, file=trim(folder)//"\"//trim(data_set_name(counter))//".mop")
!
      write(jmol,"(a)")'</td><td>  &nbsp;  &nbsp; </td><td width=450>'
!
! Table of nearby molecules goes in right box
!
      j = 0
      do i = counter - 10, counter + 10
        j = max(j, len_trim(system_name(i)))
      end do
      write(jmol,"(a,i4.4,a)")'<font face="Courier New"><P ALIGN="LEFT"><LEFT><TABLE CELLSPACING=0 BORDER=0 CELLPADDING=0 WIDTH=', &
      300 + j*10, '>  <TR>  <TD>  &nbsp; &nbsp; #&nbsp; </TD><TD> Species </TD><TD>  Formula  </TR>'
!
! Entries in table
!
    do i = counter - 10, counter + 10
      if (i > 0) then
        if (i > nmols) then
          j = i - nmols
        else
          j = i
        end if
      else 
       j = nmols + i
      end if
       if (j < 10) then
        k = 4
      else if (j < 100) then
        k = 3
      else if (j < 1000) then
        k = 2
      else
        k = 1
      end if
      if (i == counter) then
        write(jmol,"(a,i5,1x,9a)")'<TR>  <TD>'//spaces(:k*6),j, '</TD><TD> '//trim(system_name(i)), &
    '</TD><TD> ',formula(i)(:20),' </TD><TD></TR>'
      else
        write(jmol,"(a,i5,1x,9a)")'<TR>  <TD>'//spaces(:k*6),j, '</TD><TD><a href="', &
        trim(clean_data_set_name(i)),'_jmol.html">', &
        trim(system_name(i)),'</a></TD><TD>'//formula(i)(:20)//'</TD><TD></TR>'
      end if
    end do
    write(jmol,"(a)")'</TABLE></LEFT><p></td></tr></table><BR></EMBED>'
    do i = 1,4
      if (ref_data(counter,i) == " ") exit
      write(jmol,"(a)")trim(ref_data(counter,i))//"<BR>"
    end do
  !
  !  Read in mopac data set
  !
      open(unit=PM7_file, file=trim(location)//trim(working)//".xyz")
      rewind (file_name)
      read(file_name,"(a120)", iostat=io_stat)line
      read(file_name,"(a120)", iostat=io_stat)line
      read(file_name,"(a120)", iostat=io_stat)line
      
      do i = 1, len_trim(line) + 1
        if (line(i:i) >= "a" .and. line(i:i) <= "z") line(i:i) = char(ichar(line(i:i)) + ichar("A") - ichar("a")) 
      end do
      prt =  (index(line,"CCDC") == 0) 
      rewind (file_name)
      if (prt) write(jmol,"(a)")'<PRE>  '
      do keys = 1,3
        read(file_name,"(a120)", iostat=io_stat)line
        do
          if (line(1:1) /= "*") exit
          read(file_name,"(a120)", iostat=io_stat)line
        end do
        if (keys == 1) then
          call tidy_map(line, len_trim(line))
          keywrd = line
          call upcase(keywrd, len_trim(keywrd))
        end if
        if (io_stat /= 0) exit
        if (keys == 3) then
          line1 = line
          j = 0
           do i = 1, len_trim(line1) + 1
            j = j + 1
            if (line1(i:i) == "<") then
              line(j:) = "&lt"//line1(i + 1:)
              j = j + 2
            else if (line1(i:i) == ">") then
              line(j:) = "&gt"//line1(i + 1:)
              j = j + 2
            else
              line(j:) = line1(i:)
            end if
          end do
        end if
        if (prt) write(jmol,"(a)")trim(line)
      end do
      nvar = 0
      if (working(1:5) == "b2cl4") then
        continue
      end if
      call getgeo (file_name, labels, geo, coord, lopt, na, nb, nc, int)      
      if (natoms > 0 .and. numat > 0) then
        if (index(keywrd, " SYM") /= 0) then
          call getsym(locpar, idepfn, locdep, depmul) 
        else
          ndep = 0
        end if
        call gmetry(geo,coord)
         do i = 1, natoms 
           do j = 1, 3 
            if (lopt(j,i) > 0) then
              nvar = nvar + 1 
              loc(1,nvar) = i 
              loc(2,nvar) = j 
            end if
          end do 
        end do 
        nat(1) = 0
        if (prt) call geout(jmol)
        write(PM7_file,"(i5)") numat
        write(PM7_file,"(a)") " "
        j = 0
        do i = 1, natoms
          if (labels(i) < 98) then
            j = j + 1
            nat(j) = labels(i)
          end if
          if (labels(i) == 107) then
            j = j + 1
            nat(j) = labels(i)
          end if
        end do
        do i = 1, numat
          write(PM7_file,"(a, 3f12.5)")elemnt(nat(i)),(coord(j,i),j=1,3)
        end do
      end if
      if (prt) then
        write(jmol,"(a)")'</PRE>'
      else
        write(jmol,"(a)")' For X-Ray structure, contact the CCDC: '// &
          '<a href="http://www.ccdc.cam.ac.uk/">http://www.ccdc.cam.ac.uk/</a><BR>'
      end if
      close(PM7_file, STATUS = "keep")
      write(jmol,"(a)")'<BR></HTML>'
      close(jmol,  STATUS = "keep")
    end do
    write(all_mols,"(a)")'</TABLE></LEFT>'
    write(all_mols,*)" <HTML><BODY>"
end program SetUpMap
subroutine tidy_map(line, len)
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
  end subroutine tidy_map
  