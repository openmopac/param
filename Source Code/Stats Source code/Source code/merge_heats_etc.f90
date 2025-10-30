subroutine merge(nmethods, methods)
use common_texts, only: type_of_ref, heading_htm, unit, l_set, title, output_file
  character (len=30), dimension (20) :: methods
  integer :: nmethods
  integer :: ihtml = 22
!
  type_of_ref = "heats"
  heading_htm = "Heat of Formation"
  unit = "kcal/mol"
  call print_table(ihtml, nmethods, methods)
  if (l_set) return
!
  type_of_ref = "geos"
  heading_htm = "Geometries"
  unit = "Angstroms and Degrees"
  call print_table(ihtml, nmethods, methods)
!
  type_of_ref = "dips"
  heading_htm = "Dipole"
  unit = "D"
  call print_table(ihtml, nmethods, methods)
!
  type_of_ref = "ips"
  heading_htm = "I.P."
  unit = "eV"
  call print_table(ihtml, nmethods, methods)
end subroutine merge 
subroutine print_table(ihtml, nmethods, methods)
!**********************************************************************
!
!  MERGE WILL MERGE PARAMETERIZATION OUTPUT FILES SO AS TO GENERATE
!  PUBLICATION-QUALITY TABLES
!
!**********************************************************************
use common_texts, only: type_of_ref, heading_tex, heading_htm, unit, all, solid, name, name1, &
  use_HDI_system, folder, l_set, n_files, set_of_filenames, title, output_path
use elements, only: must_have, can_have, element
implicit none
integer, parameter ::maxmol = 20000
integer :: nmethods, ihtml, i, j, k, l, m, ierr, nmols, loop, ndata, ir=25, n_ref, &
  prt_ref(maxmol), naves(20), io_stat
character :: ref(maxmol)*100, line*300, formla*12, methods(20)*30, one_file_name(maxmol)*60, &
  symbol2*1, one_element*2, line1*300, line2*300, date*24
character (len=285) :: mol_name(maxmol), new_name, old_name, fmt
logical :: is_geo, okay(maxmol), ispresent(107), opend, save_file, exists
double precision :: ref_par(maxmol), calc_par(maxmol,20), dummy, calcd, aves(20)
!
!  MERGE assumes that all input files have exactly the same layout.
!
  save_file = .false.
  n_ref = 0
  do i = 1, nmethods
    line = trim(methods(i))//"."//type_of_ref
    call add_path(line, "INPUT ")
    inquire (file=trim(line), exist = exists)
    if (exists) open(unit = 6 + i, status="unknown", file=trim(line))
  end do
  inquire(unit=ihtml, opened=opend) 
  if (opend) close(ihtml)
  line = name1
  call upcase(line, len_trim(line))  
  all = (index(" "//line, " ALL") /= 0 .or. len_trim(line) == 0)
  if (solid) then
    if (all) then
      open(unit=ihtml, status="unknown", file="table_of_"//trim(heading_tex)//"_for_Solids.html") 
    else
      open(unit=ihtml, status="unknown", file="table_of_"//trim(type_of_ref)//"_for_Solids_containing_"//trim(name)//".html") 
    end if        
  else
    if (all) then
      if (title /= " " .and. type_of_ref == "heats") then
        open(unit=ihtml, status="unknown", file=trim(output_path)//trim(title)//".html")
      else if (len_trim(name) == 0) then
        open(unit=ihtml, status="unknown", file=trim(output_path)//"table_of_"//trim(type_of_ref)//".html")
      else
        open(unit=ihtml, status="unknown", file=trim(output_path)//"table_of_"//trim(type_of_ref)//"_for_"//trim(name)//".html")
      end if
    else
      if (type_of_ref == "geos") return
      open(unit=ihtml, status="unknown", file=trim(output_path)//"table_of_"//trim(type_of_ref)//"_for_"//trim(name)//".html")
    end if     
  end if    
!
! Write out the HTML header
!
  line = " "
  call fdate(date)
  if (unit == "kcal/mol") then
    line = name1
    call upcase(line, len_trim(line))  
    all = (index(line, "ALL") /= 0 .or. len_trim(line) == 0 .or. l_set)
    solid = (line == "solids") 
    write(ihtml,'(a)')" Date: "//trim(date)
    write(ihtml,'(a)')"<HTML><HEAD><TITLE>  Errors in "//trim(heading_htm)//"  </TITLE></HEAD>"
    if (solid) then 
      if (all) then
        write(ihtml,'(a)')"<H2 ALIGN=""CENTER""> Errors in Calculated HoF for Solids ("//unit(1:len_trim(unit))//")", &
        '(<a href="index.html">Back</a> <a href="Accuracy of PM7 and PM6-D3H4.html">Accuracy</a> )'// &
        '<BR> <a href="../OpenMOPAC.NET/PM6/List_of_solids.html">(Individual Species)</a></H3>'
      else
        write(ihtml,'(a)')"<H2 ALIGN=""CENTER""> Errors in Calculated HoF for Solids containing the elements: "// &
        trim(old_name)//" ("//unit(1:len_trim(unit))//")", '<a href="index.html">(Back)</a>'// &
        '<BR> <a href="../OpenMOPAC_NET/PM6/List_of_solids.html">(Individual Species)</a></H3>'
!
! Test if this code is ever run
!
          dummy = -1.d0
          dummy = sqrt(dummy)
          write(ihtml,'(a)')dummy
      end if
    else
      if (heading_tex /= " ") then
         write(ihtml,'(a)')"<H2 ALIGN=""CENTER""> Errors in Calculated "//trim(heading_tex)// &
           " for "//trim(name)//" set <BR>"
      else if (all) then
        if (title == " ") title = "Errors in Calculated Heats of Formation"          
        write(ihtml,'(a)')"<H2 ALIGN=""CENTER"">"//trim(title)//"  (kcal/mol)<BR>"
      else
        write(ihtml,'(a)')"<H2 ALIGN=""CENTER""> Errors in Calculated Heats of Formation for compounds containing"// &
        " the elements "//trim(name1)//"  ("//trim(unit)//")<BR>"
      end if  
      write(ihtml,'(a)') '</H3> <p align="center">(<a href="../index.html">Home</a> '// &
      '<a href="Accuracy of PM7 and PM6-D3H4.html">Accuracy</a> '// &
      '<a href="..\Manual\index.html">Manual</a>)'
    end if  
  else
    write(ihtml,'(a)')" Date: "//trim(date)
    write(ihtml,'(a)')"<HTML><HEAD><TITLE>  Errors in "//heading_htm//"  </TITLE></HEAD>"
    write(ihtml,'(a)')"<H2 ALIGN=""CENTER""> Errors in "//heading_htm//" ("//unit(1:len_trim(unit))//")</H2>"
  end if  
  write(ihtml,'(''<P ALIGN="CENTER"><CENTER><TABLE CELLSPACING=0 BORDER=1 ''// &
  &'' CELLPADDING=7 WIDTH=900>'')')
  is_geo = (heading_htm == "Geometries")
  write(ihtml,"(a)")'<TR>  <TD>Empirical </TD><TD> <p align="center">  Name  ' 
  if (is_geo) write(ihtml,"(a)")'<TD> Variable </TD>' 
  write(ihtml,"(a)") '</TD><TD><P ALIGN="CENTER">Expt.   </TD>'
  do i = 1, nmethods
    write(ihtml,"(4a)")"<TD><P ALIGN=""CENTER"">   ",trim(methods(i)),'  </TD>' 
  end do
  if (all) write(ihtml,"(4a)")'<TD><P ALIGN="CENTER">  Ref. </TD></TR>'
!
! Dummy read over first two lines
!
  do i = 7, 7 + nmethods
    read(i,'(a)', iostat=ierr)(new_name,j=1,2)
  end do
!
! Read in data for first method
!
  
  if (is_geo) then
! 1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890123
!  H2            Hydrogen                        H-H            0.741   0.760     a  0.018          1  1
!  BH2           BH2                             B-H            1.180   1.141     e -0.039          1  5
!  LiC5H5        Lithium cyclopentadienide       Li-C           2.110   2.374    ee  0.264          3  6
!
    fmt = "(a59,f8.2,f8.2,a6,17x,a12,a40)"
    all = .true.
  else
    fmt = "(a63,f8.2,f10.2,9x,a6, 4x, a12, a60)"
  end if
  ndata = 0
  okay = .true.
  do i = 1, maxmol
    read(7,'(a)', iostat = j)line
    if (j /= 0) exit
    read(line(1:),fmt,iostat=ierr)mol_name(i),ref_par(i), calc_par(i,1), ref(i), formla, one_file_name(i)
    if (is_geo) then
      continue
    end if
    if (mol_name(i)(15:17) == " ") then
      if (i > 1) then
        ndata = ndata + 1
        mol_name(i)(1:47) = mol_name(i - 1)(1:47)
        mol_name(i)(60:60) = char(ndata + ichar("0")) 
      end if
    else
      ndata = 0
    end if
    if (ierr /= 0) exit
    if (formla == " ") exit
    if (is_geo) then
      line = " "//formla//" 0 0 0"
      read(line,*)j, k, l, m
      write(formla,'(4a)')must_have(j), must_have(k), must_have(l), must_have(m)
      if (formla(1:1) == " ") formla = formla(2:)
    end if
    if (l_set) then
      line = trim(one_file_name(i))
      call upcase(line, len_trim(line))       
      do j = 1, len_trim(line)
        if (line(j:j) /= " ") exit
      end do
      do j = 1, n_files
        if (line == set_of_filenames(j)) exit
      end do
      okay(i) = (j <= n_files)
      cycle
    end if
    ispresent = .false.
    if (all) then
      okay(i) = .true.
    else
      j = 1
      do
        if(formla(j:j) == " ") exit
        symbol2 = formla(j + 1:j + 1)
        if(symbol2 >= "A" .and. symbol2 <= "Z" .or. symbol2 == " ") then
!
! Element has one symbol, e.g. H, C, N, O, etc.
!    
          one_element = " "//formla(j:j)
          j = j + 1
        else
!
! Element has two symbols, e.g. Cl, Br, Sb, etc.
!
          one_element = formla(j:j+1)
          j = j + 2       
        end if
!
!  Now identify the elements that are present
!
        do k = 1, 107
          if (element(k) == one_element) exit
        end do
        if ( k < 108) ispresent(k) = .true.
      end do
!
!  Identify the elements that MUST be present
!
      do k = 1, 107
        if (must_have(k) /= "**" .and. .not. ispresent(k)) exit
      end do
      okay(i) = (k == 108) 
!
! Identify the elements that MUST NOT be present, indicated by "**"
!
      if (okay(i)) then
        do k = 1, 107
          if (can_have(k) == "**" .and. ispresent(k)) exit
        end do
        okay(i) = (k == 108) 
      end if
    end if 
  end do
  nmols = i - 1
  calc_par(:,2:) = -1.d6
!
! Read in data for the other methods
!
  do loop = 2, nmethods
    do i = 1, maxmol
      read(loop + 6, fmt, iostat=ierr) new_name, dummy, calcd
      if (abs(dummy) < -10.d0) stop
      if (new_name(15:17) == " ") then
      ndata = ndata + 1
      new_name(1:47) = old_name(1:47)
      new_name(60:60) = char(ndata + ichar("0")) 
    else
      ndata = 0
    end if
      old_name = new_name
      if (ierr /= 0) exit
      do j = 1, nmols
        if (new_name(:60) == mol_name(j)(:60)) exit
      end do
      if (j <= nmols) calc_par(j,loop) = calcd
    end do
  end do
  do i = 2, nmethods
    close(6 + i)
  end do
!
!   Write out the table, formatted
!
  line2 = " "
  naves = 0
  aves = 0.d0
  do i = 1, nmols
    if (okay(i)) then
      line = mol_name(i)
!
!   Print statements
!
!
!  Start of line
!

      save_file = .true.
      line1 = "<a href=""data_molecules/"
      new_name = trim(one_file_name(i))
      j = len_trim(line1) 
      do k = 1, len_trim(new_name)
        if (new_name(k:k) == "+") then
          j = j + 4
          line1(j - 3:j) = "plus"
        else if (new_name(k:k) == "[") then
          j = j + 1
          line1(j:j) = "_"
          else if (new_name(k:k) == "]") then
          j = j + 1
          line1(j:j) = "_"
        else
          j = j + 1
          line1(j:j) = new_name(k:k)
        end if
      end do
      line1 = trim(line1)//"_jmol.html"">"//trim(line(16:47))//"</a>"
      
      if (is_geo) then
        if (line(:47) == line2(:47)) then
          write(ihtml,"(a,f8.2,a)")"<TR>  <TD> </TD><TD> </TD><TD>"//line(48:59)//"</TD><TD><p style=""text-align:center"">", ref_par(i),"</TD>"
        else
          write(ihtml,"(a,f8.2,a)")"<TR>  <TD>"//line(:15)//"</TD><TD>"//&
      & trim(line1)//"</TD><TD>"//line(48:59)//"</TD><TD><p align=""center"">", &
      & ref_par(i),"</TD>"
        end if       
        line2 = line
      else
        line1 = line(16:)
        line1 = "<a href=""data_molecules/"
        new_name = trim(one_file_name(i))
        if (is_geo) then
          continue
        end if
        j = len_trim(line1) 
        do k = 1, len_trim(new_name)
          if (new_name(k:k) == "+") then
            j = j + 4
            line1(j - 3:j) = "plus"
          else if (new_name(k:k) == "[") then
            j = j + 1
            line1(j:j) = "_"
            else if (new_name(k:k) == "]") then
            j = j + 1
            line1(j:j) = "_"
          else
            j = j + 1
            line1(j:j) = new_name(k:k)
          end if
        end do
        line1 = trim(line1)//"_jmol.html"">"//trim(line(16:))//"</a>"
        if (solid) then
          write(ihtml,"(a,f8.2,a)")"<TR>  <TD>"//line(:15)//"</TD><TD>"//'<a href="../OpenMOPAC_NET/PM6/data_solids/'// &
          & trim(line1)//'_jmol.html">'//trim(line1)//'</A>'//"</TD><TD><p align=""right"">", ref_par(i),"</TD>"
!
! Test if this code is ever run
!
          dummy = -1.d0
          dummy = sqrt(dummy)
          write(ihtml,'(a)')dummy
        else
          if (all) then
              write(ihtml,"(a,f8.2,a)")"<TR>  <TD>"//line(:15)//"</TD><TD>"//trim(line1)// &
          "</TD><TD><p align=""right"">", ref_par(i),"</TD>"
          else
            write(ihtml,"(a,f8.2,a)")"<TR>  <TD>"//line(:15)//"</TD><TD>"//trim(line1)// &
          "</TD><TD><p align=""right"">", ref_par(i),"</TD>"
          end if 
        end if
      endif      
!
!   methods
!
      do k = 1, nmethods
        if (calc_par(i,k) > -1.d5) then
          naves(k) = naves(k) + 1
          aves(k) = aves(k) + abs(calc_par(i,k) - ref_par(i))
          write(ihtml,"(a,f8.2,a)")'<TD><p align="right">', calc_par(i,k) - ref_par(i),"</TD>"
        else
          write(ihtml,"(a)") '<TD><p align="center"> - </TD>'
        end if
      end do
!
!   End of line
!
      if (all) then
        read(ref(i),'(i6)') j
        do k = 1, n_ref
          if (prt_ref(k) >= j) exit
        end do
        if (k > n_ref) then
          n_ref = n_ref + 1
          prt_ref(n_ref) = j
        end if
        write(ihtml,"(a,i6,a)")'<TD><p align="center">', k, "</TD></TR>"
      end if
    end if
  end do 
  if (l_set) then
    write(ihtml,"(a,f8.2,a)")"<TR>  <TD> </TD><TD>Average unsigned errors:</TD><TD><p align=""right""> </TD>"
    do k = 1, nmethods
      write(ihtml,"(a,f8.2,a)")'<TD><p align="right">', aves(k)/max(1,naves(k)),"</TD>"
    end do
    write(ihtml,"(a)")"</TR>"
  end if
  write(ihtml,'(a)')"</TABLE> <BR> "
  if (all) then
    rewind(7)
    do 
      read(7,'(a5)', iostat = ierr)line
      if (line(1:5) == "*****")exit
      if (ierr /= 0) exit
    end do  
    k = 1
    do i = 1, 1000
      read(7,'(i5, a300)',iostat=ierr)j, line
      if (ierr /= 0) exit
      if (j == prt_ref(k)) then
        line = line(1:1)//" "//trim(line(2:))
        write(ihtml,"(a, i6, a)")"<P align=""LEFT"">", k, trim(line)//"</P>"
        k = k + 1
        if (k > n_ref) exit
      end if      
    end do
  end if
  line = "Footnote.txt"
  call add_path(line, "INPUT ")
  inquire (file=trim(line), exist = exists)
    if (exists) then
      open(unit = 8, status="unknown", file=trim(line))
      do
        read(8,'(a)', iostat = io_stat )line
        if (io_stat /= 0) exit
        write(ihtml,"(a)") trim(line)
      end do
    end if

  write(ihtml,"(a)")"</HTML>"
  if (.not. save_file) close (ihtml, status = 'delete')    
end subroutine print_table
