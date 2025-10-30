Program Sublimation
!
! Generate a list of heats of sublimation 
!
!  Input data consists of a list of substances, one per line.
!  Each line contains:
!    Name of the chemical, set off with quotation marks
!    Name of the solid (this often has extra detail in the name), set off with quotation marks
!    Reference heat of sublimation (in kcal/mol)
!    Reference, abbreviated
!  
  use Common_C, only : ir, iw, location_of_mols, location_of_solids, list_of_compounds, &
    mol_names, sol_names, heat_of_sub, refs, ir_mol, ir_sol, all_mol_names, all_sol_names, &
    heat_of_sub, all_PM6_D3H4_g, all_PM7_g, all_PM6_D3H4_s, all_PM7_s, calc_sub_PM7, &
    calc_sub_PM6_D3H4, ref_heat_sub, sorted_mol_names, header, footer, mol_temp_names, &
    sol_temp_names, sorted_ref_heat
  implicit none
  integer :: i, j, k, l, ii, jj,loop, io_stat, n_mols, n_solids, n_subs, n_refs, max_txt
  character :: line*1000, line1*1000, timestamp*24
  double precision :: ref, PM6_D3H4_error, PM7_error 
  logical :: exists(3000)
  double precision, external :: reada
!
! Read in all the data
!
  open(ir, file = trim(list_of_compounds), status='UNKNOWN', form='FORMATTED', action='READ')
  mol_names = " "
  do n_refs = 1, 2000
    read(ir,'(a)', iostat=io_stat) line
    if (io_stat /= 0) exit
    if (line == " ") exit
    i = index(line, '"')
    j = index(line(i + 1:), '"') + i
    mol_names(n_refs) = line(i + 1: j - 1)
    i = index(line(j + 1:), '"') + j
    if (i /= j) then
      j = index(line(i + 1:), '"') + i
      sol_names(n_refs) = line(i + 1: j - 1)
    else
      sol_names(n_refs) = trim(mol_names(n_refs))
    end if      
    heat_of_sub(n_refs) = reada(line, j)
  !  j = index(line(j:), ".") + j
  !  do
  !    j = j + 1
  !    if (line(j:j) == " ") exit
  !  end do
  !  do
  !    j = j + 1
  !    if (line(j:j) /= " ") exit
  !  end do    
  !  refs(n_refs) = trim(line(j:))   
  end do
!
! Open the files containing calculated values
!
  open(ir_mol, file = trim(location_of_mols), status='UNKNOWN', form='FORMATTED', action='READ')
  open(ir_sol, file = trim(location_of_solids), status='UNKNOWN', form='FORMATTED', action='READ')
!
!  Read in all the data
!
!            Molecules
!
  do
    read(ir_mol,'(a)') line
    if (index(line, "</TD></TR>") /= 0) exit
  end do
  do i =1, n_refs
    mol_temp_names(i) = mol_names(i)
    call upcase(mol_temp_names(i), len_trim(mol_temp_names(i)))
  end do
  exists = .false.
  n_mols = 0
  do k = 1, 10000
    read(ir_mol,'(a)', iostat = io_stat) line
    if (io_stat /= 0) exit
    i = index(line, """>")
    if (i == 0) exit
    j = index(line(i + 1:), "<") + i - 1
    if (j < 0) exit
    n_mols = n_mols + 1     
    all_mol_names(n_mols) = line(i + 2: j)
    call upcase(all_mol_names(n_mols), len_trim(all_mol_names(n_mols)))  
    i = index(line(j:), "right") + j
    ref = reada(line, i)
    read(ir_mol,'(a)', iostat = io_stat) line
    if (io_stat /= 0) exit
    all_PM7_g(n_mols) = ref + reada(line, 20)
    read(ir_mol,'(a)', iostat = io_stat) line
    if (io_stat /= 0) exit
    all_PM6_D3H4_g(n_mols) = ref + reada(line, 20)    
    do i = 1, 5
      read(ir_mol,'(a)', iostat = io_stat) line
      if (io_stat /= 0) exit
    end do 
    do i = 1, n_refs
      if (trim(mol_temp_names(i)) == trim(all_mol_names(n_mols))) exit
    end do
    if (i > n_refs) then
      n_mols = n_mols - 1
    else
      exists(i) = .true.
    end if
  end do
!
!            Solids
!
  i = 0
  do
    read(ir_sol,'(a)') line
    if (index(line, "</PRE>") /= 0) then
      i = i + 1
      if (i == 2) exit
    end if
  end do
  do i =1, n_refs
    sol_temp_names(i) = sol_names(i)
    call upcase(sol_temp_names(i), len_trim(sol_temp_names(i)))
  end do
  exists = .false.
  n_solids = 0
  do 
    read(ir_sol,'(a)', iostat = io_stat) line
    if (io_stat /= 0) exit
    i = index(line, "JSmol")
    if (i == 0) exit
    j = index(line(i:i + 50), ".")
    if (j == 0) i = 150
    n_solids = n_solids + 1
    j = i + j
    do k = 1, 4
      i = j
      j = index(line(i:), ".") + i
      if (k == 1) then
        all_PM7_s(n_solids) = reada(line, j - 5)
      end if
      if (k == 2) all_PM6_D3H4_s(n_solids) = reada(line, j - 5)
    end do
!
!  Jump over three spaces
!
    i = index(line(j + 19:), "<") + j + 17
    all_sol_names(n_solids) = line(j + 19:i)
    call upcase(all_sol_names(n_solids), len_trim(all_sol_names(n_solids)))
    do i = 1, n_refs
      if (trim(sol_temp_names(i)) == trim(all_sol_names(n_solids))) exit
    end do
    if (i > n_refs) then
      n_solids = n_solids - 1
    else
      exists(i) = .true.
    end if
  end do
!
!  Merge solids and molecules to make a master list
!
  open(iw,file="M:\PARAM\HTML Files\Heats of Sublimation.html")
  n_subs = 0
  do j = 1, n_solids
    do i = 1, n_refs
      if (sol_temp_names(i) == all_sol_names(j)) exit
    end do
    if (j > n_solids) cycle
    do k = 1, n_mols
      if (mol_temp_names(i) == all_mol_names(k)) exit
    end do
    if (k > n_mols) cycle
    n_subs = n_subs + 1
    if (mol_temp_names(i)(1:4) == "IODI") then
      continue
    end if
    all_mol_names(n_subs) = mol_names(i)
    all_sol_names(n_subs) = all_sol_names(j)
    sorted_ref_heat(n_subs) = heat_of_sub(i)
    calc_sub_PM7(n_subs) = all_PM7_g(k) - all_PM7_s(j)
    calc_sub_PM6_D3H4(n_subs) = all_PM6_D3H4_g(k) - all_PM6_D3H4_s(j)
    continue
  end do
!
! n_refs:             Number of compounds for which there are reference heats of sublimation
! mol_names:          Names of reference systems
! heat_of_sub:        Reference heats of sublimation
!
!
! n_subs:             Number of compounds for which there are calculated heats of sublimation
! all_mol_names:      Names of calculated systems
! calc_sub_PM7:       PM7 calculated heats of sublimation
! calc_sub_PM6_D3H4:  PM6_D3H4 calculated heats of sublimation
! 
  close (ir)
  open (ir, file=trim(header))
  call fdate(timestamp)
  do i = 1, 3
    read(ir, '(a)', iostat = io_stat) line
    if (io_stat /= 0) exit
    write(iw,'(a)')trim(line)
  end do
  write(iw,*)"Time stamp: "//timestamp
  do
    read(ir, '(a)', iostat = io_stat) line
    if (io_stat /= 0) exit
    write(iw,'(a)')trim(line)
  end do
  max_txt = 0
  do i = 1, n_refs
    max_txt = max(max_txt, len_trim(mol_names(i)))
  end do
  k = 0
  do i = 1, n_subs
    if (ref_heat_sub(i) > 0.d0) k = k + 1
  end do
  write(iw,'(a)') '<P ALIGN="CENTER"><CENTER><TABLE CELLSPACING=0 BORDER=1'
  write(iw,'(a)') 'CELLPADDING=7 WIDTH=600>'
  write(iw,'(a)') '<TR>  <TD ALIGN="CENTER">No. </TD><TD width="141"> <p align="center">  Compound  </TD>'
  write(iw,'(a)') '<TD colspan="2" width="110"><P ALIGN="CENTER">Ref   </TD>'
  write(iw,'(a)') '<TD colspan="2" width="100"><P ALIGN="CENTER">   PM7  </TD>'
  write(iw,'(a)') '<TD colspan="2" width="131"><P ALIGN="CENTER">   PM6-D3H4  </TD>'
  PM7_error = 0.d0
  PM6_D3H4_error = 0.d0
  j = 0
  do i = 1, n_subs
    if (ref_heat_sub(i) > -330.d0) then
      j = j + 1
      line = trim(all_sol_names(i))
      line1 = " "
      jj = 0
      do ii = 1, len_trim(line)
        if (line(ii:ii) == " ") then
          jj = jj + 1
          line1(jj:jj) = "_"
        else if (line(ii:ii) == "(") then
          jj = jj + 1
          line1(jj:jj) = "_"
         else if (line(ii:ii) == ")") then
          jj = jj + 1
          line1(jj:jj) = "_"
        else
          jj = jj + 1
          line1(jj:jj) = line(ii:ii)
        end if
      end do 
      write(iw,'(a, i4, a)') '<TR>  <TD ALIGN="CENTER" style="border-bottom-style: none; '// &
        'border-bottom-width: medium; border-top-style:none; border-top-width:medium">', j, '</TD>'
      write(iw,'(a, i4, a)') '<TD style="border-bottom-style: none; border-bottom-width: medium; '// &
        'border-top-style:none; border-top-width:medium"><a href=><a href="./data_solids/'//trim(line1) &
        &//'_Jmol.html" target="_blank">'//trim(all_mol_names(i))//"</a></TD>"
      
      write(iw,'(a, f10.1, a)')'<TD align="right" style="border-right-style: none; border-right-width: '// &
        'medium; border-bottom-style: none; border-bottom-width: medium; border-top-style:none; border-top-width:medium">', &
        sorted_ref_heat(i), '</TD><TD style="border-left-style: none; border-left-width: medium; border-bottom-style: none; '// &
        'border-bottom-width: medium; border-top-style:none; border-top-width:medium">&nbsp; </TD>'
      write(iw,'(a, f10.1, a)')'<TD align="right" style="border-right-style: none; border-right-width: medium; '// &
        'border-bottom-style: none; border-bottom-width: medium; ; border-top-style:none; border-top-width:medium">', &
        calc_sub_PM7(i),  '</TD><TD style="border-left-style: none; border-left-width: medium; border-bottom-style: none; '// &
        'border-bottom-width: medium; ; border-top-style:none; border-top-width:medium">&nbsp; </TD>'
      write(iw,'(a, f10.1, a)')'<TD align="right" style="border-right-style: none; border-right-width: medium; '// &
        'border-bottom-style: none; border-bottom-width: medium; ; border-top-style:none; border-top-width:medium">', &
        calc_sub_PM6_D3H4(i), '</TD><TD style="border-left-style: none; border-left-width: medium; border-bottom-style: '// &
        'none; border-bottom-width: medium; ; border-top-style:none; border-top-width:medium"> &nbsp; </TD>'
      PM7_error = PM7_error + abs(sorted_ref_heat(i) - calc_sub_PM7(i))
      PM6_D3H4_error = PM6_D3H4_error + abs(sorted_ref_heat(i) - calc_sub_PM6_D3H4(i))
    end if
  end do   
  write(iw,'(a)') '</TABLE><BR>'
  close (ir)
  open (ir, file=trim(footer))
  read(ir, '(a)', iostat = io_stat) line
  write(iw,'(a)')trim(line)
  write(iw,'(a,f6.3,a,f6.3)')"<p align=""center"">Average Unsigned Error for PM7:", PM7_error/j, &
    ", &nbsp; PM6-D3H4:", PM6_D3H4_error/j
  do
    read(ir, '(a)', iostat = io_stat) line
    if (io_stat /= 0) exit
    write(iw,'(a)')trim(line)
  end do
end program Sublimation

