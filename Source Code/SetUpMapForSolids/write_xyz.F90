subroutine write_xyz(arc_file, xyz_file, cell_a, cell_b, cell_c, PM7_date, hof_ref1, hof_calc)
  use chanel_C, only : f_len, P_PM6_arc_files, PM6_input
  use molkst_C, only: natoms, nvar, numat, line, current_data_set_name
  USE symmetry_C, ONLY: idepfn, locdep, depmul, locpar 
  use common_arrays_C, only : labels, coord, geo, na, nb, nc, tvec, nat, loc
  use elemts_C, only : elemnt
  implicit none
  character*8 :: arc_file
  character(len=f_len) :: xyz_file 
  character :: PM7_date*26, line_1*240
 
  double precision, intent(out) :: cell_a, cell_b, cell_c, hof_ref1, hof_calc
!
  integer :: io_stat, i, jmol = 17, keys, j, lopt(3,2000)
  integer :: xyz_channel = 44, arc_channel = 45, structure_type
  logical :: prt, opend, l_sym
  double precision :: sum
  double precision, external :: reada
!
!  Write MOP file for molecular structure
!
    structure_type = 0
    if (index(arc_file,"X_R")   /= 0) structure_type = 1
    if (index(arc_file,"PM7") /= 0) structure_type = 2
    if (index(arc_file,"PM6_D3H4") /= 0) structure_type = 3
    open(unit=xyz_channel, file=trim(xyz_file), iostat = io_stat)
    if (io_stat /= 0) then
      write(*,'(a)')" Problem opening: '"//trim(xyz_file)//"' - file skipped"
      write(23,'(a)')" Problem opening: '"//trim(xyz_file)//"' - file skipped"
      return
    end if
    arc_channel = 45
    if (structure_type == 3) then  
      arc_channel = arc_channel + 3  
      inquire (unit=arc_channel, opened=opend) 
      if (opend) close(arc_channel)
      open(arc_channel, file=trim(P_PM6_arc_files)//trim(PM6_input))
    else if (structure_type == 2) then
      arc_channel = arc_channel + 1
    else
      arc_channel = arc_channel + 2
    end if
    rewind (arc_channel)
    do 
      read(arc_channel,"(a120)", iostat=io_stat)line
      if (index(line, "MOPAC20") /= 0) exit
      if (io_stat /= 0) exit
    end do
    i = index(line, "MOPAC20")
    if (structure_type == 2) then
!
! Read in date and convert to fractional date
!
      read(arc_channel,"(a80)", iostat=io_stat)PM7_date
      i = index(PM7_date, ":")
      if (PM7_date(i - 5:i - 5) == " ") PM7_date(i - 5:i - 5) = "0"
!
!  Find the day of the month
!
      read(PM7_date(i - 4:i - 3),'(f4.0)') sum
!
! Convert hour into fraction of a day and add on to the day
!
      sum = sum + reada(PM7_date, i - 2)/24.d0
!
! Convert minute into fraction of a day and add on to the day
! 
      sum = sum + reada(PM7_date, i + 1)/1440.d0
 
      write(line,"(a5, f3.1, a5)")PM7_date(i - 9:i - 5),sum, PM7_date(i + 7:i + 10)
      PM7_date = trim(line)
    end if
    hof_ref1 = -1.d10
    hof_calc = -1.d10
    do 
      read(arc_channel,"(a120)", iostat=io_stat)line
      call upcase(line, len_trim(line))
      if (index(line,"H=") /= 0 .and. index(line,"H= ") == 0 .and. index(line, "GUESS") == 0) then
        i = index(line," H=")
        if (line(i + 3:i + 3) /= " ") hof_ref1 = reada(line,i + 3)
      end if
!    if (hof_ref1 >  -1.d9) then
        if (line(10:30) == " H.O.F. PER UNIT CELL") hof_calc = reada(line,40)
!    end if
      if (index(line, "FINAL GEOMETRY OBTAINED") /= 0) exit
      if (io_stat /= 0) exit
    end do    
    prt = .true.
    do keys = 1,3
      do
        read(arc_channel,"(a120)", iostat=io_stat)line
        if (line(1:1) /= "*") exit
      end do
      if (keys == 1) then
        call tidy1(line, len_trim(line))
        line_1 = trim(line)
        call upcase(line_1, len_trim(line_1))
        l_sym = (index(line_1, " SYMM") /= 0)
        if (structure_type == 1) then
          prt = (index(current_data_set_name,"(ICSD") == 0)
          if (prt) then
            i = len_trim(current_data_set_name)
            if (current_data_set_name(i:i) /= ")") then
              write(jmol,"(a)")'<PRE>  X-Ray data set:'
              write(jmol,"(a)")trim(line)
              cycle
            end if
!
!  Test to see if the data-set name contains a CCDC identifier.
!  This would be six alphabetic characters
!
            do j = i - 1, i - 10, -1
              if (current_data_set_name(j:j) == "(") exit
            end do
            do i = j + 1, j + 6
              if (.not. (ichar(current_data_set_name(i:i)) >= ichar("A") .and. &
                ichar(current_data_set_name(i:i)) <= ichar("Z"))) exit
            end do
            prt = (i /= j + 7) 
          end if
        end if
        select case (structure_type)
          case (1)!  X-ray structure
            if (prt) write(jmol,"(a)")'<PRE>  X-Ray data set:'
          case (2)!  PM7 structure
            write(jmol,"(a)")'<PRE>  Optimized PM7 data set:'
          case (3)!  PM6_D3H4 structure
            write(jmol,"(a)")'<PRE>  Optimized PM6_D3H4 data set:'
          end select
      end if
      if (io_stat /= 0) exit
      if (prt) write(jmol,"(a)")trim(line)
    end do
    nvar = 0
    call getgeo (arc_channel, labels, geo, coord, lopt, na, nb, nc)
    call getsym(l_sym, arc_channel, locpar, idepfn, locdep, depmul) 
    if (numat > 0) then
      numat = numat - 3
      call gmetry(geo,coord)
      cell_a = sqrt(tvec(1,1)**2 + tvec(2,1)**2 + tvec(3,1)**2)
      cell_b = sqrt(tvec(1,2)**2 + tvec(2,2)**2 + tvec(3,2)**2)
      cell_c = sqrt(tvec(1,3)**2 + tvec(2,3)**2 + tvec(3,3)**2)
      do i = 1, natoms 
        do j = 1, 3 
          if (lopt(j,i) > 0)then         
            nvar = nvar + 1 
            loc(1,nvar) = i 
            loc(2,nvar) = j 
          end if
        end do 
      end do 
      if (prt) then      
        nat(1) = 0      
        call geout(-jmol)
      else
        if (index(current_data_set_name,"(ICSD") /= 0) then
          write(jmol,"(a)")' <PRE> For X-Ray structure, contact the ICSD: '// &
            '<a href="https://icsd.fiz-karlsruhe.de/search/index.xhtml/">'// &
            'https://icsd.fiz-karlsruhe.de/search/index.xhtml/</a><BR>'
        else
          write(jmol,"(a)")' <PRE> For X-Ray structure, contact the CCDC: '// &
            '<a href="http://www.ccdc.cam.ac.uk/">http://www.ccdc.cam.ac.uk/</a><BR>'
        end if
      end if
      write(xyz_channel,"(i5)") numat
      write(xyz_channel,"(a)") " "
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
        write(xyz_channel,"(a, 3f12.5)")elemnt(nat(i)),(coord(j,i) - coord(j,1),j=1,3)
      end do
    end if    
    close(xyz_channel,STATUS="keep")
    inquire (unit=arc_channel, opened=opend) 
    if (opend) close(arc_channel)
    write(jmol,"(a)")'<BR></PRE>'
    return
  end subroutine write_xyz
  