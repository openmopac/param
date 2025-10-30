subroutine bigloop(nmols, maxtxt)
!
! Bigloop generates all the files for individual solids, and the master file of all solids.
! It also makes the 
!
  use chanel_C, only : P_Xray_arc_files, P_PM7_arc_files, P_PM6_arc_files, P_output_sub_folder, &
    input_data_set, P_notes_folder, notes, pm6_input
!
  use sets_of_names_C, only : f_len, x_ray_input, PM7_in, PM6_in, all_solids, line, line1, working, line2, hof_PM7_all, &
    data_set_name, all_clean_data_set_names
!
  use molkst_C, only : keywrd, title, numat, current_data_set_name, natoms, loop, moperr
!
  use common_arrays_C, only : na, nb, nc, coord, geo, labels, numats, nat, all_accuracy_PM7, all_accuracy_PM6, &
    hof_ref, hof_PM7, hof_PM6, tvec, elements, d_ref, d_PM7, d_PM6
!
  use elemts_C, only : elemnt
!
  use molkst_C, only: timestamp
!
  implicit none
!
  integer, intent (in) :: nmols
  integer, intent (out) :: maxtxt
!
  double precision :: PM7_a, PM7_b, PM7_c, PM6_a, PM6_b, PM6_c, Xray_a, Xray_b, Xray_c, tvec_PM7(3,3), &
    tvec_PM6(3,3), tvec_ref(3,3), alpha, beta, gamma, a, b, c, volume, density, sum, &
    sum1, hof_ref1, hof_calc, ref_coord(3,2000), PM7_accuracy, PM6_accuracy, lim1, lim2, lim3, lim4, lim5, &
    PM7_acc_prt, PM6_acc_prt
!
  integer :: counter, i, j, k, l, ii, kk, io_stat, Z_ref, Z_PM7, Z_PM6, lopt(3,2000), &
    iels(2,107), niels, mers(3), primes(10), nat_store(2000), score_PM7(6)=0, score_PM6(6)=0, &
    days_in_month(12), txt_counter
!
!  All the channels used here
!
  integer :: jmol = 17, jmol_XRay_fs = 18, jmol_PM7_fs = 19, jmol_PM6_fs = 20, arc_channel = 45, &
    notes1 = 22, file_name = 10
!
  logical :: exists, pre, H_ok, Be_ok, C_ok, N_ok, O_ok, F_ok, Mg_ok, P_ok, &
    S_ok, Cl_ok, Ca_ok, Zn_ok, Br_ok, Sr_ok, Mo_ok, I_ok, first = .true.
!
  character(len=26) :: PM7_date, Xray_date, PM6_date
!
  character :: chr, clean_data_set_name*(f_len), clean_data_set_name_m1*(f_len), clean_data_set_name_p1*(f_len), &
    relative_folder*(20), els(10000)*(8), formla(10000)*40, months(12)*(3), pad*(2000), txt_ref*(100), &
    txt_PM7*(100), txt_PM6*(100), txt_start*(1000), nbsp*(1000), PM7_geo*(100), PM6_geo*(100)
!
  double precision, external :: reada, accuracy
  integer, external :: numb
!
  data primes /2, 3, 5, 7, 11, 13, 17, 19, 23, 29/
  data days_in_month/31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
  data months /"Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"/
  relative_folder = "../../jsmol/"
  do i = 1, 166
    nbsp((i - 1)*6 + 1:) = "&nbsp;"
  end do
  counter = 0
  pad = " "
  maxtxt = 0
  do loop = 1, nmols
    if (mod(nmols - loop, 10) == 0) write(*,"('Left:',i7)") nmols - loop
    x_ray_input  = trim(P_Xray_arc_files)//trim(data_set_name(loop))//" (ReF).arc"
    inquire (file = trim(x_ray_input), exist = exists)    
    PM7_in  = trim(P_PM7_arc_files)//trim(data_set_name(loop))//".arc" 
    inquire (file = trim(PM7_in), exist = exists)
    if ( .not. exists) then
      write(*,"(a)") "PM7 arc   file:  "//trim(PM7_in)//" does not exist"
      write(23,"(a)") "PM7 arc   file:  "//trim(PM7_in)//" does not exist"
      cycle
    end if
    counter = counter + 1
    data_set_name(counter) = data_set_name(loop)
    PM7_a = 0.d0
    PM7_b = 0.d0
    PM7_c = 0.d0
    Xray_a = 0.d0
    Xray_b = 0.d0
    Xray_c = 0.d0
    PM7_date = " "
    PM6_in  = trim(P_PM6_arc_files)//trim(data_set_name(loop))//" (PM6-D3H4).arc"
!
!  Write map for each molecule
!
!
!  Read in all the data, once
!
    open(arc_channel, file=trim(PM7_in), action='READ')
    open(arc_channel + 1, status = "SCRATCH", form="FORMATTED", iostat = io_stat)
    if (io_stat /= 0) cycle
    rewind (arc_channel)
    rewind (arc_channel + 1)
    do i = 1,2000
     read(arc_channel,"(a120)", iostat=io_stat)line
      if (io_stat /= 0) exit
      write(arc_channel + 1,"(a)")trim(line)
    end do
!
!
    open(arc_channel, file=trim(x_ray_input), action='READ')
    open(arc_channel + 2, status = "SCRATCH", form="FORMATTED")
    rewind (arc_channel)
    rewind (arc_channel + 2)
    do i = 1,2000
     read(arc_channel,"(a120)", iostat=io_stat)line
      if (io_stat /= 0) exit
      write(arc_channel + 2,"(a)")trim(line)
    end do
    if (i == 1) then
      write(arc_channel + 2,"(a)")" error in X-ray data set"
      write(*,"(a)") "X-ray arc file:  "//trim(x_ray_input)//" exists, but is empty! " 
      write(23,"(a)") "X-ray arc file:  "//trim(x_ray_input)//" exists, but is empty! " 
    end if
!
!  Work out the X-ray geometry
!
    arc_channel = arc_channel + 2
    rewind (arc_channel)
    do 
      read(arc_channel,"(a120)", iostat=io_stat)line
      if (index(line, "MOPAC20") /= 0) exit
      if (io_stat /= 0) exit
    end do
    i = index(line, "MOPAC20")
    read(arc_channel,"(a90)", iostat=io_stat)line
    i = index(line, ":")
    Xray_date = line(i - 9:i - 4)//line(i + 6:i + 10)//" @ "//line(i - 2:i +5)
    if (Xray_date(5:5) == " ") Xray_date = Xray_date(:4)//"0"//Xray_date(6:) 
    do 
      read(arc_channel ,"(a120)", iostat=io_stat)line
      if (index(line, "FINAL GEOMETRY OBTAINED") /= 0) exit
      if (io_stat /= 0) exit
    end do   
    do
      read(arc_channel ,"(a200)", iostat=io_stat)keywrd
      if (keywrd(1:1) /= "*") exit
    end do
    call upcase(keywrd,len_trim(keywrd))
    i = index(keywrd, "Z=") 
    if (i /= 0) then
      Z_ref = nint(reada(keywrd,i+2))
    else
      Z_ref = 0
    end if
    do i = 1,2
      read(arc_channel ,"(a120)", iostat=io_stat)line 
    end do
    call getgeo (arc_channel, labels, geo, coord, lopt, na, nb, nc)
    if (moperr) stop
    arc_channel = arc_channel - 2
!
!   Get empirical formula
!
    call sub_rab(iels, niels)
!
!  Work out the PM7 geometry
! 
    arc_channel = arc_channel + 1
    rewind (arc_channel)
    do 
      read(arc_channel ,"(a120)", iostat=io_stat)line
      if (index(line, "FINAL GEOMETRY OBTAINED") /= 0) exit
      if (io_stat /= 0) exit
    end do 
    if (io_stat /= 0) then
      write(*,'(a)')" Fault detected in PM7 '.arc' file: '"//trim(data_set_name(counter))//"'"
      write(*,'(a)')" This is a severe error that must be corrected before continuing"
      write(23,'(a)')" Fault detected in PM7 '.arc' file: '"//trim(data_set_name(counter))//"'"
      write(23,'(a)')" This is a severe error that must be corrected before continuing"
      stop    
    end if
    do
      read(arc_channel ,"(a120)", iostat=io_stat)keywrd! keywrd
      if (keywrd(1:1) /= "*") exit 
    end do
    read(arc_channel ,"(a120)", iostat=io_stat)line! title
    read(arc_channel ,"(a120)", iostat=io_stat)line! koment
    call getgeo (arc_channel, labels, geo, coord, lopt, na, nb, nc)
    tvec_PM7 = tvec
    if (natoms == 0) then
      write(*,*)" Problem with '"//trim(data_set_name(counter))//"' - no atoms - skipping this entry"
      cycle
    end if
    rewind (arc_channel)
    line1 = trim(P_output_sub_folder)//trim(data_set_name(counter))//'_(PM7).html' 
    open(48,file=trim(line1))
    rewind (arc_channel) 
    write(48,'(a)')"<HTML><PRE>"
    do
      read(arc_channel,"(a120)", iostat=io_stat)line
      if (io_stat /= 0) exit
      write(48,'(a)')trim(line)
    end do  
    write(48,'(a)')"</HTML></PRE>"
    arc_channel = arc_channel - 1
!
!  Make a very "clean" filename: no characters other than
!  a-z, A-Z, 0-9, ., and -.
!
    line = data_set_name(counter)
    line1 = " "
    j = 0
    do i = 1, len_trim(line)
      chr = line(i:i)
      k = ichar(chr)
      if (line(i:i) == "+") then
        j = j + 4
        line1(j - 3:j) = "plus"
      else if (k >= ichar("a") .and. k <= ichar("z") .or. &
               k >= ichar("A") .and. k <= ichar("Z") .or. &
               k >= ichar("0") .and. k <= ichar("9") .or. &
               k == ichar(".") .or.  k == ichar("-")) then
        j = j + 1
        line1(j:j) = chr
      else
        j = j + 1
        line1(j:j) = "_"
      end if
    end do
    clean_data_set_name = line1
    all_clean_data_set_names(counter) = trim(clean_data_set_name)
    working =  line1
    do i = 1, len_trim(working)
      if (working(i:i) == " ")working(i:i) = "_"
      if (working(i:i) == ",")working(i:i) = "_"
      if (working(i:i) == "(")working(i:i) = "_"
      if (working(i:i) == ")")working(i:i) = "_"
      if (working(i:i) == "'")working(i:i) = "_"
    end do
    line = data_set_name(counter - 1)
    line1 = " "
    j = 0
    do i = 1, len_trim(line)
      chr = line(i:i)
      k = ichar(chr)
      if (line(i:i) == "+") then
        j = j + 4
        line1(j - 3:j) = "plus"
      else if (k >= ichar("a") .and. k <= ichar("z") .or. &
               k >= ichar("A") .and. k <= ichar("Z") .or. &
               k >= ichar("0") .and. k <= ichar("9") .or. & 
               k == ichar(".") .or.  k == ichar("-")) then
        j = j + 1
        line1(j:j) = chr
      else
        j = j + 1
        line1(j:j) = "_"
      end if
    end do
    clean_data_set_name_m1 = line1
    line = data_set_name(counter + 1)
    line1 = " "
    j = 0
    do i = 1, len_trim(line)
      chr = line(i:i)
      k = ichar(chr)
      if (line(i:i) == "+") then
        j = j + 4
        line1(j - 3:j) = "plus"
      else if (k >= ichar("a") .and. k <= ichar("z") .or. &
               k >= ichar("A") .and. k <= ichar("Z") .or. &
               k >= ichar("0") .and. k <= ichar("9") .or. &
               k == ichar(".") .or.  k == ichar("-")) then
        j = j + 1
        line1(j:j) = chr
      else
        j = j + 1
        line1(j:j) = "_"
      end if
    end do
    clean_data_set_name_p1 = line1
    j = 0
    open(unit=jmol, file=trim(P_output_sub_folder)//trim(clean_data_set_name)//"_Jmol.html", iostat = i)
    if (abs(i) > 0) then
      write(*,*)" Problem with '"//trim(clean_data_set_name)//"' - skipping this entry"
      write(*,*)" (Could not open '"//trim(P_output_sub_folder)//trim(clean_data_set_name)//"_Jmol.html')"
      cycle
    end if
    j = j + abs(i)
    open(unit=jmol_XRay_fs, file=trim(P_output_sub_folder)//trim(clean_data_set_name)//"_jmol_XRay_fs.html", iostat = i)
    if (abs(i) > 0) then
      write(*,*)" Problem with '"//trim(clean_data_set_name)//"' - skipping this entry"
      write(*,*)" (Could not open '"//trim(P_output_sub_folder)//trim(clean_data_set_name)// &
      trim(clean_data_set_name)//"_jmol_XRay_fs.html')"
      cycle
    end if
    j = j + abs(i)
    open(unit=jmol_PM7_fs, file=trim(P_output_sub_folder)//trim(clean_data_set_name)//"_jmol_PM7_fs.html", iostat = i)
    if (abs(i) > 0) then
      write(*,*)" Problem with '"//trim(clean_data_set_name)//"' - skipping this entry"
      write(*,*)" (Could not open '"//trim(P_output_sub_folder)//trim(clean_data_set_name)//"_jmol_PM7_fs.html')"
      cycle
    end if
    j = j + abs(i)
    open(unit=jmol_PM6_fs, file=trim(P_output_sub_folder)//trim(clean_data_set_name)//"_jmol_PM6_fs.html", iostat = i)
    if (abs(i) > 0) then
      write(*,*)" Problem with '"//trim(clean_data_set_name)//"' - skipping this entry"
      write(*,*)" (Could not open '"//trim(P_output_sub_folder)//trim(clean_data_set_name)//"_jmol_PM6_fs.html')"
      cycle
    end if
    j = j + abs(i)
    if (j > 0) then
      j = 0
      cycle
    end if
    write(jmol,'(a)')"<!DOCTYPE html>", &
    '<HTML>', &
    '<meta charset="utf-8"> <script type="text/javascript" src="'//trim(relative_folder)//'JSmol.min.js"></script>'
    write(jmol,*)"Time stamp: "//timestamp
    write(jmol,'(a)')'<TITLE>'//trim(data_set_name(counter))//'</TITLE>', ' '
    write(jmol,'(a)') &!  First script
    '<script type="text/javascript">', &
    '$(document).ready(function() {InfoRef = {width: 560, height: 560, color: "0x000000",', & 
    'disableInitialConsole: true, addSelectionOptions: false,', & 
    'j2sPath: "'// trim(relative_folder)//'j2s",',&
    'jarPath: "'// trim(relative_folder)//'java",',&
    'use: "HTML5", script: "load '//trim(clean_data_set_name)// &
    '_(Ref).xyz; set measurementUnits ANGSTROMS; hBonds calculate; connect 0.8  1.5 (hydrogen) (phosphorus) create;"}', & 
    '$("#Ref").html(Jmol.getAppletHtml("jmolAppletRef",InfoRef))});', & 
    '</script>'
     write(jmol,'(a)') ' ', &!  Second script
    '<script type="text/javascript">', &
    '$(document).ready(function() {InfoPM7 = {width: 560, height: 560, color: "0x000000",', & 
    'disableInitialConsole: true, addSelectionOptions: false,', & 
    'j2sPath: "'// trim(relative_folder)//'j2s",',&
    'jarPath: "'// trim(relative_folder)//'java",',&
    'use: "HTML5", script: "load '//trim(clean_data_set_name)// &
    '_(PM7).xyz; set measurementUnits ANGSTROMS; hBonds calculate; connect 0.8  1.5 (hydrogen) (phosphorus) create;"}', & 
    '$("#PM7").html(Jmol.getAppletHtml("jmolAppletPM7",InfoPM7))});', & 
    '</script>'
    write(jmol,'(a)') ' ', &!  Third script
    '<script type="text/javascript">', &
    '$(document).ready(function() {InfoPM6 = {width: 560, height: 560, color: "0x000000",', & 
    'disableInitialConsole: true, addSelectionOptions: false,', & 
    'j2sPath: "'// trim(relative_folder)//'j2s",',&
    'jarPath: "'// trim(relative_folder)//'java",',&
    'use: "HTML5", script: "load '//trim(clean_data_set_name)// &
    '_(PM6_D3H4).xyz; set measurementUnits ANGSTROMS; hBonds calculate; connect 0.8  1.5 (hydrogen) (phosphorus) create;"}', & 
    '$("#PM6_D3H4").html(Jmol.getAppletHtml("jmolAppletPM6",InfoPM6))});', & 
    '</script>'
    write(jmol,'(a)') ' ','<BODY>'
    i = len_trim(input_data_set) - 1
    write(jmol,"(a,i5,1x,30(a))") "<H3>", &
      counter, trim(data_set_name(counter)), &
   '  </H3>', &
   '  <a href="'//trim(clean_data_set_name_m1)//'_Jmol.html"> ', &
   ' (Previous)</a> &nbsp;&nbsp;&nbsp;'//trim(data_set_name(counter - 1)), '<BR>', &
   ' <a href="../'//input_data_set(:i)//'.html#',trim(clean_data_set_name)//'">', &
   ' (Back)</a>  &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;Elements:'
     call upcase(keywrd, len_trim(keywrd))
     if (keywrd(1:1) /= " ") keywrd = " "//trim(keywrd)
     i = index(keywrd, "Z=") 
     if (i /= 0) then
       Z_PM7 = nint(reada(keywrd,i+2))
     else
       Z_PM7 = 0
     end if
     j = index(keywrd," MERS")
     mers = 0
     k = 0
     i = Index (keywrd(j + 1:), " ") + j
     do l = 1, 3
       j = j + k
       if (l > 1 .and. k == 0) exit
       mers(l) = Nint (reada (keywrd(j:), 1))
       k = Index (keywrd(j:i), ",")
     end do
     if (index(keywrd," Z=") /= 0) then
       i = index(keywrd," Z=")
       kk = nint(reada(keywrd,i))      
       i = index(keywrd," MERS")
       if (i /= 0) then
         i = mers(1)
         do l = 2, 3
           if (mers(l) == 0) exit
           i = i*mers(l)
         end do
         kk = kk*i
         iels(2,:niels) = iels(2,:niels)/kk
       end if
     else
!
!   Work out number of empirical formula units, kk, in the cluster
!   Allow for every prime number up to 29 (Mn)
!
       kk = 1
       do k = 1, 10
         l = primes(k)
         do
           j = 0
           do i = 1, niels
             if (mod(iels(2,i),l) /= 0 .or. iels(2,i) < l) j = 1 
           end do
           if (j == 1 ) exit
           iels(2,:niels) = iels(2,:niels)/l
           kk = kk * l
         end do
       end do
     end if
     do i = 1, niels
       if (elemnt(iels(1,i))(1:1) == " ") then
         write(jmol,"(2a,i4)")'<a href=".././'//elemnt(iels(1,i))(2:2)//'.html#',trim(clean_data_set_name)// &
           '">'//trim(elemnt(iels(1,i)))//"</A>", iels(2,i)
       else
         write(jmol,"(2a,i4)")'<a href=".././'//elemnt(iels(1,i))//'.html#',trim(clean_data_set_name)// &
           '">'//trim(elemnt(iels(1,i)))//"</A>", iels(2,i)
       end if
     end do    
     if (kk > 1) then
       write(jmol,"(a,i4,a)")' (Z = ',kk,')'
       numats(counter) = kk
     else
        numats(counter) = 1
     end if
     write(jmol,"(a)")'&nbsp; &nbsp; (<a href=".././Periodic_table_solids.html">Periodic Table</a>) <BR>'
     write(jmol,"(3a)")'  <a href="'//trim(clean_data_set_name_p1)//'_Jmol.html"> ', &
      ' (Next)</a>  &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;'//trim(data_set_name(counter + 1))
     write(jmol,"(a)")'<BR>'
     line1 = trim(data_set_name(loop))
    do i = 1, len_trim(line1)
      if (line1(i:i) ==  " ") line1(i:i) = "_"
    end do
    inquire (file = trim(P_notes_folder)//trim(line1)//".txt", exist = exists)
    if (exists) then
      close(notes)
!
!  Write out notes on this system
!
      pre = .false.
      open(unit=notes, file=trim(P_notes_folder)//trim(line1)//".txt", action='READ')
      do
        read(notes,"(a600)", iostat=io_stat) line       
        if (io_stat /= 0) exit
        if (index(line,"<PRE>") /= 0) pre = .true.
        if (index(line,"</PRE>") /= 0) pre = .false.
        if (len_trim(line) < 240) then
          line1 = trim(line)
          call upcase(line1, len_trim(line1))
        else
          line1 = " "
        end if
        if (index(line1,"./") /= 0) then 
          i = index(line1,"./")
          line = trim(P_notes_folder)//trim(line1(i + 2:))
          line1 = trim(line)          
          open(unit=notes1, file=trim(line1), action='READ', iostat = i)
          if (i /= 0) then
!
!  Current system: data_set_name(loop)
            write(jmol,'(a)')"DEADLY ERROR DETECTED IN """//trim(data_set_name(loop))//""""

            stop
          end if
          do 
           read(notes1,"(a600)", iostat=io_stat) line       
           if (io_stat /= 0) exit
           write(jmol,"(3a)")trim(line),'<BR>'
          end do  
          close(notes1)
        else
          if (pre) then
             write(jmol,"(3a)")trim(line)
          else
             write(jmol,"(3a)")trim(line),'<BR>'
          end if
        end if     
      end do
    end if
!
!   Write out unit cell parameters
!    
    call get_unit_cell_etc(arc_channel, sum, a, b, c, alpha, beta, gamma, volume, density)
    close (arc_channel)
    write(jmol,"(a)")"<PRE>                      Unit Cell Parameters:       a      b      c   "// &
      "alpha   beta  gamma   Volume  Density       Heat of Formation (Kcal/mol)"
    d_ref(counter) = density
    write(jmol,"(37x,a,3f7.2,3f7.2,f10.2,f7.3, f15.1, a)")"X-ray:  ", a, b, c, alpha, beta, gamma, &
      volume, density, sum, " calc'd using PM7"
    open(unit=file_name, file=trim(PM7_in), action='READ')
    call get_unit_cell_etc(file_name, sum, a, b, c, alpha, beta, gamma, volume, density)
    d_PM7(counter) = density
    if (sum1 > -1.d8) then
    write(jmol,"(37x,a,3f7.2,3f7.2,f10.2,f7.3, f15.1,a,f8.1,a)")"  PM7:  ", &
      a, b, c, alpha, beta, gamma, volume, density, sum, " calc'd using PM7 (ref:",sum1,")"
    else
    write(jmol,"(37x,a,3f7.2,3f7.2,f10.2,f7.3, f15.1,a)")"  PM7:  ", &
      a, b, c, alpha, beta, gamma, volume, density, sum, " calc'd using PM7"
    end if
!
    inquire (file = trim(PM6_in), exist = exists)
    if (exists) then
      open(unit=file_name, file=trim(PM6_in), action='READ')
      do 
      read(file_name,"(a120)", iostat=io_stat)line
        if (index(line, "MOPAC20") /= 0) exit
        if (io_stat /= 0) exit
      end do
      read(file_name,"(a90)", iostat=io_stat)line
      i = index(line, ":")
      PM6_date = line(i - 9:i - 4)//line(i + 6:i + 10)//" @ "//line(i - 2:i +5)
      if (PM6_date(5:5) == " ") PM6_date = PM6_date(:4)//"0"//PM6_date(6:)  
      call get_unit_cell_etc(file_name, sum, a, b, c, alpha, beta, gamma, volume, density)
      do 
        read(file_name ,"(a90)", iostat=io_stat)line
        if (index(line, "FINAL GEOMETRY OBTAINED") /= 0) exit
        if (io_stat /= 0) exit
      end do 
      if (io_stat /= 0) then
        write(*,'(a)')" Fault detected in PM6_D3H4 '.arc' file: '"//trim(data_set_name(counter))//"'"
        write(*,'(a)')" This is a severe error that must be corrected before continuing"
        write(23,'(a)')" Fault detected in PM6_D3H4 '.arc' file: '"//trim(data_set_name(counter))//"'"
        write(23,'(a)')" This is a severe error that must be corrected before continuing"
        stop    
      end if
      read(file_name ,"(a120)", iostat=io_stat)keywrd! keywrd
      read(file_name ,"(a120)", iostat=io_stat)title! keywrd
      read(file_name ,"(a120)", iostat=io_stat)line ! keywrd
      close(file_name)
      call upcase(keywrd, len_trim(keywrd))
      i = index(keywrd, "Z=") 
      if (i /= 0) then
        Z_PM6 = nint(reada(keywrd, i+2))
      else
        Z_PM6 = 0
      end if
      if (Z_ref /= Z_PM7 .or. Z_PM7 /= Z_PM6) then
        if (first) then
          first = .false.
          write(*, '(a)')"               Values of Z for Ref, PM7, and PM6_D3H4, are different"
          write(23,'(a)')"               Values of Z for Ref, PM7, and PM6_D3H4, are different"
        end if
        write(*,'(a,i4, i5, i9)')"  '"//trim(data_set_name(counter))//"':",Z_ref, Z_PM7, Z_PM6
        write(23,'(a,3i4)')"  '"//trim(data_set_name(counter))//"':",Z_ref, Z_PM7, Z_PM6
      end if
      d_PM6(counter) = density
      write(jmol,"(32x,a,3f7.2,3f7.2,f10.2,f7.3, f15.1,a)")"  PM6_D3H4:  ", &
        a, b, c, alpha, beta, gamma, volume, density, sum, " calc'd using PM6_D3H4"
    end if
    write(jmol,"(a)")'</PRE>' 
!
!  Jmol script
!
     line  = trim(P_PM6_arc_files)//trim(data_set_name(counter))//" (PM6-D3H4).arc"
    inquire (file = trim(line), exist = exists)
    if (exists) then
      i = len_trim(data_set_name(counter))
      if (i < 147) then
        write(jmol,"(30(a,/))") & 
        '<PRE>                                 <font size="5"><b> <a href="'//trim(clean_data_set_name)// &
          '_jmol_XRay_fs.html"  target="_blank">X-Ray</a> '// &
        '                                            <a href="'//trim(clean_data_set_name)// &
        '_jmol_PM7_fs.html"  target="_blank">PM7</a>                                            <a href="' &
          //trim(clean_data_set_name)//'_jmol_PM6_fs.html"  target="_blank">PM6_D3H4</a></PRE>'
      else
        write(jmol,"(30(a,/))") & 
        '<PRE>                                 <font size="5"><b> <a href="'//trim(clean_data_set_name)// &
          '_jmol_XRay_fs.html"  target="_blank">X-Ray</a> '// &
        '                                            <a href="'//trim(clean_data_set_name)// &
        '_jmol_PM7_fs.html"  target="_blank">PM7</a>                                            PM6_D3H4</PRE>'
      end if
    else
        write(jmol,"(30(a,/))") &
        '<PRE>                                 <font size="5"><b> <a href="'//trim(clean_data_set_name)// &
          '_jmol_XRay_fs.html"  target="_blank">X-Ray</a> '// &
        '                                            <a href="'//trim(clean_data_set_name)// &
        '_jmol_PM7_fs.html"  target="_blank">PM7</a>                          </PRE>'
    end if
    if (counter == 16) then
      P_ok = .false.
    end if
    line2 = " "
    
    H_ok  = .false. 
    Be_ok = .false.
    C_ok  = .false.
    N_ok  = .false.
    O_ok  = .false.
    F_ok  = .false.
    Mg_ok = .false.
    P_ok  = .false.
    S_ok  = .false.
    Cl_ok = .false.
    Ca_ok = .false.
    Zn_ok = .false.
    Br_ok = .false.
    Sr_ok = .false.
    Mo_ok = .false.
    I_ok  = .false.
    do i = 1, numat - 3
      if (.not. H_ok)  H_ok  = (nat(i) == 1)
      if (.not. Be_ok) Be_ok = (nat(i) == 4)
      if (.not. C_ok)  C_ok  = (nat(i) == 6)
      if (.not. N_ok)  N_ok  = (nat(i) == 7)
      if (.not. O_ok)  O_ok  = (nat(i) == 8)
      if (.not. F_ok)  F_ok  = (nat(i) == 9)
      if (.not. Mg_ok) Mg_ok = (nat(i) == 12)
      if (.not. P_ok)  P_ok  = (nat(i) == 15)
      if (.not. S_ok)  S_ok  = (nat(i) == 16)
      if (.not. Cl_ok) Cl_ok = (nat(i) == 17)
      if (.not. Ca_ok) Ca_ok = (nat(i) == 20)
      if (.not. Zn_ok) Zn_ok = (nat(i) == 30)
      if (.not. Br_ok) Br_ok = (nat(i) == 35)
      if (.not. Sr_ok) Sr_ok = (nat(i) == 38)
      if (.not. Mo_ok) Mo_ok = (nat(i) == 42)
      if (.not. I_ok)  I_ok  = (nat(i) == 53)
    end do
    if (O_ok  .and. H_ok)  write(line2(len_trim(line2) + 1:), "(a)")"hBonds calculate; "
    if (Be_ok .and. F_ok)  write(line2(len_trim(line2) + 1:), "(a)")"connect 0.8  1.7 (beryllium) (fluorine) create; "
    if (Mg_ok .and. N_ok)  write(line2(len_trim(line2) + 1:), "(a)")"connect 0.8  2.3 (magnesium) (nitrogen) create; "
    if (Mg_ok .and. O_ok)  write(line2(len_trim(line2) + 1:), "(a)")"connect 0.8  2.3 (magnesium) (oxygen) create; "
    if (P_ok  .and. H_ok)  write(line2(len_trim(line2) + 1:), "(a)")"connect 0.8  1.5 (hydrogen) (phosphorus) create; "
    if (P_ok  .and. C_ok)  write(line2(len_trim(line2) + 1:), "(a)")"connect 0.8  2.1 (carbon) (phosphorus) create; "
    if (P_ok  .and. S_ok)  write(line2(len_trim(line2) + 1:), "(a)")"connect 0.8  2.4 (sulfur) (phosphorus) create; "
    if (Be_ok .and. Cl_ok) write(line2(len_trim(line2) + 1:), "(a)")"connect 0.8  2.1 (beryllium) (chlorine) create; "
    if (Ca_ok .and. Br_ok) write(line2(len_trim(line2) + 1:), "(a)")"connect 0.8  3.0 (calcium) (bromine) create; "
    if (Zn_ok .and. Zn_ok) write(line2(len_trim(line2) + 1:), "(a)")"connect 3.2  3.6 (zinc) (zinc) delete; "
    if (Be_ok .and. Br_ok) write(line2(len_trim(line2) + 1:), "(a)")"connect 0.8  2.2 (beryllium) (bromine) create; "
    if (Sr_ok .and. I_ok)  write(line2(len_trim(line2) + 1:), "(a)")"connect 0.8  3.45 (strontium) (iodine) create; "
    if (Mo_ok .and. C_ok)  write(line2(len_trim(line2) + 1:), "(a)")"connect 0.8  2.8 (molybdenum) (carbon) create; "
    write(line2(len_trim(line2) + 1:), "(a,9f9.3,a)")'set measurementUnits ANGSTROMS; set zoomLarge false; ")</script>'
    line = trim(P_PM7_arc_files)//trim(data_set_name(counter))//".arc"
    open(arc_channel, file=trim(line), action='READ', iostat = io_stat)
    if (io_stat /= 0) then
      write(*,'(a)')" File '"//trim(line)//"' is inaccessible - skipping to next file."
      cycle
    end if
    rewind (arc_channel)
    do 
      read(arc_channel,"(a120)", iostat=io_stat)line
      if (index(line, "MOPAC20") /= 0) exit
      if (io_stat /= 0) exit
    end do
    read(arc_channel,"(a90)", iostat=io_stat)line
    i = index(line, ":")
    PM7_date = line(i - 9:i - 4)//line(i + 6:i + 10)//" @ "//line(i - 2:i +5)
    if (PM7_date(5:5) == " ") PM7_date = PM7_date(:4)//"0"//PM7_date(6:) 
!
!  Write table 
! 
    write(jmol,"(a)") '<TABLE>'
    write(jmol,"(a)") '<td><span id=Ref></span><a href="javascript:Jmol.script(jmolAppletRef)"></a>' 
    write(jmol,"(a)") '<font face="Courier New" style="font-size: 14pt"> Calc''d on: '//Xray_date//'</font></td>'
    write(jmol,"(a)") '<td><span id=PM7></span><a href="javascript:Jmol.script(jmolAppletPM7)"></a>'
    write(jmol,"(a)") '<font face="Courier New" style="font-size: 14pt"> '//PM7_date
    line = "./"//trim(data_set_name(counter))//'_(PM7).html'
    write(jmol,"(a)") '<a href="'//trim(line)//'"  target="_blank">ARC file</a></font></td>'
    if (exists) &
    write(jmol,"(a)") '<td><span id=PM6_D3H4></span><a href="javascript:Jmol.script(jmolAppletPM6)"></a>'
    write(jmol,"(a)") '<font face="Courier New" style="font-size: 14pt"> '//PM6_date//'</font></td>'
    write(jmol,"(a)") '</TABLE>'
!
!  Write out the X-ray structure to an xyz file
!     
    nat_store(:numat) = nat(:numat)!  Store X-ray atom sequence
    j = numat - 3
    line1 = trim(P_output_sub_folder)//trim(working)//'_(Ref).xyz'
    current_data_set_name = trim(data_set_name(counter))
    call write_xyz("X_R     ", line1, Xray_a, Xray_b, Xray_c, line, hof_ref1, hof_calc) 
    if (j /= numat) then
      write(*,'(a,i3,i4,a)')" Number of atoms in PM7 and X-ray structure '.arc' file: '"// &
        trim(data_set_name(counter))//"' are different. (",j, numat,")" 
      write(23,'(a,i3,i4,a)')" Number of atoms in PM7 and X-ray structure '.arc' file: '"// &
        trim(data_set_name(counter))//"' are different. (",j, numat,")" 
    end if    
    j = 0
    do i = 1, numat
      j = j + abs(nat_store(i) -  nat(i))
    end do
    if (j /= 0) then
      write(*,'(a)')" Fault detected in PM7 or X-ray structure '.arc' file: '"//trim(data_set_name(counter))//"'"  
      write(23,'(a)')" Fault detected in PM7 or X-ray structure '.arc' file: '"//trim(data_set_name(counter))//"'"  
    end if
    tvec_ref(:,1) = tvec(:,1)/mers(1)
    tvec_ref(:,2) = tvec(:,2)/mers(2)
    tvec_ref(:,3) = tvec(:,3)/mers(3)
    tvec_ref = tvec
    write(jmol_XRay_fs,'(a)')"<!DOCTYPE html>", &
    '<HTML>', &
    '<meta charset="utf-8"> <script type="text/javascript" src="'//trim(relative_folder)//'JSmol.min.js"></script>', &
    '<TITLE>'//trim(data_set_name(counter))//' (Ref.)</TITLE>', ' '
    write(jmol_XRay_fs,'(a)') &
    '<script type="text/javascript">', &
    '$(document).ready(function() {InfoRef = { color: "0x000000", height: "800", width: "1600",', &
    'disableInitialConsole: true, addSelectionOptions: false, ', &
    'j2sPath: "'// trim(relative_folder)//'j2s",',&
    'jarPath: "'// trim(relative_folder)//'java",'
    write(jmol_XRay_fs,"(a,9f9.3,30a)")'use: "HTML5", script: "load 	'//trim(clean_data_set_name)// &
    '_(Ref).xyz {1,1,1} packed unitcell {',tvec_Ref,'};  \'
    write(jmol_XRay_fs,'(a)') &
    'set measurementUnits ANGSTROMS; hBonds calculate; connect 0.8  1.5 (hydrogen) (phosphorus) create; set zoomLarge false;"} ', &
    '$("#Ref").html(Jmol.getAppletHtml("jmolAppletRef",InfoRef))}); ', &
    '</script>', &   
    '</head>', &
    '<body>', &
    '<span id=Ref></span><a href="javascript:Jmol.script(jmolAppletRef)"></a>', &
    '</body>', &
    '</html>'
    close (jmol_XRay_fs, STATUS = "keep")
    ref_coord(:,:numat) = coord(:,:numat)
!
! set numats(counter) equal to the number of atoms in a fundamental unit, i.e., the empirical formula
!
    i = numats(counter)! delete
    call empire (formla(counter), nat, els(counter))
    numats(counter) = numat/numats(counter)
!
! Write out the PM7 structure to an xyz file
!
    line1 = trim(P_output_sub_folder)//trim(working)//'_(PM7).xyz' 
    call write_xyz("PM7     ", line1, PM7_a, PM7_b, PM7_c, PM7_date, hof_ref1, hof_calc)
    if (hof_ref1 > -1.d9) then
      hof_ref(counter) = hof_ref1
    end if
    if (hof_calc > -1.d9 ) then
      hof_PM7(counter) = hof_calc
    end if
    tvec_PM7(:,1) = tvec(:,1)/mers(1)
    tvec_PM7(:,2) = tvec(:,2)/mers(2)
    tvec_PM7(:,3) = tvec(:,3)/mers(3)
    tvec_PM7 = tvec
    write(jmol_PM7_fs,'(a)')"<!DOCTYPE html>", &
    '<HTML>', &
    '<meta charset="utf-8"> <script type="text/javascript" src="'//trim(relative_folder)//'JSmol.min.js"></script>', &
    '<TITLE>'//trim(data_set_name(counter))//' (PM7) </TITLE>', ' '
    write(jmol_PM7_fs,'(a)') &
    '<script type="text/javascript">', &
    '$(document).ready(function() {InfoRef = { color: "0x000000", height: "800", width: "1600",', &
    'disableInitialConsole: true, addSelectionOptions: false, ', &
    'j2sPath: "'// trim(relative_folder)//'j2s",',&
    'jarPath: "'// trim(relative_folder)//'java",'
    write(jmol_PM7_fs,"(a,9f9.3,30a)")'use: "HTML5", script: "load 	'//trim(clean_data_set_name)// &
    '_(PM7).xyz {1,1,1} packed unitcell {',tvec_PM7,'};  \'
    write(jmol_PM7_fs,'(a)') &
    'set measurementUnits ANGSTROMS; hBonds calculate; connect 0.8  1.5 (hydrogen) (phosphorus) create; set zoomLarge false;"} ', &
    '$("#Ref").html(Jmol.getAppletHtml("jmolAppletRef",InfoRef))}); ', &
    '</script>', &   
    '</head>', &
    '<body>', &
    '<span id=Ref></span><a href="javascript:Jmol.script(jmolAppletRef)"></a>', &
    '</body>', &
    '</html>'
    close (jmol_PM7_fs, STATUS = "keep")
    if (PM7_a > 0.1d0 .and. Xray_a > 0.1d0) then
      PM7_accuracy = accuracy(ref_coord, coord, Xray_a, Xray_b, Xray_c, PM7_a, PM7_b, PM7_c)
    else
      PM7_accuracy = 0.d0
    end if
!
!  If exists, write out the PM6-D3H4 structure to an xyz file
!
    line = trim(P_PM6_arc_files)//trim(data_set_name(counter))//" (PM6-D3H4).arc"
    inquire (file = trim(line), exist = exists)
    PM6_accuracy = 00.d0
    if (exists) then
      line1 = trim(P_output_sub_folder)//trim(working)//'_(PM6_D3H4).xyz'
      pm6_input = trim(data_set_name(counter))//' (PM6-D3H4).arc'
      call write_xyz("PM6_D3H4", line1, PM6_a, PM6_b, PM6_c, line, hof_ref1, hof_calc)
      j = 0
      do i = 1, numat
        j = j + abs(nat_store(i) -  nat(i))
      end do
      if (j /= 0) then
        write(*,'(a)')" Fault detected in PM6_D3H4 '.arc' file: '"//trim(data_set_name(counter))//"': sequence of atoms is different"
        write(23,'(a)')" Fault detected in PM6_D3H4 '.arc' file: '"//trim(data_set_name(counter))//"': sequence of atoms is different"   
      end if
      tvec_PM6(:,1) = tvec(:,1)/mers(1)
      tvec_PM6(:,2) = tvec(:,2)/mers(2)
      tvec_PM6(:,3) = tvec(:,3)/mers(3)
      tvec_PM6 = tvec
     write(jmol_PM6_fs,'(a)')"<!DOCTYPE html>", &
    '<HTML>', &
    '<meta charset="utf-8"> <script type="text/javascript" src="'//trim(relative_folder)//'JSmol.min.js"></script>', &
    '<TITLE>'//trim(data_set_name(counter))//' (PM6_D3H4) </TITLE>', ' '
    write(jmol_PM6_fs,'(a)') &
    '<script type="text/javascript">', &
    '$(document).ready(function() {InfoRef = { color: "0x000000", height: "800", width: "1600",', &
    'disableInitialConsole: true, addSelectionOptions: false, ', &
    'j2sPath: "'// trim(relative_folder)//'j2s",',&
    'jarPath: "'// trim(relative_folder)//'java",'
     write(jmol_PM6_fs,"(a,9f9.3,30a)")'use: "HTML5", script: "load 	'//trim(clean_data_set_name)// &
    '_(PM6_D3H4).xyz {1,1,1} packed unitcell {',tvec_PM6,'};  \'
    write(jmol_PM6_fs,'(a)') &
    'set measurementUnits ANGSTROMS; hBonds calculate; connect 0.8  1.5 (hydrogen) (phosphorus) create; set zoomLarge false;"} ', &
    '$("#Ref").html(Jmol.getAppletHtml("jmolAppletRef",InfoRef))}); ', &
    '</script>', &   
    '</head>', &
    '<body>', &
    '<span id=Ref></span><a href="javascript:Jmol.script(jmolAppletRef)"></a>', &
    '</body>', &
    '</html>'
      close (jmol_PM6_fs, STATUS = "keep")      
      if (hof_calc > -1.d9) then 
        hof_PM6(counter) = hof_calc
      end if     
      if (PM6_a > 0.1d0 .and. Xray_a > 0.1d0) then
        PM6_accuracy = accuracy(ref_coord, coord, Xray_a, Xray_b, Xray_c, PM6_a, PM6_b, PM6_c)
      else
        PM6_accuracy = 0.d0
      end if
    else
      write(*,"(a)") "  PM6_D3H4 arc   file:  "//trim(line)//" does not exist"
      write(23,"(a)") "  PM6_D3H4 arc   file:  "//trim(line)//" does not exist" 
    end if
    write(jmol,"(a)")'</BODY></HTML>'
    close(jmol, STATUS = "keep")

!
!  Write master list of all molecules
!
    lim1 = 2.d0
    lim2 = 5.d0
    lim3 = 10.d0
    lim4 = 20.d0
    lim5 = 50.d0
    if (counter < 10) then
      txt_counter = 4
    else if (counter < 100) then
      txt_counter = 3
    else if (counter < 1000) then
      txt_counter = 2
    else
      txt_counter = 1
    end if
!
!  Write text for start of line
!
    txt_start = trim(line)
    write(txt_start,'(a,i5, a)') nbsp(:txt_counter*6), counter, &
        nbsp(:6)//' <a href="./data_solids/'//trim(clean_data_set_name)//'_Jmol.html">(JSmol)</a>'
!
! Write text for Ref. HoF
!
    j = numb(hof_PM7(counter), 1)
    j = 8 - j
    if (hof_ref(counter) > -1.d9) then
       i = numb(hof_ref(counter), 1)
       i = 4 - i
      write(txt_ref,'(a, f8.1)')nbsp(:i*6), hof_ref(counter)
      sum = (Abs(hof_ref(counter) - hof_PM7(counter)))
      if (sum < lim1) then
        line1 = '<B><FONT color="green">'
        line2 = '</FONT></B>'
        score_PM7(1) = score_PM7(1) + 1
      else if (sum < lim2) then
        line1 = '<FONT color="green">'
        line2 = '</FONT>'
        score_PM7(2) = score_PM7(2) + 1
      else if ((sum < lim3)) then
        line1 = '<B> '
        line2 = '</B>'
        score_PM7(3) = score_PM7(3) + 1
      else if ((sum < lim4)) then
        line1 = '<FONT color="red">'
        line2 = '</FONT>'
        score_PM7(4) = score_PM7(4) + 1
      else if ((sum < lim5)) then
        line1 = '<B><FONT color="red">'
        line2 = '</FONT></B>'
        score_PM7(5) = score_PM7(5) + 1
      else 
        line1 = '<B><blink><FONT color="red">'
        line2 = '</FONT></B>'
        score_PM7(6) = score_PM7(6) + 1
      end if
      write(txt_PM7,'(a, f9.1, a)') nbsp(:j*6)//trim(line1), hof_PM7(counter), trim(line2)
    else
      write(txt_ref,'(a, f9.1)')nbsp(:8*6)
      write(txt_PM7,'(a, f9.1)')nbsp(:j*6), hof_PM7(counter)
    end if      
    if (hof_PM6(counter) > -1.d9) then
      j = numb(hof_PM6(counter), 1)
      j = 7 - j
      if (hof_ref(counter) > -1.d9) then
        sum = (Abs(hof_ref(counter) - hof_PM6(counter)))
        if (sum < lim1) then
          line1 = '<B><FONT color="green">'
          line2 = '</FONT></B>'
          score_PM6(1) = score_PM6(1) + 1
        else if (sum < lim2) then
          line1 = '<FONT color="green">'
          line2 = '</FONT>'
          score_PM6(2) = score_PM6(2) + 1
        else if ((sum < lim3)) then
          line1 = '<B> '
          line2 = '</B>'
          score_PM6(3) = score_PM6(3) + 1
        else if ((sum < lim4)) then
          line1 = '<FONT color="red">'
          line2 = '</FONT>'
          score_PM6(4) = score_PM6(4) + 1
        else if ((sum < lim5)) then
          line1 = '<B><FONT color="red">'
          line2 = '</FONT></B>'
          score_PM6(5) = score_PM6(5) + 1
        else 
          line1 = '<B><FONT color="red">'
          line2 = '</FONT></B>'
          score_PM6(6) = score_PM6(6) + 1
        end if
      else
        line1 = " "
        line2 = " "
      end if
      write(txt_PM6,'(a, f10.1, a)')nbsp(:j*6)//trim(line1), hof_PM6(counter), trim(line2)
    else
      write(txt_PM6,'(a, a, a)')nbsp(:8*6), "-.-"    
    end if
    PM6_acc_prt = 100.d0 - max(0.1d0, PM6_accuracy)
    if (PM6_acc_prt < 9.95d0) then
      j = 3
    else
      j = 2
    end if
    if (PM6_acc_prt > 99.d0) then
      write(PM6_geo,'(a, a)')nbsp(:4*6), "-.-"
    else      
    write(PM6_geo,'(a, f9.1)')nbsp(:j*6), PM6_acc_prt
    end if
    
    PM7_acc_prt = 100.d0 - max(0.1d0, PM7_accuracy)
    if (PM7_acc_prt < 9.95d0) then
      k = 7
    else
      k = 6
    end if
    write(PM7_geo,'(a, f9.1)')nbsp(:k*6), PM7_acc_prt
    line = trim(txt_start)//trim(txt_ref)//trim(txt_PM7)//trim(txt_PM6)// &
      trim(PM7_geo)//trim(PM6_geo)//nbsp(:3*6)// &
      trim(data_set_name(counter))//'<a name="'//trim(clean_data_set_name)//'"></a> <BR>'
!
!  Set elements that are present in this solid
!
    do l = 1, natoms
      elements(labels(l),counter) = .true.
      elements(labels(l),0) = .true.
    end do
    do l = 1, 12
      if (PM7_date(:3) == months(l)) exit
    end do
    sum = days_in_month(l) + reada(PM7_date, 4) + 0.0005d0
    sum1 = 0.d0
    if (sum > sum1 + 183.d0) sum = sum - 365.d0
      write(all_solids(counter),"(a)") trim(line)
    all_accuracy_PM7(counter) = PM7_accuracy
    all_accuracy_PM6(counter) = PM6_accuracy
     if (hof_PM7(counter) < -1.d9 .or. hof_ref(counter) < -1.d9) cycle
        line1 = " "
!
!  Write formula and name, and an "@" sign
!
        i = len_trim(data_set_name(counter))     
        write(line,'(a15,a,a)') &
        formla(counter), '<a href="./data_solids/'//trim(clean_data_set_name)//'_jmol.html">'// &
        trim(data_set_name(counter))//"</a>", pad(:70 - i)//"@"
!
! Write the numerical data and element symbols
!
        j = len_trim(line)
        write(line(j + 1:),'(f9.2,f10.2,f9.2,9x,a,a)') Hof_ref(counter), hof_PM7(counter), &
        hof_PM7(counter) - hof_ref(counter), els(counter)
        k = index(line, "</a>") - i
        if (k > maxtxt) then
          maxtxt = k
        end if
        write(hof_PM7_all(counter),'(a)')trim(line)
        if (hof_PM6(counter) < -1.d9 .or. hof_ref(counter) < -1.d9) cycle
    end do 
  return
  end subroutine bigloop
  