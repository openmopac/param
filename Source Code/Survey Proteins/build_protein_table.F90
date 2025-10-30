  program build_protein_table
!
! Given folders of .arc and .pdb files, make a table in HTML
!
    use folder_names_C, only : location, relative_location, input_data, arc_files
    implicit none
    character :: line*600, nice_line*600, one_line(100)*600, pre_font*200, post_font*200, file_name*120, &
      method1*8, method2*3, line1*600, timestamp*24
    integer :: i, j, k, ir = 5, iw = 6, io_stat, idata = 7, nlines, loop, &
      natoms, molecule, datum, n_3, n_10, n_opt
    logical :: exists
    double precision :: g_lim, sum, &
      PM_gnorm_pdb,   PM_gnorm_3,  PM_gnorm_10,  PM_gnorm_opt,  &
      PM_hof_pdb,     PM_hof_3,    PM_hof_10,    PM_hof_opt,    &
      PM_rms_pdb,     PM_rms_3,    PM_rms_10,    PM_rms_opt,    &
      PM_cutoff_pdb,  PM_cutoff_3, PM_cutoff_10, PM_cutoff_opt, &
      PM_charge_pdb,  PM_charge_3, PM_charge_10, PM_charge_opt, &
      sum_3, sum_10, sum_opt
!
!  Location and relative location of various folders
!
    input_data = "M:/PARAM/Source Code/Survey Proteins/" ! Location of data-files used by the program to control events
    location = "M:/PARAM/HTML Files/"                    ! Location of folder that will hold the HTML files
    relative_location = "../../jsmol/"                   ! Location of the folder that holds the JSmol files, relative to "location"
    arc_files = "M:\PARAM\Analysis\Proteins\"            ! Location of the ARC files used in this analysis 
    call fdate(timestamp)
    g_lim = 22.d0 ! 0.6d0
    method1 = "PM7"
    do loop = 1, 2
      open(unit = idata, file = trim(location)//trim(method1)//"_PDB", iostat = io_stat)
      open(unit = iw, file = trim(location)//trim(method1)//"_Proteins.html", iostat = io_stat)
!
!  Write HTML header
!
      open(unit = ir, file = trim(input_data)//"header_for_table.txt", iostat = io_stat)
      write(iw,'(a)')"<HTML>"
      write(iw,*)"Time stamp: "//timestamp
      do k = 1, 3
      read(ir,'(a)', iostat = io_stat)line
      if (io_stat /= 0) exit
      i = index(line, "PM7")
      if (i /= 0) then
        j = i
        do
          i = index(line(j:), "PM7")
          if (i == 0) exit
          i = i + j - 1
          line(i:) = trim(method1)//line(i+3:)
          j = j + i - 1
        end do
      end if
      write(iw,'(a)')trim(line)
      end do
      do
        read(ir,'(a)', iostat = io_stat)line
        if (io_stat /= 0) exit
        write(iw,'(a)')trim(line)
      end do
!
!  Read in "one_line"  This formats one line of the HTML file
!
      open(unit = ir, file = trim(input_data)//"one_line.txt", iostat = io_stat)
      nlines = 0
      do
        read(ir,'(a)', iostat = io_stat)line
        if (io_stat /= 0) exit
        nlines = nlines + 1
        one_line(nlines) = line
      end do
!
! write body of table
!
      open(unit = ir, file = trim(input_data)//"List of ARC files to be used.txt", iostat = io_stat)
      sum_3 = 0.d0
      sum_10 = 0.d0
      sum_opt = 0.d0
      n_3 = 0
      n_10 = 0
      n_opt = 10
      do molecule = 1, 10000
    !    if (debugging .and. molecule == 2) exit
        read(ir,'(a)', iostat = io_stat)line
        if (io_stat /= 0) exit
        j = index(line,'.arc') 
        nice_line = trim(line(:j - 1))
        file_name = trim(nice_line)
        do i = 1, j
          if (nice_line(i:i) == "_") nice_line(i:i) = " "
        end do
        PM_cutoff_pdb = 0.d0
        PM_cutoff_3   = 0.d0
        PM_cutoff_10  = 0.d0
        PM_cutoff_opt = 0.d0
        PM_charge_pdb = 0.d0
        PM_charge_3   = 0.d0
        PM_charge_10  = 0.d0
        PM_charge_opt = 0.d0
! 
!  Get data on this system
!
        line1 = trim(location)//trim(method1)//"_PDB\"//trim(line)
        inquire (file = trim(line1), exist = exists)
        if (exists) then
          open(unit = idata, file = trim(line1))
          call getdata(natoms, PM_gnorm_pdb, PM_hof_pdb, PM_rms_pdb, idata, PM_cutoff_pdb, PM_charge_pdb)
        else
          PM_hof_PDB = 0.d0
        end if
!
        line1 = trim(arc_files)//trim(method1)//"_3\"//trim(line)
        inquire (file = trim(line1), exist = exists)
        if (exists) then
          open(unit = idata, file = trim(line1))
          call getdata(natoms, PM_gnorm_3, PM_hof_3, PM_rms_3, idata, PM_cutoff_3, PM_charge_3)
        else
          PM_gnorm_3  = 0.d0
          PM_hof_3    = 0.d0
          PM_rms_3    = 0.d0
        end if
!
        line1 = trim(arc_files)//trim(method1)//"_10\"//trim(line)
        inquire (file = trim(line1), exist = exists)
        if (exists) then
          open(unit = idata, file = trim(line1))
          call getdata(natoms, PM_gnorm_10, PM_hof_10, PM_rms_10, idata, PM_cutoff_10, PM_charge_10)
        else
          PM_gnorm_10  = 0.d0
          PM_hof_10    = 0.d0
          PM_rms_10    = 0.d0
        end if
!
        line1 = trim(arc_files)//trim(method1)//"_"//trim(method1)//"\"//trim(line)
        inquire (file = trim(line1), exist = exists)
        if (exists) then
          open(unit = idata, file = trim(line1))
          call getdata(natoms, PM_gnorm_opt, PM_hof_opt, PM_rms_opt, idata, PM_cutoff_opt, PM_charge_opt)
        else
          PM_gnorm_opt  = 0.d0
          PM_hof_opt    = 0.d0
          PM_rms_opt    = 0.d0
        end if
        i = 0
        j = 0
        k = 0
        if (PM_cutoff_pdb > 0.d0 .or. PM_cutoff_3 > 0.d0 .or. PM_cutoff_10 > 0.d0 .or. PM_cutoff_opt > 0.d0) then
          if (abs(PM_cutoff_pdb - PM_cutoff_3) > 0.01d0 .or. abs(PM_cutoff_pdb - PM_cutoff_10) > 0.01d0 .or. &
              abs(PM_cutoff_pdb - PM_cutoff_10) > 0.01d0 .or. abs(PM_cutoff_pdb - PM_cutoff_opt) > 0.01d0) i = 1
        end if
        sum = PM_hof_3 - PM_hof_opt
        if (sum < 0.d0) j = j + 1
        sum = PM_hof_10 - PM_hof_3
        if (sum < 0.d0) j = j + 1
        sum = PM_hof_PDB - PM_hof_10
        if (sum < 0.d0) j = j + 1
        sum = PM_charge_3 - PM_charge_opt
        if (abs(sum) > 0.1d0) k = k + 1
        sum = PM_charge_10 - PM_charge_opt
        if (abs(sum) > 0.1d0) k = k + 1
        sum = PM_charge_PDB - PM_charge_opt
        if (abs(sum) > 0.1d0) k = k + 1
        if (.false. .and. i + j + k > 0) then
          inquire(unit=74, opened = exists) 
          if (.not. exists) then
            open(unit = 74, file = trim(location)//'Protein_table_bugs.txt')
            write(74,'(/,a)') &
"                   Protein                                PDB           10             3 "// &
            "Opt           Cutoff           Charge"
          end if
          write(line,'(a50,4f14.3, 8i4)')trim(method1)//"  "//line(:40), PM_hof_PDB, PM_hof_10, PM_hof_3, PM_hof_opt, &
            nint(PM_cutoff_PDB), nint(PM_cutoff_10), nint(PM_cutoff_3), nint(PM_cutoff_opt), &
            nint(PM_charge_PDB), nint(PM_charge_3), nint(PM_charge_10), nint(PM_charge_opt)
          if (i == 0) line(108:123) = " "
          if (j == 0) line(50:108)  = " "
          if (k == 0) line(123:)    = " "
          write(74,'(a)')trim(line)
        end if
        j = 0
        do i = 1, nlines
          line = one_line(i)
          k = index(line,"%") 
          if (k /= 0) then
            j = j + 1
            datum = 1
            if (j == 1) then
              write(line(k:),'(i3)')molecule
              line(len_trim(line) + 2:) = one_line(i)(k + 1:)
            end if
            datum = datum + 1
            if (j == datum) then
              line1 = 'PM7_PDB/Notes on '//trim(nice_line)//'.html'
              inquire (file=trim(location)//trim(line1), exist = exists)
              if (exists) then
                line = '<p align="left"><a href="'//trim(line1)//'" target="_blank">'//trim(nice_line)//'</a></td>'
              else
                line(k:) = trim(nice_line)
                line(len_trim(line) + 1:) = one_line(i)(k + 1:)
              end if
            end if
            if (abs(PM_hof_PDB) < 0.1d0) then
              write(iw,'(a)')trim(line)
              cycle
            end if
            datum = datum + 1
            if (j == datum) then
              write(line(k:),'(i5)')natoms
              line(len_trim(line) + 1:) = one_line(i)(k + 1:)
            end if

!
!   Start of <method>_PDB
!
            datum = datum + 1
            if (j == datum) then
              if (PM_gnorm_pdb > g_lim*sqrt(1.d0*natoms) .or. PM_hof_pdb > -1.d0) then
                pre_font = '<FONT color="red">'
                post_font = '</FONT>'
              else
                pre_font = ' '
                post_font = ' '
              end if       
               write(line(k:),'(a,f8.1,a)')trim(pre_font),PM_hof_pdb,trim(post_font)
              line(len_trim(line) + 1:) = one_line(i)(k + 1:)
            end if   
            
            datum = datum + 1
            if (j == datum) then
              line(k:) = &
                '<a href="'//trim(method1)//'_PDB/'//trim(file_name)//'.html" target="_blank">'//'Jmol</a> '// &
                '<a href="'//trim(method1)//'_PDB/'//trim(file_name)//'.arc">'//'ARC</a> '// &
                '<a href="'//trim(method1)//'_PDB/'//trim(file_name)//'.pdb">'//'PDB</a> '
              line(len_trim(line) + 1:) = one_line(i)(k + 1:)
            end if
!
!   End of <method>_PDB
!

!
!   Start of <method>_10
!
            datum = datum + 1
            if (j == datum) then 
              n_10 = n_10 + 1
              sum_10 = sum_10 + PM_rms_10
               write(line(k:),'(f6.3)')PM_rms_10
              line(len_trim(line) + 1:) = one_line(i)(k + 1:)
            end if

            datum = datum + 1
            if (j == datum) then
              if (PM_gnorm_10 > g_lim*sqrt(1.d0*natoms) .or. PM_hof_10 > -1.d0) then
                pre_font = '<FONT color="red">'
                post_font = '</FONT>'
              else
                pre_font = ' '
                post_font = ' '
              end if       
              write(line(k:),'(a,f8.1,a)')trim(pre_font),PM_hof_10,trim(post_font)
              line(len_trim(line) + 1:) = one_line(i)(k + 1:)
            end if
            
            datum = datum + 1
            if (j == datum) then
              if (abs(PM_hof_10) < 0.1d0) then
                line = "&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; "
              else
                line(k:) = &
                  '<a href="'//trim(method1)//'_10/'//trim(file_name)//'.html" target="_blank">'//'Jmol</a> '// &
                  '<a href="'//trim(method1)//'_10/'//trim(file_name)//'.arc">'//'ARC</a> '// &
                  '<a href="'//trim(method1)//'_10/'//trim(file_name)//'.pdb">'//'PDB</a> '
                line(len_trim(line) + 1:) = one_line(i)(k + 1:)
              end if
            end if      
!
!   End of <method>_10
!
       
!
!   Start of <method>_3
!
           

            datum = datum + 1
            if (j == datum) then
              n_3 = n_3 + 1
              sum_3 = sum_3 + PM_rms_3
              write(line(k:),'(f6.3)')PM_rms_3
              line(len_trim(line) + 1:) = one_line(i)(k + 1:)
            end if 

            datum = datum + 1
            if (j == datum) then
            if (PM_gnorm_3 > g_lim*sqrt(1.d0*natoms) .or. PM_hof_3 > -1.d0) then
                pre_font = '<FONT color="red">'
                post_font = '</FONT>'
              else
                pre_font = ' '
                post_font = ' '
              end if       
              write(line(k:),'(a,f8.1,a)')trim(pre_font),PM_hof_3,trim(post_font)
                  line(len_trim(line) + 1:) = one_line(i)(k + 1:)
            end if
            
            datum = datum + 1
            if (j == datum) then
              if (abs(PM_hof_3) < 0.1d0) then
                line = "&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; "
              else
                line(k:) = &
                  '<a href="'//trim(method1)//'_3/'//trim(file_name)//'.html" target="_blank">'//'Jmol</a> '// &
                  '<a href="'//trim(method1)//'_3/'//trim(file_name)//'.arc">'//'ARC</a> '// &
                  '<a href="'//trim(method1)//'_3/'//trim(file_name)//'.pdb">'//'PDB</a> '
                line(len_trim(line) + 1:) = one_line(i)(k + 1:)
              end if
            end if
!
!   End of <method>_3
!

!
!   Start of <method>_<method>
!                   
            datum = datum + 1
            if (j == datum) then
              n_opt = n_opt + 1
              sum_opt = sum_opt + PM_rms_opt
              write(line(k:),'(f6.3)')PM_rms_opt
              line(len_trim(line) + 1:) = one_line(i)(k + 1:)
            end if

            datum = datum + 1
            if (j == datum) then
              if (PM_gnorm_opt > g_lim*sqrt(1.d0*natoms) .or. PM_hof_opt > -1.d0) then
                pre_font = '<FONT color="red">'
                post_font = '</FONT>'
              else
                pre_font = ' '
                post_font = ' '
              end if       
              write(line(k:),'(a,f8.1,a)')trim(pre_font),PM_hof_opt,trim(post_font)
              line(len_trim(line) + 1:) = one_line(i)(k + 1:)
            end if
            
            datum = datum + 1
            if (j == datum) then
              if (abs(PM_hof_opt) < 0.1d0) then
                line = "&nbsp; &nbsp; &nbsp; &nbsp; &nbsp; "
              else
                line(k:) = &
                '<a href="'//trim(method1)//'_'//trim(method1)//'/'//trim(file_name)//'.html" target="_blank">'//'Jmol</a> '// &
                '<a href="'//trim(method1)//'_'//trim(method1)//'/'//trim(file_name)//'.arc">'//'ARC</a> '// &
                '<a href="'//trim(method1)//'_'//trim(method1)//'/'//trim(file_name)//'.pdb">'//'PDB</a> '
                line(len_trim(line) + 1:) = one_line(i)(k + 1:)
              end if
            end if
!
!   End of <method>_<method>
!
          end if 
          write(iw,'(a)')trim(line)
        end do
      end do
!
! Write averages for RMS deviations
!
      write(iw,'(a)')'<tr><td></td></tr><tr><td></td><td></td><td></td><td></td><td></td>'
      write(iw,'(a,f6.3,a)')'<td colspan="2"><p align="right">Average RMS Error:</td><td><p align="center">', &
        sum_10/n_10, '</td>'
      write(iw,'(a,f6.3,a)')'<td></td><td></td><td><p align="center">', &
        sum_3/n_3, '</td>'
      write(iw,'(a,f6.3,a)')'<td></td><td></td><td><p align="center">', &
        sum_opt/n_opt, '</td></tr>'
!
!  Write footer for the HTML file
!
      open(unit = ir, file = trim(input_data)//"footer_for_table.txt", iostat = io_stat)
      do
        read(ir,'(a)', iostat = io_stat)line
        if (io_stat /= 0) exit
        write(iw,'(a)')trim(line)
      end do
      method1 = "PM6-D3H4"
    end do
  end program build_protein_table
  subroutine write_script(script_name, idata, file_location, file_name)
!
! Write one script for a protein
!
   use folder_names_C, only : location, relative_location
   implicit none
   integer :: idata
   character:: file_name*120, file_location*11, script_name*7, line*(300)
   logical :: exists
   write(idata,'(a)') ' ' 
   write(idata,'(a)') '<script type="text/javascript">'
   write(idata,'(a)') '$(document).ready(function() {Info'//script_name//' = {width: 570, height: 570, color: "0x000000",'
   write(idata,'(a)') 'disableInitialConsole: true, addSelectionOptions: false,'
   write(idata,'(a)') 'j2sPath: "'//trim(relative_location)//'j2s",'
   write(idata,'(a)') 'jarPath: "'//trim(relative_location)//'java",'
   write(idata,'(a)') 'use: "HTML5", script: "load '''//trim(file_location)//trim(file_name)// &
     '.pdb''; connect 0.3  3.6 (all) (all) Delete; \'
   write(idata,'(a)') 'connect; hBonds calculate; set measurementUnits ANGSTROMS; calculate structure; '// &
     'set hermiteLevel 4; set ribbonAspectRatio 12; trace; "}'
   write(idata,'(a)') '$("#'//trim(script_name)//'").html(Jmol.getAppletHtml("jmolApplet' &
     //trim(script_name)//'",Info'//trim(script_name)//'))});'
   write(idata,'(a)') '</script>'
   line = trim(location)//trim(file_location(4:))//trim(file_name)//'.pdb'
   inquire (file=trim(line), exist = exists)
   if (.not. exists) then
      inquire(unit=74, opened=exists) 
      if (.not. exists) open(unit = 74, file = trim(location)//'Protein_table_bugs.txt')
      write(74,'(a)')" File '"//trim(line)//"' does not exist"   
    end if      
   return
  end subroutine write_script

