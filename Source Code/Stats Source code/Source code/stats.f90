Program stats
!
! Statistics generates statistical information on one to 20 methods.
! The results are put into the folder "../HTML files"
!
!
   use common_texts, only: method, type_of_ref, heading, unit, name, name1, all, &
     folder, n_files, set_of_filenames, l_set, title, output_path, input_path, output_path, &
     heading_tex, output_file
   use elements, only: must_have, element, can_have, UP_elemnt
   implicit none
   integer, parameter :: max_methods = 20
   character (len=30), dimension (max_methods) :: methods
   character (len=30), dimension (3) :: quantities, headings, units
   character (len=300) :: line, store_line, external_file, file_name, line1
   character :: delete_old_files(9)*48
   double precision, dimension (-1:107,max_methods*2):: results = 0.d0

   integer, dimension (-1:107,max_methods*2):: no_in_set = 0
   integer :: imethod, nmethods, nquantity, i, j, k, ihtml, ir = 5, &
    mquantities
   logical :: first, exists
   !
   !.. External Functions ..
   integer, external :: iargc
   data quantities /"heats","dips","ips"/
   data units /"kcal/mol","Debye","eV"/
   data headings /"Heat of Formation","Dipole Moment","Ionization Potential"/
   data delete_old_files / &
   "Stats_for_angles.html                           ", &
   "Stats_for_bond_lengths.html                     ", &
   "Stats_for_Dipole_Moments_per_Element.html       ", &
   "Stats_for_Heats_of_Formation_per_Element.html   ", &
   "Stats_for_Ionization_Potentials_per_Element.html", &
   "table_of_dips_for_all_elements.html             ", &
   "table_of_geos_for_all_elements.html             ", &
   "table_of_heats_for_all_elements.html            ", &
   "table_of_ips_for_all_elements.html              " /
   output_path = "..\..\HTML files\"
!
! Delete all results files
!
    do i = 1, 9
      inquire (file=trim(output_path)//trim(delete_old_files(i)), exist = exists) 
      if (exists) then
        open(unit=30, file=trim(output_path)//trim(delete_old_files(i)), status='UNKNOWN', position='asis', iostat=j)
        close(30, status = 'delete', iostat=j)
      end if
    end do
!
!  Read in methods to be used
!
    heading = "not used"
    i = iargc()
    if (i > 0) then
      call getarg(1, line)
      do
        i = index(line, "/")
        if (i == 0) exit
        line(i:i) = "\"
      end do
      name = trim(line)
      file_name = trim(name)
      do i = len_trim(file_name), 1, -1
        if (file_name(i:i) == "\") exit
      end do
      input_path = file_name(:i)
      file_name = trim(file_name(i + 1:))
      i = max(1,len_trim(file_name) - 3)
      if (file_name(i:i) == ".") file_name = file_name(:i - 1)
      call upcase(line, len_trim(line))
      if (index(line, ".TXT") == 0) then
        line = trim(line)//".txt"
      else
        name = name(:len_trim(name) - 4)
      end if     
      do
        i = index(name, " ")
        if (i == 0) exit
        if (name(i:) == " ") exit
        name(i:i) = "_"
      end do
      name1 = name
      inquire (file=trim(line), exist = exists) 
      if (.not. exists) then
        write(*,'(//10x,a)') "File """//trim(line)//""" not found in this folder"
        stop
      end if
      open (unit = ir, file = trim(line), iostat = i) 
      rewind (ir)
    else
      write(*,'(10x,a)') " "
      write(*,'(10x,a)') " Program ""Stats"" analyzes the results of surveys and generates HTML files"
      write(*,'(10x,a)') " that allow average errors of the various methods to be compared easily."
      write(*,'(10x,a)') " "
      write(*,'(10x,a)') " To run the program, move to a folder that contains one or more surveys, and run"
      write(*,'(10x,a)') " the command 'stats ""<type>.txt""' Replace ""<type>.txt"" with the name of a local file,"
      write(*,'(10x,a)') " e.g., '""CHNO only.txt""' that contains the type of statistics to be used."
      write(*,'(10x,a)') " The file-name will be used in the name of an HTML file that will be made."
      write(*,'(10x,a)') " "
      write(*,'(10x,a)') " "
      write(*,'(10x,a)') " Contents of local file"
      write(*,'(10x,a)') " "
      write(*,'(10x,a)') " The following keywords define and control the type of statistical analysis: "
      write(*,'(10x,a)') " "
      write(*,'(10x,a)') " OUTPUT, METHODS, HEADING, ELEMENTS, SET "
      write(*,'(10x,a)') " "
      write(*,'(10x,a)') " Each word should be at the start of a line, and is followed by options."
      write(*,'(10x,a)') " The control keywords can be in any order, but ""SET"", if present, must be the last keyword"
      write(*,'(10x,a)') " "
      write(*,'(10x,a)') " Comment lines, lines that start with ""*"", can be added anywhere."
      write(*,'(10x,a)') " "
      write(*,'(10x,a)') "  ""OUTPUT:"" Gives the name of the folder where the results are to be put"
      write(*,'(10x,a)') "  ""METHODS:"" Defines the methods to be used, e.g., PM7 PM6 AM1 PM6-D3H4"
      write(*,'(10x,a)') "  ""HEADING:"" (An optional keyword) Defines the text of the heading for the web-page"
      write(*,'(10x,a)') "  ""ELEMENTS:"" Definition of the elements to be used:"
      write(*,'(10x,a)') "         The word ""MUST=(elements)"" to indicate which elements MUST be present"
      write(*,'(10x,a)') "         The word ""AND=(elements)"" to indicate which elements can also be present"
      write(*,'(10x,a)') "         The word ""NOT=(elements)"" to indicate which other elements MUST NOT be present"
      write(*,'(10x,a)') " "      
      write(*,'(10x,a)') "  ""SET"" On its own, this must be the last keyword, and the filenames to be used in the"
      write(*,'(10x,a)') "         analysis should be listed, one per line, starting on the next line."
      write(*,'(10x,a)') "         This option is an alternative to ""ELEMENTS"""
      write(*,'(10x,a)') " "
      write(*,'(10x,a)') "  ""SET=file"" The filenames to be used are in a separate file, ""file.txt"""
      write(*,'(10x,a)') " Any filename in the list of filenames that is not present in one or more of the methods"
      write(*,'(10x,a)') " will be ignored in that method."
      write(*,'(10x,a)') " "
      write(*,'(20x,a)') " Examples of command"
      write(*,'(10x,a)') " "
      write(*,'(10x,a)') " stats ""CHNO only.txt"""
      write(*,'(10x,a)') " stats ""H-bonds in water set.txt"""
      write(*,'(10x,a)') " "
      write(*,'(10x,a)') " "
      write(*,'(20x,a)') " Examples of local files"
      write(*,'(10x,a)') " "
      write(*,'(10x,a)') " Line 1: ""OUTPUT: ../../HTML files"""
      write(*,'(10x,a)') " Line 2: ""METHODS: PM7  PM6  PM5  PM3  PM6-D3H4  PM6-DH2X  PM6-DH+ AM1 MNDO RM1"""
      write(*,'(10x,a)') " Line 3: ""ELEMENTS: MUST=(H C N O)"" - Selects compounds of H C N and O only."
      write(*,'(10x,a)') " "
      write(*,'(10x,a)') " Line 1: ""OUTPUT: ../../Sandbox"""
      write(*,'(10x,a)') " Line 2: ""METHODS: PM7  PM6  PM5  PM3  PM6-D3H4  PM6-DH2X  PM6-DH+ AM1 MNDO RM1"""
      write(*,'(10x,a)') " Line 3: ""ELEMENTS: MUST=(C) AND=(H N O F CL BR I)"" - Selects simple organic compounds "
      write(*,'(10x,a)') "   such as CH3OH and CF4, but excludes, e.g., F2 and NH3"
      write(*,'(10x,a)') "   Note that C, C3, C4, CCl, CCl2, etc. are also selected"
      write(*,'(10x,a)') " "
      write(*,'(10x,a)') " Line 1: ""OUTPUT: ../../HTML files"""
      write(*,'(10x,a)') " Line 2: ""METHODS: PM7  PM6  PM5  PM3  PM6-D3H4  PM6-DH2X  PM6-DH+ AM1 MNDO RM1"""
      write(*,'(10x,a)') " Line 3: ""HEADING: Errors in Calculated Interaction Energies between two Water molecules"""
      write(*,'(10x,a)') " Line 4: ""SET"""
      write(*,'(10x,a)') " Line 5: ""water dimer structure 0 (non-planar open cs)"""
      write(*,'(10x,a)') " Line 6: ""water dimer structure 10 (planar bifurcated c2v)"""
      write(*,'(10x,a)') " Line 7: ""water dimer structure 2 (open ci)"""
      write(*,'(10x,a)') " "
      write(*,'(10x,a)') " Line 1: ""OUTPUT: ../../HTML files"""
      write(*,'(10x,a)') " Line 2: ""METHODS: PM7  PM6  PM5  PM3  PM6-D3H4  PM6-DH2X  PM6-DH+ AM1 MNDO RM1"""
      write(*,'(10x,a)') " Line 3: ""HEADING: Errors in Calculated Interaction Energies between two Water molecules"""
      write(*,'(10x,a)') " Line 4: ""SET=S66 set for survey"""
      stop
    end if
    title = " "
    all = .true.
    name = " "
    name1 = " "
    l_set = .false.
    output_file = " "
    heading_tex = " "
    k = 0
    do 
      read(ir,'(a)', iostat = i)line
      if (i /= 0) exit
      k = k + 1
      if (line(1:1) == "*") cycle
      store_line = trim(line)
      call upcase(line, len_trim(line))
      i = index(line(:8), "OUTPUT") 
      if (i /= 0) then
        do
          j = index(line, "/")
          if (j == 0) exit
          line(j:j) = "\"
        end do
        do
          j = index(store_line, "/")
          if (j == 0) exit
          store_line(j:j) = "\"
        end do
        do i = i + 7 , 100
          if (line(i:i) /= " ") exit
        end do
        line1 = trim(store_line(i:))
        if (line1(1:1) =='"') then
          line1 = line1(2:)
          i = len_trim(line1)
          line1(i:i) = " "
        end if
        i = len_trim(line1)
        if (i == 0) then
          output_path = "../../HTML Files/"
        else
          if (line1(i-4:i-4) == ".") then
!
!  OUTPUT contains the file and path name, not just the path name,
!  so split it into output_path and output_file
!
            output_file = trim(line1)
            do i = len_trim(line1), 1, -1
              if (line1(i:i) == "\") exit
            end do
            if (i > 0) then
              line1(i + 1:) = " "
              output_file = trim(output_file(i + 1:))
            end if
          else        
            if (line1(i:i) /= "\") line1(i + 1:i + 1) = "\"
          end if
          if (index(line1, "\") /= 0) then
            output_path = trim(line1)
          else
            output_path = "../../HTML Files/"
          end if
        end if        
        inquire (directory=trim(output_path), exist = exists) 
        if (.not. exists) then
          if (output_path(:2) == "..") then
            output_path = trim(output_path(4:))
            inquire (directory=trim(output_path), exist = exists) 
          end if
          if (.not. exists) then
            if (output_path(:2) == "..") then
              output_path = trim(output_path(4:))
              inquire (directory=trim(output_path), exist = exists) 
            end if
            if (.not. exists) then
              write(*,'(10x,a)') " "
              write(*,'(10x,a)') " OUTPUT folder does not exist"
              write(*,'(10x,a)') " "
              call sleep(40)
              stop
            end if
          end if
        end if
      end if
      i = index(line(:8), "METHOD") 
      if (i > 0) then
        do i = i , 100
          if (line(i:i) == " ") exit
        end do
        line = trim(line(i:))
        do 
          i = index(line, ",")
          if (i == 0) exit
          line(i:i) = " "
        end do
        i = 1
        nmethods = 0
        line = " "//trim(line)
        do
          do j = i, len_trim(line)
            if (line(j:j) /= " ") exit
          end do
          if (i == len_trim(line)) exit
          i = index(line(j + 1:), " ") + j
          if (i - j < 3) exit
          method = line(j:i)
          line1 = trim(method)//".heats"
          call add_path(line1, "INPUT ")
          inquire (file=trim(line1), exist = exists)
          if (exists) then
            nmethods = nmethods + 1
            methods(nmethods) = method
          end if
        end do  
        continue
      end if 
      i = index(line(:10), "HEADING") 
      if (i /= 0) then
        i = index(line, " ") + 1
        heading_tex = trim(store_line(i:))
      end if
      i = index(line(:10), "TITLE") 
      if (i /= 0) then
        i = index(line, " ") + 1
        title = trim(store_line(i:))
      end if
      i = index(line(:6), "SET") 
      if (i > 0) then
        l_set = .true.
        external_file = trim(line)
        if (index(line, "=") == 0) exit
      end if
      i = index(line(:7), "ELEME") 
      if (i > 0) then
        do i = i , 100
          if (line(i:i) == " ") exit
        end do
        line = trim(line(i:))
        all = .false.
        must_have = "**"
        if (.not. all) then
          if (line(1:1) /= " ") line=" "//line(1:len_trim(line))
          UP_elemnt = element//" "
          do i = 1,83
            call upcase(UP_elemnt(i),2) 
          end do
          do 
            i = index(line, ",")
            if (i == 0) exit
            line(i:i) = " "
          end do
          store_line = line
    !
    ! List of elements that MUST be present
    !
          i = index(line, "MUST")
          if (i > 0) then
            j = index(line(i:), ")") + i - 2
            i = index(line(i:), "(") + i
            line = " "//line(i:j)
            do i = 1, 83
              if (index(line,UP_elemnt(i)) /= 0) then
                must_have(i) = element(i)
              end if
            end do
          end if
          line = store_line
    !
    ! List of elements that CAN be present
    !
          i = index(line, "AND")
          can_have = must_have
          can_have(0) = "ZZ"
          if (i > 0) then
            j = index(line(i:), ")") + i - 2
            i = index(line(i:), "(") + i
            line = " "//line(i:j)
            do i = 1, 83
              if (index(line,UP_elemnt(i)) /= 0) then
                can_have(i) = element(i)
              end if
            end do
          end if
          if (index(line, "NOT") > 0) then
    !
    ! List of elements that must not be present
    !
            can_have = element
            can_have(0) = "ZZ"
            i = index(line, "NOT")
            if (i > 0) then
              j = index(line(i:), ")") + i - 2
              i = index(line(i:), "(") + i
              line = " "//line(i:j)
              do i = 1, 83
                if (index(line,UP_elemnt(i)) /= 0) then
                  can_have(i) = "**"
                end if
              end do
            end if
          end if
        end if
        if (.not. all) then  
          j = 1
          k = 1
          do i = 1, 83
            if (must_have(i) /= "**") then
              if (must_have(i)(1:1) == " ") then
                name(j:j) = must_have(i)(2:2)
                name1(k:k) = must_have(i)(2:2)
                j = j + 1
                k = k + 2
              else
                name(j:j + 1) = must_have(i)(1:2)
                name1(k:k + 1) = must_have(i)(1:2)
                j = j + 2
                k = k + 3
              end if
            end if
          end do
        end if
        if (name == " ") then
          name = trim(file_name)
          name1 = trim(name)
          can_have(1:107) = element(1:107)
        end if
      end if
    end do
    if (k == 0) then
      write(*,'(//10x,a)')"File """//trim(line)//""" exists, but is empty"
      stop
    end if   
    rewind(ir)

    if (l_set) then
      line = trim(external_file)
      i = index(line, "=")
      if (i > 0) then
!
! Use a defined file for the list of filenames
!
        i = i + 1
        j = len_trim(line)
        if (line(i:i) == '"') i = i + 1
        if (line(j:j) == '"') j = j - 1
        if (line(j - 3:j) /= ".TXT") then
          line = line(:j)//".TXT"
          j = j + 4
        end if
 
        close (ir)
        inquire (file=line(i:j), exist = exists)
        if (.not. exists) then
          write(*,'(//10x,a,/)')"File """//line(i:j)//""" does not exist."
          stop
        end if
        open (unit = ir, file = line(i:j), iostat = i) 
        name1 = " "
        k = 0
      end if        
!
!   A set of reference data file-names has been supplied.
!
      do i = 1, k
        read(ir,'(a)', iostat = j) line
      end do
      do n_files = 1, 1000
        do
          read(ir,'(a)', iostat = i)line
          if (line(1:1) /= "*") exit
        end do
        call upcase(line, len_trim(line))       
        if (i /= 0) exit
!
!  Delete ".mop" if present
!
        i = len_trim(line)
        if (i < 4) exit
        if (line(i - 3:i) == ".MOP") line(i - 3:) = " "
!
! Remove leading spaces
!
        do i = 1, len_trim(line)
          if (line(i:i) /= " ") exit
        end do
        if (line == " ") exit
        set_of_filenames(n_files) = trim(line(i:))
      end do
      n_files = n_files - 1  
      all = .true.
      goto 99
    end if
!
!  Generate statistics for geometries
!
    do imethod = 1, nmethods
      method = methods(imethod)
      call geos(results(-1, imethod), no_in_set(-1, imethod), &
      results(-1, imethod + max_methods), no_in_set(-1, imethod + max_methods))
    end do
    ihtml = 8
    unit = "Angstroms"
    if (name == " ") then
      line = trim(output_path)//"Stats_for_bond_lengths.html"
    else
      line = trim(output_path)//"Stats_for_bond_lengths_for_sets_of_elements.html"
    end if
    inquire (file=trim(line), exist = first)             
    first = (.not. first)
    open (unit = ihtml, file = trim(line), iostat = i)
    if (i /= 0) then
      write(*,'(/10x,a)')" Cannot write to file """//trim(line)//""""
      stop
    end if
    call print_nicely_html("Bond Lengths",results(-1,1), no_in_set(-1,1), ihtml, methods, nmethods, first)
    close (ihtml)
    unit = "Degrees"
    if (name == " ") then
      line = trim(output_path)//"Stats_for_angles.html"
    else
      line = trim(output_path)//"Stats_for_angles_for_sets_of_elements.html"
    end if
    inquire (file=trim(line), exist = first)             
    first = (.not. first)
    open (unit = ihtml, file = trim(line))
    call print_nicely_html("Bond Angles",results(-1,max_methods + 1), no_in_set(-1,max_methods + 1), &
    ihtml, methods, nmethods, first)
    close (ihtml)
 99 continue 
    if (l_set) then
      mquantities = 1
    else
      mquantities = 3
    end if 
!
!    Outer loop - over different quantities
!
    do nquantity = 1, mquantities
      if (line == "solids" .and. nquantity /= 1) cycle
      type_of_ref = quantities(nquantity)
      heading = headings(nquantity)
      unit = units(nquantity)
!
!  Inner loop - over different methods
!
      do imethod = 1, nmethods
        method = methods(imethod)
      call heats_dipoles_and_ips(results(-1,imethod), no_in_set(-1,imethod))
      end do
      ihtml = 8   
      if (line == "solids") then
        if (all) then
          open (unit = ihtml, file = trim(output_path)//"Stats_for_Heats_of_Formation_for_solids_per_Element.html")
        else
          open (unit = ihtml, file = trim(output_path)//"Stats_for_Heats_of_Formation_for_solids.html")
        end if  
      else        
        if (all .and. .not. l_set) then
          if (nquantity == 1) then
            if (output_file == " ") then
              open (unit = ihtml, file = trim(output_path)//"Stats_for_Heats_of_Formation_per_Element.html")
            else
              open (unit = ihtml, file = trim(output_path)//trim(output_file))
            end if
          end if
          if (nquantity == 2) open (unit = ihtml, file = trim(output_path)//"Stats_for_Dipole_Moments_per_Element.html")
          if (nquantity == 3) open (unit = ihtml, file = trim(output_path)//"Stats_for_Ionization_Potentials_per_Element.html")
        else
          if (nquantity == 1 .and. .not. l_set) then            
            if (output_file == " ") then
              line = trim(output_path)//"Stats_for_Heats_of_Formation.html"
            else
              line = trim(output_path)//trim(output_file)
            end if
            inquire (file=trim(line), exist = first)             
            first = (.not. first)
            open (unit = ihtml, file = trim(line))
          end if
          if (nquantity == 2) then
            line = trim(output_path)//"Stats_for_Dipole_Moments.html"
            inquire (file=trim(line), exist = first)             
            first = (.not. first)
            open (unit = ihtml, file = trim(line))
          end if
          if (nquantity == 3) then
            line = trim(output_path)//"Stats_for_Ionization_Potentials.html"
            inquire (file=trim(line), exist = first)             
            first = (.not. first)
            open (unit = ihtml, file = trim(line))
          end if
        end if     
      end if
      if (.not. l_set) call print_nicely_html(heading,results(-1,1), no_in_set(-1,1), ihtml, methods, nmethods, first)
    end do
    call merge(nmethods, methods)
end program stats

    !
  subroutine print_nicely_html(name_of_property, results, no_in_set, ihtml, methods, nmethods, first)
  use elements, only: atom_names
  use common_texts, only:  unit, name, name1, all, solids, output_file, heading_tex
  implicit none
  character(len=*) :: name_of_property
  double precision, dimension (-1:107,20):: results
  integer :: no_in_set(-1:107,20)
  logical :: first
!
!  Local
!
  character (len=30) :: methods(20), starts1(20), starts2
  character :: line*300, table(300)*300
  character :: no*2, nd, np
  double precision :: const = 0.999d0, errors(20)
  integer :: nmethods, nn(20), nlines
  integer ::  z(109)
  integer :: ihtml, i, j, k, max_i, ii
  logical :: save_file
    save_file = .false.
    do i = 1,107
      z(i) = i
    end do
    z(108) = -1
    z(109) = 0
    solids = (line == "solids")
    write(no,'(i2)')nmethods
!
! Now that all the results for the geometries are available, print them "nicely"
!
    if (first) then
      write(ihtml,'(a)')"<HTML><HEAD><TITLE> Average Errors in "&
      & //trim(name_of_property)//", per Element  </TITLE>"
      write(ihtml,'(a)')"<style>", &
        "table, th, td {", &
        "border-collapse:collapse;", &
        "border: 1px solid black;", &
        "padding: 5px;", &
        "}", &
        ".auto-style1 {", &
	      "border-style: solid;", &
        "}", &
        "</style></HEAD>"
      line = " "
      call fdate(line) 
      write(ihtml,'(a)')" Date:"//trim(line)
      if (solids) then
        if (all) then
           write(ihtml,'(a)')"<H2 ALIGN=""CENTER""> Average Unsigned Errors in Predicted Heat of Formation of Solids (" &
           //trim(unit)//") (<a href=""index.html"">Back</a>)<BR>"
           write(ihtml,'(a)') '<a href="PM6\List_of_solids.html">(Individual Species)</a></H2>'
        else
           write(ihtml,'(a)')"<H2 ALIGN=""CENTER""> Average Unsigned Errors in Predicted "//trim(name_of_property)// &
           " for sets of Elements  ("//trim(unit)//") <BR>"
           write(ihtml,'(a)') '</H2>'
        end if
      else
        if (all) then
          write(ihtml,'(a)')"<H2 ALIGN=""CENTER""> Average Unsigned Errors in Predicted "//trim(name_of_property)//" ("//trim(unit)// &
          ") <BR>"
          if (name_of_property(:4) == "Bond") then
            write(ihtml,'(a)') '<a href="table_of_geos.html">(Individual Species)</a></H2>'
          else if (name_of_property(:4) == "Dipo") then
            write(ihtml,'(a)') '<a href="table_of_dips.html">(Individual Species)</a></H2>'
          else if (name_of_property(:4) == "Heat") then
            write(ihtml,'(a)') '<a href="table_of_heats.html">(Individual Species)</a></H2>'
          else if (name_of_property(:4) == "Ioni") then
            write(ihtml,'(a)') '<a href="table_of_ips.html">(Individual Species)</a></H2>'
          end if          
        else
          if (heading_tex == " ") then
            write(ihtml,'(a)')"<H2 ALIGN=""CENTER""> Average Unsigned Errors in Predicted "//trim(name_of_property)// &
          " for sets of Elements  ("//trim(unit)//") <BR>"
          else
            write(ihtml,'(a)')"<H2 ALIGN=""CENTER""> Average Unsigned Errors in Predicted "//trim(heading_tex)//" <BR>"
          end if
          write(ihtml,'(a)') '</H2>'
        end if
      end if   
      write(ihtml,'(''<P ALIGN="CENTER"><div align="center"><center><TABLE>'')')
      write(ihtml,"(3a)")"<TR>  <TD>  <p  align=""CENTER"">  Set  </TD>"
      do i = 1, nmethods
        write(ihtml,"(a)")"<TD><p  align=""CENTER""> "//Trim(methods(i))//" </TD><TD>No. in set</TD>"
      end do
    else
      open (unit=9, form="FORMATTED", status="SCRATCH")
      rewind (ihtml)
      do 
        read(ihtml,'(a)', iostat = j)line
        if (line == " ") exit
        if (j /= 0) exit
        if (index(line, "</TABLE>") /= 0) exit
        write(9,'(a)')trim(line)
      end do
      rewind (9)
      rewind (ihtml)
      do 
        read(9,'(a)', iostat = j)line
        if (line == " ") exit
        if (j /= 0) exit
        write(ihtml,'(a)')trim(line)
      end do
      close (9)
    end if
    
    max_i = 0
    starts1(1) = '<font color="#000000">'
!
!  Print results for each element
!
    if (unit == "Angstroms") then
      np = "6"
      nd = "3"
    else
      np = "7"
      nd = "2"
    end if
    do ii = 1, 109
      i = z(ii)
      if (no_in_set(i,1) > 0 .or. .not. all) then
!
! Work out average errors for element i for each method
!
        j = i
        if (.not. all) j = -1
        do k = 1, nmethods
          nn(k) = no_in_set(j,k)
          errors(k) = results(j,k)/nn(k)
        end do
!
!  Set up colors for each average.
!  Colors are relative to the first method.
!  The first method is black if it is the best, red otherwise.
!
        j = 1
        do k = 2, nmethods
          if(errors(k) > const*errors(1) .or. errors(1) < 1.d-3) then
            starts1(k) = '<font color="#FF0000">'
          else if(errors(k) > errors(1)/const) then
            starts1(k) = '<font color="#00FF00">'
          else
            starts1(k) = '<font color="#0000FF">'
            j = 2
          end if
        end do
        if (j /= 1) starts1(1) = '<font color="#FF0000">'
        max_i = max_i + 1
!
!   Convert atomic number into a text string
!
        k = i/10
        if (k == 0) then
          line = "    &nbsp; "//char(i + ichar("0"))
        else
          line = char(k + ichar("0"))//char(i - 10*k + ichar("0"))
        end if
        if (unit == "Angstroms" .or. unit == "Degrees") then
          if (i == -1) then
            write(ihtml,"(a)")"<TR>  <TD> All elements, <br>avoiding double counting</TD>"
          else if (i == 0) then
            if (unit == "Degrees") then
              write(ihtml,"(a)")"<TR>  <TD>All elements </TD>"
            else
              write(ihtml,"(a)")"<TR>  <TD>Same as above, but<br>with double counting </TD>"
            end if 
          else
            if (name1 == " ") then
              write(ihtml,"(a,a,"//no//"(a,f8.3,a), i20,a)") &
          & "<TR>  <TD><b><font color=""#000000""><p style=""text-align:left"">"//atom_names(i), "</font></b></TD>"
            else
              write(ihtml,"(a,a,"//no//"(a,f8.3,a), i20,a)") &
          & "<TR>  <TD><b><font color=""#000000""><p style=""text-align:left"">  &nbsp;  &nbsp; " &
            //trim(name1), "</font></b></TD>"              
            end if
          end if
          do k = 1, nmethods
             j = nn(k)
             if (j > 0) then
               if (j == nn(1)) then
                 starts2 = '<font color="#000000">'
               else
                 starts2 = '<font color="#FF0000">'
               end if
               write(ihtml,"(a,f"//np//"."//nd//",a,i5,a)") "<TD><B>"//trim(starts1(k))//"<p style=""text-align:center"">", errors(k), &
              "</font></b></TD><TD>"//trim(starts2)//"<p style=""text-align:center""> ",j,"</font></TD>"
               save_file = .true.
             else
               write(ihtml,"(a,8x,a,4x,a)") "<TD><p style=""text-align:center""> ", "- </TD><TD><p align=""center"">  - ","</TD>"
             end if
          end do
          write(ihtml,"(a)")"</TR>"          
        else
!
!   Print statements
!
!
!  Start of line
!
          if (i == -1) then
            write(ihtml,"(a)")"<TR>  <TD> All elements, <br>avoiding double counting</TD>"
          else if (i == 0) then
            if (unit == "Degrees") then
              write(ihtml,"(a)")"<TR>  <TD>All elements </TD>"
            else
              write(ihtml,"(a)")"<TR>  <TD>Same as above, but<br>with double counting </TD>"
            end if            
          else
            if (.not. all) then
              if (solids) then
                write(ihtml,"(a)") "<TR>  <TD><A href=""table_of_heats_for_solids_containing_" &
                //trim(name)//".html"">"//trim(name1)//"</A></TD>"
              else
                if (unit == "kcal/mol") then
                  write(ihtml,"(a)") "<TR>  <TD><A href=""table_of_heats_for_"//trim(name)//".html"">" &
                  //trim(name1)//"</A></TD>"
                else if (unit == "Debye") then
                  write(ihtml,"(a)") "<TR>  <TD><A href=""table_of_dips_for_"//trim(name)//".html"">" &
                  //trim(name1)//"</A></TD>"
                 else if (unit == "eV") then
                  write(ihtml,"(a)") "<TR>  <TD><A href=""table_of_ips_for_"//trim(name)//".html"">" &
                  //trim(name1)//"</A></TD>"
                end if
              end if              
            else
              write(ihtml,"(a)")"<TR>  <TD><b><font color=""#000000""><p style=""text-align:left"">" &
              //trim(line)//atom_names(i)//" </P></font></b></TD>"
            end if
          end if
          do k = 1, nmethods
            if (nn(k) > 0) then
              if (nn(k) == nn(1)) then
                 starts2 = '<font color="#000000">'
               else
                 starts2 = '<font color="#FF0000">'
               end if
              write(ihtml,"(a,f"//np//"."//nd//",a,i5,a)") "<TD><B>"//trim(starts1(k))// &
              "<p style=""text-align:center"">", errors(k), "</P></font></b></TD><TD>"//trim(starts2)// &
              "<p style=""text-align:center""> ",nn(k)," </P></font></TD>"
               save_file = .true.
            else
              write(ihtml,"(a,8x,a,4x,a)") "<TD><p style=""text-align:center""> ", &
              "- </P></TD><TD><p style=""text-align:center"">  - ","</P></TD>"
            end if              
          end do
!
!   End of line
!
          write(ihtml,"(a)")" </TR>"
        end if  
        if (.not. all) exit      
      end if
    end do
    if(results(0,2)/no_in_set(0,2) > results(0,1)/no_in_set(0,1)) then
      starts1(2) = '<TD><b><font color="#FF0000">'
    else
      starts1(2) = '<TD><b><font color="#0000FF">'
    end if
    if(results(0,3)/no_in_set(0,3) > results(0,1)/no_in_set(0,1)) then
      starts1(3) = '<TD><b><font color="#FF0000">'
    else
      starts1(3) = '<TD><b><font color="#0000FF">'
    end if
    if(results(0,4)/no_in_set(0,4) > results(0,1)/no_in_set(0,1)) then
      starts1(4) = '<TD><b><font color="#FF0000">'
    else
      starts1(4) = '<TD><b><font color="#0000FF">'
    end if
    write(ihtml,'(a)')"</TABLE>  </center></div>" 
    write(ihtml,'(a)')"</HTML>"
    if (.not. save_file) then
      close (ihtml, status = 'delete')  
    end if
    return
  end subroutine print_nicely_html


