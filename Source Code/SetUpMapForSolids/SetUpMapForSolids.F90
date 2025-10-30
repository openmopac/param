Program SetUpMapForSolids
!
! Construct a HTML list of a set of crystals, and print the list and related files
!
  use chanel_C, only : set, input, input_folder, input_data_set, P_Xray_arc_files, &
    P_PM7_arc_files
!
!
  use common_arrays_C, only : hof_ref, hof_PM7, hof_PM6, elements
!
  use sets_of_names_C, only : hof_PM7_all, data_set_name, line
!
  use molkst_C, only: timestamp
!
!  Build a simple HTM document to list all the molecules
!
  implicit none
  character :: x_ray_input*240, PM7_in*300
  integer :: io_stat, counter, i, nmols, j, nheader, loop, maxtxt
  logical :: exists
  call fdate(timestamp)
  call get_paths_and_files  
  set = trim(input_folder)//trim(input_data_set)
  i = len_trim(input_data_set)
  do j = i, i - 5, -1
    if (input_data_set(j:j) == ".") exit
  end do
  input_data_set(j:) = "/"  
  open(unit=input, file=trim(set), status='UNKNOWN', form='FORMATTED', action='READ') 
!
!  Header for set of molecules file
!  
!  This is a "dummy read" - to make sure that the read done later on will work.
!
  do i = 1, 1000
    read(input,"(a120)", iostat=io_stat)line
    if (index(line,"++++") /= 0) exit
    if (io_stat /= 0) exit
  end do
  if (i > 999) then
    write(*,*)"The first set of '++++' in file '"//trim(set)//".txt' is missing"
    write(23,*)"The first set of '++++' in file '"//trim(set)//".txt' is missing"
    call sleep(2)
    stop
  end if
!
!  Header for each molecule
!
  do nheader = 1, 200
    read(input,"(a180)", iostat=io_stat)line
    if (index(line,"++++") /= 0) exit
    if (io_stat /= 0) exit
  end do
  if (io_stat /= 0) then
    write(*,*)" Deadly error during read between first and second set of '++++' in file '"//trim(set)//".txt'"
    write(23,*)" Deadly error during read between first and second set of '++++' in file '"//trim(set)//".txt'"
    call sleep(2)
    stop
  end if
  if (i > 199) then
    write(*,*)"The second set of '++++' in file '"//trim(set)//".txt' is missing"
    write(23,*)"The second set of '++++' in file '"//trim(set)//".txt' is missing"
    call sleep(2)
    stop
  end if
!
!  End of "Dummy read"
!
  nheader = nheader - 1
  maxtxt = 0
  counter = 0
  do i = 1, 100000
    counter = counter + 1
    data_set_name(counter) = " "
    read(input,"(a180)", iostat=io_stat)data_set_name(counter)
    if (len_trim(data_set_name(counter)) == 0) then
      counter = counter - 1
      exit
    end if
    do
      if (data_set_name(counter)(1:1) /= " ") exit
      data_set_name(counter) = data_set_name(counter)(2:)
    end do
    j = index(data_set_name(counter), ".arc")
    if (j /= 0) data_set_name(counter) = data_set_name(counter)(:j - 1)
    if (data_set_name(counter)(1:1) == "*") counter = counter - 1
    if (io_stat /= 0) exit
  end do 
  
    
  nmols = counter 
  if (nmols == 0) then
    write(*,*) " No names of solids found!"
    write(23,*) " No names of solids found!"
    call sleep(2)
    stop
  end if
  hof_ref(:nmols) = -1.d10
  hof_PM7(:nmols) = -1.d10
  hof_PM6(:nmols) = -1.d10
  j = len_trim(input_data_set)
!
!   Restrict set to the system which have both PM7 and REF
!
  counter = 0
  do loop = 1, nmols
    x_ray_input  = trim(P_Xray_arc_files)//trim(data_set_name(loop))//" (ReF).arc"
    inquire (file = trim(x_ray_input), exist = exists)
    if (exists) then
         PM7_in  = trim(P_PM7_arc_files)//trim(data_set_name(loop))//".arc"
      inquire (file = trim(PM7_in), exist = exists)
      if (exists) then
        counter = counter + 1
      else
        write(*,"(a)") "PM7 arc file:  "//trim(PM7_in)//" does not exist " 
        write(23,"(a)") "~/utilities/run_mopac.csh """//trim(data_set_name(loop))//".mop""" 
      end if
    else
      if (index(x_ray_input, "test") /= 0) cycle
      if (index(x_ray_input, "try") /= 0) cycle
      write(*,"(a)") "X-ray arc file:  "//trim(x_ray_input)//" does not exist " 
      write(23,"(a)") "X-ray arc file:  "//trim(x_ray_input)//" does not exist " 
    end if
    if (exists) data_set_name(counter) = data_set_name(loop)
  end do
  close (input)
  nmols = counter
!
! Create "Round Robin"
  do i = 1,11
    data_set_name(nmols + i) = data_set_name(i)
    data_set_name(-11 + i) = data_set_name(nmols - 11 + i)
  end do
!
  elements(:,:nmols) = .false.
  counter = 0 
  hof_PM7_all(1:nmols) = " "  
  call bigloop(nmols, maxtxt)
!
! Write out all densities, in order of there error
!
  call write_big_table(nmols)
  call write_elements_table(nmols)
  call write_densities(nmols)
  call write_heats(nmols, maxtxt)
end program SetUpMapForSolids

