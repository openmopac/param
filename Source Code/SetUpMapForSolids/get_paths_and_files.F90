subroutine get_paths_and_files
!
! Read in data for assigning paths and files, and
! use these to construct all paths and files that will be used.
!
  use chanel_C, only : f_len, input,  &
    input_folder, input_data_set, &
    P_Xray_arc_files, P_PM6_arc_files, P_PM7_arc_files, P_output_folder, P_notes_folder, &
    P_output_sub_folder, iw
  use dfport, only : system
!
  implicit none
  character :: line*(f_len)
  integer :: i
  logical :: exists
  i = iargc()
  call getarg(i, line)
!
! Clean up data-set file path, then identify path and input data-set file
!
  do i = 1, len_trim(line)
    if (line(i:i) == "\") line(i:i) = "/"
  end do
  do i = len_trim(line), 1, -1
    if (line(i:i) == "/") exit
  end do
  input_folder = line(:i)
  input_data_set = trim(line(i + 1:))
!
! Verify that the data-set exists
!
  inquire (file = trim(line), exist = exists)
  open(unit=23, file=trim(input_folder)//"Log-file.txt", status='UNKNOWN', form='FORMATTED')
  if (.not. exists) then
    write(*,*) "solids file '",trim(line)//"' does not exist"
    write(23,*) "solids file '",trim(line)//"' does not exist"
    call sleep(3)
    stop
  end if
!
!  Identify folders where the various types of files are found
!
  open(unit=input, file=trim(line), status='UNKNOWN', form='FORMATTED', action='READ') 
  line = "XRAY"
  call get_folder(line, input_folder, P_Xray_arc_files )
  line = "PM7"
  call get_folder(line, input_folder, P_PM7_arc_files )
  line = "PM6_D3H4"
  call get_folder(line, input_folder, P_PM6_arc_files )
   line = "OUT"
  call get_folder(line, input_folder, P_output_folder )
   line = "NOTES"
  call get_folder(line, input_folder, P_notes_folder )
!
! Make folders, if they don't already exist
!
  inquire (DIRECTORY = trim(P_output_folder), exist = exists)
  if ( .not. exists) then
    line = "mkdir """//trim(P_output_folder)//""""
    i = SYSTEM(trim(line))
    if (i /= 0) then
    write(*,'(///10x,a,///)') "   Folder """//trim(P_output_folder)//""" does not exist"
    call sleep(60)
    stop
    end if
  end if
  P_output_sub_folder = trim(P_output_folder)//"data_solids/"
  inquire (DIRECTORY = trim(P_output_sub_folder), exist = exists)
  if ( .not. exists) then
    line = "mkdir """//trim(P_output_sub_folder)//""""
    i = SYSTEM(trim(line))
    if (i /= 0) then
    write(*,'(///10x,a,///)') "   Folder """//trim(P_output_sub_folder)//""" does not exist"
    call sleep(60)
    stop
    end if
  end if
  open(unit=23, file=trim(input_folder)//"Log-file.txt", status='UNKNOWN', form='FORMATTED')
  return
end subroutine get_paths_and_files
