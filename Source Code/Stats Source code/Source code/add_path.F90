 subroutine add_path(file, txt)
    use common_texts, only : input_path, output_path
    implicit none
    character, intent (inout) :: file*(*)
    character, intent(in) :: txt*6
    integer :: i, j, n
    character :: path*241, path_to_file*241
    logical :: need_path, windows
!
!  Test to see if the file "file" already has an absolute path.
!
    windows = .true.
    need_path =  ((file(2:2) /= ":") .and. &
    (file(1:1) /= "\") .and. file(2:2) /= "\") 
    if (.not. need_path) return
    if (txt == "INPUT ") then
      path_to_file = trim(input_path)
    else
      path_to_file = trim(output_path)
    end if      
!
!  Test to see if the job has an absolute path.
!    
    need_path = (path_to_file(2:2) == ":" .or. index(path_to_file, "\") + index(path_to_file, "/") > 0)
    if (.not. need_path) return      
!
!  file "file" does not include the path and the path is present in path_to_file
!  therefore add path to the file "file"
!
!  First check if a relative path is used
!
    n = 0
    do 
      if (file(1:3) /= "..\") exit
      n = n + 1
      file = trim(file(4:))
    end do
!
!  Is the file defined relative to the current folder?
!  If so, delete the definition.
!
    if (file(1:2) == ".\") file = trim(file(3:))
!
!  Delete the data-set-name from the path, and delete the folders up to the relative folder.
!
    path = trim(path_to_file)
    i = len_trim(path)
    do j = 1, n
      do i = len_trim(path) - 1, 1, -1
        if (path(i:i) == "\") exit
      end do
      path = path(:i)
    end do
!
!  Finally, join the path to the start of the file
!
    file = path(:i)//trim(file)   
    return
  end subroutine add_path