 subroutine get_folder(folder, input_folder, path)
  use chanel_C, only : f_len, input
  implicit none
  character*(f_len), intent (in) :: input_folder
  character*(f_len), intent (inout) :: folder
  character, intent (out) :: path*(f_len)
!
!  Local
!
  logical :: first = .true.
  integer, parameter :: maxfolders=10
  integer :: i, j, n_folders
  character :: folders(maxfolders)*80, keyword*(f_len)
  save
  if (first) then
    first = .false.
    do i = 1, maxfolders
      read(input,'(a80)') folders(i)
      if (folders(i)(1:4) == "++++") exit
      call upcase(folders(i), len_trim(folders(i)))
    end do
    n_folders = i - 1    
  end if
!
!  Read in the names of important folders 
!
  keyword = folder
  call upcase(keyword, len_trim(keyword))
  do i = 1, n_folders
    j = index(folders(i), ":")
    if (index(folders(i)(1:j), trim(keyword)) > 0) exit
  end do
  folder = trim(folders(i)(j + 1:))
  if (folder == " ") then
    write(*,*) "Folder corresponding to """//trim(keyword)//""" does not exist"
    call sleep (40)
    stop
  end if   
  do
    if (folder(1:1) /= " ") exit
    folder = trim(folder(2:))
  end do
  do i = 1, len_trim(folder)
    if (folder(i:i) == "\") folder(i:i) = "/"
  end do
  i = len_trim(folder)
  if (folder(i:i) /= "/") folder(i + 1:i + 1) = "/"
!
! Construct the path.
!
  if (folder(2:2) /= ":") then
    path = trim(input_folder)//trim(folder)
  else
    path = trim(folder)
  end if
  return
 end subroutine get_folder