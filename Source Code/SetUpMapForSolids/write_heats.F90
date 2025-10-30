subroutine write_heats(nmols, maxtxt)
  use common_arrays_C, only :  hof_ref, hof_PM7
  use sets_of_names_C, only : line, hof_PM7_all
  use chanel_C, only : input_data_set, P_output_folder
  use molkst_C, only: timestamp
  implicit none
  integer, intent (in) :: nmols, maxtxt
  double precision :: sum, hof_pm7_err(nmols)
  integer :: counter, i, j, k, loop, offset1 = -58, offset2 = -59
  character :: pad*100
  logical :: l_heats = .false.
  pad = " "
!
!  Identify and compact all entries that have a reference heat of formation
!  Remove some spaces at the end of the system name.
!  The desired length to the end of the "</a>" is one or two spaces
!
  loop = 0
  do counter = 1, nmols
    if (hof_PM7_all(counter) /= " ") then
      loop = loop + 1
      i = index(hof_PM7_all(counter), "@")
      line = hof_PM7_all(counter)(1:i - 122 + maxtxt)//trim(hof_PM7_all(counter)(i + 1:))
      hof_PM7_all(loop) = trim(line)
      hof_pm7_err(loop) = hof_PM7(counter) - hof_ref(counter)
    end if
  end do
  do i = 1, loop
    sum = 1.d9
    k = 0
    do j = 1, loop
      if (hof_pm7_err(j) < sum) then
        k = j
        sum = hof_pm7_err(j)
      end if
    end do
    if ( .not. l_heats) then
      l_heats = .true.
      j = len_trim(input_data_set)
      open(34,file=trim(P_output_folder)//input_data_set(:j - 1)//"_"//"PM7_heats.html")
      write(34,*)"Time stamp: "//timestamp
      call write_header(34)
      write(34,'(a)')"<h2>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"// &
      "Comparison of Heats of Formation in order of error in PM7 </h2><PRE>"
      write(34,'(a)')"<HTML><HEAD><TITLE> Here</TITLE> </HEAD><PRE>"
      write(34,"(a)")" Formula        Molecule"//pad(:maxtxt + offset1)//"Heat of Formation    Diff.", &
        & pad(:maxtxt + offset2)//"                            Exp.     Calc.   "
    end if
    write(34,'(a)')trim(hof_PM7_all(k))
    hof_pm7_err(k) = 1.d9
  end do
  counter = nmols
  if (l_heats) then
    write(34,'(a)')"</PRE></HTML>"
    close (34)
  end if
  return  
end subroutine write_heats
  