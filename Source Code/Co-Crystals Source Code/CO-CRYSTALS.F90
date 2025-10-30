Program Co_Crystals
!
! Given a list of co-crystals and their pre-cursors,
! calculate the enery of the co-crystal relative to its pre-cursors.
!
  implicit none
  integer, parameter :: max_crystals = 1000
  integer :: i, j, k, input = 11, ir = 5, iw = 6, n_p(max_crystals), io_stat,  &
    loop, n_method, No = 0
  character :: co_crystal(max_crystals)*200, pre_cursor(10, max_crystals)*200, timestamp*24, &
    line*200, chr*1, method(2)*10, suffix(2)*12
  double precision :: HoF_COC(2, max_crystals), HoF_PRE(2, max_crystals), heat(max_crystals), &
    sum, n_pre_cursor(10, max_crystals)
  logical :: exists, salt(max_crystals)
  double precision, external :: reada
  method(1) = "PM7"
  method(2) = "PM6-D3H4"
  suffix(1) = " "
  suffix(2) = " (PM6-D3H4)"
  call fdate(timestamp)
  i = iargc()
  call getarg(i, line)
  open(unit=input, file=trim(line), status='UNKNOWN', form='FORMATTED', action='READ') 
  line = "M:\PARAM\HTML Files\Heats of Formation of Co-Crystals.html"
  open(unit=iw, file=trim(line), status='UNKNOWN', form='FORMATTED', action='WRITE') 
  HoF_PRE = 0.d0
  do loop = 1, max_crystals
    read(input,'(a)', iostat = io_stat) line
    if (io_stat /= 0) exit
      salt(loop) = (index(line, "SALT") /= 0)
!
! Identify the co-crystal and its pre-cursors
!
!  First, the co-crystal
!
    i = index(line, '"') + 1
    j = index(line(i + 1:), '"') + i - 1
    co_crystal(loop) = line(i:j)
!
! Now the pre-cursors
!
    n_p(loop) = 0
    n_pre_cursor = 1
    do
      if (index(line(j + 2:), '"') == 0) exit
      n_p(loop) = n_p(loop) + 1
      i = index(line(j + 2:), '"') + j + 2
      j = index(line(i:), '"') + i - 2
      pre_cursor(n_p(loop), loop) = line(i:j)
      if (line(i - 2:i - 2) /= " ") then
        do i = i - 2, 1, -1
          if (line(i:i) == " ") exit
        end do
        n_pre_cursor(n_p(loop), loop) = reada(line, i)
      end if
    end do
    if (n_p(loop) == 0) exit
    do n_method = 1,2
!
!  Get the HoF of the co-crystal
!
      line = "./"//trim(method(n_method))//"/"//trim(co_crystal(loop))//trim(suffix(n_method))//".arc"
      inquire (file=trim(line), exist = exists)
      if (.not. exists) then
        write(iw,'(a)')"<align=""left"">File """//trim(line)//""" does not exist<br>"
        goto 99
      end if
      open(unit=ir, file=trim(line), status='UNKNOWN', form='FORMATTED', action='READ') 
      do i = 1,100
        read(ir,'(a)') line
        if (index(line, "H.o.F. per unit") /= 0) exit
      end do
      HoF_COC(n_method, loop) = reada(line,index(line,"="))
      close(ir)
!
!  Get the HoF of the pre-cursors
!
99     do i = 1, n_p(loop)
        line = "./"//trim(method(n_method))//"/"//trim(pre_cursor(i, loop))//trim(suffix(n_method))//".arc"
        inquire (file=trim(line), exist = exists)
        if (.not. exists) then
          write(iw,'(a)')"<align=""left"">File """//trim(pre_cursor(i, loop))//trim(suffix(n_method))//".arc"" does not exist<br>"
          cycle
        end if
        open(unit=ir, file=trim(line), status='UNKNOWN', form='FORMATTED', action='READ') 
        do j = 1,100
          read(ir,'(a)') line
          if (index(line, "H.o.F. per unit") /= 0) exit
        end do
        HoF_PRE(n_method, loop) = HoF_PRE(n_method, loop) + reada(line,index(line,"="))*n_pre_cursor(i, loop)
        close(ir)
      end do
    end do
!
!  Sort using method 2
!
    heat(loop) = HoF_COC(2, loop) - HoF_PRE(2, loop)
  end do
  call header(iw)
    loop = loop - 1
    do i = 1, loop
      sum = -1.d4
      do j = 1, loop
        if (heat(j) > sum) then
          k = j
          sum = heat(j)
        end if
      end do
      No = No + 1
      call one_line(iw, No, HoF_COC(1, k), HoF_PRE(1, k), HoF_COC(2, k), HoF_PRE(2, k), co_crystal(k), pre_cursor(1,k), salt(k))
      heat(k) = -1.d5
    end do
    write(iw,'(a)')"    </table>"
    write(iw,'(a)')"  </center>"
    write(iw,'(a)')"</div>"
    write(iw,'(a)')"<p align=""center"">  All heats in Kcal/mol</p>"
    write(iw,'(a)')" <HTML><BODY>"  
    stop
  end program Co_Crystals