subroutine one_line(iw, No, HoF_COC_1, HoF_PRE_1, HoF_COC_2, HoF_PRE_2, co_crystal, pre_cursor, salt)
!
!  Write one line of the table
!
  implicit none
  integer, intent (in) :: iw, No
  double precision, intent (in) :: HoF_COC_1, HoF_PRE_1, HoF_COC_2, HoF_PRE_2
  character, intent (in) :: co_crystal*200, pre_cursor(2)*200
  logical, intent (in) :: salt
!
  integer :: i, j, k
  character :: line*200, chr*1, spacer*100, number*100, line1*100
  spacer = '<td width="70" style="border-style: none; border-width: medium" align="right">&nbsp;</td>'
  number = '<td width="85" style="border-style: none; border-width: medium" align="right">'
!
    write(iw,'(a)')"    <tr>"
    write(iw,'(a,i3, a)')'<td width="40" style="border-style: none; border-width: medium" align="right">', No, ' &nbsp; </td>'
    write(iw,'(a)')'<td width="1700" style="border-left: 2px none #008000; border-right-style: none; border-right-width: medium;'//&
      'border-top-style: none; border-top-width: medium; border-bottom-style: none; border-bottom-width: medium" height="20">' 
    call clean_line(line, co_crystal) 
!
!  Name of co-crystal
!
    write(iw,'(a)')'<a href="./data_solids/'//trim(line)//'_Jmol.html" target="_blank">'
    write(iw,'(a)')'<span style="font-size:12.0pt">'//trim(co_crystal)//'</span></a></td>'  
!
! Names of pre-cursors
!
    if (salt) then
      write(iw,'(a)')'<td width="76" style="border-style: none; border-width: medium" align="center"> <b><font color="#008000">SALT</font></b></td>'
    else
      write(iw,'(a)')'<td width="76" style="border-style: none; border-width: medium" align="center"> </td>'
    end if
    call clean_line(line, pre_cursor(1)) 
    write(iw,'(a)')'<td width="76" style="border-style: none; border-width: medium" align="center">'
    write(iw,'(a)')'<a href="./data_solids/'//trim(line)//'_Jmol.html"  target="_blank">A-Jmol</td>'
    call clean_line(line, pre_cursor(2)) 
    write(iw,'(a)')'<td width="76" style="border-style: none; border-width: medium" align="center">'
    write(iw,'(a)')'<a href="./data_solids/'//trim(line)//'_Jmol.html"  target="_blank">B-Jmol</td>'
!
!  First three numbers
!
    write(iw,'(a,f8.1, a)')trim(number), HoF_COC_1, '</td>'
    write(iw,'(a)')       trim(spacer)
    write(iw,'(a,f8.1, a)')trim(number), HoF_PRE_1, '</td>'
    write(iw,'(a)')       trim(spacer)
    write(iw,'(a,f8.1, a)')trim(number), HoF_COC_1 - HoF_PRE_1, '</td>'
    write(iw,'(a)')       trim(spacer)
!
! Empty column
!
    write(iw,'(a)')'<td width="65" style="border-style: none; border-width: medium" align="right">&nbsp;</td>'
    write(iw,'(a)')       trim(spacer)
!
!  Second three numbers
!
    write(iw,'(a,f8.1, a)')trim(number), HoF_COC_2, '</td>'
    write(iw,'(a)')       trim(spacer)
    write(iw,'(a,f8.1, a)')trim(number), HoF_PRE_2, '</td>'
    write(iw,'(a)')       trim(spacer)
    write(iw,'(a,f8.1, a)')trim(number), HoF_COC_2 - HoF_PRE_2, '</td>'
    write(iw,'(a)')       trim(spacer)
    write(iw,'(a)')'    </tr>'
    return
  end subroutine one_line
  subroutine clean_line(new_line, old_line)
    implicit none
    character, intent (in) :: old_line*200
    character, intent (out) :: new_line*200
    integer :: i, j, k
    character :: chr*1
    new_line = " "
    j = 0
    do i = 1, len_trim(old_line)
      chr = old_line(i:i)
      k = ichar(chr)
      if (old_line(i:i) == "+") then
        j = j + 4
        new_line(j - 3:j) = "plus"
      else if (k >= ichar("a") .and. k <= ichar("z") .or. &
               k >= ichar("A") .and. k <= ichar("Z") .or. &
               k >= ichar("0") .and. k <= ichar("9") .or. &
               k == ichar(".") .or.  k == ichar("-")) then
        j = j + 1
        new_line(j:j) = chr
      else
        j = j + 1
        new_line(j:j) = "_"
      end if
    end do
    return
    end subroutine clean_line
   