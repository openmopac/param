subroutine write_big_table(nmols)
  use common_arrays_C, only :  hof_ref, hof_PM7, all_accuracy_PM7, all_accuracy_PM6, &
    hof_PM6
  use sets_of_names_C, only : line, all_solids, line1
  use chanel_C, only : input_data_set, P_output_folder, input, set
  use molkst_C, only: timestamp
  implicit none
  integer, intent (in) :: nmols
  double precision :: sum, sum1, work(10000), sumg7, sumg6
  integer :: all_mols=21, i, j, k, batch_file = 50
  character :: chr*1, chr2*2
  logical :: opend
  i = len_trim(input_data_set)
  line = input_data_set(:i - 1)
  line1 = trim(P_output_folder)//trim(line)//".html"
  open(unit=all_mols, file=trim(line1))
  line1 = trim(line)//".html"
  write(all_mols,"(a)")"<!-- saved from url=(0020)http://OpenMOPAC.net -->"
  write(all_mols,*)" <HTML><BODY>"
  write(all_mols,*)"Time stamp: "//timestamp
  call write_header(all_mols)
  if (i > 999) then
    write(*,*)" First set of '++++' in data set '"//trim(set)//".txt' not found"
    write(23,*)" First set of '++++' in data set '"//trim(set)//".txt' not found"
  end if
  call fdate(line)
  write(all_mols,*)'HoF Errors: <B><FONT color="green"> < 2 kcal/mol </B>&nbsp;</FONT>'
  write(all_mols,*)'<FONT color="green"> < 5 kcal/mol     &nbsp;</FONT>'
  write(all_mols,*)'<B><FONT color="black"> < 10 kcal/mol &nbsp;</FONT></B>'
  write(all_mols,*)'<FONT color="red"> < 20 kcal/mol </B>&nbsp;</FONT>'
  write(all_mols,*)'<B><FONT color="red"> < 50 kcal/mol</B></FONT>. &nbsp;'
  write(all_mols,*)'Geometries: 0: good, 100: bad <P>'
  do i = 1,1000
    read(input,"(a300)")line
    if (index(line,"++++") /= 0) exit
    write(all_mols,"(a)")line(:len_trim(line))
  end do  
  if (i > 999) then
    write(*,*)" Second set of '++++' in data set '"//trim(set)//".txt' not found"
    write(23,*)" Second set of '++++' in data set '"//trim(set)//".txt' not found"
  end if  
  inquire (unit=batch_file, opened=opend)
  if (opend) close (unit=batch_file, status="KEEP")
!   write(*,*)" writing "//trim(line)
!   write(23,*)" writing "//trim(line)
  work(:nmols) = all_accuracy_PM7(:nmols)
  do i = 1, nmols
    do j = i + 1, nmols
      if (work(j) < work(i)) then
        sum = work(i)
        work(i) = work(j)
        work(j) = sum
      endif
    end do
  end do
  sum = work(max(1,nmols/2))
  do i = 1, nmols
    do j = i + 1, nmols
      if (work(j) < work(i)) then
        sum = work(i)
        work(i) = work(j)
        work(j) = sum
      endif
    end do
  end do
  sum = 0.d0
  k = 0
  do i = 1, nmols
    sum = sum + all_accuracy_PM7(i)
    k = k + 1
  end do
   
!
!  Allow for non-existent PM6_D3H4 results
!
  j = 0
  sum1 = 0.d0
  do i = 1, nmols
    if (all_accuracy_PM6(i) > 0.002d0) then
      sum1 = sum1 + all_accuracy_PM6(i)
      j = j + 1
    end if
  end do
  sumg7 = 100.d0 - sum/k
  sumg6 = 100.d0 - sum1/j
  sum = 0.d0
  sum1 = 0.d0
  j = 0
  k = 0
  do i = 1, nmols
    if (hof_PM7(i) > -1.d9 .and. hof_ref(i) > -1.d9) then
      j = j + 1
      sum = sum + abs(hof_PM7(i) - hof_ref(i))! abs(hof_ref(i) - hof_PM7(i))
    end if
    if (hof_PM6(i) > -1.d9 .and. hof_ref(i) > -1.d9) then
      k = k + 1
      sum1 = sum1+ abs(hof_ref(i) - hof_PM6(i))
    end if
  end do

    if (k > 0) then
      chr  = char(ichar("5") + max(0, int(log10(sum1/k))))
    else
      chr = "5"
    end if
    chr2 = char(ichar("2") + int(log10(nmols*1.0)))
    if (j /= 0) write(all_mols,'(a,f'//chr//'.2,a,i4,a,f6.2,a,i'//chr2//',a)') &
      '<PRE> Average unsigned error in PM7      HoF:',sum/j,' kcal/mol, for',j,' solids; geometries:',sumg7,' for',nmols,' solids '
    if (k /= 0) then
      write(all_mols,'(a,f'//chr//'.2,a,f6.2,a)') &
      '                           PM6_D3H4 HoF:',sum1/k,' kcal/mol;                 geometries:',sumg6,' </PRE>'
    else
      write(all_mols,'(a)') ' </PRE>'
    end if      
  write(all_mols,'(a)') '<font face="Courier New">'
  write(all_mols,'(a)')'<P ALIGN="LEFT"> <PRE STYLE="font-family:''Courier New''">' 
  write(all_mols,'(a)')'     #     GUI   Heat of formation (Kcal/formula unit) Geometry        Species'
  write(all_mols,'(a)')'                    Ref.        PM7    PM6-D3H4       PM7 PM6-D3H4 </PRE>'
!
!  Write out list of solids, and list of solids for each element
!
  do i = 1, nmols  
    write(all_mols,"(a)")trim(all_solids(i))
  end do
  write(all_mols,*)" <HTML><BODY>"
  return  
  end subroutine write_big_table
  subroutine write_header(file)
  use chanel_C, only : input, set
  implicit none
  integer, intent (in) :: file
  integer :: i
  character :: line*300
  open(unit=input, file=trim(set), status='UNKNOWN', form='FORMATTED') 
  rewind(input)
  line = " "   
  do i =  1, 1000
    read(input,"(a300)")line
    if (index(line,"++++") /= 0) exit
  end do  
  do i =  1, 1000
    read(input,"(a300)")line
    if (index(line,"++++") /= 0) exit
    write(file,"(a)")trim(line)
  end do 
  return
  end subroutine write_header
  