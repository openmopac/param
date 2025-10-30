subroutine write_elements_table(nmols)
  use common_arrays_C, only :  hof_ref, hof_PM7, elements, all_accuracy_PM7, all_accuracy_PM6, &
    hof_PM6
  use sets_of_names_C, only : line, all_solids, line1
  
  use chanel_C, only : P_output_folder, input, set, input_data_set
  use molkst_C, only: timestamp
  use elemts_C, only : elemnt, atom_names
  implicit none
  integer, intent (in) :: nmols
  double precision :: sum, sum1, work1(10000), sumg7, sumg6
  integer :: all_mols=21, i, j, k, loop, batch_file = 50, ihtml = 52
  character :: num*1, all_solids_file*100
  logical :: opend
  work1(:nmols) = all_accuracy_PM7(:nmols)
  do
    open(unit=input, file=trim(set), status='UNKNOWN', form='FORMATTED', iostat = i) 
    if (i == 0) exit
    call sleep(1)
  end do
  open(unit = ihtml, file = trim(P_output_folder)//"/Errors_per_element.html")
  write(ihtml,'(a)')"<HTML><HEAD><TITLE> Average Errors in solids per Element  </TITLE>"
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
  write(ihtml,'(a)')" Date:"//trim(timestamp)
        write(ihtml,'(a)')"<H2 ALIGN=""CENTER""> Average Errors in Heats of Formation and Geometries of Solids, per element"  
  write(ihtml,'(''<P ALIGN="CENTER"><TABLE>'')')
  write(ihtml,"(3a)")"<TD>  &nbsp;</TD> <TD colspan=""4""><p  align=""CENTER"">  Heats of Formation (Kcal/mol)</TD>"
  write(ihtml,"(3a)")" <TD colspan=""4""><p  align=""CENTER"">  Geometric Error (Arb. Units)</TD></TR>"

  write(ihtml,"(3a)")"<TR>  <TD>  <p  align=""CENTER"">  Element  </TD>"
  write(ihtml,"(3a)")"<TD>  <p  align=""CENTER"">  PM7  </TD>"
  write(ihtml,"(3a)")"<TD>  <p  align=""CENTER"">  No. in set  </TD>"
  write(ihtml,"(3a)")"<TD>  <p  align=""CENTER"">  PM6-D3H4  </TD>"
  write(ihtml,"(3a)")"<TD>  <p  align=""CENTER"">  No. in set </TD>"
  write(ihtml,"(3a)")"<TD>  <p  align=""CENTER"">  PM7  </TD>"
  write(ihtml,"(3a)")"<TD>  <p  align=""CENTER"">  PM6-D3H4  </TD></TR>"  
  open(unit=input, file=trim(set), status='UNKNOWN', form='FORMATTED') 
  i = len_trim(input_data_set)
  all_solids_file = "./"//input_data_set(:i - 1)//".html"
  do loop = 1, 83
    if (.not. elements(loop,0)) cycle
    if (nmols < 10) cycle!  Do not generate element-specific files unless a large number of solids are used 
    rewind(input)
    all_accuracy_PM7(:nmols) = work1(:nmols)
    line = elemnt(loop)
    if (line(1:1) == " ") line = line(2:)
    line1 = trim(P_output_folder)//trim(line)//".html"
    open(unit=all_mols, file=trim(line1))
    line1 = trim(line)//".html"
    write(all_mols,"(a)")"<!-- saved from url=(0020)http://OpenMOPAC.net -->"
    write(all_mols,*)" <HTML><BODY>"
    write(all_mols,*)"Time stamp: "//timestamp
    line = " "   
    do i =  1, 1000
      read(input,"(a300)")line
      if (index(line,"++++") /= 0) exit
    end do  
    do i =  1, 1000
      read(input,"(a300)")line
      if (index(line,"++++") /= 0) exit
    end do  
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
    end do  
    if (i > 999) then
      write(*,*)" Second set of '++++' in data set '"//trim(set)//".txt' not found"
      write(23,*)" Second set of '++++' in data set '"//trim(set)//".txt' not found"
    end if  
    j = 0
    do i = 1, nmols
      if (elements(loop,i)) j = j + 1
    end do
    num = char(ichar("1") + int(log10(j + 0.05)))
    line = atom_names(loop)
    do
      if (line(1:1) /= " ") exit
      line = line(2:)
    end do
    write(all_mols,'(a,i'//num//',8a)')'<h2 align="center">Comparison of Structures of ',j, &
      &  ' '//trim(line)//'-containing Solids predicted using PM7 and PM6-D3H4 with X-Ray </h2>'// &
      &  ' <p align="center"><a href="'//trim(all_solids_file)//'"> (All solids</a> - ' // &
      &  '<a href="Periodic_table_solids.html">Periodic Table</a> - ', & 
      &  '<a href="../index.html">Home</a> - <a href="Accuracy of PM7 and PM6-D3H4.html">Accuracy</a>', &
      & ' - <a href="../manual/index.html">Manual</a>)'     
    inquire (unit=batch_file, opened=opend)
    if (opend) close (unit=batch_file, status="KEEP")
    sum = 0.d0
    k = 0
    do i = 1, nmols
      if (elements(loop,i)) then
        sum = sum + all_accuracy_PM7(i)
        k = k + 1
      end if
    end do
   
!
!  Allow for non-existent PM6-D3H4 results
!
    j = 0
    sum1 = 0.d0
    do i = 1, nmols
      if (elements(loop,i) .and. all_accuracy_PM6(i) > 0.002d0) then
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
      if ( .not. elements(loop,i)) cycle
      if (hof_PM7(i) > -1.d9 .and. hof_ref(i) > -1.d9) then
        j = j + 1
        sum = sum + abs(hof_PM7(i) - hof_ref(i))! abs(hof_ref(i) - hof_PM7(i))
      end if
      if (hof_PM6(i) > -1.d9 .and. hof_ref(i) > -1.d9) then
        k = k + 1
        sum1 = sum1+ abs(hof_ref(i) - hof_PM6(i))
      end if
    end do
    if (j /= 0) write(all_mols,'(a,f9.2,a,i5,a)') &
      '<P> Average unsigned error in PM7      HoF:',sum/j,' kcal/mol, for',j,' solids<br> '
    if (k /= 0) write(all_mols,'(a,f9.2,a,i5,a,f6.2,a)')   ' Average unsigned error in PM6-D3H4 HoF:', &
      sum1/k,' kcal/mol, for',k,' solids </P>'
    write(ihtml,"(a, a, a)")"<TR>  <TD>  <p  align=""LEFT""> <a href="""//trim(line1)//""" </a>"// atom_names(loop)//"</TD>"
    if (j > 0) then
      write(ihtml,"(a, f8.2, a)")"<TD>  <p  align=""CENTER"">",  sum/j, "</TD>"
      write(ihtml,"(a, i3, a)")"<TD>  <p  align=""CENTER"">",  j, "</TD>"
    else
      write(ihtml,"(a, a, a)")"<TD><p style=""text-align:center""> ", &
            "- </P></TD><TD><p style=""text-align:center"">  - ","</P></TD>"
    end if      
    if (k > 0) then
      write(ihtml,"(a, f8.2, a)")"<TD>  <p  align=""CENTER"">",  sum1/k, "</TD>"
      write(ihtml,"(a, i3, a)")"<TD>  <p  align=""CENTER"">",  k, "</TD>"
    else
      write(ihtml,"(a, a, a)")"<TD><p style=""text-align:center""> ", &
            "- </P></TD><TD><p style=""text-align:center"">  - ","</P></TD>"
    end if      
    write(ihtml,"(a, f6.2, a)")"<TD>  <p  align=""CENTER"">",  sumg7, "  </TD>"
    write(ihtml,"(a, f6.2, a)")"<TD>  <p  align=""CENTER"">",  sumg6, "  </TD>"
    write(ihtml,"(a, i3, a)")"</TR>"
    write(all_mols,'(a)') '<font face="Courier New">'
    write(all_mols,'(a)')'<P ALIGN="LEFT"> <PRE STYLE="font-family:''Courier New''">' 
    write(all_mols,'(a)')'     #     GUI   Heat of formation (Kcal/formula unit) Geometry        Species'
    write(all_mols,'(a)')'                    Ref.        PM7     PM6-D3H4      PM7 PM6-D3H4 </PRE>'
!
!  Write out list of solids, and list of solids for each element
!
    do i = 1, nmols          
      if (elements(loop,i)) then
        if (len_trim(all_solids(i)) > 10) write(all_mols,"(a)")trim(all_solids(i)) 
      end if
    end do
  end do
    write(ihtml,'(a)')"</TABLE>" 
  write(ihtml,'(a)')"</HTML>"
  return  
end subroutine write_elements_table
  