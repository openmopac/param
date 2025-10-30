subroutine header(iw)
  implicit none
  integer, intent (in) :: iw
  character :: timestamp*24
  call fdate(timestamp)
  write(iw,'(a)')"  <!-- saved from url=_0020_http://OpenMOPAC.net -->"
  write(iw,'(a)')"  <HTML>"
  write(iw,'(a)')"  <style>"
  write(iw,'(a)')"  <!--"
  write(iw,'(a)')".font5"
  write(iw,'(a)')"	{color:black;"
  write(iw,'(a)')"	font-size:11.0pt;"
  write(iw,'(a)')"	font-weight:400;"
  write(iw,'(a)')"	font-style:normal;"
  write(iw,'(a)')"	text-decoration:none;"
  write(iw,'(a)')"	font-family:Symbol, serif;"
  write(iw,'(a)')"	}"
  write(iw,'(a)')".font0"
  write(iw,'(a)')"	{color:black;"
  write(iw,'(a)')"	font-size:11.0pt;"
  write(iw,'(a)')"	font-weight:400;"
  write(iw,'(a)')"	font-style:normal;"
  write(iw,'(a)')"	text-decoration:none;"
  write(iw,'(a)')"	font-family:Calibri, sans-serif;"
  write(iw,'(a)')"	}"
  write(iw,'(a)')".font6"
  write(iw,'(a)')"	{color:black;"
  write(iw,'(a)')"	font-size:10.0pt;"
  write(iw,'(a)')"	font-weight:400;"
  write(iw,'(a)')"	font-style:normal;"
  write(iw,'(a)')"	text-decoration:none;"
  write(iw,'(a)')"	font-family:Calibri, sans-serif;"
  write(iw,'(a)')"	}"
  write(iw,'(a)')"-->"
  write(iw,'(a)')"  </style>"
  write(iw,'(a)')"  <BODY>"
  write(iw,'(a)')"  Time-stamp: "//timestamp
  write(iw,'(a)')"<h1 align=""center"">Heats of Formation of Co-Crystals </h1>"
  write(iw,'(a)')"<p align=""center"">  <a href=""SURVEY_OF_SOLIDS.html"">(All solids</a> - <a href=""Periodic_table_solids.html"">Periodic Table</a>"
  write(iw,'(a)')" - <a href=""../index.html"">Home</a> - <a href=""Accuracy of PM7 and PM6-D3H4.html"">Accuracy</a> - <a href=""../manual/index.html"">Manual</a>)"
  write(iw,'(a)')"</p>"
  write(iw,'(a)')"<div align=""center"">"
  write(iw,'(a)')"  <center>"
   write(iw,'(a)')'   <table border="1" cellpadding="0" cellspacing="0" style="border-collapse: collapse; border-left-width:0; border-right-width:0; border-top-width:2; border-bottom-width:2" bordercolor="#008000" width="1113">'
    write(iw,'(a)')'    <tr>'
     write(iw,'(a)')'   <td  align="center" style="border-left-style: none; border-left-width: medium; border-right-style: none; border-right-width: medium; border-top-style: solid; border-top-width: 1; border-bottom-style: solid; border-bottom-width: 1">'
      write(iw,'(a)')'No.</td>'
    write(iw,'(a)')'      <td  align="center" style="border-left-style: none; border-left-width: medium; border-right-style: none; border-right-width: medium; border-top-style: solid; border-top-width: 1; border-bottom-style: solid; border-bottom-width: 1">'
    write(iw,'(a)')'      <p align="center">Co-Crystal</td>'
    write(iw,'(a)')'<td   style="border-left-style: none; border-left-width: medium; border-right-style: none; border-right-width: medium; border-top-style: solid; border-top-width: 1; border-bottom-style: solid; border-bottom-width: 1">'
     write(iw,'(a)')'  <p align="center">   </td>'

    write(iw,'(a)')' <td  colspan="2" style="border-left-style: none; border-left-width: medium; border-right-style: none;'
    write(iw,'(a)')'border-right-width: medium; border-top-style: solid; border-top-width: 1; border-bottom-style: solid; border-bottom-width: 1"><p align="center">Precursors</td>'

    write(iw,'(a)')'      <td colspan="2" style="border-left-style: none; border-left-width: medium; border-right-style: none; border-right-width: medium; border-top-style: solid; border-top-width: 1; border-bottom-style: solid; border-bottom-width: 1">'
   write(iw,'(a)')'       <p align="center">PM7 &Delta;H<sub>f </sub> <br>'
    write(iw,'(a)')'      of <br> Cocrystal</td>'
    write(iw,'(a)')'      <td colspan="2" style="border-left-style: none; border-left-width: medium; border-right-style: none; border-right-width: medium; border-top-style: solid; border-top-width: 1; border-bottom-style: solid; border-bottom-width: 1">'
    write(iw,'(a)')'      <p align="center">PM7 &Delta;H<sub>f </sub> <br> of <br> Precusors</td>'
    write(iw,'(a)')'      <td  colspan="2" style="border-left-style: none; border-left-width: medium; border-right-style: none; border-right-width: medium; border-top-style: solid; border-top-width: 1; border-bottom-style: solid; border-bottom-width: 1">'
    write(iw,'(a)')'      <p align="center">PM7 Energy of Cocrystal</td>'
    write(iw,'(a)')'      <td  colspan="2" style="border-left-style: none; border-left-width: medium; border-right-style: none; border-right-width: medium; border-top-style: solid; border-top-width: 1; border-bottom-style: solid; border-bottom-width: 1">'
    write(iw,'(a)')'      <p align="center">&nbsp;</td>'
     write(iw,'(a)')'      <td  colspan="2" style="border-left-style: none; border-left-width: medium; border-right-style: none; border-right-width: medium; border-top-style: solid; border-top-width: 1; border-bottom-style: solid; border-bottom-width: 1">'
    write(iw,'(a)')'     <p align="center">PM6-D3H4<br> &Delta;H<sub>f </sub> &nbsp; of <br>  Cocrystal</td>'
    write(iw,'(a)')'      <td  colspan="2" style="border-left-style: none; border-left-width: medium; border-right-style: none; border-right-width: medium; border-top-style: solid; border-top-width: 1; border-bottom-style: solid; border-bottom-width: 1">'
    write(iw,'(a)')'     <p align="center">PM6-D3H4 <br> &Delta;H<sub>f </sub> &nbsp; of <br>  Precusors</td>'
   write(iw,'(a)')'      <td   colspan="2" style="border-left-style: none; border-left-width: medium; border-right-style: none; border-right-width: medium; border-top-style: solid; border-top-width: 1; border-bottom-style: solid; border-bottom-width: 1">'
    write(iw,'(a)')'     <p align="center">PM6-D3H4 Energy of Cocrystal</td>'
   write(iw,'(a)')'    </tr>'
   return
  end subroutine header
