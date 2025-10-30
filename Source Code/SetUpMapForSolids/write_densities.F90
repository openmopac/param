subroutine write_densities(nmols)
  use common_arrays_C, only : d_ref, d_PM7, d_PM6
  use sets_of_names_C, only : data_set_name, all_clean_data_set_names, line
  use chanel_C, only : input_data_set, P_output_folder
  use molkst_C, only: timestamp
  implicit none
  integer, intent (in) :: nmols
  double precision :: sum, den_error(nmols)
  integer :: iw_dens=21, counter, j, k
    j = len_trim(input_data_set)
    open(iw_dens, file=trim(P_output_folder)//input_data_set(:j - 1)//"_"//"densities.html")
    write(iw_dens,'(a)')"<HTML>Time stamp: "//timestamp
    call write_header(iw_dens)
    write(iw_dens,'(a)')"<h2>&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"// &
      "Comparison of Densities in order of error in PM7 </h2><PRE>"
    write(iw_dens,'(a)')"    Ref     PM7   PM6-D3H4   PM7%err   PM6-D3H4%err      Data-set name"  
!
!  Print all entries, sorted in order of the error in the PM7 density
!
    den_error(1:nmols) = (d_PM7(1:nmols) - d_ref(1:nmols))/d_ref(1:nmols)      
    do counter = 1, nmols
      sum = -1.d6
      do j = 1, nmols
        if (sum < den_error(j)) then
          k = j
          sum = den_error(k)
        end if
      end do
      den_error(k) = -1.d7          
      line = "<a href=""./data_solids/"//trim(all_clean_data_set_names(k))//"_jmol.html"">"//trim(data_set_name(k))//"</a>"
      write(iw_dens,'(3f8.3,2f12.3, 3x, a)') d_ref(k), d_PM7(k), d_PM6(k), &
        100.d0*(d_PM7(k) - d_ref(k))/d_ref(k), &
        100.d0*(d_PM6(k) - d_ref(k))/d_ref(k), "    "//trim(line)
    end do 
    write(iw_dens,'(a)')"</PRE></HTML>"
  return  
end subroutine write_densities
  