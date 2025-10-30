subroutine heats_dipoles_and_ips(results, no_in_set)   
   use elements, only: must_have, can_have, atom_names, element
   use common_texts, only: type_of_ref, method, all, use_HDI_system, folder, &
     n_files, set_of_filenames, l_set
   implicit none
   double precision, dimension (-1:107) :: results
   integer, dimension (-1:107) :: no_in_set
!
!  Local
!
   character :: line*80, formla*14, symbol2*1, one_element*2, compound_type(10000)*10, &
     line1*80, one_file_name*60
   double precision, dimension(0:107) :: unsigned_error, signed_error, rms_error
   double precision, dimension(10000) ::  unsigned_type, signed_type, rms_type
   double precision :: calc
   integer :: i, j, ir, no_of_types, max_i, loop, number_of_type(10000), iw = 26
   double precision, dimension(6000) :: calc_all
   logical :: valid, ispresent(107), exists 
!
! ... Executable Statements ...
!
    results(-1) = 0.d0
    unsigned_error = 0.d0
    signed_error =  0.d0
    rms_error = 0.d0
    no_in_set = 0
    unsigned_type = 0.d0
    signed_type = 0.d0
    number_of_type = 0
    rms_type = 0.d0
    no_of_types = 0
    ir = 5
    line = trim(method)//"."//trim(type_of_ref)
    call add_path(line, "INPUT ")
    inquire (file=trim(line), exist = exists)
    if (.not. exists) return
    open (unit = ir, file = trim(line))  
! 
!  Now read all the data, one datum per line, but first
!  dummy read over the first two lines.  
!
    rewind (ir)
    read(ir,"(a)", iostat = i) line
    if (i == 0) read(ir,"(a)", iostat = i) line
    do loop = 1,30000
!
!  Read in all data on one system
!
      read(ir,"(15x,a68,f8.3,8x,a13,a60)",iostat=i)line, calc_all(loop), formla, one_file_name
      if (i /= 0) exit
      calc = calc_all(loop)
      if (formla == " " .or. formla(1:1) /= " ") exit
      formla = formla(2:)
!
!  Decompose formula into elements
!
      j = 1
      valid = .true.
      ispresent = .false.
      do 
        if(formla(j:j) == " ") exit
        symbol2 = formla(j + 1:j + 1)
        if(symbol2 >= "A" .and. symbol2 <= "Z" .or. symbol2 == " ") then
!
! Element has one symbol, e.g. H, C, N, O, etc.
!    
          one_element = " "//formla(j:j)
          j = j + 1
        else
!
! Element has two symbols, e.g. Cl, Br, Sb, etc.
!
          one_element = formla(j:j+1)
          j = j + 2       
        end if
!
!  Now identify the elements that are present
!
        do i = 1, 107
          if (element(i) == one_element) exit
        end do
        if ( i < 108) ispresent(i) = .true.
      end do
      if (l_set) then
!
!  Is the name of the system one of the names in the set set_of_filenames?
!
        line1 = trim(one_file_name) 
        call upcase(line1, len_trim(line1))
        do i = 1, n_files
          if (line1 == set_of_filenames(i)) exit
        end do
        if (i > n_files) cycle
      else if (.not. all) then
!
!  Identify the elements that MUST be present
!
        do i = 1, 107
          if (must_have(i) /= "**" .and. .not. ispresent(i)) exit
        end do
        if (i < 108) cycle
!
! Identify the elements that MUST NOT be present, indicated by "**"
!
        do i = 1, 107
          if (can_have(i) == "**" .and. ispresent(i)) exit
        end do
        if (i < 108) cycle
      end if
      do i = 1, 83
        if( ispresent(i)) then
          unsigned_error(i) = unsigned_error(i) + Abs(calc)
          signed_error(i) =  signed_error(i) + calc
          rms_error(i) = rms_error(i) + calc**2
          no_in_set(i) = no_in_set(i) + 1
      !    write(66,'(a)')trim(line)
        end if
      end do        
      results(-1) = results(-1) + Abs(calc)! For use in Unsigned Mean Error avoiding double-counting
      no_in_set(-1) = no_in_set(-1) + 1
!
!  Group compounds according to the elements they contain
!
      do i = 1, no_of_types
        if (compound_type(i) == formla) then
          unsigned_type(i) = unsigned_type(i) + Abs(calc)
          signed_type(i) = signed_type(i) + calc
          rms_type(i) = rms_type(i) + calc**2
          number_of_type(i) = number_of_type(i) + 1
          exit
        end if
      end do
      if (i > no_of_types) then
        no_of_types = no_of_types + 1
        compound_type(i) = formla
        unsigned_type(i) =  Abs(calc)
        signed_type(i) = calc
        rms_type(i) = calc**2
        number_of_type(i) = 1
      end if
    end do
!
!  Statistics for elements for each method
!
    max_i = 0
    do i = 1, 107
      j = no_in_set(i)
      if(j > 0) then
        max_i = max_i + 1
        unsigned_error(0) = unsigned_error(0) + unsigned_error(i)
        signed_error(0) = signed_error(0) + signed_error(i)
        rms_error(0) = rms_error(0) + rms_error(i)
        no_in_set(0) = no_in_set(0) + no_in_set(i)
        results(i) = unsigned_error(i)
      end if
    end do
    results(0) = unsigned_error(0)
  return
end subroutine heats_dipoles_and_ips