subroutine geos(error_in_element_bl, no_of_element_bl, &
   & error_in_element_angle, no_of_element_angle)
   use elements, only: must_have, can_have
   use common_texts, only: method, folder
   implicit none
   double precision, dimension(-1:107) :: error_in_element_bl, error_in_element_angle
   integer, dimension (-1:107) :: no_of_element_bl, no_of_element_angle
!
!  Local quantities
!
   character :: line(20000)*58, filename*200
   integer :: i, j, k, ir, no_of_types, loop, no_of_bonds, new_ranking(20000)
   integer, dimension(5,20000) ::  no_of_type, nabc, nsort   
   double precision :: sum
   double precision, dimension(20000) :: unsigned_error, signed_error, calc, ref
   logical :: okay, ref_okay, exists
   !
   ! ... Executable Statements ...
   !
    ir = 5  
    filename = trim(method)
    call add_path(filename, "INPUT ") 
    inquire (file=trim(filename)//".geos", exist = exists)
    if (.not. exists) return
    open (unit = ir, file = trim(filename)//".geos")
    rewind (ir)    
!
!  Dummy read over the first two lines.  These are always title lines.
!
    read(ir,"(a)", end=99) line(1), line(1)
    no_of_types = 0
    unsigned_error = 0.d0
    no_of_bonds = 0
    no_of_element_bl = 0
    error_in_element_bl = 0.d0
    do i = 1, 83
      if (must_have(i) /= "**") exit
    end do
    ref_okay = (i > 83)
    do_loop: do loop = 1,20000
!
!  Read in all data on one system
!
      read(ir,'(a60,f7.3,f8.3,22x,4i3)',end=99,err=99)line(loop),ref(loop),calc(loop),&
      & (nabc(i,loop),i=4,1,-1)
      okay = ref_okay
      do i = 1, 4
        if (nabc(i,loop) == 99 .or. (can_have(nabc(i,loop)) == "**")) cycle do_loop
        if (must_have(nabc(i,loop)) /= "**") okay = .true.
      end do
      if (.not. okay) cycle
      calc(loop) = calc(loop) - ref(loop)
      if (line(loop)(1:45) == " ") line(loop)(1:45) = line(loop-1)(1:45)
!
!  Tidy up Bond lengths
!
      if (nabc(2,loop) == 0 .and. nabc(3,loop) /= 0) then
        if (calc(loop) > 1.d0) cycle ! Exclude bond-lengths of decomposed systems.
        if (nabc(3,loop) < nabc(4,loop)) then
          i = nabc(4,loop)
          nabc(4,loop) = nabc(3,loop)
          nabc(3,loop) = i
        end if
! For use in Unsigned Mean Error avoiding double-counting
        error_in_element_bl(-1) = error_in_element_bl(-1) + abs(calc(loop))
        no_of_element_bl(-1) = no_of_element_bl(-1) + 1
        i = i + 1
      else if (nabc(1,loop) == 0 .and. nabc(2,loop) /= 0) then
!
!  Tidy up Bond angles
!
        if (nabc(2,loop) < nabc(4,loop)) then
          i = nabc(4,loop)
          nabc(4,loop) = nabc(2,loop)
          nabc(2,loop) = i
        end if
! For use in Unsigned Mean Error avoiding double-counting
        error_in_element_angle(-1) = error_in_element_angle(-1) + abs(calc(loop))
        no_of_element_angle(-1) = no_of_element_angle(-1) + 1
      else if (nabc(1,loop) /= 0) then
!
!  Tidy up dihedrals
!
        if (nabc(1,loop) < nabc(4,loop)) then
          i = nabc(4,loop)
          nabc(4,loop) = nabc(1,loop)
          nabc(1,loop) = i
        end if
      end if
!
!    Sort all geometric variables into classes
!
      do i = 1, no_of_types
        do j = 1,4
          if (no_of_type(j,i) /= nabc(j,loop)) exit
        end do
        if (j > 4) then
          no_of_type(5,i) = no_of_type(5,i) + 1
          unsigned_error(i) =unsigned_error(i) + Abs(calc(loop))
          signed_error(i) = signed_error(i) + calc(loop)
          exit
        end if
      end do
      if (i > no_of_types) then
         no_of_types = no_of_types + 1
         do j = 1,4
           no_of_type(j,i) = nabc(j,loop)
         end do
           unsigned_error(i) =unsigned_error(i) + Abs(calc(loop))
           signed_error(i) = signed_error(i) + calc(loop)
           no_of_type(5,i) = 1
      end if
    end do do_loop
 99 continue
    loop = loop - 1
    if (.true.) then
!
!  Contents of "no_of_type"
!  (1,i) - torsion atom (if a torsion, otherwise zero)
!  (2,i) - angle atom (if a torsion or an angle, otherwise zero)
!  (3,i) - bond-length (if a bond-length, angle, or torsion, otherwise zero)
!  (4,i) - the atom itself (unconditional)
!  (5,i) - number of occurrences of that type
!
! Analysis of bond-lengths
!
    
    do i = 1, no_of_types
      if (no_of_type(2,i) > 0) cycle ! exclude angles and dihedrals
      if (no_of_type(3,i) == 0) cycle ! exclude Cartesian coordinates
      no_of_bonds = no_of_bonds + 1
        nsort(1,no_of_bonds) = no_of_type(3,i)
        nsort(2,no_of_bonds) = no_of_type(4,i)
      j = nsort(1,no_of_bonds) 
      no_of_element_bl(j) = no_of_element_bl(j) + no_of_type(5,i)
      error_in_element_bl(j) = error_in_element_bl(j) + unsigned_error(i)
      if (j == nsort(2,no_of_bonds)) cycle
      j = nsort(2,no_of_bonds) 
      no_of_element_bl(j) = no_of_element_bl(j) + no_of_type(5,i)
      error_in_element_bl(j) = error_in_element_bl(j) + unsigned_error(i)
    end do
!
    sum = 0.d0
    k = 0
    do i = 1,92
      if (no_of_element_bl(i) > 0) then
                k = k + no_of_element_bl(i)
        sum = sum + error_in_element_bl(i)
      end if
    end do
!
      
    error_in_element_bl(0) = sum
    no_of_element_bl(0) = k
!
    ref = 0.d0
    do i = 1, loop
      if (nabc(2,i) == 0 .and. nabc(3,i) /= 0) ref(i) = Abs(calc(i))
    end do
    call sort(ref, loop, new_ranking)
!
! Analysis of bond-angles
!
    no_of_element_angle(0:) = 0
    error_in_element_angle(0:) = 0.d0
    do i = 1, no_of_types
      if (no_of_type(1,i) > 0) cycle ! exclude dihedrals
      if (no_of_type(2,i) == 0) cycle ! exclude bond lengths
      j = no_of_type(3,i)
      no_of_element_angle(j) = no_of_element_angle(j) + no_of_type(5,i)
      error_in_element_angle(j) = error_in_element_angle(j) + unsigned_error(i)
    end do
!
    sum = 0.d0
    k = 0
    do i = 1,92
      if (no_of_element_angle(i) > 0) then        
        k = k + no_of_element_angle(i)
        sum = sum + error_in_element_angle(i)
      end if
    end do
!
    error_in_element_angle(0) = sum
    no_of_element_angle(0) = k
!
!
! Analysis of dihedral angles
!
end if
  return
end subroutine geos