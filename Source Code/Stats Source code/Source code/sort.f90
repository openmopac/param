subroutine sort(ranking, no_of_types,new_ranking)
   implicit none
!
!   sort "sorts" the entries in 'ranking' and puts the sorted mapping into new_ranking
!
   integer :: no_of_types, new_ranking(no_of_types)
   double precision, dimension(no_of_types) :: ranking
!
!  Local variables
!
   integer :: i, j, k
   double precision minimum_value
   do i = 1, no_of_types
     minimum_value = 1.d10
     do j = 1, no_of_types
       if (ranking(j) < minimum_value) then
         k = j
         minimum_value = ranking(j)
       end if
     end do
     ranking(k) = 1.d20
     new_ranking(i) = k
    end do
  return
end subroutine sort