 double precision function accuracy(ref, calcd, Xray_a, Xray_b, Xray_c, calcd_a, calcd_b, calcd_c) 
!
!  Calculate accuracy based on interatomic distances
!  and unit-cell translation vectors
!
  use molkst_C, only : numat
  use mod_atomradii, only: atom_radius_covalent
  use common_arrays_C, only : nat
  implicit none
  double precision :: ref(3,numat), calcd(3,numat), Xray_a, Xray_b, Xray_c, calcd_a, calcd_b, calcd_c
  double precision :: sum, acc, Rab(2000,2000), R12
  integer :: i, j, k, n_ab
    acc = 0.d0
    do i = 1, numat
      do j = 1, i - 1
        R12 = 0.d0
        do k = 1,3
          R12 = R12 + (ref(k,i) - ref(k,j))**2
        end do
        R12 = 1.3d0*sqrt(R12)/(atom_radius_covalent(nat(i)) + atom_radius_covalent(nat(j)))
        Rab(i,j) = R12
      end do
    end do
!
!  Rab holds the bonded interatomic atom distances
!
    n_ab = 1
    do i = 1, numat
      do j = 1, i - 1            
        R12 = 0.d0
        do k = 1,3
          R12 = R12 + (calcd(k,i) - calcd(k,j))**2
        end do
        R12 = 1.3d0*sqrt(R12)/(atom_radius_covalent(nat(i)) + atom_radius_covalent(nat(j)))
        if (Rab(i,j) < 2.d0) then
          sum = R12/Rab(i,j)
          if (sum > 1.d0) sum = 1.d0/sum
          acc = acc + 1.d0 - sum
          n_ab = n_ab + 1
        end if
        if (R12 < 2.d0) then
          sum = R12/Rab(i,j)
          if (sum > 1.d0) sum = 1.d0/sum
          acc = acc + 1.d0 - sum
          n_ab = n_ab + 1
        end if
      end do
    end do
!
!  set accuracy to a simple function of error in bond-lengths
!
    acc = 100.d0*(1.d0 - 2.d0*acc/n_ab)
!
!  modify accuracy to include error in Tv
!
    if (calcd_a > Xray_a) then
      acc = acc*Xray_a/calcd_a
    else
      acc = acc*calcd_a/Xray_a
    end if
    if (calcd_b > Xray_b) then
      acc = acc*Xray_b/calcd_b
    else
      acc = acc*calcd_b/Xray_b
    end if
    if (calcd_c > Xray_c) then
      acc = acc*Xray_c/calcd_c
    else
      acc = acc*calcd_c/Xray_c
    end if
    accuracy = max(0.01d0, acc) - 1.d-10
    return
  end function accuracy