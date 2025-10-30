  subroutine setup_mopac_arrays(n) 
  use common_arrays_C, only : geo, coord, na, nb, nc, simbol, atmass, &
     & labels, loc, xparam, nat,   txtatm
!
!
  use symmetry_C, only : depmul, locpar, idepfn, locdep
!
  implicit none
    integer :: n
!
! Create essential arrays only at this point
!
        allocate(geo(3,n), coord(3,n), na(n), nb(n), nc(n))
        allocate(simbol(n*3), atmass(n), labels(n), loc(2, 3*n))
        allocate(xparam(3*n), nat(n))
        allocate(depmul(3*n))
        allocate(txtatm(n), locpar(3*n), idepfn(3*n), locdep(3*n))
   
  end subroutine setup_mopac_arrays
