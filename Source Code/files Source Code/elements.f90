module common_elements
    implicit none
    character (len=3), dimension (2000) :: numbrs 
    character (len=2), dimension (107) :: elemnt
    character (len=2), dimension (107, 3) :: lowel
    integer :: inorg, maxele


  !.. Data Declarations ..
    data elemnt / "H ", "He", "Li", "Be", "B ", "C ", "N ", "O ", "F ", "Ne", &
   &"Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar", "K ", "Ca", "Sc", "Ti", &
   &"V ", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", &
   &"Br", "Kr", "Rb", "Sr", "Y ", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", &
   &"Ag", "Cd", "In", "Sn", "Sb", "Te", "I ", "Xe", "Cs", "Ba", "La", "Ce", &
   &"Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", &
   &"Lu", "Hf", "Ta", "W ", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", &
   &"Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U ", "Np", "Pu", &
   &"Am", "Cm", "Bk", "Cf", "XX", "Fm", "Md", "Cb", "++", " +", "--", " -", &
   &"Tv" /
!
!  JANAF sequence
!
    data lowel(1:107,1) / "Al", "As", "B ", "Ba", "Be", "Bi", "Br", "C ", "Ca", "Cd", &
    &"Cl", "Cs", "F ", "Ga", "Ge", "H ", "Hg", "I ", "In", "K ", "Li", &
    &"Mg", "N ", "Na", "O ", "P ", "Pb", "Rb", "S ", "Sb", "Se", "Si", &
    &"Sn", "Sr", "Te", "Tl", "Zn", "Ti", "Fe", "Cu", 67 * "XX"/
!
!  Cox & Pilcher sequence
!
    data lowel(1:107,2) / "H ", "C ", "O ", &
    &"N ", "S ", "F ", "Cl", "Br", "I ", "Li", "Na", "K ", "Rb", "Cs", &
    &"Be", "Mg", "Ca", "Sr", "Ba", "B " , "Al", "Ga", "In", "Tl", "Si", &
    &"Ge", "Sn", "Pb", "P ", "As", "Zn", "Cd", "Hg", "Sb", "Se", "Te", &
    &"Bi", "Ti", "Fe", "Cu", 67 * "XX" /
!
!  Periodic table sequence
!
    data lowel(1:107,3) / "H ", "He", "Li", "Be", "B ", "C ", "N ", "O ", "F ", "Ne", &
   &"Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar", "K ", "Ca", "Sc", "Ti", &
   &"V ", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", &
   &"Br", "Kr", "Rb", "Sr", "Y ", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", &
   &"Ag", "Cd", "In", "Sn", "Sb", "Te", "I ", "Xe", "Cs", "Ba", "La", "Ce", &
   &"Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", &
   &"Lu", "Hf", "Ta", "W ", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", &
   &"Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U ", "Np", "Pu", &
   &"Am", "Cm", "Bk", "Cf", "XX", "Fm", "Md", "Cb", "++", " +", "--", " -", &
   &"Tv" /
!
! Hard-wire periodic table sequence
!
    data inorg, maxele / 3, 83 /
    end module common_elements
