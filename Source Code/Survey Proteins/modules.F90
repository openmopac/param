  module molkst_C
    integer :: natoms=1000, numat, ndep, nvar, maxtxt, moperr, id
    character :: keywrd*2000, title*200, line*1000, errtxt*1000, current_data_set_name*360, num3*(3), &
      timestamp*24
    logical :: gui
  end module molkst_C
  
  module chanel_C
    integer, parameter :: f_len=1000
    integer :: ir = 5, iw = 6, input = 11, notes = 16
    character(len=f_len) ::  &
      PM6_input,       & !  The name of a PM6 arc file
      set,             & !
      input_folder,    & !  Path to data for the program
      input_data_set     !  Name of program dat-set file
!
!  Paths to the various folders
!
  character(len=f_len) :: &
      P_Xray_arc_files,     & !
      P_PM6_arc_files,      & !
      P_PM7_arc_files,      & !
      P_output_folder,      & !
      P_notes_folder,       & !
      P_output_sub_folder     !
  end module chanel_C
  
  module folder_names_C
    integer, parameter :: f_len=1000
    character (LEN=f_len) :: &
      input_data,            & ! Location of data-files used by the program to control events
      relative_location,     & ! Location of JSmol, relative to the HTML file
      arc_files,             & ! Location of the ARC files used in this analysis 
      location                 ! Location of folder that will hold the HTML files

      
    character (LEN=f_len) :: &
      X_ray_input,    & !
      PM7_in,         & !
      PM6_in            !
    
    character (LEN=f_len) :: line, line1, line2, working
    
  end module folder_names_C
  
  
  module symmetry_C
    integer, dimension(6000)            :: locpar, idepfn, locdep
    double precision, dimension(6000)   :: depmul
  end module symmetry_C

  module Common_arrays_C
    integer, parameter :: ilen=10000
    integer, dimension (2000)            :: labels, nat, na, nb, nc
    integer, dimension (2,6000)          :: loc
    double precision                               :: all_accuracy_PM7(ilen), all_accuracy_PM6(ilen), tvec(3,3)
    double precision, dimension (3,2000) :: coord, geo
    double precision, dimension (ilen)             :: xparam, hof_ref, hof_PM7, hof_PM6, d_ref, d_PM7, d_PM6
    integer :: numats(10000)
    logical :: elements(0:107,0:10000)
  end module Common_arrays_C
  
  module funcon_C 
      double precision, parameter :: pi = 3.14159265358979323846d0
      double precision, parameter :: twopi = 2.0d0 * pi
  end module funcon_C 
  
  module parameters_C 
   double precision, dimension(107) :: ams
  end module parameters_C
  
  module elemts_C 
!...Created by Pacific-Sierra Research 77to90  4.4G  08:50:09  03/09/06  
      character, dimension(107) :: elemnt*2
      character (len=12), dimension(107) :: atom_names
      data elemnt/ ' H', 'He', 'Li', 'Be', ' B', ' C', ' N', ' O', ' F', 'Ne', &
        'Na', 'Mg', 'Al', 'Si', ' P', ' S', 'Cl', 'Ar', ' K', 'Ca', 'Sc', 'Ti'&
        , ' V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', &
        'Se', 'Br', 'Kr', 'Rb', 'Sr', ' Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh'&
        , 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', ' I', 'Xe', 'Cs', 'Ba', &
        'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er'&
        , 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', ' W', 'Re', 'Os', 'Ir', 'Pt', 'Au', &
        'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa'&
        , ' U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Mi', 'XX', '+3', '-3', 'Cb', &
        '++', ' +', '--', ' -', 'Tv'/  
               data atom_names / &
   & "    Hydrogen", "      Helium", "     Lithium", "   Beryllium", &
   & "       Boron", "      Carbon", "    Nitrogen", "      Oxygen", &
   & "    Fluorine", "        Neon", "      Sodium", "   Magnesium", &
   & "    Aluminum", "     Silicon", "  Phosphorus", "      Sulfur", &
   & "    Chlorine", "       Argon", "   Potassium", "     Calcium", &
   & "    Scandium", "    Titanium", "    Vanadium", "    Chromium", &
   & "   Manganese", "        Iron", "      Cobalt", "      Nickel", &
   & "      Copper", "        Zinc", "     Gallium", "   Germanium", &
   & "     Arsenic", "    Selenium", "     Bromine", "     Krypton", &   
   & "    Rubidium", "   Strontium", "     Yttrium", "   Zirconium", &
   & "     Niobium", "  Molybdenum", "  Technetium", "   Ruthenium", &
   & "     Rhodium", "   Palladium", "      Silver", "     Cadmium", &
   & "      Indium", "         Tin", "    Antimony", "   Tellurium", &
   & "      Iodine", "       Xenon", "      Cesium", "      Barium", &
   & "   Lanthanum", "      Cerium", "Praseodymium", "   Neodymium", &
   & "  Promethium", "    Samarium", "    Europium", "  Gadolinium", &
   & "     Terbium", "  Dysprosium", "     Holmium", "      Erbium", &
   & "     Thulium", "   Ytterbium", "    Lutetium", &
   & "     Hafnium", "    Tantalum", "    Tungsten", &
   & "     Rhenium", "      Osmium", "     Iridium", "    Platinum", &
   & "        Gold", "     Mercury", "    Thallium", "        Lead", &
   & "     Bismuth", "    Polonium", "    Astatine", "       Radon", &
   & "    Francium", "      Radium", "    Actinium", "     Thorium", &
   & "Protactinium", "     Uranium", "   Neptunium", "   Plutonium", &
   & "   Americium", "      Curium", "   Berkelium", "     Mithril", &
   & "  Dummy atom", &
   & "  3+ Sparkle", "  3- Sparkle", " Capped bond", "  ++ Sparkle", &
   & "   + Sparkle", "  -- Sparkle", "   - Sparkle", "     Tv     "/
  end module elemts_C 

  
    
