  program BEGDB_to_MOPAC_converter
!
! Convert a set of files from Hobza's database into MOPAC format.
! One file represents the geometry of the interacting molecules, the other the geometry of the separated molecules
!
    implicit none
    character :: filenames(1000)*100, energies(1000)*12, file_filenames*100, file_energies*100, file*100, &
      file_separated*100, elem(1000)*3     
    integer :: i, j, loop, nmols, ir = 5, iw1 = 6, iw2 = 7, io_stat, natoms, mol_no
    double precision :: Z_matrix(3,1000), Rab, coords(3,1000), radii(1000)
    logical :: mol_1(1000), mol_n(1000), growing
    i = iargc()
    if (i /= 2) then
      write(*,'(a)') " "
      write(*,'(a)') " Program RtoM converts reference data files from Hobza's BEGDB database into MOPAC files."
      write(*,'(a)') " Two MOPAC files are made from each reference data file: One is the un-modified data set,"
      write(*,'(a)') " the other is the same set but with the molecules well-separated"
      write(*,'(a)') " "
      write(*,'(a)') " RtoM is run in the folder that contains the files from Hobza's BEGDB database."
      write(*,'(a)') " "
      write(*,'(a)') " INPUT DATA"
      write(*,'(a)') " "
      write(*,'(a)') " RtoM uses two files.  The first file consists of the filenames for the geometries,"
      write(*,'(a)') " one filename per line, the other file contains the interaction energies from the BEGDB database,"
      write(*,'(a)') " again one per line."
      write(*,'(a)') " "
      write(*,'(a)') " IMPORTANT: Each interaction energy must match up with the filename for the corresponding geometry "
      write(*,'(a)') " "
      write(*,'(a)') " Tip: Use EXCEL to convert the table in the BEGDB database into the interaction energy file."
      write(*,'(a)') " "
      write(*,'(a)') " "
      stop
    end if
    call getarg(1, file_filenames)
    i =  len_trim(file_filenames)
    if (file_filenames(i - 3:i - 3) /= ".") file_filenames = file_filenames(:i)//".txt"
    call getarg(2, file_energies)
    i =  len_trim(file_energies)
    if (file_energies(i - 3:i - 3) /= ".") file_energies = file_energies(:i)//".txt"
    open(ir, file = trim(file_filenames))
!
!  Count the number of files to be processed
!
    do nmols = 1, 100
      read(ir,'(a)', iostat = io_stat)filenames(nmols)
      if (io_stat /= 0) exit
      do
        if (filenames(nmols)(1:1) /= " ") exit
        filenames(nmols) = filenames(nmols)(2:)
      end do
    end do
    nmols = nmols - 1
    close(ir)
    open(ir, file = trim(file_energies))
    do i = 1, 100
      read(ir,'(a)', iostat = io_stat)energies(i)
      if (io_stat /= 0) exit
      if (energies(i) == " ") exit
      do
        if (energies(i)(1:1) /= " ") exit
        energies(i) = energies(i)(2:)
      end do
    end do
    if (nmols /= i - 1) then
      continue
    end if
!
! All file-names and energies for "nmols" systems have been read in.
! Now process the files, one at a time
!
    do loop = 1, nmols
      close (ir)
      file = trim(filenames(loop))
      open(ir, file = trim(file))
      i = len_trim(file) - 3
      if (file(i:i) == ".")file(i:) = " "
      file_separated = trim(file)//" separated.mop"
      file = trim(file)//".mop"
      read(ir,*) natoms
      read(ir,*)
      open(iw1, file = trim(file))
      open(iw2, file = trim(file_separated))
      write(iw1,'(a)') " 0SCF HTML"
      write(iw2,'(a)') " 0scf HTML"
      write(iw1,'(a)') file(:i - 1)
      write(iw2,'(a)') file(:i - 1)//" separated"
      write(iw1,'(a)') " H="//trim(energies(loop))//"+"""//trim(file(:i - 1))//" separated.mop"" HR=CCSDT HWT=5"
      write(iw2,'(a)') " H=0 HR=REF"
      radii(:natoms) = 0.d0
      do i = 1, natoms
        read(ir,*)elem(i), coords(:,i)
        if (index(elem(i), "I") + index(elem(i), "Br") > 0) radii(i) = 0.5d0
        if (index(elem(i), "H")                        > 0) radii(i) = -0.5d0
        write(iw1,'(2x, a3, 3x, 3(f14.9, " 0 "))') elem(i), coords(:,i)
      end do
      close (iw1) 
!
!  Now find all atoms in the first molecule
!
      mol_1(:natoms) = .false.
      mol_no = 0
      do
        mol_n(:natoms) = .false.
        j = 0
        do i = 1, natoms
          if (.not. mol_1(i)) then
            mol_1(i) = .true.
            mol_n(i) = .true.
            j = 1
            exit
          end if
        end do
        if (j == 0) exit
!
! Start with one atom, and find all atoms in that molecule
!
        do
          growing = .false.
          do i = 1, natoms
            if (mol_1(i)) cycle
            do j = 1, natoms
              if (j == i) cycle
              if (mol_1(j)) then
                Rab = (coords(1,i) - coords(1,j))**2 + &
                      (coords(2,i) - coords(2,j))**2 + &
                      (coords(3,i) - coords(3,j))**2
                if (Rab < (2.d0 + radii(i) + radii(j))**2) then
                  mol_1(i) = .true.
                  mol_n(i) = .true.
                  growing = .true.
                  exit
                end if
              end if
            end do    
          end do
          if (.not. growing) exit
        end do
        do i = 1, natoms
          if (mol_n(i)) write(iw2,'(2x, a3, 3x, 3(f14.9, " 0 "))') elem(i), coords(:2,i), coords(3,i) + mol_no*50.d0
        end do
        mol_no = 1
      end do
      close (iw2)
    end do
    stop 
  end program BEGDB_to_MOPAC_converter