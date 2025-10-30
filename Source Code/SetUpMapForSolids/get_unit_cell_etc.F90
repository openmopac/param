  subroutine get_unit_cell_etc(channel, HoF, a, b, c, alpha, beta, gamma, volume, density)
!
! Read in unit cell quantities from an arc file, and put then in the argument list.
!
    implicit none
    integer, intent (in) :: channel
    double precision, intent (out) :: HoF, a, b, c, alpha, beta, gamma, volume, density
!
!  Local variables
!
    integer :: i, j
    integer, parameter :: nflags = 9
    character :: line*150, flag(nflags)*6
    double precision :: cell(nflags)
    double precision, external :: reada
    data flag/"H.o.F.", "VOLUME", "DENSIT", " A   =", " B   =", " C   =", "ALPHA ", " BETA ", "GAMMA "/
!
!
!
    rewind (channel)
    do i = 1, nflags
      do 
        read(channel,'(a)', iostat = j) line
        if (j /= 0 .or. index(line, flag(i)) /= 0) exit
      end do
      if (j /= 0) exit
      cell(i) = reada(line, 30)
    end do
    if (j == 0) then 
      HoF     = cell(1)
      volume  = cell(2)
      density = cell(3)
      a       = cell(4)
      b       = cell(5)
      c       = cell(6)
      alpha   = cell(7)
      beta    = cell(8)
      gamma   = cell(9)
      return
    end if
!
!  Must be old format for ARC file
!
    rewind (channel)
    do i = 1,100
      read(channel,'(a)') line
      if (index(line, "H.o.F. per unit") /= 0) exit
    end do
    HoF = reada(line,index(line,"="))
    do i = 1,100
      read(channel,'(a)') line
      if (index(line, "alpha, beta, gamma") /= 0) exit
    end do
    do i = 1, len_trim(line)
      if (line(i:i) ==  " ") exit
    end do
    a       = reada(line, 35)
    b       = reada(line, 42)
    c       = reada(line, 49)
    alpha   = reada(line, 55)
    beta    = reada(line, 62)
    gamma   = reada(line, 69)
    volume  = reada(line, 81) 
    density = reada(line, 98) 
    return    
  end subroutine get_unit_cell_etc