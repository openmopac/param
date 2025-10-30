      subroutine gmetry(geo, coord) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      use molkst_C, only : natoms, id, nvar
      use common_arrays_C, only : labels, na, nb, nc, loc, &
      xparam, txtatm
      use elemts_C, only : elemnt
      use chanel_C, only : iw
      use funcon_C, only : pi
!***********************************************************************
!DECK MOPAC
!...Translated by Pacific-Sierra Research 77to90  4.4G  10:47:18  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      double precision:: geo(3,natoms) 
      double precision, intent(out) :: coord(3,natoms) 
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, i, mb, mc, ma, k
      double precision :: sum, error, ccos, cosa, xb, yb, zb, rbc, xa, ya, za, xyb&
        , xpa, xpb, costh, sinth, ypa, sinph, cosph, zqa, yza, coskh, sinkh, &
        sina, sind, cosd, xd, yd, zd, ypd, zpd, xpd, zqd, xqd, yqd, xrd 
      double precision, dimension(3, natoms) :: geovec
!-----------------------------------------------
!***********************************************************************
!
!    GMETRY  COMPUTES COORDINATES FROM BOND-ANGLES AND LENGTHS.
! *** IT IS ADAPTED FROM THE PROGRAM WRITTEN BY M.J.S. DEWAR.
!
!  ON INPUT:
!         GEO    = ARRAY OF INTERNAL COORDINATES.
!         NATOMS = NUMBER OF ATOMS, INCLUDING DUMMIES.
!         NA     = ARRAY OF ATOM LABELS FOR BOND LENGTHS.
!
!  ON OUTPUT:
!         COORD  = ARRAY OF CARTESIAN COORDINATES
!
!***********************************************************************
!                                     OPTION (A)
      error = 0.D0 
      geovec = 0.d0
      geo = geo - error*geovec(:,:natoms) 
!                                     OPTION (B)
      coord(:,1) = geo(:,1)
      if (natoms == 1) return
      if (na(2) == 1) then
        coord(1,2)   = coord(1,1) + geo(1,2)
        coord(2:3,2) = coord(2:3,1)
      else
        coord(:,2) = geo(:,2)  
      end if      
      if (natoms /= 2) then 
        if (na(3) == 0) then
          coord(:,3) = geo(:,3) 
        else
          ccos = cos(geo(2,3)) 
          if (na(3) == 1) then 
            coord(1,3) = coord(1,1) + geo(1,3)*ccos 
          else 
            coord(1,3) = coord(1,2) - geo(1,3)*ccos 
          endif 
          coord(2,3) = coord(2,2) + geo(1,3)*sin(geo(2,3)) 
          coord(3,3) = coord(3,2) 
        end if
          
        do i = 4, natoms 
          if (na(i) == 0) then
            coord(:,i) = geo(:,i)  ! Coordinate is already Cartesian
          end if
        end do  
        do i = 4, natoms 
          if (na(i) /= 0) then
            cosa = cos(geo(2,i)) 
            mb = nb(i) 
            mc = na(i) 
            xb = coord(1,mb) - coord(1,mc) 
            yb = coord(2,mb) - coord(2,mc) 
            zb = coord(3,mb) - coord(3,mc) 
            rbc = xb*xb + yb*yb + zb*zb 
            if (rbc < 1.D-16) then 
  !
  !     TWO ATOMS ARE COINCIDENT.  A FATAL ERROR.
  !
              write (iw, '(/10x, A,I4,A,I4,A, 18x, a)') 'ATOMS', mb, ' AND', mc, ' ARE COINCIDENT', 'CARTESIAN COORDINATES' 
              write(iw,'(58x, a, 2(11x,a))')"X", "Y", "Z" 
              write(iw, '(10x, a, i5, a, a, 3f12.5)')"Atom", mb, ": ", elemnt(labels(mb))//"("//trim(txtatm(mb))//")",coord(:,mb)
              write(iw, '(10x, a, i5, a, a, 3f12.5,/)')"Atom", mc, ": ", elemnt(labels(mc))//"("//trim(txtatm(mc))//")",coord(:,mc) 
              write (iw,'(a)')" Geometry at point of failure:"
              call geout(1) 
              call mopend ('TWO ATOMS ARE COINCIDENT.  A FATAL ERROR.') 
              return  
            else 
              rbc = 1.0D00/sqrt(rbc) 
            endif 
            ma = nc(i) 
            xa = coord(1,ma) - coord(1,mc) 
            ya = coord(2,ma) - coord(2,mc) 
            za = coord(3,ma) - coord(3,mc) 
  !
  !     ROTATE ABOUT THE Z-AXIS TO MAKE YB=0, AND XB POSITIVE.  IF XYB IS
  !     TOO SMALL, FIRST ROTATE THE Y-AXIS BY 90 DEGREES.
  !
            xyb = sqrt(xb*xb + yb*yb) 
            k = -1 
            if (xyb <= 0.009d0) then 
              xpa = za 
              za = -xa 
              xa = xpa 
              xpb = zb 
              zb = -xb 
              xb = xpb 
              xyb = sqrt(xb*xb + yb*yb) 
              if (xyb < 0.009d0) then
                write (iw, '(/10x, A,I4,A,I4,A, 18x, a)') 'ATOMS', ma, ' AND', mc, ' ARE COINCIDENT', 'CARTESIAN COORDINATES' 
                write(iw,'(58x, a, 2(11x,a))')"X", "Y", "Z" 
                write(iw, '(10x, a, i5, a, a, 3f12.5)')"Atom", ma, ": ", elemnt(labels(ma))//"("//trim(txtatm(ma))//")",coord(:,ma)
                write(iw, '(10x, a, i5, a, a, 3f12.5,/)')"Atom", mc, ": ", &
                  elemnt(labels(mc))//"("//trim(txtatm(mc))//")",coord(:,mc)
                write (iw,'(a)')" Geometry at point of failure:"
                call geout(1) 
                call mopend ('TWO ATOMS ARE COINCIDENT.  A FATAL ERROR.') 
                return   
              end if
              k = 1 
            endif 
  !
  !     ROTATE ABOUT THE Y-AXIS TO MAKE ZB VANISH
  !
            costh = xb/xyb 
            sinth = yb/xyb 
            xpa = xa*costh + ya*sinth 
            ypa = ya*costh - xa*sinth 
            sinph = zb*rbc 
            cosph = sqrt(abs(1.D00 - sinph*sinph)) 
            zqa = za*cosph - xpa*sinph 
  !
  !     ROTATE ABOUT THE X-AXIS TO MAKE ZA=0, AND YA POSITIVE.
  !
            yza = sqrt(ypa**2 + zqa**2) 
            if (yza >= 1.D-4) then 
              coskh = ypa/yza 
              sinkh = zqa/yza 
            else 
  !
  !   ANGLE TOO SMALL TO BE IMPORTANT
  !
              coskh = 1.D0 
              sinkh = 0.D0 
            endif 
  !
  !     COORDINATES :-   A=(???,YZA,0),   B=(RBC,0,0),  C=(0,0,0)
  !     NONE ARE NEGATIVE.
  !     THE COORDINATES OF I ARE EVALUATED IN THE NEW FRAME.
  !
            sina = sin(geo(2,i)) 
            sind = -sin(geo(3,i)) 
            cosd = cos(geo(3,i)) 
            xd = geo(1,i)*cosa 
            yd = geo(1,i)*sina*cosd 
            zd = geo(1,i)*sina*sind 
  !
  !     TRANSFORM THE COORDINATES BACK TO THE ORIGINAL SYSTEM.
  !
            ypd = yd*coskh - zd*sinkh 
            zpd = zd*coskh + yd*sinkh 
            xpd = xd*cosph - zpd*sinph 
            zqd = zpd*cosph + xd*sinph 
            xqd = xpd*costh - ypd*sinth 
            yqd = ypd*costh + xpd*sinth 
            if (k >= 1) then 
              xrd = -zqd 
              zqd = xqd 
              xqd = xrd 
            endif 
            coord(1,i) = xqd + coord(1,mc) 
            coord(2,i) = yqd + coord(2,mc) 
            coord(3,i) = zqd + coord(3,mc) 
!
!  If the bond-angle cosine is near zero or 180 degrees, then re-define the torsion or dihedral connectivity
!  unless it is a solid.  In solids, the re-definition would need to be restricted to the same unit cell.
!
            if (Abs (cosa) > 0.9998d0 .and. i > 4 .and. id == 0) then
              call bangle (coord, na(i), nb(i), nc(i), sum)
                if (sum > 3.05d0 .or. sum < 0.09d0) then
                  j = nc(i)
                  call renum (coord, na, nb, nc, i, natoms)
                  
                  if (nc(i) == 0) nc(i) = j
                  !
                  !  Correct GEO to show the new dihedral
                  !
                  if (nc(i) /= j) then
                    call dihed (coord, i, na(i), nb(i), nc(i), geo(3, i))
                    if (sum >= pi/2) sum = pi - sum
!
!  If angle is not in the domain 0 - 180 degrees, increase dihedral by 180 degrees,
!  because the dihedral is for the wrong angle.
!
                    if (Mod (Int(geo(2,i)/pi + 100.d0), 2) == 1) &
                 &  geo(3,i) = geo(3,i) + pi
!
! Update xparam, if necessary
!
                  do j = 1, nvar
                    if(loc(2,j) == 3 .and. loc(1,j) == i) then                    
                      xparam(j) = geo(3,i)
                      exit
                    end if
                  end do
                end if
              end if
            end if
          end if
        end do 
!
! *** NOW REMOVE THE TRANSLATION VECTORS, IF ANY, FROM THE ARRAY COORD
!
      endif 
      k = natoms 
      do while(labels(k) == 107) 
        k = k - 1 
        if (k == 0) then
           call mopend ('SOLIDS WITHOUT ATOMS ARE NOT ALLOWED.') 
           return
        end if
      end do 
      k = k + 1 
      j = 0 
      do i = 1, natoms 
        if (labels(i)==99) cycle  
        j = j + 1 
        coord(:,j) = coord(:,i) 
      end do 
      return  
      end subroutine gmetry 

     subroutine dihed(xyz, i, j, k, l, angle) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
   use funcon_C, only : pi
!...Translated by Pacific-Sierra Research 77to90  4.4G  10:47:09  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
   
   implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
   integer , intent(in) :: i 
   integer , intent(in) :: j 
   integer , intent(in) :: k 
   integer , intent(in) :: l 
   double precision, intent (out)  :: angle 
   double precision, intent(in) :: xyz(3,*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
   double precision :: xi1, xj1, xl1, yi1, yj1, yl1, zi1, zj1, zl1, dist, cosa, &
     ddd, yxdist, xi2, xl2, yi2, yl2, costh, sinth, cosph, sinph, yj2, yi3, &
     yl3
!-----------------------------------------------
!********************************************************************
!
!      DIHED CALCULATES THE DIHEDRAL ANGLE BETWEEN ATOMS I, J, K,
!            AND L.  THE CARTESIAN COORDINATES OF THESE ATOMS
!            ARE IN ARRAY XYZ.
!
!     DIHED IS A MODIFIED VERSION OF A SUBROUTINE OF THE SAME NAME
!           WHICH WAS WRITTEN BY DR. W. THEIL IN 1973.
!
!********************************************************************
     xi1 = xyz(1,i) - xyz(1,k) 
     xj1 = xyz(1,j) - xyz(1,k) 
     xl1 = xyz(1,l) - xyz(1,k) 
     yi1 = xyz(2,i) - xyz(2,k) 
     yj1 = xyz(2,j) - xyz(2,k) 
     yl1 = xyz(2,l) - xyz(2,k) 
     zi1 = xyz(3,i) - xyz(3,k) 
     zj1 = xyz(3,j) - xyz(3,k) 
     zl1 = xyz(3,l) - xyz(3,k) 
  !      ROTATE AROUND Z AXIS TO PUT KJ ALONG Y AXIS
   dist = sqrt(xj1*xj1 + yj1*yj1 + zj1*zj1) 
   cosa = zj1/dist 
   cosa = min(1.0D0,cosa) 
   cosa = dmax1(-1.0D0,cosa) 
   ddd = 1.0D0 - cosa**2 
   if (ddd <= 0.0D0) go to 10 
   yxdist = dist*sqrt(ddd) 
   if (yxdist > 1.0D-6) go to 20 
    10 continue 
   xi2 = xi1 
   xl2 = xl1 
   yi2 = yi1 
   yl2 = yl1 
   costh = cosa 
   sinth = 0.D0 
   go to 30 
    20 continue 
   cosph = yj1/yxdist 
   sinph = xj1/yxdist 
   xi2 = xi1*cosph - yi1*sinph 
   xl2 = xl1*cosph - yl1*sinph 
   yi2 = xi1*sinph + yi1*cosph 
   yj2 = xj1*sinph + yj1*cosph 
   yl2 = xl1*sinph + yl1*cosph 
!      ROTATE KJ AROUND THE X AXIS SO KJ LIES ALONG THE Z AXIS
   costh = cosa 
   sinth = yj2/dist 
    30 continue 
   yi3 = yi2*costh - zi1*sinth 
   yl3 = yl2*costh - zl1*sinth 
   call dang (xl2, yl3, xi2, yi3, angle) 
!     6.2831853  IS 2 * 3.1415926535 = 180 DEGREE
   if (angle < 0.) angle = pi*2.d0 + angle 
   if (angle >= 6.28318530717959D0) angle = 0.D0 
   return  
  end subroutine dihed 
        subroutine dang(a1, a2, b1, b2, rcos) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
!...Translated by Pacific-Sierra Research 77to90  4.4G  08:33:27  03/09/06  
!...Switches: -rl INDDO=2 INDIF=2 
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      double precision, intent(inout) :: a1 
      double precision, intent(inout) :: a2 
      double precision, intent(inout) :: b1 
      double precision, intent(inout) :: b2 
      double precision, intent(out) :: rcos 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      double precision :: zero, anorm, bnorm, sinth, costh 
!-----------------------------------------------
!*********************************************************************
!
!    DANG  DETERMINES THE ANGLE BETWEEN THE POINTS (A1,A2), (0,0),
!          AND (B1,B2).  THE RESULT IS PUT IN RCOS.
!
!*********************************************************************
      zero = 1.0D-6 
      if (abs(a1) >= zero .or. abs(a2) >= zero) then 
        if (abs(b1) >= zero .or. abs(b2) >= zero) then  
          anorm = 1.0D0/sqrt(a1**2 + a2**2) 
          bnorm = 1.0D0/sqrt(b1**2 + b2**2) 
          a1 = a1*anorm 
          a2 = a2*anorm 
          b1 = b1*bnorm 
          b2 = b2*bnorm 
          sinth = a1*b2 - a2*b1 
          costh = a1*b1 + a2*b2 
          costh = min(1.0D0,costh) 
          costh = dmax1(-1.0D0,costh) 
          rcos = acos(costh) 
          if (abs(rcos) >= 4.0D-5) then 
            if (sinth > 0.D0) rcos = 6.28318530717959D0 - rcos 
            rcos = -rcos 
            return  
          endif 
        endif 
      endif 
      rcos = 0.0D0 
      return  
  end subroutine dang 
        subroutine bangle(xyz, i, j, k, angle) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      use molkst_C, only : numat
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: i 
      integer , intent(in) :: j 
      integer , intent(in) :: k 
      double precision, intent(out) :: angle 
      double precision, intent(in) :: xyz(3,numat + 3) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      double precision :: d2ij, d2jk, d2ik, xy, temp
!-----------------------------------------------
!********************************************************************
!
! BANGLE CALCULATES THE ANGLE BETWEEN ATOMS I,J, AND K. THE
!        CARTESIAN COORDINATES ARE IN XYZ.
!
!********************************************************************
        d2ij = (xyz(1,i)-xyz(1,j))**2 + (xyz(2,i)-xyz(2,j))**2 + (xyz(3,i)-xyz(3,j))**2 
        d2jk = (xyz(1,j)-xyz(1,k))**2 + (xyz(2,j)-xyz(2,k))**2 + (xyz(3,j)-xyz(3,k))**2 
        d2ik = (xyz(1,i)-xyz(1,k))**2 + (xyz(2,i)-xyz(2,k))**2 + (xyz(3,i)-xyz(3,k))**2 
      xy = sqrt(d2ij*d2jk)
      if (xy < 1.d-20) then
        angle = 0.d0
        return
      end if
      temp = 0.5D0*(d2ij + d2jk - d2ik)/xy 
      temp = min(1.0D0,temp) 
      temp = dmax1(-1.0D0,temp) 
      angle = acos(temp) 
      return  
  end subroutine bangle 
    subroutine renum (coord, na, nb, nc, ii, natoms)
      implicit none
      integer, intent (in) :: ii, natoms
      integer, dimension (natoms), intent (in) :: na, nb
      integer, dimension (natoms), intent (inout) :: nc
      double precision, dimension (3, natoms), intent (in) :: coord
      integer :: i, jj, nai, nbi
      double precision :: angle, rab, rmin, theta
      intrinsic Asin
 !***********************************************************************
 !
 !  Renumber the NC of atom II.  On input, the angle NA(II)-NB(II)-NC(II)
 !  is too near to 0 or 180 degrees.  Find a new atom for NC(II), so that
 !  the angle will be acceptable (as large as possible)
 !***********************************************************************
      nai = na(ii)
      nbi = nb(ii)
 !
 !   Theta = 45 degrees
 !
      theta = 0.7853d0
      jj = 0
      rmin = 1.d10
      do
        do i = 1, ii - 1
          if (i /= nai .and. i /= nbi) then
            call bangle (coord, nai, nbi, i, angle)
            if (angle > 1.5707963d0) then
              angle = 2.d0 * Asin (1.d0) - angle
            end if
            if (angle >= theta) then
             !
             !   Angle is OK.  Now find atom of lowest distance
             !
              rab = (coord(1, nbi)-coord(1, i)) ** 2 + (coord(2, &
             & nbi)-coord(2, i)) ** 2 + (coord(3, nbi)-coord(3, i)) ** 2
              if (rab < rmin) then
                jj = i
                rmin = rab
              end if
            end if
          end if
        end do
        if (jj /= 0) then
       !
       !  Best NC is JJ; best angle is THMIN
       !
          nc(ii) = jj
          exit
        end if
    !
    !   No atom inside the allowed angle - reduce the angle
    !
        theta = theta * 0.5d0
        if (theta < 0.0174533d0) theta = 0.d0            
      end do
   end subroutine renum






