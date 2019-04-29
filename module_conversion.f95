! Standard functions for the conversion from intrinsic to apparent galaxy properties
   
module module_conversion

   use module_constants
   use module_system
   use module_linalg
   use module_cosmology
   
   public
   
   real*4   :: Rvector(3,3), Rpseudo(3,3) ! argument used by function rotate, defined as global variable to
                                          ! avoid passing them through the use module

contains

function rotate(rotationmatrix,x) result(y)

   implicit none
   real*4,intent(in)    :: x(3)
   real*4,intent(in)    :: rotationmatrix(3,3)
   real*4               :: y(3)
   
   y = matmul(rotationmatrix,x)
   
end function rotate

subroutine make_redshift(x,v,zobs,zcmb,zcos)

   implicit none
   real*4,intent(in)             :: x(3)     ! [Mpc] position vector
   real*4,intent(in)             :: v(3)     ! [km/s] velocity of galaxy
   real*4,intent(out),optional   :: zobs     ! redshift in observer-frame
   real*4,intent(out),optional   :: zcmb     ! redshift in CMB frame
   real*4,intent(out),optional   :: zcos     ! cosmological redshift without peculiar motions
   real*4                        :: zv1      ! redshift due to the peculiar motion of object relative to Hubble flow
   real*4                        :: zv2      ! redshift due to the peculiar motion of the observer relative to the Hubble flow
   real*4                        :: elos(3)  ! unit vector pointing along the line of sight
   real*4                        :: dc       ! [box side-length] comoving distance
   real*4                        :: zobs_,zcmb_,zcos_
   
   dc = norm(x)
   if (dc<=epsilon(dc)) then
      zobs_ = 0
      zcmb_ = 0
      zcos_ = 0
   else   
      elos = x/dc ! unit-vector along the line of slight
      zcos_ = dc_to_redshift(dc)
      zv1 = min(0.1,max(-0.1,sum(v*elos)/c*1e3)) ! limited to 0.1 to avoid relativistic regime, ok for all practical purposes
      zv2 = min(0.1,max(-0.1,-sum(para%velocity_car*elos)/c*1e3)) ! limited to 0.1 to avoid relativistic regime, ok for all practical purposes
      zcmb_ = (1+zcos_)*(1+zv1)-1 ! following Davis & Scrimgeour 2014
      zobs_ = (1+zcos_)*(1+zv1)*(1+zv2)-1 ! following Davis & Scrimgeour 2014
   end if
   
   if (present(zobs)) zobs = zobs_
   if (present(zcmb)) zcmb = zcmb_
   if (present(zcos)) zcos = zcos_
   
end subroutine make_redshift

subroutine make_inclination_and_pa(x,J,inclination,pa)

   implicit none
   real*4,intent(in)                :: x(3)           ! [arbitrary units] position in cone
   real*4,intent(in)                :: J(3)           ! [arbitrary units] angular momentum in cone
   real*4,intent(out),optional      :: inclination    ! [rad] inclination (angle between LOS and galaxy axis)
   real*4,intent(out),optional      :: pa             ! [rad] position angle from North to East
   real*4                           :: eLOS(3)        ! unit vector along x = line-of-sight (LOS)
   real*4                           :: eJ(3)          ! unit vector pointing along galaxy axis
   real*4                           :: eMajor(3)      ! unit vector pointing along the major axis (orthogonal to LOS)
   real*4                           :: eNorth(3)      ! unit vector pointing north (or south) (orthoconal to LOS)
   real*4                           :: eEast(3)       ! unit vector pointing east (or west) (orthoconal to LOS)
   real*4                           :: normx,normJ
   real*4                           :: rand(2)
   
   normx = norm(x)
   if (normx<=epsilon(normx)) call error('make_inclination_and_pa: norm of x is zero.')
   normJ = norm(J)
   
   if (normJ<=epsilon(normJ)) then
      
      ! assign random inclination
      call random_number(rand)
      pa = rand(1)*2*pi
      inclination = acos(rand(2))
      
   else
   
      eLOS = x/normx
      eJ = J/normJ
   
      inclination = acos(sum(eLOS*eJ))
      if (inclination>pi/2.0) inclination = pi-inclination
   
      eMajor = cross_product(eLOS,eJ)
      eMajor = eMajor/norm(eMajor)
   
      eNorth = (/0,0,1/)-eLOS(3)*eLOS
      eNorth = eNorth/norm(eNorth)
      eEast = cross_product(eNorth,eLOS)
   
      pa = atan2(sum(eMajor*eEast),sum(eMajor*eNorth))
      if (pa<0) pa = 2*pi+pa
      
   end if
   
end subroutine make_inclination_and_pa

function convert_stellarmass2absmag(M,M2L,magSun) result(mag)
      
   implicit none
   real*4,intent(in)          :: M              ! [Msun] stellar mass
   real*4,intent(in)          :: M2L            ! [Msun/Lsun] stellar mass-to-light ratio
   real*4,intent(in),optional :: magSun         ! absolute magnitude of the sun
   real*4                     :: luminosity     ! [Lsun] luminosity
   real*4                     :: mag            ! absolute magnitude
   real*4,parameter           :: magSunV = 4.83 ! absolute visual magnitude of sun
   
   if (M==0) then
      mag = 99
   else
      luminosity = M/M2L ! [Lsun]
      if (present(magSun)) then
         mag = magSun-2.5*log10(luminosity)
      else
         mag = magSunV-2.5*log10(luminosity)
      end if
   end if

end function convert_stellarmass2absmag

function convert_absmag2appmag(absmag,dl) result(appmag)

   implicit none
   real*4,intent(in) :: absmag   ! absolute magnitude
   real*4,intent(in) :: dl       ! [Mpc] luminosity distance
   real*4            :: appmag   ! apparent magnitude

   if (absmag>=99) then
      appmag = 99
   else
      appmag = absmag+5*log10(dl)+25
   end if

end function convert_absmag2appmag

function convert_luminosity2flux(L,dl) result(S)

   implicit none
   real*8,intent(in) :: L     ! [W] Luminosity
   real*4,intent(in) :: dl    ! [Mpc] luminosity distance
   real*8            :: S     ! [W/m^2] Flux

   S = L/real(dl,8)**2/ASphereMpc

end function convert_luminosity2flux

function convert_vector(x,rotation) result(y)

   implicit none
   real*4,intent(in)    :: x(3)
   integer*4,intent(in) :: rotation
   real*4               :: y(3)

   select case (abs(rotation))
      case (1) ! identity
         y = x
      case(2) ! invert x-axis, while permuting y and z
         y = (/-x(1),x(3),x(2)/)
      case(3) ! invert y-axis, while permuting z and x
         y = (/x(3),-x(2),x(1)/)
      case(4) ! invert z-axis, while permuting x and y
         y = (/x(2),x(1),-x(3)/)
      case(5) ! permute (x,y,z) -> (y,z,x)
         y = (/x(2),x(3),x(1)/)
      case(6) ! permute (x,y,z) -> (z,x,y)
         y = (/x(3),x(1),x(2)/)
      case default
         call error('Unknown rotation')
   end select

   if (rotation<0) y = -y ! inversion

end function convert_vector

function convert_pseudovector(x,rotation) result(y)

   implicit none
   real*4,intent(in)    :: x(3)
   integer*4,intent(in) :: rotation
   real*4               :: y(3)

   y = convert_vector(x,abs(rotation))

end function convert_pseudovector

end module module_conversion