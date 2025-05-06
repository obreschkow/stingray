! Standard functions for the conversion from intrinsic to apparent galaxy properties
   
module module_conversion

   use shared_module_core
   use shared_module_cosmology
   use shared_module_maths
   use shared_module_vectors
   use shared_module_constants
   use module_global
   use module_parameters
   
   public

contains

function rotate_vector(x,base,ispseudovector) result(y)

   implicit none
   real*4,intent(in)          :: x(3)        ! vector to be rotated
   type(type_base),intent(in) :: base
   logical,intent(in)         :: ispseudovector
   real*4                     :: sgn
   real*4                     :: y(3)
   
   if (ispseudovector) then
      sgn = 1
   else
      sgn = 1-2*log2int(base%transformation%inverted)
   end if
   y = sgn*matmul(base%transformation%rotation,x)
   
end function rotate_vector

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
      zv1 = min(0.1,max(-0.1,sum(v*elos)/const%c*1e3)) ! limited to 0.1 to avoid relativistic regime, ok for all practical purposes
      zv2 = min(0.1,max(-0.1,-sum(para%velocity_car*elos)/const%c*1e3)) ! limited to 0.1 to avoid relativistic regime, ok for all practical purposes
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
   
   normx = norm(x)
   normJ = norm(J)
   
   if ((normx<=epsilon(normx)).or.(normJ<=epsilon(normJ))) then
      
      ! assign random inclination
      pa = get_random_uniform_number(0.0,2.0*pi,modern=para%modern_prng)
      inclination = acos(get_random_uniform_number(0.0,1.0,modern=para%modern_prng))
      
   else
   
      eLOS = unitvector(x)
      eJ = unitvector(J)
   
      inclination = acos(min(1.0,max(-1.0,eLOS.dot.eJ)))
      if (inclination>pi/2.0) inclination = pi-inclination
   
      eMajor = unitvector(eLOS.cross.eJ)
      eNorth = unitvector((/0,0,1/)-eLOS(3)*eLOS)
      eEast = eNorth.cross.eLOS
   
      pa = atan2(eMajor.dot.eEast,eMajor.dot.eNorth)
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
   real*8            :: S     ! [W/m^2] Integrated flux density
   real*8,parameter  :: ASphereMpc = 1.1964952e+46_8 ! [m^2] surface area of a sphere with 1 Mpc radius (=4*pi*(Mpc/m)^2)
   
   S = L/real(dl,8)**2/ASphereMpc

end function convert_luminosity2flux

function convert_intflux2velintflux(S,lambda,z) result(S_V)

   implicit none
   real*4,intent(in) :: S           ! [W/m^2] integrated flux density
   real*4,intent(in) :: lambda      ! [m] rest-frame wave length
   real*4,intent(in) :: z           ! [-] observed redshift
   real*4            :: S_V         ! [Jy km/s] velocity-integrated flux density
   
   S_V = 1e23*S*lambda*(1+z)

end function convert_intflux2velintflux

end module module_conversion