module module_conversion

! standard functions for the conversion from intrinsic to apparent galaxy properties

   use module_constants
   use module_system

contains

function convert_stellarmass2absmag(M,M2L) result(mag)
      
   implicit none
   real*4,intent(in)    :: M     ! [Msun] stellar mass
   real*4,intent(in)    :: M2L   ! [Msun/Lsun] stellar mass-to-light ratio
   real*4               :: mag   ! absolute magnitude
   real*4,parameter     :: magSun = 4.83 ! absolute magnitude of sun
   
   mag = magSun-2.5*real(log10(M/M2L),4)

end function convert_stellarmass2absmag

function convert_luminosity2flux(L,dl) result(S)

   implicit none
   real*8,intent(in) :: L     ! [W] Luminosity
   real*4,intent(in) :: dl    ! [Mpc] luminosity distance
   real*8            :: S     ! [W/m^2] Flux

   S = L/real(dl,8)**2/ASphereMpc

end function

function convert_abs2appmag(M,dl) result(mag)

   implicit none
   real*4,intent(in) :: M     ! absolute magnitude
   real*4,intent(in) :: dl    ! [Mpc] luminosity distance
   real*4            :: mag   ! apparent magnitude

   mag = M+5*log10(dl)+25

end function convert_abs2appmag

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
         call out('ERROR: Unknown rotation.')
         stop
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