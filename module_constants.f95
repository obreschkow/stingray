! This module initializes and defines global constants shared by different modules

module module_constants

   public
   
   character(*),parameter  :: version = '0.11'
   
   ! Universal constants
   real*4,parameter        :: pi = 3.14159265359
   real*4,parameter        :: G = 6.67408e-11 ! [m^3/kg/s^2] gravitational constant
   real*4,parameter        :: c = 299792458.0 ! [m/s] light speed
   real*8,parameter        :: kb = 1.38064852e-23 ! [m^2/kg/s^2/K] Boltzmann constant
   real*8,parameter        :: Msun = 1.989e30 ! [kg] solar mass
   real*8,parameter        :: Lsun = 3.828e26 ! [W] solar luminosity
   real*4,parameter        :: pc = 3.0856776e+16 ! [m]
   real*4,parameter        :: kpc = 3.0856776e+19 ! [m]
   real*4,parameter        :: Mpc = 3.0856776e+22 ! [m]
   real*4,parameter        :: Gpc = 3.0856776e+25 ! [m]
   real*4,parameter        :: degree = pi/180.0 ! [rad]
   real*8,parameter        :: ASphereMpc = 1.1964952e+46_8 ! [m^2] surface area of a sphere with 1 Mpc radius (=4*pi*(Mpc/m)^2)
   real*4,parameter        :: L2MHI = 6.27e-9 ! (LHI/Lsun)/(MHI/Msun)
   
   ! 90-degree rotation matrices
   real*4                  :: rot(3,3,-6:6)
   
contains

   subroutine initialize_constants
   
      implicit none
      integer*4 :: i
      
      rot(:,:,0) = 0
      
      ! proper rotations
      rot(:,:,1) = reshape((/+1,+0,+0,+0,+1,+0,+0,+0,+1/),(/3,3/))   ! identity
      rot(:,:,2) = reshape((/-1,+0,+0,+0,+0,+1,+0,+1,+0/),(/3,3/))   ! invert x-axis, while permuting y and z
      rot(:,:,3) = reshape((/+0,+0,+1,+0,-1,+0,+1,+0,+0/),(/3,3/))   ! invert y-axis, while permuting z and x
      rot(:,:,4) = reshape((/+0,+1,+0,+1,+0,+0,+0,+0,-1/),(/3,3/))   ! invert z-axis, while permuting x and y
      rot(:,:,5) = reshape((/+0,+1,+0,+0,+0,+1,+1,+0,+0/),(/3,3/))   ! permute (x,y,z) -> (y,z,x)
      rot(:,:,6) = reshape((/+0,+0,+1,+1,+0,+0,+0,+1,+0/),(/3,3/))   ! permute (x,y,z) -> (z,x,y)
      
      ! improper rotations (with axis flip)
      do i = 1,6
         rot(:,:,-i) = -rot(:,:,i)
      end do
      
   end subroutine initialize_constants

end module module_constants