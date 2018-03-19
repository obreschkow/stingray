! This module initializes and defines global constants shared by different modules

module module_constants

   use module_types

   public
   
   ! Universal constants
   
   real*4,parameter        :: pi = 3.14159265359
   real*8,parameter        :: Msun = 1.989e30 ! [kg] solar mass
   real*4,parameter        :: Lsun = 3.828e26 ! [W] solar luminosity
   real*4,parameter        :: pc = 3.0856776e+16 ! [m]
   real*4,parameter        :: kpc = 3.0856776e+19 ! [m]
   real*4,parameter        :: Mpc = 3.0856776e+22 ! [m]
   real*4,parameter        :: Gpc = 3.0856776e+25 ! [m]
   real*4,parameter        :: c = 299792458 ! [m/s] light speed
   real*8,parameter        :: ASphereMpc = 1.1964952e+46_8 ! [m^2] surface area of a sphere with 1 Mpc radius (=4*pi*(Mpc/m)^2)
   real*4,parameter        :: LMratioHI = 6.27e-9 ! (LHI/Lsun)/(MHI/Msun)
   
   ! Constant parameters, specific to a mock cone
   
   type(type_para)         :: para

end module module_constants