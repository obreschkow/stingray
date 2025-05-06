! **********************************************************************************************************************************
! Shared Fortran module defining physical and mathematical constants
! Developed by Danail Obreschkow
! **********************************************************************************************************************************

module shared_module_constants

   private
   
   public   :: const    ! mathematical and physical constants
   public   :: unit     ! non-SI units expressed in SI units
   public   :: planck   ! planck units expressed in SI units
   public   :: pi,pi_8
   public   :: example_constants ! subroutine giving some examples of using this module
   
   type type_const
   
      ! dimensionless mathematical constants
      real*4   :: pi = 3.14159265359
      real*4   :: e = 2.7182818284
      real*4   :: golden = 1.618033988749894 ! golden ratio
      real*4   :: identity2(2,2) = reshape((/1.0,0.0,0.0,1.0/),(/2,2/))
      real*4   :: identity3(3,3) = reshape((/1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0/),(/3,3/))
      
      ! physical constants in SI units
      real*4   :: G = 6.67430e-11 ! [m^3/kg/s^2] gravitational constant
      real*4   :: c = 299792458.0 ! [m/s] vacuum light speed
      real*4   :: q = 1.602176634e-19 ! [C] elementary charge (=charge of the electron)
      real*4   :: h = 6.62607015e-34 ! [Js] Planck constant
      real*4   :: hbar = 1.0545718e-34 ! [Js] reduced Planck constant = h/(2*pi)
      real*4   :: kb = 1.38064852e-23 ! [m^2/kg/s^2/K] Boltzmann constant
      real*4   :: avogadro = 6.02214076e23 ! number of particles in a mol
      real*4   :: epsilon0 = 8.8541878128e-12 ! [F/m] vacuum electric permittivity
      real*4   :: mu0 = 1.25663706212e-6 ! [N/A^2] vacuum magnetic permeability
      real*4   :: me = 9.1093837015-31 ! [kg] electron mass
      real*4   :: mp = 1.6726219e-27 ! [kg] proton mass
      real*4   :: alpha = 7.2973525693e-3 ! [-] fine structure constant
      real*4   :: R = 8.31446261815324 ! [J/K/mol] ideal gas constant
      real*4   :: Rydberg = 10973731.568160 ! [1/m] Rydberg constant
      
   end type type_const
   
   type type_planck
      
      ! planck units
      real*4   :: length = 1.616255e-35 ! [m] Planck length
      real*4   :: mass = 2.176435-8 ! [kg] Planck mass
      real*4   :: time = 5.391247-44 ! [s] Planck time
      real*4   :: charge = 1.87554603778e-18 ! [C] Planck charge
      real*4   :: temperature = 1.416785e32 ! [K] Planck temperature
   
   end type type_planck
   
   type type_unit
   
      ! time
      real*4   :: year = 31557600 ! [s] Julian year, defined as exactly 365.25 days
      real*4   :: day = 86400 ! [s]
      real*4   :: hour = 3600 ! [s]
      real*4   :: min = 60 ! [s]
      
      ! length
      real*4   :: inch = 0.0254 ! [m] exact definition of one inch
      real*4   :: ft = 0.3044 ! [m] exact definition of one foot (=12inch)
      real*4   :: yd = 0.9144 ! [m] exact definition of one yard (=3ft)
      real*4   :: mile = 1609.344 ! [m] exact international definition of a mile
      real*4   :: nauticalmile = 1852 ! [m] exact international definiton of a nautical mile (~ 1 arc second on Earth)
      real*4   :: AU = 1.495978707e11 ! [m] astronomical unit
      real*4   :: pc = 3.0856776e+16 ! [m]
      real*4   :: kpc = 3.0856776e+19 ! [m]
      real*4   :: Mpc = 3.0856776e+22 ! [m]
      real*4   :: Gpc = 3.0856776e+25 ! [m]
      
      ! mass
      real*4   :: pound = 0.45359237 ! [kg] exact definition
      real*4   :: Msun = 1.98847e30 ! [kg] solar mass (approximate)
      real*4   :: Mearth = 5.97237e24 ! [kg] Earth mass (approximate)
      
      ! power/luminosity
      real*4   :: Lsun = 3.828e26 ! [W] solar luminosity (approximate)
      
      ! velocity
      real*4   :: kmh = 1.0/3.6 ! [m/s] exact value of one kilometer per hour
      real*4   :: mph = 0.44704 ! [m/s] exact value of one mile per hour
      
      ! energy
      real*4   :: eV = 1.602176634e-19 ! [J] electron volt
      real*4   :: erg = 1e-7 ! [J]
      
      ! angles
      real*4   :: degree = 3.14159265359/180.0 ! [rad] one degree
      real*4   :: arcmin = 3.14159265359/180.0/60.0 ! [rad] arc minute
      real*4   :: arcsec = 3.14159265359/180.0/3600.0 ! [rad] arc second
   
   end type type_unit
   
   type(type_const),parameter    :: const = type_const()
   type(type_unit),parameter     :: unit = type_unit()
   type(type_planck),parameter   :: planck = type_planck()
   real*4,parameter              :: pi = 3.14159265359
   real*8,parameter              :: pi_8 = 3.1415926535897932384626433
   
contains

   subroutine example_constants
      
      write(*,*) 'Planck length = ',planck%length
      write(*,*) 'Speed of light = ',const%c
      write(*,*) 'Astronomical unit = ',unit%au
      
   end subroutine example_constants

end module shared_module_constants