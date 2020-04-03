! Cosmology module
! NB: only functions needed by stringray are included

module module_cosmology

   use module_constants
   use module_types
   
   public

contains

   pure elemental function redshift_to_dc(zmax)
      
      implicit none
      real,intent(in)      :: zmax
      real                 :: redshift_to_dc
      real                 :: dz,integral
      real,allocatable     :: z(:) ! fct which the integrand depends on
      real,allocatable     :: f(:) ! integrand
      integer,allocatable  :: w(:) ! weight
      integer              :: i,n
      real                 :: DH
      real                 :: OmegaK
      
      OmegaK = 1-para%omega_m-para%omega_l
      DH = c/para%h*1e-5   ! [Mpc] Hubble distance = c/H0
      dz = 0.01
      
      ! Simpson Integration
      n = int(zmax/dz/2.0)*2+2    ! even number of subintervals
      dz = zmax/real(n)            ! length of subinterval
      allocate(z(0:n),w(0:n),f(0:n))
      z = (/(i*dz,i=0,n)/)
      w(0)=1; w(1)=4; w(n)=1
      do i=2,n-2,2 ; w(i)=2 ; w(i+1)=4 ; end do
      f = 1/sqrt(para%omega_m*(1+z)**3+OmegaK*(1+z)**2+para%omega_l) ! = 1/E(z)
      integral = sum(f*w)*dz/3.0
      ! End Simpson Integration
      
      redshift_to_dc = DH*integral
      
   end function redshift_to_dc

   pure elemental function dc_to_redshift(dc)
  
      ! redshift from comoving distance in Mpc

      implicit none
      real,intent(in)  :: dc
      real             :: dc_to_redshift
      real             :: zmin,zmax,zmean,dcmean,diff
      integer          :: i
      real,parameter   :: precision = 0.00001

      zmin = 0
      zmax = 100.0

      do i=1,100
         zmean  = (zmin+zmax)/2.0
         dcmean = redshift_to_dc(zmean)
         diff   = dcmean-dc
         if (abs(diff)<precision) then
            exit
         else
            if (diff<0.0) then
               zmin = zmean
            else
               zmax = zmean
            end if
         end if
      end do
      
      dc_to_redshift = zmean
     
  end function dc_to_redshift

end module module_cosmology
