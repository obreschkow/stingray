! **********************************************************************************************************************************
! Shared Fortran module with cosmology functions
! Developed by Danail Obreschkow
! **********************************************************************************************************************************

module shared_module_cosmology

use shared_module_core

private

public   :: cosmology      ! structure with cosmological parameters
public   :: set_cosmology  ! set cosmological parameters, either to a default set or to any custom values
public   :: redshift_to_dc
public   :: dc_to_redshift

! cosmological parameters used in this module with Planck defaults (Planck 2015, p. 32, table 4, last column)
type type_cosmology

   character(255) :: name
   real*4         :: h = 0.6774        ! dimensionless Hubble parameter at z=0, H0=h*100km/s/Mpc
   real*4         :: omega_m = 0.3089  ! matter energy density relative to closure density at z=0
   real*4         :: omega_l = 0.6911  ! dark energy density relative to closure density at z=0
   
end type type_cosmology

type(type_cosmology),protected :: cosmology

contains

subroutine set_cosmology(name,h,omega_m,omega_l)
   
   implicit none
   character(*),intent(in)    :: name
   real*4,intent(in),optional :: h,omega_m,omega_l
   
   select case (trim(name))
   case ('737')
      cosmology%h = 0.7
      cosmology%omega_m = 0.3
      cosmology%omega_l = 0.7
   case ('planck')
      cosmology%h = 0.6774
      cosmology%omega_m = 0.3089
      cosmology%omega_l = 0.6911
   case ('surfs')
      cosmology%h = 0.6751
      cosmology%omega_m = 0.3121
      cosmology%omega_l = 0.6879
   case default
      if (present(h).and.present(omega_m).and.present(omega_l)) then
         cosmology%h = h
         cosmology%omega_m = omega_m
         cosmology%omega_l = omega_l
      else
         call error(trim(name)//' is an unknown cosmology.')
      end if
   end select
   
   cosmology%name = trim(name)

end subroutine set_cosmology

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
   real*4,parameter     :: c = 299792458.0 ! [m/s]
   
   OmegaK = 1-cosmology%omega_m-cosmology%omega_l
   DH = c/cosmology%h*1e-5 ! [Mpc] Hubble distance = c/H0
   dz = 0.01
   
   ! Simpson Integration
   n = int(zmax/dz/2.0)*2+2    ! even number of subintervals
   dz = zmax/real(n)            ! length of subinterval
   allocate(z(0:n),w(0:n),f(0:n))
   z = (/(i*dz,i=0,n)/)
   w(0)=1; w(1)=4; w(n)=1
   do i=2,n-2,2 ; w(i)=2 ; w(i+1)=4 ; end do
   f = 1/sqrt(cosmology%omega_m*(1+z)**3+OmegaK*(1+z)**2+cosmology%omega_l) ! = 1/E(z)
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

end module shared_module_cosmology
