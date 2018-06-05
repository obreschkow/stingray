! This module is old and needs to be checked and partially rewritten

module module_cosmology

   use module_constants
   use module_types
   
   public

contains

   pure elemental function redshift_to_dc(zmax)
      ! accumulated comoving distance in Mpc
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
      OmegaK = 1-para%OmegaM-para%OmegaL
      DH = 2.9979e8/para%h*1e-5    ! Hubble distance [Mpc] = c/H0
      dz = 0.01
      ! Simpson Integration
      n   = int(zmax/dz/2.0)*2+2  ! even number of subintervals
      dz  = zmax/real(n)              ! length of subinterval
      allocate(z(0:n),w(0:n),f(0:n))
      z    = (/(i*dz,i=0,n)/)
      w(0)=1; w(1)=4; w(n)=1
      do i=2,n-2,2 ;  w(i+1)=4 ; w(i)=2 ; end do
      f = 1/sqrt(para%OmegaM*(1+z)**3+OmegaK*(1+z)**2+para%OmegaL) ! = 1/E(z)
      integral = sum(f*w)*dz/3.0
      ! End Simpson Integration
      redshift_to_dc = DH*integral
   end function redshift_to_dc

  pure elemental function redshift_to_tl(zmax)
    ! lookback time in yrs
    implicit none
    real,intent(in)     :: zmax
    real                :: redshift_to_tl
    real                :: dz,integral
    real,allocatable    :: z(:) ! fct which the integrand depends on
    real,allocatable    :: f(:) ! integrand
    integer,allocatable :: w(:) ! weight
    integer             :: i,n
    real                :: OmegaK
    real                :: tH
    tH = 9.77793e9/para%h ! Hubble time [yrs]
    OmegaK = 1-para%OmegaM-para%OmegaL
    dz = 0.01
    ! Simpson Integration
    n   = int(zmax/dz/2)*2+2  ! even number of subintervals
    dz  = zmax/n              ! length of subinterval
    allocate(z(0:n),w(0:n),f(0:n))
    z    = (/(i*dz,i=0,n)/)
    w(0)=1; w(1)=4; w(n)=1
    do i=2,n-2,2 ;  w(i+1)=4 ; w(i)=2 ; end do
    f        = 1/(1+z)/sqrt(para%OmegaM*(1+z)**3+OmegaK*(1+z)**2+para%OmegaL) ! = 1/(1+z)/E(z)
    integral = sum(f*w)*dz/3.0
    ! End Simpson Integration
    redshift_to_tl = tH*integral
  end function redshift_to_tl

  pure elemental function redshift_to_da(z)
    ! angular diameter distance in Mpc
    implicit none
    real,intent(in) :: z
    real            :: redshift_to_da
    redshift_to_da = dc_to_dm(redshift_to_dc(z))/(1.0+z)
  end function redshift_to_da

  pure elemental function dc_to_redshift(dc)
    ! redshift from accumulated comoving distance in Mpc
    implicit none
    real,intent(in) :: dc
    real            :: dc_to_redshift
    real            :: zmin,zmax,zmean,precision,dcmean,diff
    integer         :: i
    precision = 0.0001 ! used 0.001 for online s-cubed catalog, which is just a bit too low
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

  pure elemental function dc_to_dm(dc)
    ! comoving distance in Mpc to traverse comoving distance in Mpc
    implicit none
    real,intent(in)  :: dc
    real             :: dc_to_dm
    real             :: DH
    real             :: OmegaK
    OmegaK = 1-para%OmegaM-para%OmegaL
    DH = 2.9979e8/para%h*1e-5    ! Hubble distance [Mpc] = c/H0
    if (OmegaK==0) then
       dc_to_dm = dc
    else
       dc_to_dm = DH/sqrt(abs(OmegaK))*sinh(sqrt(abs(OmegaK))*dc/DH)
    end if
  end function dc_to_dm
  
  function universe_age(z)
    implicit none
    real,intent(in) :: z
    real :: HoverH0
    real :: universe_age
    real :: OmegaK
    OmegaK = 1-para%OmegaM-para%OmegaL
    HoverH0 = sqrt(para%OmegaM*(1.0+z)**3.0+OmegaK*(1.0+z)**2.0+para%OmegaL) 
    universe_age = 9.739e9/para%h/HoverH0 ! [years]
  end function universe_age

end module module_cosmology
