module module_user_selection

! **********************************************************************************************************************************
! This module defines four types of selections functions, which are called sequentially:
! 1) pos_selection: select purely on the comoving position in the sky
! 2) sam_selection: select purely on sam properties
! 3) pre_selection: select purely on sam properties AND position in sky
! 4) sky_selection: final selection using all information, including apparent properties
! Each consecutive selection function should be as restrictive as possible to avoid unnecessary computations
! **********************************************************************************************************************************

use shared_module_core
use shared_module_maths
use shared_module_cosmology
use module_global
use module_conversion
use module_user_routines

public   :: pos_selection
public   :: sam_selection
public   :: pre_selection
public   :: sky_selection

private
real*4,parameter  :: wallaby_channel = 4.0 ! [km/s]
real*4,parameter  :: wallaby_beam = 30.0 ! [arcsec] beam diameter (at z=0, scales as 1+z)
real*4,parameter  :: wallaby_sn = 5.0 ! integrated S/N for a detection
real*4,parameter  :: wallaby_noise_jy = 1.6e-3 ! [Jy] channel noise per beam
real*4,parameter  :: wallaby_noise_wm2 = wallaby_noise_jy*(wallaby_channel*1e3/0.21)*1e-26 ! [W/m^2] channel noise per beam
real*4,parameter  :: wallaby_zmax = 0.26

contains

logical function pos_selection(dc,ra,dec) result(selected)
   
   implicit none
   real*4,intent(in) :: dc    ! [simulation units] comoving distance
   real*4,intent(in) :: ra    ! [deg] right ascension
   real*4,intent(in) :: dec   ! [deg] declination
   
   call nil(dc,ra,dec) ! dummy to avoid compiler warnings for unused arguments
   
   select case (trim(para%survey))
   case ('test')
      selected = ((ra>=165.0).and.(ra<=195.0).and.(dec>=-10.0).and.(dec<=10.0).and.(dc<200.0))
   case ('devils')
      selected = ((ra>= 34.000).and.(ra<= 37.050).and.(dec>= -5.200).and.(dec<= -4.200)).or. &
               & ((ra>= 52.263).and.(ra<= 53.963).and.(dec>=-28.500).and.(dec<=-27.500)).or. &
               & ((ra>=149.380).and.(ra<=150.700).and.(dec>= +1.650).and.(dec<= +2.790))
   case ('gama')
      selected = ((ra>= 30.200).and.(ra<= 38.800).and.(dec>=-10.250).and.(dec<= -3.720)).or. &
               & ((ra> 129.000).and.(ra<=141.000).and.(dec>= -2.000).and.(dec<= +3.000)).or. &
               & ((ra>=174.000).and.(ra<=186.000).and.(dec>= -3.000).and.(dec<= +2.000)).or. &
               & ((ra>=211.500).and.(ra<=223.500).and.(dec>= -2.000).and.(dec<= +3.000)).or. &
               & ((ra>=339.000).and.(ra<=351.000).and.(dec>=-35.000).and.(dec<=-30.000))
   case ('alfalfa')
      selected = ((dec>= 0.000).and.(dec<= 36.000)).and. &
               & (((ra>= 112.500).and.(ra<= 247.500)).or.((ra>= 330.000).or.(ra<= 45.000))).and. &
               & (dc<260.0)
   case ('wallaby_micro')
      selected = (dec<=30.000).and.(dc<=60.0)
   case ('wallaby_medi')
      selected = (dec<=30.000).and.(dc>60.0).and.(dc<=750.0)
   case default
      selected = .false. ! to avoid compiler warning of uninitialised variable
      call error('unknown survey name: ',trim(para%survey))
   end select

end function pos_selection

logical function sam_selection(sam) result(selected)

   implicit none
   class(type_sam),intent(in) :: sam
   
   call nil(sam)  ! dummy to avoid compiler warnings for unused arguments
   
   select case (trim(para%survey))
   case ('test')
      selected = (sam%mstars_disk>1e8)
   case ('devils')
      selected = (sam%mstars_disk>1e8)
   case ('gama')
      selected = (sam%mstars_disk>1e8)
   case ('alfalfa')
      selected = (sam%mgas_disk>1e6).or.((sam%matom_disk>1e6))
   case ('wallaby_micro','wallaby_medi')
      selected = .true.
   case default
      selected = .false. ! to avoid compiler warning of uninitialised variable
      call error('unknown survey name: ',trim(para%survey))
   end select
   
end function sam_selection

logical function pre_selection(sam,dc,ra,dec) result(selected)

   implicit none

   type(type_sam),intent(in)           :: sam
   real*4,intent(in)                   :: dc    ! [simulation units] comoving distance
   real*4,intent(in)                   :: ra    ! [deg] right ascension
   real*4,intent(in)                   :: dec   ! [deg] declination
   real*4,parameter                    :: wallaby_fmin = 4.98e27*wallaby_sn*wallaby_noise_wm2*sqrt(23.5/wallaby_channel)
   real*4                              :: mhi
   
   call nil(sam,dc,ra,dec) ! dummy to avoid compiler warnings for unused arguments
   
   select case (trim(para%survey))
   case ('test')
      selected = .true.
   case ('devils')
      selected = .true.
   case ('gama')
      selected = .true.
   case ('alfalfa')
      selected = .true.
   case ('wallaby_micro','wallaby_medi')
      mhi = (sam%matom_disk+sam%matom_bulge)/para%h/1.35 ! [Msun] HI mass
      selected = mhi>wallaby_fmin*(dc/para%h)**2
   case default
      selected = .false. ! to avoid compiler warning of uninitialised variable
      call error('unknown survey name: ',trim(para%survey))
   end select
   
end function pre_selection

logical function sky_selection(sky,sam) result(selected)

   implicit none

   class(type_sky_galaxy),intent(in)   :: sky
   type(type_sam),intent(in)           :: sam
   real*4,parameter                    :: dmag = 2.0
   real*4                              :: wallaby_noise_integrated ! [W/m^2]
   real*4                              :: n_channels,n_beams,dhi,mhi,da,n_beams_major,n_beams_minor,d_beam
   
   call nil(sky,sam) ! dummy to avoid compiler warnings for unused arguments
   
   select case (trim(para%survey))
   case ('test')
      selected = .true.
   case ('devils')
      selected = sky%mag<=21.2+dmag
   case ('gama')
      selected = ((sky%mag<=19.8+dmag).and.(sky%ra<330.0*unit%degree)).or.(sky%mag<=19.2+dmag)
   case ('alfalfa')
      selected = sam%matom_disk>10000.0*sky%dc**2
   case ('wallaby_micro','wallaby_medi')
      da = sky%dc/(1+sky%zobs)/para%h*1e3 ! [kpc] angular diameter distance
      mhi = (sam%matom_disk+sam%matom_bulge)/para%h/1.35 ! [Msun] HI mass
      dhi = 10.0**(0.506*log10(mhi)-3.293)/da/unit%arcsec ! [arcsec] apparent HI diameter to fixed surface density of 1 Msun/pc^2, empirical relation from J. Wang
      d_beam = wallaby_beam*(1+sky%zobs) ! [arcsec] wallaby beam at observed frequency
      n_beams_major = max(1.0,dhi/d_beam) ! number of spherical beams along major axis
      n_beams_minor = max(1.0,dhi*cos(sky%inclination)/d_beam) ! number of spherical beams along minor axis
      n_beams = nint(0.25*pi*n_beams_major*n_beams_minor) ! number of beams
      n_channels = max(1.0,sky%hiline_shape%w50/wallaby_channel) ! number of channels
      wallaby_noise_integrated = wallaby_noise_wm2*sqrt(n_channels*n_beams) ! [W/m^2] noise threshold level
      selected = (sky%hiline_flux_int>wallaby_noise_integrated*wallaby_sn).and.(sky%zobs<=wallaby_zmax)
   case default
      selected = .false. ! to avoid compiler warning of uninitialised variable
      call error('unknown survey name: ',trim(para%survey))
   end select
   
end function sky_selection

end module module_user_selection