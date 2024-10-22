! **********************************************************************************************************************************
! This module specifies selection functions of different galaxy surveys
!
! For an illustration with comments of how a selection function is coded, please look at the subroutine selection_example().
!
! To add a new selection function, follow these steps:
! 1) Copy the function selection_example() to the bottom of the code
! 2) Change the function name to the name your survey: selection_[survey name]()
! 3) Code up the selection function by completing the five 'case' clauses; use selected = .true. for clauses that are not required
! 4) Link the selection function to a survey name by adding a new 'case' clause in the subroutine assign_selection_function
!
! Notes:
! + see the example of the WALLABY survey for an illustration of multiple surveys with very similar selection functions
! **********************************************************************************************************************************


module module_user_selection

! **********************************************************************************************************************************
! MODULE INTERFACE (DO NOT EDIT)
! **********************************************************************************************************************************

use shared_module_core
use shared_module_maths
use shared_module_cosmology
use module_global
use module_parameters
use module_conversion
use module_user_routines
use module_selection_tools

private

public   :: assign_selection_function

contains

! **********************************************************************************************************************************
! LINKING SURVEY NAMES TO SELECTION FUNCTIONS
! **********************************************************************************************************************************

! All custom selection functions must be linked to a survey name, i.e. the parameter "survey" in the parameter file, via a
! case clause in the following subroutine.

subroutine assign_selection_function

   select case (trim(para%survey))
      case('example');              selection_function => selection_example
      case('gama');                 selection_function => selection_gama
      case('devils');               selection_function => selection_devils
      case('waves-g23');            selection_function => selection_waves_g23
      case('deep-optical');         selection_function => selection_deep_optical
      case('deep-optical_narrow');  selection_function => selection_deep_optical_narrow
      case('wallaby-micro');        selection_function => selection_wallaby_micro
      case('wallaby-medi');         selection_function => selection_wallaby_medi
      case('dsa-wide');             selection_function => selection_dsa_wide
      case('dsa-pulsar');           selection_function => selection_dsa_pulsar
      case('dsa-deep');             selection_function => selection_dsa_deep
   case default
      call selection_function_unknown
   end select   

end subroutine assign_selection_function


! **********************************************************************************************************************************
! CUSTOM SELECTION FUNCTIONS (private, only accessible via the public pointer "selection_function")
! **********************************************************************************************************************************

! Default example (do not edit this example, but use it as a template for new selection functions)

subroutine selection_example(pos,sam,sky,range,selected)

   ! do not edit
   implicit none
   type(type_spherical),intent(in),optional     :: pos      ! has components dc [length unit of parameterfile], ra [deg], dec [deg]
   type(type_sam),intent(in),optional           :: sam      ! has components as defined in the "module_user_routines_..."
   type(type_sky_galaxy),intent(in),optional    :: sky      ! has components as defined in the "module_user_routines_..."
   type(type_fov),intent(inout),optional        :: range    ! has components dc(2) [length unit], ra(2) [deg], dec(2) [deg]
   logical,intent(inout),optional               :: selected
   ! end do not edit
   
   select case (selection_type(pos,sam,sky,range,selected))
   case (return_position_range)
   
      ! here enter the individual maximal ranges of comoving distance, right ascension and declination covered by the survey,
      ! as restrictive as possible; these ranges are mandatory
      ! Note: it is possible to have the range%ra go from a larger number (e.g. 330) to a smaller one (e.g. 30). In this case,
      !       the sky defined by the wedge from 330 deg (=-30deg) to 30 deg is considered.
      range%dc = (/0.0,300.0/)      ! [simulation length units, here Mpc/h] comoving distance range
      range%ra = (/150.0,210.0/)    ! [deg] range of right ascensions, bound to 0 to 360
      range%dec = (/-10.0,10.0/)    ! [deg] range of declinations, bound to -90 to +90
      
   case (select_by_pos)
   
      ! here enter additional restrictions on the position (pos)
      ! if all positional restrictions are already covered by the ranges above, leave this clause empty
      selected = pos%ra>=190.0 .or. pos%ra<=170.0
      
   case (select_by_sam)
      
      ! here add additional, maximal restrictions that only use the SAM-properties in type type_sam;
      ! if no such restrictions exist, leave this clause empty
      selected = sam%mstars_disk+sam%mstars_bulge>1e8
      
   case (select_by_pos_and_sam)
   
      ! here add additional, maximal restrictions that require both the position (pos) *and* SAM-properties (sam);
      ! if no such restrictions exist, leave this clause empty
      selected = (sam%mstars_disk+sam%mstars_bulge)/pos%dc**2>1e4 ! rough preselection, only for acceleration
      
   case (select_by_all)
   
      ! here add additional, maximal restrictions that require apparent sky properties (sky), as defined in module_user_routines,
      ! possibly combined with position and SAM properties;
      ! if no such restrictions exist, leave this clause empty
      selected = sky%mag<19.23 .and. sky%zobs<0.1 ! select by apparent magnitude and redshift
      
   end select
   
end subroutine

! **********************************************************************************************************************************

! GAMA survey
! This function defines a pre-selection of the GAMA galaxy survey, before applying the final cut by apparent magnitude. To apply
! the final cuts the stingray output must first be post-processed by Viperfish (by A. Robotham), which generates observer-frame
! SEDs as a function of the star formation history and position of each source.

subroutine selection_gama(pos,sam,sky,range,selected)

   implicit none
   type(type_spherical),intent(in),optional     :: pos
   type(type_sam),intent(in),optional           :: sam
   type(type_sky_galaxy),intent(in),optional    :: sky
   type(type_fov),intent(inout),optional        :: range
   logical,intent(inout),optional               :: selected
   
   ! computation variables
   real*4            :: mstars ! [Msun] stellar mass
   real*4            :: dl ! [Mpc] comoving distance
   real*4            :: mag ! generic apparent magnitude assuming M/L=1
   real*4,parameter  :: dmag = 4.0 ! magnitude tolerance
   
   ! selection function
   select case (selection_type(pos,sam,sky,range,selected))
   case (return_position_range)
      range%dc = (/0.0,1700.0/)     ! [simulation length units, here Mpc/h] distance range (to z~0.7)
      range%ra = (/30.2,351.0/)     ! [deg] range of right ascensions, bound to 0 to 360
      range%dec = (/-35.0,3.0/)     ! [deg] range of declinations, bound to -90 to +90
   case (select_by_pos)
      selected = ((pos%ra>= 30.2).and.(pos%ra<= 38.8).and.(pos%dec>=-10.25).and.(pos%dec<= -3.72)).or. & ! field G02
               & ((pos%ra> 129.0).and.(pos%ra<=141.0).and.(pos%dec>= -2.00).and.(pos%dec<= +3.00)).or. & ! field G09
               & ((pos%ra>=174.0).and.(pos%ra<=186.0).and.(pos%dec>= -3.00).and.(pos%dec<= +2.00)).or. & ! field G12
               & ((pos%ra>=211.5).and.(pos%ra<=223.5).and.(pos%dec>= -2.00).and.(pos%dec<= +3.00)).or. & ! field G15
               & ((pos%ra>=339.0).and.(pos%ra<=351.0).and.(pos%dec>=-35.00).and.(pos%dec<=-30.00))       ! field G23
   case (select_by_sam)
      selected = (sam%mstars_disk+sam%mstars_bulge)/para%h>1e6
   case (select_by_pos_and_sam)
      ! note: The selection below is based on a rough estimate of a generic apparent magnitude mag.
      !       This same magnitude is computed later and stored in sky%mag, which is used onder 'selected_by_all'.
      !       The point of pre-computing the magnitude here, using the comoving distance as an
      !       inferior limit for the luminosity distance, is to accelerate the computation by avoiding time-consuming
      !       computations of many apparent galaxy properties
      mstars = (sam%mstars_disk+sam%mstars_bulge)/para%h ! [Msun]
      dl = pos%dc/para%h ! [Mph] comoving distance as an inferior limit for the luminosity distance, which would require sky%zobs
      mag = convert_absmag2appmag(convert_stellarmass2absmag(mstars,1.0),dl)
      selected = ((mag<=19.8+dmag).and.(pos%ra<330.0)).or.(mag<=19.2+dmag)
   case (select_by_all)
      selected = sky%mag<=21.2+dmag
   end select
   
end subroutine

! **********************************************************************************************************************************

! DEVILS survey
! This function defines a pre-selection of the DEVILS galaxy survey, before applying the final cut by apparent magnitude. To apply
! the final cuts the stingray output must first be post-processed by Viperfish (by A. Robotham), which generates observer-frame
! SEDs as a function of the star formation history and position of each source.

subroutine selection_devils(pos,sam,sky,range,selected)

   implicit none
   type(type_spherical),intent(in),optional     :: pos
   type(type_sam),intent(in),optional           :: sam
   type(type_sky_galaxy),intent(in),optional    :: sky
   type(type_fov),intent(inout),optional        :: range
   logical,intent(inout),optional               :: selected
   
   ! computation variables
   real*4            :: mag ! rough estimate of a apparent magnitude assuming M/L=1
   real*4,parameter  :: dmag = 4.0 ! magnitude tolerance
   
   ! selection function
   select case (selection_type(pos,sam,sky,range,selected))
   case (return_position_range)
      range%dc = (/0.0,2350.00/) ! [simulation length units, here Mpc/h] distance range (to z~1.0)
      range%ra = (/34.0,150.70/) ! [deg] range of right ascensions, bound to 0 to 360
      range%dec = (/-28.5,2.79/) ! [deg] range of declinations, bound to -90 to +90
   case (select_by_pos)
      selected = ((pos%ra>= 34.000).and.(pos%ra<= 37.050).and.(pos%dec>= -5.200).and.(pos%dec<= -4.200)).or. & ! D02 (XMM-LSS)
               & ((pos%ra>= 52.263).and.(pos%ra<= 53.963).and.(pos%dec>=-28.500).and.(pos%dec<=-27.500)).or. & ! D03	(ECDFS)
               & ((pos%ra>=149.380).and.(pos%ra<=150.700).and.(pos%dec>= +1.650).and.(pos%dec<= +2.790))       ! D10	(COSMOS)
   case (select_by_sam)
      selected = (sam%mstars_disk+sam%mstars_bulge)/para%h>1e6
   case (select_by_pos_and_sam)
      ! note: see comments in selection_gama
      mag = convert_absmag2appmag(convert_stellarmass2absmag((sam%mstars_disk+sam%mstars_bulge)/para%h,1.0),pos%dc/para%h)
      selected = mag<=21.2+dmag
   case (select_by_all)
      selected = sky%mag<=21.2+dmag
   end select
   
end subroutine

! **********************************************************************************************************************************

! WAVES G23 survey

subroutine selection_waves_g23(pos,sam,sky,range,selected)

   implicit none
   type(type_spherical),intent(in),optional     :: pos
   type(type_sam),intent(in),optional           :: sam
   type(type_sky_galaxy),intent(in),optional    :: sky
   type(type_fov),intent(inout),optional        :: range
   logical,intent(inout),optional               :: selected
   
   ! computation variables
   real*4            :: mag ! rough estimate of a apparent magnitude assuming M/L=1
   real*4,parameter  :: dmag = 4.0 ! magnitude tolerance
   
   ! selection function
   select case (selection_type(pos,sam,sky,range,selected))
   case (return_position_range)
      range%dc = (/0.0,2350.0/) ! [simulation length units, here Mpc/h]
      range%ra = (/339.0,351.0/) ! [deg] range of right ascensions, bound to 0 to 360
      range%dec = (/-35.0,-30.0/) ! [deg] range of declinations, bound to -90 to +90
   case (select_by_pos)
   case (select_by_sam)
      selected = (sam%mstars_disk+sam%mstars_bulge)/para%h>1e6
   case (select_by_pos_and_sam)
      ! note: see comments in selection_gama
      mag = convert_absmag2appmag(convert_stellarmass2absmag((sam%mstars_disk+sam%mstars_bulge)/para%h,1.0),pos%dc/para%h)
      selected = mag<=24.0+dmag
   case (select_by_all)
      selected = sky%mag<=24.0+dmag
   end select
   
end subroutine

! **********************************************************************************************************************************

! deep-optical survey

subroutine selection_deep_optical(pos,sam,sky,range,selected)

   implicit none
   type(type_spherical),intent(in),optional     :: pos
   type(type_sam),intent(in),optional           :: sam
   type(type_sky_galaxy),intent(in),optional    :: sky
   type(type_fov),intent(inout),optional        :: range
   logical,intent(inout),optional               :: selected
   
   ! computation variables
   real*4,parameter  :: dmag = 4.0 ! magnitude tolerance
   
   ! selection function
   select case (selection_type(pos,sam,sky,range,selected))
   case (return_position_range)
      range%dc = (/0.0,5695.765/) ! [simulation length units, here Mpc/h] (out to z=6)
      range%ra = (/211.500,223.500/) ! [deg]
      range%dec = (/-4.5,4.5/) ! [deg]
   case (select_by_pos)
   case (select_by_sam)
      selected = (sam%mstars_disk+sam%mstars_bulge)/para%h>1e8
   case (select_by_pos_and_sam)
   case (select_by_all)
      selected = sky%mag<=28.0+dmag
   end select
   
end subroutine

! **********************************************************************************************************************************

! deep-optical-narrow survey

subroutine selection_deep_optical_narrow(pos,sam,sky,range,selected)

   implicit none
   type(type_spherical),intent(in),optional     :: pos
   type(type_sam),intent(in),optional           :: sam
   type(type_sky_galaxy),intent(in),optional    :: sky
   type(type_fov),intent(inout),optional        :: range
   logical,intent(inout),optional               :: selected
   
   ! computation variables
   real*4,parameter  :: dmag = 4.0 ! magnitude tolerance
   
   ! selection function
   select case (selection_type(pos,sam,sky,range,selected))
   case (return_position_range)
      range%dc = (/0.0,6173.688/) ! [simulation length units, here Mpc/h] (out to z=8)
      range%ra = (/211.500,223.500/) ! [deg]
      range%dec = (/-2.5,2.5/) ! [deg]
   case (select_by_pos)
   case (select_by_sam)
      selected = (sam%mstars_disk+sam%mstars_bulge)/para%h>1e6
   case (select_by_pos_and_sam)
   case (select_by_all)
      selected = sky%mag<=38.0+dmag
   end select
   
end subroutine

! **********************************************************************************************************************************

! WALLABY survey on the ASKAP telescope

subroutine selection_wallaby_micro(pos,sam,sky,range,selected)

   implicit none
   type(type_spherical),intent(in),optional     :: pos
   type(type_sam),intent(in),optional           :: sam
   type(type_sky_galaxy),intent(in),optional    :: sky
   type(type_fov),intent(inout),optional        :: range
   logical,intent(inout),optional               :: selected
   
   call selection_wallaby(dcmin=0.0,dcmax=60.0,pos=pos,sam=sam,sky=sky,range=range,selected=selected)

end subroutine

subroutine selection_wallaby_medi(pos,sam,sky,range,selected)

   implicit none
   type(type_spherical),intent(in),optional     :: pos
   type(type_sam),intent(in),optional           :: sam
   type(type_sky_galaxy),intent(in),optional    :: sky
   type(type_fov),intent(inout),optional        :: range
   logical,intent(inout),optional               :: selected
   
   call selection_wallaby(dcmin=60.0,dcmax=1000.0,pos=pos,sam=sam,sky=sky,range=range,selected=selected)

end subroutine

subroutine selection_wallaby(dcmin,dcmax,pos,sam,sky,range,selected)

   implicit none
   real*4,intent(in)                            :: dcmin ! [sim length unit = Mpc/h] minimum comoving distance
   real*4,intent(in)                            :: dcmax ! [sim length unit = Mpc/h] maximum comoving distance
   type(type_spherical),intent(in),optional     :: pos
   type(type_sam),intent(in),optional           :: sam
   type(type_sky_galaxy),intent(in),optional    :: sky
   type(type_fov),intent(inout),optional        :: range
   logical,intent(inout),optional               :: selected
   
   ! wallaby survey parameters
   real*4,parameter  :: wallaby_channel = 4.0 ! [km/s]
   real*4,parameter  :: wallaby_beam = 30.0 ! [arcsec] FWHM synthesised beam (at z=0, scales as 1+z)
   real*4,parameter  :: wallaby_sn = 5.0 ! integrated S/N for a detection
   real*4,parameter  :: wallaby_noise_jy = 1.6e-3 ! [Jy] channel noise per beam
   real*4,parameter  :: wallaby_zmax = 0.26
   
   ! derived survey parameters
   real*4,parameter  :: wallaby_noise_wm2 = wallaby_noise_jy*(wallaby_channel*1e3/0.21)*1e-26 ! [W/m^2] channel noise per beam
   real*4,parameter  :: wallaby_fmin = 4.98e27*wallaby_sn*wallaby_noise_wm2*sqrt(23.5/wallaby_channel) ! min. flux in sim units
   
   ! computation variables
   real*4   :: mhi            ! [Msun] HI mass, excluding helium
   real*4   :: da             ! [kpc] angular diameter distance
   real*4   :: dhi            ! [arcsec] apparent HI diameter at surface density of 1 Msun/pc^2, empirical relation from J. Wang
   real*4   :: d_beam         ! [arcsec] wallaby beam FWHM at observed frequency
   real*4   :: n_beams_major  ! number of spherical beams along major axis
   real*4   :: n_beams_minor  ! number of spherical beams along minor axis
   real*4   :: n_beams        ! total number of beams
   real*4   :: n_channels     ! number of channels at W50
   real*4   :: noise          ! [W/m^2] integrated noise
   
   ! selection function
   select case (selection_type(pos,sam,sky,range,selected))
   case (return_position_range)
      range%dc = (/dcmin,dcmax/)
      range%ra = (/0.0,360.0/)
      range%dec = (/-90.0,30.0/)
   case (select_by_pos)
   case (select_by_sam)
   case (select_by_pos_and_sam)
      mhi = (sam%matom_disk+sam%matom_bulge)/para%h/1.35 ! [Msun] HI mass
      selected = mhi>wallaby_fmin*(pos%dc/para%h)**2 ! rough preselection to accelerate computation
   case (select_by_all)
      mhi = (sam%matom_disk+sam%matom_bulge)/para%h/1.35
      da = sky%dc/(1+sky%zobs)/para%h*1e3 ! [kpc]
      dhi = 10.0**(0.506*log10(mhi)-3.293)/da/unit%arcsec ! [arcsec] 
      d_beam = wallaby_beam*(1+sky%zobs) ! [arcsec]
      n_beams_major = max(1.0,dhi/d_beam)
      n_beams_minor = max(1.0,dhi*cos(sky%inclination)/d_beam)
      n_beams = nint(n_beams_major*n_beams_minor)
      n_channels = max(1.0,sky%hiline_shape%w50/wallaby_channel)
      noise = wallaby_noise_wm2*sqrt(n_channels*n_beams) ! [W/m^2] noise threshold level
      selected = sky%hiline_flux_int>noise*wallaby_sn .and. sky%zobs<=wallaby_zmax
   end select
   
end subroutine

! **********************************************************************************************************************************

! Hypothetical DSA-2000 survey (for Fabian Walter)

subroutine selection_dsa_wide(pos,sam,sky,range,selected)

   implicit none
   type(type_spherical),intent(in),optional     :: pos
   type(type_sam),intent(in),optional           :: sam
   type(type_sky_galaxy),intent(in),optional    :: sky
   type(type_fov),intent(inout),optional        :: range
   logical,intent(inout),optional               :: selected
   
   !call selection_dsa(decmin=-30.00000,shimin=45e-3,pos=pos,sam=sam,sky=sky,range=range,selected=selected) ! 3pi survey (old)
   call selection_dsa(decmin=-30.00000,shimin=5e-3,pos=pos,sam=sam,sky=sky,range=range,selected=selected) ! 3pi survey

end subroutine

subroutine selection_dsa_pulsar(pos,sam,sky,range,selected)

   implicit none
   type(type_spherical),intent(in),optional     :: pos
   type(type_sam),intent(in),optional           :: sam
   type(type_sky_galaxy),intent(in),optional    :: sky
   type(type_fov),intent(inout),optional        :: range
   logical,intent(inout),optional               :: selected
   
   call selection_dsa(decmin=64.56022,shimin=18e-3,pos=pos,sam=sam,sky=sky,range=range,selected=selected) ! 2000 sqdeg

end subroutine

subroutine selection_dsa_deep(pos,sam,sky,range,selected)

   implicit none
   type(type_spherical),intent(in),optional     :: pos
   type(type_sam),intent(in),optional           :: sam
   type(type_sky_galaxy),intent(in),optional    :: sky
   type(type_fov),intent(inout),optional        :: range
   logical,intent(inout),optional               :: selected
   
   call selection_dsa(decmin=86.90943,shimin=9e-3,pos=pos,sam=sam,sky=sky,range=range,selected=selected) ! 30 sqdeg

end subroutine

subroutine selection_dsa(decmin,shimin,pos,sam,sky,range,selected)

   implicit none
   real*4,intent(in)                            :: decmin ! [deg] minimum declination
   real*4,intent(in)                            :: shimin ! [Jy km/s] minimum HI flux
   type(type_spherical),intent(in),optional     :: pos
   type(type_sam),intent(in),optional           :: sam
   type(type_sky_galaxy),intent(in),optional    :: sky
   type(type_fov),intent(inout),optional        :: range
   logical,intent(inout),optional               :: selected
   
   real*4   :: mhi            ! [Msun] HI mass, excluding helium
   
   select case (selection_type(pos,sam,sky,range,selected))
   case (return_position_range)
      range%dc = (/0.0,2400.0/) ! [Mpc/h] corresponds to about zcos=0.0-1.05
      range%ra = (/0.0,360.0/)
      range%dec = (/decmin,90.0/)
   case (select_by_pos)
   case (select_by_sam)
   case (select_by_pos_and_sam)
      mhi = (sam%matom_disk+sam%matom_bulge)/para%h/1.35 ! [Msun] HI mass
      selected = mhi>2.35e5*(pos%dc/para%h)**2*shimin/2.1 ! rough preselection to accelerate computation (2.1 replaces 1+zobs)
   case (select_by_all)
      selected = sky%hiline_flux_int_vel>shimin .and. sky%zobs<=1.0
   end select
   
end subroutine

! **********************************************************************************************************************************
   
end module module_user_selection
