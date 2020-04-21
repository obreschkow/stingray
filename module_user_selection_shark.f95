! **********************************************************************************************************************************
! This module specifies selection functions of different galaxy surveys
!
! To add a new selection function, follow these steps:
! 1) Copy the function selection_example() to the bottom of the code
! 2) Change the function name to the name your survey: selection_[survey name]()
! 3) Code up the selection function by completing the four 'case' clauses; use selected = .true. for clauses that are not required
! 4) Link the selection function to a survey name by adding a new 'case' clause in the subroutine assign_selection_function
!
! Notes:
! + see the example of the WALLABY survey for an example of concatenating two simulations
! **********************************************************************************************************************************


module module_user_selection

! **********************************************************************************************************************************
! INTERFACE (DO NOT EDIT)
! **********************************************************************************************************************************

use shared_module_core
use shared_module_maths
use shared_module_cosmology
use module_global
use module_interface
use module_conversion
use module_user_routines
   
public

procedure(selection_all),pointer,protected :: selection_function => NULL()

contains

! Default selection function, that always returns .true.
! Do not remove or edit this function

logical function selection_all(pos,sam,sky) result(selected)

   implicit none
   type(type_pos),intent(in),optional           :: pos
   type(type_sam),intent(in),optional           :: sam
   class(type_sky_galaxy),intent(in),optional   :: sky
   
   call nil(pos,sam,sky)   ! avoids compiler warnings for unused arguments
   selected = .true.       ! all galaxies are selected, which satisfy the position selection in the parameter file
   
end function


! **********************************************************************************************************************************
! LINKING SURVEY NAMES TO SELECTION FUNCTIONS
! **********************************************************************************************************************************

! All custom selection functions must be linked to a survey name (i.e. the parameter "survey" in the parameter file) via a
! case clause in the following subroutine

subroutine assign_selection_function

   select case (trim(para%survey))
   case('all');               selection_function => selection_all
   case('example');           selection_function => selection_example
   case('gama');              selection_function => selection_gama
   case('wallaby-micro');     selection_function => selection_wallaby_micro
   case('wallaby-medi');      selection_function => selection_wallaby_medi
   case default
      call error('unknown survey name: ',trim(para%survey))
   end select   

end subroutine assign_selection_function


! **********************************************************************************************************************************
! CUSTOM SELECTION FUNCTIONS
! **********************************************************************************************************************************

! Default example (do not edit)

logical function selection_example(pos,sam,sky) result(selected)

   ! do not edit
   implicit none
   type(type_pos),intent(in),optional           :: pos
   type(type_sam),intent(in),optional           :: sam
   class(type_sky_galaxy),intent(in),optional   :: sky
   ! end do not edit
   
   select case (selection_type(pos,sam,sky))
   case (select_by_pos)
   
      ! here enter a selection function that is a restrictive as possible, using only pos%ra, pos%dec, pos%dc
      ! selection is applied in addition to the selection already specified in the parameter file
      selected = pos%ra>=165.0 .and. pos%ra<=195.0 .and. pos%dec>=-10.0 .and. pos%dec<=10.0 .and. pos%dc<200.0
      
   case (select_by_sam)
      
      ! here enter a selection function that is a restrictive as possible, using only the SAM-properties in type type_sam
      selected = sam%mstars_disk>1e8
      
   case (select_by_pos_and_sam)
   
      ! here enter any additional selections that require both the position *and* SAM-properties
      selected = sam%mstars_disk/pos%dc**2>1e4 ! rough preselection in the stellar mass-distance plane, only for acceleration
      
   case (select_by_all)
   
      ! here enter any additional selections that require apparent properties in the class type_sky_galaxy (see user_routines)
      selected = sky%mag<19.0 .and. sky%zobs<0.1 ! select by apparent magnitude and redshift
      
   end select
   
end function

! **********************************************************************************************************************************

! GAMA survey
! This function defines a pre-selection of the GAMA galaxy survey, before applying the final cut by apparent magnitude. To apply
! the final cuts the stingray output must first be post-processed by Viperfish (by A. Robotham), which generates observer-frame
! SEDs as a function of the star formation history and position of each source.

logical function selection_gama(pos,sam,sky) result(selected)

   implicit none
   type(type_pos),intent(in),optional           :: pos
   type(type_sam),intent(in),optional           :: sam
   class(type_sky_galaxy),intent(in),optional   :: sky
   
   real*4            :: mstars ! [Msun] stellar mass
   real*4            :: dl ! [Mpc] comoving distance
   real*4            :: mag ! generic apparent magnitude assuming M/L=1
   real*4,parameter  :: dmag = 2.0 ! magnitude tolerance
   
   select case (selection_type(pos,sam,sky))
   case (select_by_pos)
      selected = ((pos%ra>= 30.200).and.(pos%ra<= 38.800).and.(pos%dec>=-10.250).and.(pos%dec<= -3.720)).or. & ! field G02
               & ((pos%ra> 129.000).and.(pos%ra<=141.000).and.(pos%dec>= -2.000).and.(pos%dec<= +3.000)).or. & ! field G09
               & ((pos%ra>=174.000).and.(pos%ra<=186.000).and.(pos%dec>= -3.000).and.(pos%dec<= +2.000)).or. & ! field G12
               & ((pos%ra>=211.500).and.(pos%ra<=223.500).and.(pos%dec>= -2.000).and.(pos%dec<= +3.000)).or. & ! field G15
               & ((pos%ra>=339.000).and.(pos%ra<=351.000).and.(pos%dec>=-35.000).and.(pos%dec<=-30.000))       ! field G23
   case (select_by_sam)
      selected = sam%mstars_disk>1e8
   case (select_by_pos_and_sam)
      ! note: The selection below is based on a rough estimate of a generic apparent magnitude mag.
      !       This same magnitude is computed later and stored in sky%mag, hence one could also use sky%mag directly under the
      !       case clause 'select_by_all'. The point of pre-computing the magnitude here, using the comoving distance as an
      !       inferior limit for the luminosity distance, is just to accelerate the computation by avoiding time-consuming
      !       computations of many apparent galaxy properties
      mstars = (sam%mstars_disk+sam%mstars_bulge)/para%h ! [Msun]
      dl = pos%dc/para%h ! [Mph] comoving distance as an inferior limit for the luminosity distance, which would require sky%zobs
      mag = convert_absmag2appmag(convert_stellarmass2absmag(mstars,1.0),dl)
      selected = ((mag<=19.8+dmag).and.(pos%ra<330.0)).or.(mag<=19.2+dmag)
   case (select_by_all)
      selected = .true.
   end select
   
end function

! **********************************************************************************************************************************

! WALLABY survey on the ASKAP telescope

logical function selection_wallaby_micro(pos,sam,sky) result(selected)

   implicit none
   type(type_pos),intent(in),optional           :: pos
   type(type_sam),intent(in),optional           :: sam
   class(type_sky_galaxy),intent(in),optional   :: sky
   
   selected = selection_wallaby(dcmin=0.0,dcmax=60.0,pos=pos,sam=sam,sky=sky)

end function

logical function selection_wallaby_medi(pos,sam,sky) result(selected)

   implicit none
   type(type_pos),intent(in),optional           :: pos
   type(type_sam),intent(in),optional           :: sam
   class(type_sky_galaxy),intent(in),optional   :: sky
   
   selected = selection_wallaby(dcmin=60.0,dcmax=1000.0,pos=pos,sam=sam,sky=sky)

end function

logical function selection_wallaby(dcmin,dcmax,pos,sam,sky) result(selected)

   implicit none
   real*4,intent(in)                            :: dcmin ! [sim length unit = Mpc/h] minimum comoving distance
   real*4,intent(in)                            :: dcmax ! [sim length unit = Mpc/h] maximum comoving distance
   type(type_pos),intent(in),optional           :: pos
   type(type_sam),intent(in),optional           :: sam
   class(type_sky_galaxy),intent(in),optional   :: sky
   
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
   
   select case (selection_type(pos,sam,sky))
   case (select_by_pos)
      selected = pos%dec<=30.0 .and. pos%dc>=dcmin .and. pos%dc<dcmax      
   case (select_by_sam)
      selected = .true.
   case (select_by_pos_and_sam)
      mhi = (sam%matom_disk+sam%matom_bulge)/para%h/1.35
      selected = mhi>wallaby_fmin*(pos%dc/para%h)**2 ! rough preselection to accelerate computation
   case (select_by_all)
      mhi = (sam%matom_disk+sam%matom_bulge)/para%h/1.35
      da = sky%dc/(1+sky%zobs)/para%h*1e3 ! [kpc]
      dhi = 10.0**(0.506*log10(mhi)-3.293)/da/unit%arcsec ! [arcsec] 
      d_beam = wallaby_beam*(1+sky%zobs) ! [arcsec]
      n_beams_major = max(1.0,dhi/d_beam)
      n_beams_minor = max(1.0,dhi*cos(sky%inclination)/d_beam)
      n_beams = nint(0.25*pi*n_beams_major*n_beams_minor)
      n_channels = max(1.0,sky%hiline_shape%w50/wallaby_channel)
      noise = wallaby_noise_wm2*sqrt(n_channels*n_beams) ! [W/m^2] noise threshold level
      selected = sky%hiline_flux_int>noise*wallaby_sn .and. sky%zobs<=wallaby_zmax
   end select
   
end function

! **********************************************************************************************************************************
   
end module module_user_selection