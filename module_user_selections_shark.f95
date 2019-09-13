module module_user_selection

! ==============================================================================================================
! This module defines four types of selections functions, which are called sequentially:
! 1) pos_selection: select purely on the comoving position in the sky
! 2) sam_selection: select purely on sam properties
! 3) sky_selection: final selection using all information, including apparent properties
! Each consecutive selection function should be as restrictive as possible
! ==============================================================================================================

use module_constants
use module_system
use module_types
use module_io
use module_linalg
use module_cosmology
use module_conversion
use module_user_routines

contains

logical function pos_selection(dc,ra,dec) result(selected)
   
   implicit none
   real*4,intent(in) :: dc    ! [simulation units] comoving distance
   real*4,intent(in) :: ra    ! [deg] right ascension
   real*4,intent(in) :: dec   ! [deg] right ascension
   
   call nil(dc,ra,dec) ! dummy to avoid compiler warnings for unused arguments
   
   select case (trim(para%survey))
   case ('test')
      selected = ((ra>=0.0).and.(ra<=60.0).and.(dec>=-5.0).and.(dec<=12.0).and.(dc<1000.0))
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
   case ('wallaby')
      selected = ((dec>= -90.000).and.(dec<=30.000).and.(ra>=0.00).and.(ra<=360.0).and.&
               & dc<=780.0)
   case default
      selected = .true.
   end select

end function pos_selection

logical function sam_selection(sam) result(selected)

   implicit none
   class(type_sam),intent(in) :: sam
   
   call nil(sam)  ! dummy to avoid compiler warnings for unused arguments
   
   select case (trim(para%survey))
   case ('test')
      selected = (sam%mstars_disk>1e10)
   case ('devils')
      selected = (sam%mstars_disk>1e8)
   case ('gama')
      selected = (sam%mstars_disk>1e8)
   case ('alfalfa')
      selected = (sam%mgas_disk>1e6).or.((sam%matom_disk>1e6))
   case ('wallaby')
      selected = (sam%mstars_disk>1e5)
   case default
      selected = .true.
   end select
   
end function sam_selection

logical function sky_selection(sky,sam) result(selected)

   class(type_sky_galaxy),intent(in)   :: sky
   type(type_sam),intent(in)           :: sam
   real*4,parameter                    :: dmag = 2.0
   
   call nil(sky,sam) ! dummy to avoid compiler warnings for unused arguments
   
   select case (trim(para%survey))
   case ('test')
      selected = .true.
   case ('devils')
      selected = sky%mag<=21.2+dmag
   case ('gama')
      selected = ((sky%mag<=19.8+dmag).and.(sky%ra<330.0*degree)).or.(sky%mag<=19.2+dmag)
   case ('alfalfa')
      selected = sam%matom_disk>10000.0*sky%dc**2
   case ('wallaby')
      selected = sky%s_hi_int>4e-24
   case default
      selected = .true.
   end select
   
end function sky_selection

end module module_user_selection