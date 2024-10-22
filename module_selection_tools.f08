module module_selection_tools

   use shared_module_core
   use module_global
   use module_parameters
   use module_user_routines

   private
   public   :: selection_function
   public   :: selection_type
   public   :: return_position_range
   public   :: select_by_pos
   public   :: select_by_sam
   public   :: select_by_pos_and_sam
   public   :: select_by_all
   public   :: selection_function_unknown

   procedure(selection_all),pointer :: selection_function => NULL()

   ! selection types
   integer*4,parameter  :: return_position_range = 0
   integer*4,parameter  :: select_by_pos = 1
   integer*4,parameter  :: select_by_sam = 2
   integer*4,parameter  :: select_by_pos_and_sam = 3
   integer*4,parameter  :: select_by_all = 4
   
contains

   subroutine selection_all(pos,sam,sky,range,selected)

      implicit none
      type(type_spherical),intent(in),optional     :: pos
      type(type_sam),intent(in),optional           :: sam
      type(type_sky_galaxy),intent(in),optional    :: sky
      type(type_fov),intent(inout),optional        :: range
      logical,intent(inout),optional               :: selected
   
      call nil(pos,sam,sky,range,selected) ! avoids compiler warnings for unused arguments
   
   end subroutine selection_all
   
   integer*4 function selection_type(pos,sam,sky,range,selected) result(sel)

      implicit none
      class(*),optional :: pos,sam,sky,range,selected
      integer*4         :: q
   
      q = log2int(present(pos))+2*log2int(present(sam))+4*log2int(present(sky))+ &
          & 8*log2int(present(range))+16*log2int(present(selected))
   
      select case(q)
      case(0+0+0+8+00); sel = return_position_range
      case(1+0+0+0+16); sel = select_by_pos
      case(0+2+0+0+16); sel = select_by_sam
      case(1+2+0+0+16); sel = select_by_pos_and_sam
      case(1+2+4+0+16); sel = select_by_all
      case default
         write(*,*) q,present(selected)
         call deverror('unknown selection type')
      end select

   end function selection_type
   
   ! developer selection functions *************************************************************************************************
   
   subroutine selection_function_unknown
   
      select case (trim(para%survey))
         case('dev1'); selection_function => selection_dev1
         case('dev2'); selection_function => selection_dev2
      case default
         call error('unknown survey name: ',trim(para%survey))
      end select   
   
   end subroutine selection_function_unknown
   
   subroutine selection_dev1(pos,sam,sky,range,selected)

      ! do not edit
      implicit none
      type(type_spherical),intent(in),optional     :: pos      ! has components dc [length unit of parameterfile], ra [deg], dec [deg]
      type(type_sam),intent(in),optional           :: sam      ! has components as defined in the "module_user_routines_..."
      type(type_sky_galaxy),intent(in),optional    :: sky      ! has components as defined in the "module_user_routines_..."
      type(type_fov),intent(inout),optional        :: range    ! has components dc(2) [length unit], ra(2) [deg], dec(2) [deg]
      logical,intent(inout),optional               :: selected
      ! end do not edit
      
      real*4   :: x(3)
      real*4,parameter  :: fraction = 0.2
   
      select case (selection_type(pos,sam,sky,range,selected))
      case (return_position_range)
   
         ! here enter the individual maximal ranges of comoving distance, right ascension and declination covered by the survey,
         ! as restrictive as possible; these ranges are mandatory
         range%dc = (/0.0,840.0/)      ! [simulation length units, here Mpc/h] comoving distance range
         range%ra = (/0.0,360.0/)      ! [deg] range of right ascensions, bound to 0 to 360
         range%dec = (/-30.0,90.0/)    ! [deg] range of declinations, bound to -90 to +90
      
      case (select_by_pos)
      
         call sph2car(pos%dc,pos%ra*unit%degree,pos%dec*unit%degree,x)
         selected = abs(x(1)/para%box_side)<0.4999
      
      case (select_by_sam)
      
         x = sam%get_position()
         selected = abs(x(1)/para%box_side-0.5)<=0.5*fraction
   
      case (select_by_pos_and_sam)
         
      case (select_by_all)
      
      end select
   
   end subroutine
   
   subroutine selection_dev2(pos,sam,sky,range,selected)

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
         range%dc = (/0.0,3000.0/)      ! [simulation length units, here Mpc/h] comoving distance range
         range%ra = (/0.0,360.0/)    ! [deg] range of right ascensions, bound to 0 to 360
         range%dec = (/-30.0,90.0/)    ! [deg] range of declinations, bound to -90 to +90
      
      case (select_by_pos)
      
      case (select_by_sam)
      
         selected = sam%is_group_center().and.modulo(sam%get_groupid(),10000_8)==1_8
         
      case (select_by_pos_and_sam)
         
      case (select_by_all)
      
      end select
   
   end subroutine

end module module_selection_tools