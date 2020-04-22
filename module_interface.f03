module module_interface

   use shared_module_core
   use shared_module_parameters
   use shared_module_maths
   use shared_module_constants
   use module_global
   
   private
   
   public   :: is_in_fov
   public   :: selection_type
   
contains

logical function is_in_fov(pos)

   implicit none
   
   type(type_pos),intent(in)  :: pos
   
   is_in_fov = (pos%dc>=para%dc_min).and.(pos%dc<=para%dc_max).and. &
             & (pos%ra>=para%ra_min).and.(pos%ra<=para%ra_max).and. &
             & (pos%dec>=para%dec_min).and.(pos%dec<=para%dec_max)
             
end function is_in_fov

integer*4 function selection_type(pos,sam,sky) result(sel)

   implicit none
   class(*),optional :: pos
   class(*),optional :: sam
   class(*),optional :: sky
   
   if (present(pos).and.(.not.present(sam)).and.(.not.present(sky))) then
      sel = select_by_pos
   else if ((.not.present(pos)).and.present(sam).and.(.not.present(sky))) then
      sel = select_by_sam
   else if (present(pos).and.present(sam).and.(.not.present(sky))) then
      sel = select_by_pos_and_sam
   else if (present(pos).and.present(sam).and.present(sky)) then
      sel = select_by_all
   else
      call deverror('unknown selection type')
   end if

end function selection_type
             
end module module_interface