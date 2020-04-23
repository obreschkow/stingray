module module_interface

   use shared_module_core
   use shared_module_parameters
   use shared_module_maths
   use shared_module_constants
   use module_global
   
   private
   
   public   :: selection_type
   public   :: get_keywords
   public   :: keyword
   
   integer*4,parameter           :: n_keywords_max = 20
   integer*4,protected           :: n_keywords
   character(len=255),protected  :: keyword_name(n_keywords_max)
   
contains

integer*4 function selection_type(pos,sam,sky,range,selected) result(sel)

   implicit none
   class(*),optional :: pos,sam,sky,range,selected
   integer*4         :: q
   
   q = log2int(present(pos))+2*log2int(present(sam))+4*log2int(present(sky))+8*log2int(present(range))+16*log2int(present(selected))
   
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

subroutine get_keywords
   
   implicit none
   integer*4   :: i,imin,imax,n
   
   n_keywords = 0
   
   if (.not.isempty(para%options)) then
   
      n = len(trim(para%options))
      imin = 1
      do i = 1,n
         if ((para%options(i:i)==',').or.(i==n)) then
            if (para%options(i:i)==',') then
               imax = i-1
            else
               imax = i
            end if
            if (imax<imin) call error('cannot interpret parameter "options"')
            n_keywords = n_keywords+1
            keyword_name(n_keywords) = para%options(imin:imax)
            imin = imax+2
         end if
      end do
      
   end if
   
end subroutine get_keywords

logical function keyword(string)

   implicit none
   character(*),intent(in) :: string
   integer*4               :: i
   
   keyword = .false.
   do i = 1,n_keywords
      if (trim(keyword_name(i))==trim(string)) then
         keyword = .true.
         exit
      end if
   end do
   
end function keyword
             
end module module_interface