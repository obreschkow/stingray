module module_sky_apparent

   use module_constants
   use module_types
   use module_system
   use module_cosmology
   use module_user
   use module_parameters
   use module_tiling
   
   private
   public   :: make_sky_apparent
   
contains

subroutine make_sky_apparent

   implicit none
   character(len=255)   :: filename
   integer*8            :: n,i,m
   integer*4            :: sky_sel
   type(type_base)      :: base
   type(type_sam)       :: sam   ! intrinsic galaxy properties from SAM
   type(type_sky)       :: sky  ! apparent galaxy properties
   
   ! write user info
   call tic
   call out('CONVERT INTRINSIC TO APPARENT PROPERTIES IN MOCK SKY')
   
   ! load previous steps
   call load_parameters
   call load_box_list
   
   ! open intrinsic sky
   filename = trim(para%path_output)//'mocksky_intrinsic.bin'
   open(2,file=trim(filename),action='read',form='unformatted',status='old',access='stream')
   read(2) n ! number of galaxies
   
   ! write user info
   call out('Number of galaxies in intrinsic sky:',n)
   
   ! initialize master one
   filename = trim(para%path_output)//'mocksky.bin'
   open(1,file=trim(filename),action='write',form="unformatted",status='replace',access='stream')
   write(1) 0_8 ! place holder for number of objects in mock sky
   
   ! convert galaxy properties and write master sky
   m = 0
   do i = 1,n
      read(2) base,sam
      Rvector = tile(base%tile)%Rvector
      Rpseudo = tile(base%tile)%Rpseudo
      call rotate_vectors(sam)
      sky = convert_properties(sam,m+1,base%dc,base%ra,base%dec,base%tile)
      sky_sel = sky_selection(sky)
      if (sky_sel>0) then
         base%sky_selection = sky_sel
         m = m+1
         write(1) base,sky
      end if
   end do
   
   ! add number of objects to beginning of file
   write(1,pos=1) m
   
   ! close files
   close(1)
   close(2)
   
   if (m==0) call error('No objects in the apparent sky. Consider changing selection function.')
   
   ! user output
   call out('Number of objects in apparent sky:',m)
   call toc
   
end subroutine make_sky_apparent

end module module_sky_apparent