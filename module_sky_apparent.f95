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
   integer*8            :: n,i
   integer*8,allocatable:: m(:)
   integer*4            :: isky
   type(type_base)      :: base
   type(type_sam)       :: sam   ! intrinsic galaxy properties from SAM
   
   ! write user info
   call tic
   call out('CONVERT INTRINSIC TO APPARENT PROPERTIES IN MOCK SKY')
   
   ! load previous steps
   call load_parameters
   call load_box_list
   
   ! set random seed (in case the assignment of sky-properties involves randomness)
   call set_seed(para%seed)
   
   ! open intrinsic sky
   filename = trim(para%path_output)//'mocksky_intrinsic.bin'
   open(1,file=trim(filename),action='read',form='unformatted',status='old',access='stream')
   read(1) n ! number of galaxies
   
   ! write user info
   call out('Number of galaxies in intrinsic sky:',n)
   
   ! initialize master one
   do isky = 1,size(ptr)
      write(filename,'(A,A,I02,A)') trim(para%path_output),'mocksky_class',isky,'.bin'
      open(isky+1,file=trim(filename),action='write',form="unformatted",status='replace',access='stream')
      write(isky+1) 0_8 ! place holder for number of objects in mock sky
   end do
   
   ! convert galaxy properties and write master sky
   allocate(m(size(ptr)))
   m = 0
   do i = 1,n
      read(1) base,sam
      Rvector = tile(base%tile)%Rvector
      Rpseudo = tile(base%tile)%Rpseudo
      call rotate_vectors(sam)
      do isky = 1,size(ptr)
         call ptr(isky)%sky%convertSAM(sam,m(isky)+1,sum(m)+1,base%dc*para%L,base%ra,base%dec,base%tile)
         if (ptr(isky)%sky%selected(sam)) then
            m(isky) = m(isky)+1
            call ptr(isky)%sky%writeToFile(isky+1)
         end if
      end do
   end do
   
   ! add number of objects to beginning of file & close files
   do isky = 1,size(ptr)
      write(isky+1,pos=1) m(isky)
      close(isky+1)
   end do
   close(1)
   
   ! check number of objects
   do isky = 1,size(ptr)
      if (m(isky)==0) call error('No objects in the apparent sky. Consider changing selection function.')
   end do
   
   ! user output
   call out('Number of objects in apparent sky:',sum(m))
   call toc
   
end subroutine make_sky_apparent

end module module_sky_apparent