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
   do isky = 1,size(skyclass)
      write(filename,'(A,A,A,A)') trim(para%path_output),'mocksky_',trim(skyclass(isky)%ptr%name()),'.bin'
      open(isky+1,file=trim(filename),action='write',form="unformatted",status='replace',access='stream')
      write(isky+1) 0_8 ! place holder for number of objects in mock sky
   end do
   
   ! convert galaxy properties and write master sky
   allocate(m(size(skyclass)))
   m = 0
   do i = 1,n
      read(1) base,sam
      Rvector = tile(base%tile)%Rvector
      Rpseudo = tile(base%tile)%Rpseudo
      call rotate_vectors(sam)
      do isky = 1,size(skyclass)
         call skyclass(isky)%ptr%convertSam(sam,m(isky)+1,sum(m)+1,base%dc*para%L,base%ra,base%dec,base%tile)
         if (skyclass(isky)%ptr%selected(sam)) then
            m(isky) = m(isky)+1
            call skyclass(isky)%ptr%writeToFile(isky+1)
         end if
      end do
   end do
   
   ! add number of objects to beginning of file & close files
   do isky = 1,size(skyclass)
      write(isky+1,pos=1) m(isky)
      close(isky+1)
   end do
   close(1)
   
   ! check number of objects
   do isky = 1,size(skyclass)
      if (m(isky)==0) call error('No objects in the apparent sky. Consider changing selection function.')
   end do
   
   ! user output
   call out('Number of objects in apparent sky:',sum(m))
   do isky = 1,size(skyclass)
      call out('  '//trim(skyclass(isky)%ptr%name())//':',m(isky))
   end do
   call toc
   
end subroutine make_sky_apparent

end module module_sky_apparent