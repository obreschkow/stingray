! Module converts intrinsic into apparent properties

module module_sky

   use module_constants
   use module_system
   use module_types
   use module_io
   use module_cosmology
   use module_conversion
   use module_user
   use module_parameters
   use module_tiling
   
   private
   public   :: make_sky
   
contains

subroutine make_sky

   implicit none
   character(len=255)   :: filename
   integer*8            :: n,i,n100
   integer*8,allocatable:: m(:)
   integer*4            :: isky
   type(type_base)      :: base
   type(type_sam)       :: sam   ! intrinsic galaxy properties from SAM
   character(len=3)     :: str
   
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
      write(filename,'(A,A,A,A)') trim(para%path_output),'mocksky_',trim(skyclass(isky)%ptr%get_class_name()),'.bin'
      open(isky+1,file=trim(filename),action='write',form="unformatted",status='replace',access='stream')
      write(isky+1) 0_8 ! place holder for number of objects in mock sky
   end do
   
   ! convert galaxy properties and write master sky
   allocate(m(size(skyclass)))
   m = 0
   n100 = n/100
   do i = 1,n
      if (modulo(i-1,n100)==0) then
         write(str,'(I3)') nint(real(i-1,8)/n*100)
         call out('Progress: '//str//'%')
      end if
      read(1) base,sam
      Rvector = tile(base%tile)%Rvector
      Rpseudo = tile(base%tile)%Rpseudo
      call sam%rotate_vectors
      do isky = 1,size(skyclass)
         call skyclass(isky)%ptr%make_from_sam(sam,base,m(isky)+1,sum(m)+1)
         if (skyclass(isky)%ptr%is_selected(sam)) then
            m(isky) = m(isky)+1
            call skyclass(isky)%ptr%write_to_file(isky+1)
         end if
      end do
   end do
   call out('Progress: 100%')
   
   ! add number of objects to beginning of file & close files
   do isky = 1,size(skyclass)
      write(isky+1,pos=1) m(isky)
      close(isky+1)
   end do
   close(1)
   
   ! check number of objects
   !do isky = 1,size(skyclass)
      if (sum(m)==0) call error('No objects in the apparent sky. Consider changing selection function.')
   !end do
   
   ! user output
   call out('Number of objects in apparent sky:',sum(m))
   do isky = 1,size(skyclass)
      call out('  '//trim(skyclass(isky)%ptr%get_class_name())//':',m(isky))
   end do
   call toc
   
end subroutine make_sky

end module module_sky