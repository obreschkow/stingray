module module_developer

   use module_user
   
   type(type_galaxy_sam),allocatable   :: galaxy(:)
   
contains

subroutine make_fake_data(ngalaxies)

   implicit none
   integer*4,intent(in)    :: ngalaxies
   integer*4               :: snapshot
   integer*4               :: i
   
   call tic
   call out('MAKE FAKE DATA')

   if (allocated(galaxy)) deallocate(galaxy)
   allocate(galaxy(ngalaxies))
   
   do snapshot = para%snapshot_min,para%snapshot_max
   
      call set_seed(snapshot)
      do i = 1,ngalaxies
         call random_number(galaxy(i)%position)
         galaxy(i)%id = i+int(1e5)*snapshot
         galaxy(i)%haloid = galaxy(i)%haloid
         galaxy(i)%position = galaxy(i)%position*para%L
         galaxy(i)%mag  = get_normal_random_number(-23.0,1.0)
         galaxy(i)%MHI  = 10.0**(get_normal_random_number(9.0,1.0))
         galaxy(i)%v(1) = get_normal_random_number(0.0,1e2)
         galaxy(i)%v(2) = get_normal_random_number(0.0,1e2)
         galaxy(i)%v(3) = get_normal_random_number(0.0,1e2)
         galaxy(i)%j(1) = get_normal_random_number(0.0,1e3)
         galaxy(i)%j(2) = get_normal_random_number(0.0,1e3)
         galaxy(i)%j(3) = get_normal_random_number(0.0,1e3)
      end do
      call save_snapshot(snapshot)
      
   end do
   
   call save_redshifts
   
   call toc
   
   contains
   
   subroutine save_redshifts

      implicit none
      integer*4   :: isnapshot
   
      call out('Save redshifts')

      open(1,file=trim(para%path_input)//'redshifts.txt',action='write',form='formatted',status='replace')
      do isnapshot = para%snapshot_min,para%snapshot_max
         write(1,*) isnapshot,(para%snapshot_max-isnapshot)*0.1
      end do
      close(1)
    
   end subroutine save_redshifts

   subroutine save_snapshot(index)

      ! variable declaration
      implicit none
      integer*4,intent(in) :: index
      character(len=255)   :: fn,txt
      integer*8            :: i,n
   
      ! write user info
      write(fn,'(A,A,I0.3)') trim(para%path_input),'snapshot_',index
      call out('Save snapshot '//trim(fn))
   
      ! write header
      open(1,file=trim(fn),action='write',form='unformatted',status='replace',access='stream')
   
      ! write IDs and positions
      n = size(galaxy)
      do i = 1,n
         write(1) galaxy(i)
      end do
   
      ! read IDs
      close(1)
   
      ! output basic statistics
      call out('Number of galaxies:',n)
      write(txt,'(A,F0.4,A,F0.4)') 'Position range: ',minval((/minval(galaxy(:)%position(1)),minval(galaxy(:)%position(2)), &
      & minval(galaxy(:)%position(3))/)) &
      &,' to ',maxval((/maxval(galaxy(:)%position(1)),maxval(galaxy(:)%position(2)),maxval(galaxy(:)%position(3))/))
      call out(txt)
      write(txt,'(A,I0,A,I0)')     'Identifier range: ',minval(galaxy%id),' to ',maxval(galaxy%id)
      call out(txt)
   
   end subroutine save_snapshot
   
end subroutine make_fake_data

end module module_developer