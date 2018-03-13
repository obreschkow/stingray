module module_developer

   use module_user
   
   type(type_galaxy_sam),allocatable   :: fake(:)
   
contains

subroutine make_fake_data(ngalaxies)

   implicit none
   integer*4,intent(in)    :: ngalaxies
   integer*4               :: snapshot
   integer*4               :: i

   if (allocated(fake)) deallocate(fake)
   allocate(fake(ngalaxies))
   
   do snapshot = para%snapshot_min,para%snapshot_max
   
      call set_seed(snapshot)
      do i = 1,ngalaxies
         call random_number(fake(i)%position)
         fake(i)%position = fake(i)%position*para%L
         fake(i)%id = i+int(1e5)*snapshot
         fake(i)%haloid = fake(i)%haloid
      end do
      call save_snapshot(snapshot)
      
   end do
   
   call save_redshifts
   
end subroutine make_fake_data

subroutine save_redshifts

   implicit none
   integer*4   :: isnapshot

   open(1,file=trim(para%path_input)//trim(para%file_redshifts),action='write',form='formatted',status='replace')
   do isnapshot = para%snapshot_min,para%snapshot_max
      write(1,*) isnapshot,(para%snapshot_max-isnapshot)*0.2
   end do
   close(1)
    
end subroutine save_redshifts

subroutine save_snapshot(index)

   ! variable declaration
   implicit none
   integer*4,intent(in) :: index
   character(len=255)   :: fn,txt
   integer*8            :: i,n
   
   call tic
   
   ! write user info
   fn = snapshot_filename(index)
   call out('SAVE SNAPSHOT '//trim(fn))
   
   ! write header
   open(1,file=trim(fn),action='write',form='unformatted',status='replace',access='stream')
   
   ! write IDs and positions
   n = size(fake)
   do i = 1,n
      write(1) fake(i)
   end do
   
   ! read IDs
   close(1)
   
   ! output basic statistics
   call out('Number of galaxies:',n)
   write(txt,'(A,F0.4,A,F0.4)') 'Position range: ',minval((/minval(fake(:)%position(1)),minval(fake(:)%position(2)), &
   & minval(fake(:)%position(3))/)) &
   &,' to ',maxval((/maxval(fake(:)%position(1)),maxval(fake(:)%position(2)),maxval(fake(:)%position(3))/))
   call out(txt)
   write(txt,'(A,I0,A,I0)')     'Identifier range: ',minval(fake%id),' to ',maxval(fake%id)
   call out(txt)
   
   call toc
   
end subroutine save_snapshot
   
end module module_developer