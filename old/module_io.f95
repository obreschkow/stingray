module module_io

   use module_types
   use module_parameters
   use module_system

contains

subroutine load_redshifts

   implicit none
   integer*4          :: sn,snapshot

   allocate(redshift(para%snapshot_min:para%snapshot_max))

   open(1,file=trim(para%file_redshifts),form="formatted")
   do snapshot = para%snapshot_min,para%snapshot_max
      read(1,*) sn,redshift(snapshot)
      if (snapshot.ne.sn) then
         call out('ERROR: snapshot_min and snapshot_max inconsistent with snapshots in redshift file.')
         close(1)
         stop
      end if
   end do
   close(1)
    
end subroutine load_redshifts

subroutine load_snapshot(index)

   ! variable declaration
   implicit none
   integer*4,intent(in) :: index
   character(len=255)   :: fn,txt
   integer*8            :: i,n
   
   call tic
   
   ! write user info
   write(fn,'(A,'//trim(para%file_snapshot_extension)//')') para%file_snapshot_base,index
   call out('LOAD SNAPSHOT '//trim(fn))
   
   ! read header
   open(1,file=trim(fn),action='read',form='unformatted',status='old',access='stream')
   read(1) n
   call out('Number of galaxies:',n)
   
   ! allocate
   allocate(galaxy(n))
   
   ! read IDs and positions
   read(1) (galaxy(i),i=1,n)
   
   ! read IDs
   close(1)
   
   ! output basic statistics
   write(txt,'(A,F0.4,A,F0.4)') 'Position range: ',minval(galaxy(:)%x(1)),' to ',maxval(galaxy(:)%x(1))
   call out(txt)
   write(txt,'(A,I0,A,I0)')     'Identifier range: ',minval(galaxy%id),' to ',maxval(galaxy%id)
   call out(txt)
   
   call toc
   
end subroutine load_snapshot

subroutine save_snapshot(index)

   ! variable declaration
   implicit none
   integer*4,intent(in) :: index
   character(len=255)   :: fn,txt
   integer*8            :: i,n
   
   call tic
   
   ! write user info
   write(fn,'(A,'//trim(para%file_snapshot_extension)//')') trim(para%file_snapshot_base),index
   call out('SAVE SNAPSHOT '//trim(fn))
   
   ! write header
   open(1,file=trim(fn),action='write',form='unformatted',status='replace',access='stream')
   n = size(galaxy)
   write(1) n
   call out('Number of galaxies:',n)
   
   ! write IDs and positions
   write(1) (galaxy(i),i=1,n)
   
   ! read IDs
   close(1)
   
   ! output basic statistics
   write(txt,'(A,F0.4,A,F0.4)') 'Position range: ',minval((/minval(galaxy(:)%x(1)),minval(galaxy(:)%x(2)), &
   & minval(galaxy(:)%x(3))/)) &
   &,' to ',maxval((/maxval(galaxy(:)%x(1)),maxval(galaxy(:)%x(2)),maxval(galaxy(:)%x(3))/))
   call out(txt)
   write(txt,'(A,I0,A,I0)')     'Identifier range: ',minval(galaxy%id),' to ',maxval(galaxy%id)
   call out(txt)
   
   call toc
   
end subroutine save_snapshot

end module module_io