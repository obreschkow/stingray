module module_galaxies

   use module_parameters
   use module_system
   use module_cosmology
   use module_geometry
   
   type(type_galaxy),allocatable       :: galaxy(:)
   type(type_properties)
   type(type_snapshot),allocatable     :: snapshot(:)
   integer*8                           :: nmockgalaxies
   
contains

subroutine make_mock_galaxies

   implicit none
   integer*4            :: isnapshot,ibox,igalaxy
   logical              :: newsnapshot
   real*4,allocatable   :: x(:,:)
   
   call load_geometry
   call load_redshifts
   call assign_distance_ranges
   
   nmockgalaxies = 0
   do isnapshot = para%snapshot_min,para%snapshot_max
      if ((snapshot(isnapshot)%dmax>=minval(box%dmin)).and.(snapshot(isnapshot)%dmin<=maxval(box%dmax))) then
         call load_snapshot(isnapshot)
         newsnapshot = .true.
         do ibox = 1,size(box)
            if ((snapshot(isnapshot)%dmax>=box(ibox)%dmin).and.(snapshot(isnapshot)%dmin<=box(ibox)%dmax)) then
               if (newsnapshot) then
                  if (allocated(x)) deallocate(x)
                  allocate(x(size(galaxy),3))
                  do igalaxy = 1,size(galaxy)
                     x(igalaxy,:) = galaxy(igalaxy)%x
                  end do
                  newsnapshot = .false.
               else
                  do igalaxy = 1,size(galaxy)
                     galaxy(igalaxy)%x = x(igalaxy,:)
                  end do
               end if
               call translate_and_rotate_snapshot(ibox)
               call distribute_snapshot_into_box(isnapshot,ibox)
            end if
         end do
      end if
   end do
   
end subroutine make_mock_galaxies

subroutine translate_and_rotate_snapshot(ibox)

   implicit none
   integer*4,intent(in) :: ibox
   integer*4            :: i
   real*4,allocatable   :: tmp(:)
   
   ! normalization
   do i = 1,3
      galaxy%x(i) = galaxy%x(i)/para%L
   end do
   
   ! inversion
   if (box(ibox)%rotation<0) then
      do i = 1,3
         galaxy%x(i) = -galaxy%x(i)
      end do
   end if
   
   ! proper rotation (case 1 = identity)
   select case(abs(box(ibox)%rotation))
      case(2) ! invert x-axis, while permuting y and z
         galaxy%x(1) = -galaxy%x(1)
         tmp = galaxy%x(2)
         galaxy%x(2) = galaxy%x(3)
         galaxy%x(3) = tmp
      case(3) ! invert y-axis, while permuting z and x
         galaxy%x(2) = -galaxy%x(2)
         tmp = galaxy%x(3)
         galaxy%x(3) = galaxy%x(1)
         galaxy%x(1) = tmp
      case(4) ! invert z-axis, while permuting x and y
         galaxy%x(3) = -galaxy%x(3)
         tmp = galaxy%x(1)
         galaxy%x(1) = galaxy%x(2)
         galaxy%x(2) = tmp
      case(5) ! permute (x,y,z) -> (y,z,x)
         tmp = galaxy%x(1)
         galaxy%x(1) = galaxy%x(2)
         galaxy%x(2) = galaxy%x(3)
         galaxy%x(3) = tmp
      case(6) ! permute (x,y,z) -> (z,x,y)
         tmp = galaxy%x(1)
         galaxy%x(1) = galaxy%x(3)
         galaxy%x(3) = galaxy%x(2)
         galaxy%x(2) = tmp
   end select
      
   ! translation
   do i = 1,3
      galaxy%x(i) = modulo(galaxy%x(i)+box(ibox)%translation(i),1.0)+box(ibox)%ix(i)
   end do

end subroutine translate_and_rotate_snapshot

subroutine distribute_snapshot_into_box(isnapshot,ibox)

   implicit none
   integer*4,intent(in)                :: isnapshot
   integer*4,intent(in)                :: ibox
   integer*4                           :: igalaxy
   integer*4                           :: imockgalaxy
   real*4                              :: d
   type(type_mockgalaxy),allocatable   :: mockgalaxy(:)
   
   ! initialize
   allocate(mockgalaxy(size(galaxy)))
   imockgalaxy = 0
   
   ! make mock galaxies
   do igalaxy = 1,size(galaxy)
      d = sqrt(sum(galaxy(igalaxy)%x**2))
      if ((d>=snapshot(isnapshot)%dmin).and.(d<snapshot(isnapshot)%dmax)) then
         imockgalaxy = imockgalaxy+1
         mockgalaxy(imockgalaxy)%id = galaxy(igalaxy)%id
         mockgalaxy(imockgalaxy)%groupid = galaxy(igalaxy)%groupid
         mockgalaxy(imockgalaxy)%x = galaxy(igalaxy)%x
         mockgalaxy(imockgalaxy)%box = ibox
         mockgalaxy(imockgalaxy)%snapshot = isnapshot
         mockgalaxy(imockgalaxy)%dc = d*para%L
         mockgalaxy(imockgalaxy)%rotation = box(ibox)%rotation
      end if
   end do
   
   ! write mock galaxies into stream
   ! ...
   
   ! finalize
   deallocate(mockgalaxy)
   nmockgalaxies = nmockgalaxies+imockgalaxy

end subroutine distribute_snapshot_into_box

subroutine assign_distance_ranges

   implicit none
   integer*4            :: i
   real*4,allocatable   :: d(:)
   real*4,parameter     :: Mpc = 3.0856776e+22 ! unit conversion factor Mpc/m
   
   allocate(d(para%snapshot_min:para%snapshot_max))
   do i = para%snapshot_min,para%snapshot_max
      d(i) = redshift_to_dc(snapshot(i)%redshift)*(Mpc/para%length_unit)/para%L ! [box side-length] comoving distance to redshift of the box
   end do
   
   do i = para%snapshot_min,para%snapshot_max
      if (i==para%snapshot_max) then
         snapshot(i)%dmin = 0
      else
         snapshot(i)%dmin = 0.5*(d(i+1)+d(i))
      end if
      if (i==para%snapshot_min) then
         snapshot(i)%dmax = 1e20
      else
         snapshot(i)%dmax = 0.5*(d(i)+d(i-1))
      end if
   end do

end subroutine assign_distance_ranges

subroutine load_redshifts

   implicit none
   integer*4   :: isnapshot,i
   real*4      :: redshift

   allocate(snapshot(para%snapshot_min:para%snapshot_max))

   open(1,file=trim(para%file_redshifts),form="formatted")
   do isnapshot = para%snapshot_min,para%snapshot_max
      read(1,*) i,redshift
      if (isnapshot.ne.i) then
         call out('ERROR: snapshot_min and snapshot_max inconsistent with snapshots in redshift file.')
         close(1)
         stop
      else
         snapshot(isnapshot)%redshift = redshift
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
   write(fn,'(A,'//trim(para%file_snapshot_extension)//')') trim(para%file_snapshot_base),index
   call out('LOAD SNAPSHOT '//trim(fn))
   
   ! read header
   open(1,file=trim(fn),action='read',form='unformatted',status='old')
   read(1) n
   call out('Number of galaxies:',n)
   
   ! allocate
   if (allocated(galaxy)) deallocate(galaxy)
   allocate(galaxy(n))
   
   ! read IDs and positions
   do i = 1,n
      read(1) galaxy(i)
   end do
   
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
   open(1,file=trim(fn),action='write',form='unformatted',status='replace')
   n = size(galaxy)
   write(1) n
   call out('Number of galaxies:',n)
   
   ! write IDs and positions
   do i = 1,n
      write(1) galaxy(i),34,'test'
   end do
   
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
   
end module module_galaxies