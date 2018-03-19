module module_cone_intrinsic

   use module_constants
   use module_system
   use module_cosmology
   use module_user
   use module_parameters
   use module_geometry
   
   private
   public   :: make_cone_intrinsic
   
   type(type_box),allocatable          :: box(:)
   type(type_snapshot),allocatable     :: snapshot(:)
   type(type_galaxy_base),allocatable  :: base(:)
   type(type_galaxy_sam),allocatable   :: sam(:)
   integer*8                           :: nmockgalaxies
   character(len=255)                  :: filename_cone_intrinsic
   real*4                              :: a1(3),a2(3) ! two orthonormal basis vectors, orthogonal to the cone axis
   
contains

subroutine make_cone_intrinsic

   implicit none
   integer*4            :: isnapshot,isubsnapshot,ibox
   character(len=100)   :: snapshotname
   integer*8            :: bytes
   
   call tic
   call out('MAKE INTRINSIC CONE')
   
   ! load previous steps
   call load_parameters
   call load_geometry(box)
   
   ! make snapshot properties
   call load_redshifts(snapshot)
   call assign_distance_ranges
   
   ! create new file
   filename_cone_intrinsic = trim(para%path_output)//'cone_intrinsic.bin'
   open(1,file=trim(filename_cone_intrinsic),action='write',form="unformatted",status='replace',access='stream')
   close(1)
   
   ! fill galaxies with intrinsic properties into cone
   if (para%square_base==1) call make_cone_basis(a1,a2)
   nmockgalaxies = 0
   do isnapshot = para%snapshot_min,para%snapshot_max
      if ((snapshot(isnapshot)%dmax>=minval(box%dmin)).and.(snapshot(isnapshot)%dmin<=maxval(box%dmax))) then
         do isubsnapshot = para%subsnapshot_min,para%subsnapshot_max
            call load_sam_snapshot(isnapshot,isubsnapshot,sam,snapshotname)
            call out('Process '//trim(snapshotname))
            call initialize_base_properties(isnapshot,isubsnapshot)
            do ibox = 1,size(box)
               if ((snapshot(isnapshot)%dmax>=box(ibox)%dmin).and.(snapshot(isnapshot)%dmin<=box(ibox)%dmax)) then
                  call translate_and_rotate_snapshot(ibox)
                  call write_subsnapshot_into_box(isnapshot)
               end if
            end do
         end do
      end if
   end do
   
   ! deallocate arrays
   if (allocated(box)) deallocate(box)
   if (allocated(snapshot)) deallocate(snapshot)
   if (allocated(base)) deallocate(base)
   if (allocated(sam)) deallocate(sam)
   
   ! write info
   inquire(file=filename_cone_intrinsic, size=bytes)
   open(1,file=trim(para%path_output)//'cone_intrinsic_info.txt',action='write',form="formatted",status='replace')
   write(1,'(A,I10)') 'Number.of.galaxies.in.intrinsic.cone        ',nmockgalaxies
   write(1,'(A,I10)') 'Number.of.bytes.per.galaxy.in.intrinsic.cone',bytes/nmockgalaxies
   close(1)
   
   ! finalize output
   call out('Number of galaxies in intrinsic cone:',nmockgalaxies)
   call toc
   
end subroutine make_cone_intrinsic

subroutine translate_and_rotate_snapshot(ibox)

   implicit none
   integer*4,intent(in) :: ibox
   integer*4            :: i
   real*4,allocatable   :: tmp(:)
   
   ! save box-properties
   base%box = ibox
   base%rotation = box(ibox)%rotation
   
   ! copy xbox -> xcone and normalize to box side-length
   do i = 1,3
      base%xcone(i) = base%xbox(i)/para%L
   end do
   
   ! inversion
   if (box(ibox)%rotation<0) then
      do i = 1,3
         base%xcone(i) = 1.0-base%xcone(i)
      end do
   end if
   
   ! proper rotation (case 1 = identity)
   select case(abs(box(ibox)%rotation))
      case(2) ! invert x-axis, while permuting y and z
         base%xcone(1) = 1.0-base%xcone(1)
         tmp = base%xcone(2)
         base%xcone(2) = base%xcone(3)
         base%xcone(3) = tmp
      case(3) ! invert y-axis, while permuting z and x
         base%xcone(2) = 1.0-base%xcone(2)
         tmp = base%xcone(3)
         base%xcone(3) = base%xcone(1)
         base%xcone(1) = tmp
      case(4) ! invert z-axis, while permuting x and y
         base%xcone(3) = 1.0-base%xcone(3)
         tmp = base%xcone(1)
         base%xcone(1) = base%xcone(2)
         base%xcone(2) = tmp
      case(5) ! permute (x,y,z) -> (y,z,x)
         tmp = base%xcone(1)
         base%xcone(1) = base%xcone(2)
         base%xcone(2) = base%xcone(3)
         base%xcone(3) = tmp
      case(6) ! permute (x,y,z) -> (z,x,y)
         tmp = base%xcone(1)
         base%xcone(1) = base%xcone(3)
         base%xcone(3) = base%xcone(2)
         base%xcone(2) = tmp
   end select
      
   ! periodic translation
   do i = 1,3
      base%xcone(i) = modulo(base%xcone(i)+box(ibox)%translation(i),1.0)
   end do
   
   ! translate to apparent position
   do i = 1,3
      base%xcone(i) = base%xcone(i)+box(ibox)%ix(i)-0.5
   end do

end subroutine translate_and_rotate_snapshot

subroutine write_subsnapshot_into_box(isnapshot)

   implicit none
   integer*4,intent(in)                :: isnapshot
   integer*4                           :: igalaxy
   real*4                              :: d
   
   ! open file
   open(1,file=trim(filename_cone_intrinsic),action='write',form='unformatted',status='old',position='append',access='stream')
   
   ! write mock galaxies
   do igalaxy = 1,size(base)
      d = sqrt(sum(base(igalaxy)%xcone**2))
      if ((d>=snapshot(isnapshot)%dmin).and.(d<snapshot(isnapshot)%dmax)) then
         if (geometry_selection(base(igalaxy))) then
            if (pre_selection(sam(igalaxy))) then
               nmockgalaxies = nmockgalaxies+1
               write(1) base(igalaxy),sam(igalaxy)
            end if
         end if
      end if
   end do
   
   ! close file
   close(1)
   
   contains
   
   logical function geometry_selection(base)
   
      implicit none
      type(type_galaxy_base),intent(in)   :: base
      real*4                              :: d,dc,alpha,a(3),alpha1,alpha2
      
      ! check distance
      d = sqrt(real(sum(base%xcone**2)))  ! distance to the galaxy in box side lengths
      dc = d*para%L                       ! distance in simulation units
      if ((dc>para%dc_max).or.(dc<para%dc_min)) then
         geometry_selection = .false.
         return
      end if
      
      ! check cone
      if (para%square_base==1) then
         a = base%xcone
         alpha1 = abs(atan2(sum(a*a1),sum(a*para%axis)))
         alpha2 = abs(atan2(sum(a*a2),sum(a*para%axis)))
         alpha = max(alpha1,alpha2)
      else
         alpha = acos(sum(para%axis*base%xcone)/d) ! angle between cone axis and line-of-sight to the galaxy centre
      end if
      geometry_selection = alpha<=para%angle
      
   end function geometry_selection

end subroutine write_subsnapshot_into_box

subroutine assign_distance_ranges

   implicit none
   integer*4            :: i
   real*4,allocatable   :: d(:)
   
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

subroutine initialize_base_properties(index,subindex)

   implicit none
   integer*4,intent(in)    :: index,subindex
   integer*8   :: i,n
   n = size(sam)
   if (allocated(base)) deallocate(base)
   allocate(base(n))
   do i = 1,n
      base(i) = extract_base(sam(i))
   end do
   base%snapshot = index
   base%subsnapshot = subindex
   
end subroutine initialize_base_properties
   
end module module_cone_intrinsic