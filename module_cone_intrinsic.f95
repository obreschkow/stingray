module module_cone_intrinsic

use module_constants
use module_parameters
use module_system
use module_cosmology
use module_user
use module_geometry

type(type_galaxy_base),allocatable  :: base(:)
type(type_galaxy_sam),allocatable   :: sam(:)
type(type_snapshot),allocatable     :: snapshot(:)
integer*8                           :: nmockgalaxies
character(len=255)                  :: filename_cone_intrinsic
   
contains

subroutine make_cone_intrinsic

   implicit none
   integer*4            :: isnapshot,ibox
   
   call tic
   call out('MAKE INTRINSIC CONE')
   
   call load_geometry
   call load_redshifts
   call assign_distance_ranges
   
   ! create new file
   filename_cone_intrinsic = trim(para%path_output)//'cone_intrinsic.bin'
   open(1,file=trim(filename_cone_intrinsic),action='write',form="unformatted",status='replace',access='stream')
   close(1)
   
   ! fill galaxies with intrinsic properties into cone
   nmockgalaxies = 0
   do isnapshot = para%snapshot_min,para%snapshot_max
      if ((snapshot(isnapshot)%dmax>=minval(box%dmin)).and.(snapshot(isnapshot)%dmin<=maxval(box%dmax))) then
         call out('Process snapshot ',isnapshot*1_8)
         call load_sam_snapshot(isnapshot,sam)
         call initialize_base_properties(isnapshot)
         do ibox = 1,size(box)
            if ((snapshot(isnapshot)%dmax>=box(ibox)%dmin).and.(snapshot(isnapshot)%dmin<=box(ibox)%dmax)) then
               call translate_and_rotate_snapshot(ibox)
               call distribute_snapshot_into_box(isnapshot)
            end if
         end do
      end if
   end do
   
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

subroutine distribute_snapshot_into_box(isnapshot)

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
      real*4                              :: d,dmax,dc,alpha
      
      dmax = para%dc_max/para%L
      
      if ((base%xcone(1)*dmax>para%x_max*base%xcone(3)).or.(base%xcone(1)*dmax<para%x_min*base%xcone(3))) then
         geometry_selection = .false.
         return
      end if
      
      if ((base%xcone(2)*dmax>para%y_max*base%xcone(3)).or.(base%xcone(2)*dmax<para%y_min*base%xcone(3))) then
         geometry_selection = .false.
         return
      end if
      
      d = sqrt(real(sum(base%xcone**2))) ! distance to the galaxy in box side lengths
      dc = d*para%L
      if ((dc>para%dc_max).or.(dc<para%dc_min)) then
         geometry_selection = .false.
         return
      end if
      
      alpha = acos(sum(para%axis*base%xcone)/d) ! angle between cone axis and line-of-sight to the galaxy centre
      geometry_selection = alpha<=para%angle
      
   end function geometry_selection

end subroutine distribute_snapshot_into_box

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

subroutine load_redshifts

   implicit none
   integer*4   :: isnapshot,i
   real*4      :: redshift

   allocate(snapshot(para%snapshot_min:para%snapshot_max))

   open(1,file=trim(para%path_input)//trim(para%file_redshifts),form="formatted")
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

subroutine initialize_base_properties(index)
   implicit none
   integer*4,intent(in)    :: index
   integer*8   :: i,n
   n = size(sam)
   if (allocated(base)) deallocate(base)
   allocate(base(n))
   do i = 1,n
      base(i) = extract_base(sam(i))
   end do
   base%snapshot = index
end subroutine initialize_base_properties
   
end module module_cone_intrinsic