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
   
   ! check number of galaxies
   if (nmockgalaxies==0) then
      call error('No galaxy in cone. Consider widening the cone geometry or relaxing the selection criteria.')
   end if
   
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
   
   ! save box index
   base%box = ibox
   
   do i = 1,size(base)
   
      ! copy xbox -> xcone and normalize to box side-length
      base(i)%xcone = base(i)%xbox/para%L
   
      ! rotation
      base(i)%xcone = matmul(rot(:,:,box(ibox)%rotation),base(i)%xcone)
   
      ! periodic translation
      base(i)%xcone = modulo(base(i)%xcone+box(ibox)%translation,1.0)
   
      ! translate to apparent position
      base(i)%xcone = base(i)%xcone+box(ibox)%ix-0.5
      
   end do

end subroutine translate_and_rotate_snapshot

subroutine write_subsnapshot_into_box(isnapshot)

   implicit none
   integer*4,intent(in)                :: isnapshot
   integer*4                           :: igalaxy
   real*4                              :: d,dc
   
   ! open file
   open(1,file=trim(filename_cone_intrinsic),action='write',form='unformatted',status='old',position='append',access='stream')
   
   ! write mock galaxies
   do igalaxy = 1,size(base)
      d = norm(base(igalaxy)%xcone)
      
      ! check distance
      if ((d>=snapshot(isnapshot)%dmin).and.(d<snapshot(isnapshot)%dmax)) then
         dc = d*para%L ! comoving distance in simulation units
         if ((dc>=para%dc_min).or.(dc<=para%dc_max)) then
         
            ! check footprint on sky (ra,dec)
            base(igalaxy)%xcone = matmul(para%sky_rotation,base(igalaxy)%xcone) ! rotates cone axis onto central RA, Dec
            call make_sky_coordinates(base(igalaxy)%xcone,base(igalaxy)%dc,base(igalaxy)%ra,base(igalaxy)%dec)
            if (check_footprint(base(igalaxy))) then
            
               ! user checks
               if (intrinsic_selection(sam(igalaxy))) then
                  
                  ! write selected galaxy into intrinsic cone file
                  nmockgalaxies = nmockgalaxies+1
                  write(1) base(igalaxy),sam(igalaxy)
                  
               end if
            end if
         end if
      end if
   end do
   
   ! close file
   close(1)
   
   contains
   
   logical function check_footprint(b)
   
      implicit none
      type(type_galaxy_base),intent(in)   :: b
      real*4                              :: d
      
      if (para%dec_min>pi) then
         ! check circular cone
         d = norm(b%xcone)
         check_footprint = acos(min(1.0,sum(para%axis*b%xcone)/d))<=para%angle
      else
         ! check ra and dec range
         check_footprint = (b%ra>=para%ra_min).and.(b%ra<=para%ra_max).and. &
         & (b%dec>=para%dec_min).and.(b%dec<=para%dec_max)
      end if
                  
   end function check_footprint

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