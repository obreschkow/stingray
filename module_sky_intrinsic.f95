module module_sky_intrinsic

   use module_constants
   use module_system
   use module_types
   use module_cosmology
   use module_user
   use module_parameters
   use module_tiling
   
   private
   public   :: make_sky_intrinsic
   
   type(type_base),allocatable   :: base(:)
   type(type_sam),allocatable    :: sam(:)
   integer*8                     :: nmockgalaxies
   character(len=255)            :: filename_sky_intrinsic
   
contains

subroutine make_sky_intrinsic

   implicit none
   integer*4            :: isnapshot,isubsnapshot,itile
   character(len=100)   :: snapshotname
   
   call tic
   call out('MAKE INTRINSIC SKY')
   
   ! load previous steps
   call load_parameters
   call load_box_list
   
   ! make snapshot properties
   if (allocated(snapshot)) deallocate(snapshot)
   allocate(snapshot(para%snapshot_min:para%snapshot_max))
   call make_redshifts
   call make_distance_ranges
   call save_snapshot_list
   
   ! create new file
   filename_sky_intrinsic = trim(para%path_output)//'mocksky_intrinsic.bin'
   open(1,file=trim(filename_sky_intrinsic),action='write',form="unformatted",status='replace',access='stream')
   write(1) 0_8 ! place holder for number of galaxies
   close(1)
   
   ! fill galaxies with intrinsic properties into sky
   nmockgalaxies = 0
   do isnapshot = para%snapshot_min,para%snapshot_max
      if ((snapshot(isnapshot)%dmax>=minval(tile%dmin)).and.(snapshot(isnapshot)%dmin<=maxval(tile%dmax))) then
         do isubsnapshot = para%subsnapshot_min,para%subsnapshot_max
            call load_sam_snapshot(isnapshot,isubsnapshot,sam,snapshotname)
            call out('Process '//trim(snapshotname))
            allocate(base(size(sam)))
            do itile = 1,size(tile)
               if ((snapshot(isnapshot)%dmax>=tile(itile)%dmin).and.(snapshot(isnapshot)%dmin<=tile(itile)%dmax)) then
                  call convert_position_sam_to_sky(itile)
                  call write_subsnapshot_into_tile(isnapshot)
               end if
            end do
         end do
      end if
   end do
   
   ! deallocate arrays
   if (allocated(tile)) deallocate(tile)
   if (allocated(snapshot)) deallocate(snapshot)
   if (allocated(base)) deallocate(base)
   if (allocated(sam)) deallocate(sam)
   
   ! check number of galaxies
   if (nmockgalaxies==0) then
      call error('No galaxy in sky. Consider widening the sky geometry or relaxing the selection criteria.')
   end if
   
   ! add number of objects to beginning of file
   open(1,file=trim(filename_sky_intrinsic),action='write',form='unformatted',status='old',access='stream')
   write(1,pos=1) nmockgalaxies
   close(1) 
   
   ! finalize output
   call out('Number of galaxies in intrinsic sky:',nmockgalaxies)
   call toc
   
end subroutine make_sky_intrinsic

subroutine convert_position_sam_to_sky(itile)

   implicit none
   integer*4,intent(in) :: itile
   integer*4            :: i
   real*4               :: x(3),d
   
   ! save box index
   base%tile = itile
   
   do i = 1,size(base)
      
      x = sam(i)%getPosition()/para%L ! change length units to units of box side-length
      x = matmul(rot(:,:,tile(itile)%rotation),x) ! random 90-degree rotation/inversion
      x = modulo(x+tile(itile)%translation,1.0) ! periodic translation
      x = x+tile(itile)%ix-0.5 ! translate coordinates to tile position
      x = matmul(para%sky_rotation,x) ! convert SAM-coordinates to Sky-coordinates
      call car2sph(x,d,base(i)%ra,base(i)%dec)
      base(i)%dc = d*para%L
      
   end do

end subroutine convert_position_sam_to_sky

subroutine write_subsnapshot_into_tile(isnapshot)

   implicit none
   integer*4,intent(in) :: isnapshot
   integer*4            :: i,ps,is
   
   ! open file
   open(1,file=trim(filename_sky_intrinsic),action='write',form='unformatted',status='old',position='append',access='stream')
   
   ! write mock galaxies
   do i = 1,size(base)
      
      ! check distance relative to snapshot
      if ((base(i)%dc>=snapshot(isnapshot)%dmin*para%L).and.(base(i)%dc<snapshot(isnapshot)%dmax*para%L)) then
         
         ! check full position-selection
         if (is_in_fov(base(i)%dc,base(i)%ra,base(i)%dec)) then
            ps = position_selection(base(i)%dc,base(i)%ra/degree,base(i)%dec/degree)
            if (ps>0) then
               !base(i)%position_selection = ps
                 
               ! check intrinsic property-selection
               is = intrinsic_selection(sam(i))
               if (is>0) then
                  !base(i)%intrinsic_selection = is
                  
                  ! write selected galaxy into intrinsic sky file
                  nmockgalaxies = nmockgalaxies+1
                  write(1) base(i),sam(i)
                  
               end if
            end if
         end if
      end if
   end do
   
   ! close file
   close(1)

end subroutine write_subsnapshot_into_tile

subroutine make_distance_ranges

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

end subroutine make_distance_ranges
   
end module module_sky_intrinsic