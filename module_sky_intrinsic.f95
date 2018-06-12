module module_sky_intrinsic

   use module_constants
   use module_system
   use module_types
   use module_cosmology
   use module_sort
   use module_user
   use module_parameters
   use module_tiling
   
   private
   public   :: make_sky_intrinsic
   
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
            do itile = 1,size(tile)
               if ((snapshot(isnapshot)%dmax>=tile(itile)%dmin).and.(snapshot(isnapshot)%dmin<=tile(itile)%dmax)) then
                  call write_subsnapshot_into_tile(itile,isnapshot)
               end if
            end do
         end do
      end if
   end do
   
   ! deallocate arrays
   if (allocated(tile)) deallocate(tile)
   if (allocated(snapshot)) deallocate(snapshot)
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

subroutine convert_position_sam_to_sky(position,itile,dc,ra,dec,xbox)

   implicit none
   real*4,intent(in)    :: position(3) ! [box side-length] position in box
   integer*4,intent(in) :: itile
   real*4,intent(out)   :: dc,ra,dec ! [box side-length,rad,rad] position on sky
   real*4,intent(out)   :: xbox(3) ! [box side-length] position in box after symmetry operations
   real*4               :: x(3)
   
   xbox = modulo(matmul(rot(:,:,tile(itile)%rotation),position)+tile(itile)%translation,1.0) ! apply symmetries
   x = xbox+tile(itile)%ix-0.5 ! translate coordinates to tile position
   x = matmul(para%sky_rotation,x) ! convert SAM-coordinates to Sky-coordinates
   call car2sph(x,dc,ra,dec)
   
end subroutine convert_position_sam_to_sky

subroutine write_subsnapshot_into_tile(itile,isnapshot)

   implicit none
   integer*4,intent(in)    :: itile
   integer*4,intent(in)    :: isnapshot
   integer*4               :: i,j,k,l,n
   type(type_base)         :: base
   real*4                  :: xbox(3)
   integer*8,allocatable   :: index(:,:)
   integer*8               :: old_group
   integer*4               :: group_start
   real*4                  :: xmin(3),xmax(3)
   integer*4               :: n_in_group
   integer*4               :: n_in_snapshot
   integer*4               :: n_in_survey_volume
   logical                 :: ok
   
   ! set tile
   base%tile = itile
   
   ! sort by groupID
   n = size(sam)
   allocate(index(n,2))
   do i = 1,n
      index(i,1) = getGroupID(sam(i))
      index(i,2) = i
   end do
   call merge_sort_list(index)
   
   ! open file
   open(1,file=trim(filename_sky_intrinsic),action='write',form='unformatted',status='old',position='append',access='stream')
   
   ! initialise group flagging
   call reset_flagging(1,index(1,1))
   
   ! iterate over all mock galaxies, ordered by group
   do k = 1,n+1
   
      ! assign group flags
      if ((k>n).or.(index(min(k,n),1).ne.old_group)) then
      
         ! make group flag
         base%flag = 0
         if (n_in_survey_volume.ne.n_in_group) base%flag = base%flag+1  ! group truncated by survey volume
         if (n_in_snapshot.ne.n_in_group) base%flag = base%flag+2       ! group truncated by snapshot limit
         if (maxval(xmax-xmin)>0.5) base%flag = base%flag+4              ! group is wrapped around
         
         ! re-iterate over all group members
         do l = group_start,k-1
            i = int(index(l,2),4)
            ! check intrinsic property-selection
            if (i>0) then
               if (sam_selection(sam(i))) then
                  nmockgalaxies = nmockgalaxies+1
                  write(1) base,sam(i)
               end if
            end if
         end do
         
         ! exit do-loop
         if (k>n) exit
         
         ! reset group flagging
         call reset_flagging(k,index(k,1))
            
      end if
      
      i = int(index(k,2),4)
      ok = .true.
      
      ! compute sky position and determine range
      call convert_position_sam_to_sky(sam(i)%getPosition()/para%L,itile,base%dc,base%ra,base%dec,xbox)
      do j = 1,3
         xmin(j) = min(xmin(j),xbox(j))
         xmax(j) = max(xmax(j),xbox(j))
      end do
   
      ! check if distance is in the range covered by snapshot isnapshot
      if ((base%dc>=snapshot(isnapshot)%dmin).and.(base%dc<snapshot(isnapshot)%dmax)) then
         n_in_snapshot = n_in_snapshot+1
      else
         ok = .false.
      end if
      
      ! check full position-selection
      if (is_in_fov(base%dc*para%L,base%ra,base%dec).and. &
         & pos_selection(base%dc*para%L,base%ra/degree,base%dec/degree)) then
         n_in_survey_volume = n_in_survey_volume+1
      else
         ok = .false.
      end if
      
      n_in_group = n_in_group+1
      
      if (.not.ok) index(k,2) = 0
      
   end do
   
   ! close file
   close(1)

contains
   
   subroutine reset_flagging(start,group)
      implicit none
      integer*4,intent(in) :: start
      integer*8,intent(in) :: group
      xmin = 2.0
      xmax = -1.0
      n_in_group = 0
      n_in_snapshot = 0
      n_in_survey_volume = 0
      old_group = group
      group_start = start
   end subroutine reset_flagging

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
         snapshot(i)%dmax = 9999.0
      else
         snapshot(i)%dmax = 0.5*(d(i)+d(i-1))
      end if
   end do

end subroutine make_distance_ranges
   
end module module_sky_intrinsic