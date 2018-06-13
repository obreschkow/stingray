module module_positioning

   use module_constants
   use module_system
   use module_types
   use module_cosmology
   use module_sort
   use module_user
   use module_parameters
   use module_tiling
   
   private
   public   :: make_positioning
   
   type(type_sam),allocatable    :: sam(:)
   integer*8                     :: nmockgalaxies
   character(len=255)            :: filename_sky_intrinsic
   logical,allocatable           :: sam_sel(:)
   integer*4                     :: n_sam_selected
   
contains

subroutine make_positioning

   implicit none
   integer*4            :: isnapshot,isubsnapshot,itile
   character(len=100)   :: snapshotname
   
   call tic
   call out('POSITION OBJECTS INTO SKY')
   
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
            call preprocess_snapshot(n_sam_selected)
            if (n_sam_selected>0) then
               do itile = 1,size(tile)
                  if ((snapshot(isnapshot)%dmax>=tile(itile)%dmin).and.(snapshot(isnapshot)%dmin<=tile(itile)%dmax)) then
                     call write_subsnapshot_into_tile(itile,isnapshot)
                  end if
               end do
            end if
         end do
      end if
   end do
   
   ! deallocate arrays
   if (allocated(tile)) deallocate(tile)
   if (allocated(snapshot)) deallocate(snapshot)
   if (allocated(sam)) deallocate(sam)
   if (allocated(sam_sel)) deallocate(sam_sel)
   
   ! check number of galaxies
   if (nmockgalaxies==0) then
      !call error('No galaxy in sky. Consider widening the sky geometry or relaxing the selection criteria.')
   end if
   
   ! add number of objects to beginning of file
   open(1,file=trim(filename_sky_intrinsic),action='write',form='unformatted',status='old',access='stream')
   write(1,pos=1) nmockgalaxies
   close(1) 
   
   ! finalize output
   call out('Number of galaxies in intrinsic sky:',nmockgalaxies)
   call toc
   
end subroutine make_positioning

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

subroutine preprocess_snapshot(n_sam_selected)
   
   implicit none
   integer*4,intent(out)   :: n_sam_selected
   integer*8,allocatable   :: id(:,:)
   integer*4               :: group_n ! number of objects in current group
   integer*4               :: i,n_objects,n
   logical                 :: group_selected
   
   !call out('preprocess')
   
   ! determine if the objects pass the SAM selection
   n = size(sam)
   if (allocated(sam_sel)) deallocate(sam_sel)
   allocate(sam_sel(n))
   do i = 1,size(sam)
      sam_sel(i) = sam_selection(sam(i))
   end do
   
   ! order galaxies by group id
   allocate(id(n,2))
   do i = 1,n
      id(i,1) = getGroupID(sam(i))
      id(i,2) = i
   end do
   call merge_sort_list(id)
   
   ! only retain groups, where at least one object is selected based on its SAM properties
   n_objects = 0
   group_selected = .false.
   group_n = 0
   do i = 1,n
      group_n = group_n+1
      group_selected = group_selected.or.sam_sel(id(i,2))
      if (id(i,1).ne.id(modulo(i,n)+1,1)) then ! check if group id differs from next group id
         if (group_selected) then
            id(n_objects+1:n_objects+group_n,:) = id(i-group_n+1:i,:)
            n_objects = n_objects+group_n 
         end if
         group_selected = .false.
         group_n = 0
      end if
   end do
   
   ! apply reordering and selection to arrays sam(:) and sam_sel(:)
   sam = sam(id(1:n_objects,2))
   sam_sel = sam_sel(id(1:n_objects,2))
   n_sam_selected = n_objects
   deallocate(id)
   
end subroutine preprocess_snapshot

subroutine write_subsnapshot_into_tile(itile,isnapshot)

   implicit none
   integer*4,intent(in)    :: itile
   integer*4,intent(in)    :: isnapshot
   integer*4               :: i,j,n
   type(type_base)         :: base
   real*4                  :: xbox(3)
   real*4                  :: xmin(3),xmax(3)
   integer*4               :: n_in_snapshot
   integer*4               :: n_in_survey_volume
   real*4,allocatable      :: dc(:),ra(:),dec(:)
   logical,allocatable     :: ok(:)
   logical                 :: group_selected
   
   !call out('write')
   
   ! set tile
   base%tile = itile
   
   ! allocate arrays
   n = size(sam)
   allocate(dc(n),ra(n),dec(n),ok(n))
   ok = .true.
   
   ! open file
   open(1,file=trim(filename_sky_intrinsic),action='write',form='unformatted',status='old',position='append',access='stream')
   
   ! initialise group flagging
   call reset_flagging
   
   ! iterate over all mock galaxies, ordered by group, and check how groups have been truncated
   do i = 1,n
      
      base%group_ntot = base%group_ntot+1
      
      ! compute sky position and determine range
      call convert_position_sam_to_sky(sam(i)%getPosition()/para%L,itile,dc(i),ra(i),dec(i),xbox)
      do j = 1,3
         xmin(j) = min(xmin(j),xbox(j))
         xmax(j) = max(xmax(j),xbox(j))
      end do
   
      ! check if distance is in the range covered by snapshot isnapshot
      if ((dc(i)>=snapshot(isnapshot)%dmin).and.(dc(i)<snapshot(isnapshot)%dmax)) then
         n_in_snapshot = n_in_snapshot+1
         group_selected = .true.
      else
         ok(i) = .false. ! => this object will be rejected
      end if
      
      ! check full position-selection
      if (is_in_fov(dc(i)*para%L,ra(i),dec(i)).and. &
         & pos_selection(dc(i)*para%L,ra(i)/degree,dec(i)/degree)) then
         n_in_survey_volume = n_in_survey_volume+1
         group_selected = .true.
      else
         ok(i) = .false. ! => this object will be rejected
      end if
   
      ! assign group flags
      if (getGroupID(sam(i)).ne.getGroupID(sam(modulo(i,n)+1))) then ! check if group id differs from next group id
      
         if (group_selected) then
         
            ! re-iterate over all accepted group members to check if they also pass the SAM selection
            base%group_nsel = 0
            do j = i-base%group_ntot+1,i
               if (ok(j)) then
                  if (sam_sel(j)) then
                     base%group_nsel = base%group_nsel+1
                  else
                     ok(j) = .false.
                  end if
               end if
            end do
            
            if (base%group_nsel>0) then
      
               ! make group flag
               base%group_flag = 0
               if (n_in_survey_volume.ne.base%group_ntot) base%group_flag = base%group_flag+1   ! group truncated by survey volume
               if (n_in_snapshot.ne.base%group_ntot) base%group_flag = base%group_flag+2        ! group truncated by snapshot limit
               if (maxval(xmax-xmin)>0.5) base%group_flag = base%group_flag+4                   ! group wrapped around
               
               ! re-iterate over all accepted group members to write them in file
               do j = i-base%group_ntot+1,i
                  if (ok(j)) then
                     nmockgalaxies = nmockgalaxies+1
                     base%dc = dc(j)
                     base%ra = ra(j)
                     base%dec = dec(j)
                     write(1) base,sam(j)
                  end if
               end do
            end if
         
         end if
         
         ! reset group flagging
         call reset_flagging
            
      end if
      
   end do
   
   ! close file
   close(1)

contains
   
   subroutine reset_flagging
      implicit none
      xmin = 2.0
      xmax = -1.0
      base%group_ntot = 0
      n_in_snapshot = 0
      n_in_survey_volume = 0
      group_selected = .false.
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
   
end module module_positioning