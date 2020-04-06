module module_sky

   use module_constants
   use module_system
   use module_types
   use module_io
   use module_linalg
   use module_cosmology
   use module_sort
   use module_user_routines
   use module_user_selection
   use module_parameters
   use module_tiling
   
   private
   public   :: make_sky
   
   !character(len=255)            :: filename_sky_intrinsic
   
contains

subroutine make_sky

   implicit none
   integer*4                     :: isnapshot,isubvolume,itile,i
   integer*4                     :: nsub,nsub_max
   integer*4,allocatable         :: index(:,:)
   type(type_sam),allocatable    :: sam(:)
   logical,allocatable           :: sam_sel(:)
   integer*4,allocatable         :: sam_replica(:)
   integer*8                     :: n_distinct_galaxies
   integer*4                     :: n_replica_max
   real*4                        :: n_replica_mean
   integer*8                     :: n_galaxies,n_galaxies_tot
   integer*8                     :: n_groups,n_groups_tot
   character(len=255)            :: str,strt
   character(len=255)            :: filename
   integer,external              :: OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS
   
   call tic
   call out('POSITION OBJECTS INTO SKY AND MAKE APPARENT PROPERTIES')
   
   ! load previous steps
   call load_parameters
   call load_tile_list
   
   ! make snapshot properties
   if (allocated(snapshot)) deallocate(snapshot)
   allocate(snapshot(para%snapshot_min:para%snapshot_max))
   call make_redshifts
   call make_distance_ranges
   call count_tiles_of_snapshot
   call save_snapshot_list
   
   ! initialize apparent sky (galaxies)
   filename = trim(para%path_output)//fn_galaxies
   open(1,file=trim(filename),action='write',form="unformatted",status='replace',access='stream')
   write(1) 0_8,0_4,0_4 ! place holder for number of objects in mock sky and replication values
   
   ! initialize apparent sky (groups)
   if (para%make_groups==1) then
      filename = trim(para%path_output)//fn_groups
      open(2,file=trim(filename),action='write',form="unformatted",status='replace',access='stream')
      write(2) 0_8 ! place holder for number of objects in mock sky
   end if
   
   ! make snapshot indices
   nsub_max = (para%snapshot_max-para%snapshot_min+1)*(para%subvolume_max-para%subvolume_min+1)
   allocate(index(nsub_max,2))
   nsub = 0
   do isnapshot = para%snapshot_min,para%snapshot_max
      if ((snapshot(isnapshot)%dmax>=minval(tile%dmin)).and.(snapshot(isnapshot)%dmin<=maxval(tile%dmax))) then
         do isubvolume = para%subvolume_min,para%subvolume_max
            nsub = nsub+1
            index(nsub,:) = (/isnapshot,isubvolume/)
         end do
      end if
   end do
   
   ! fill galaxies with intrinsic properties into sky
   n_galaxies_tot = 0
   n_groups_tot = 0
   n_distinct_galaxies = 0
   n_replica_max = 0
   
   !$OMP PARALLEL PRIVATE(itile,sam,sam_sel,sam_replica,n_galaxies,n_groups)
   !$OMP DO SCHEDULE(DYNAMIC)
   do i = 1,nsub
      if (snapshot(index(i,1))%n_tiles>0) then
      
         ! Only do the following by one thread at a time
         !$OMP CRITICAL
         call load_sam_snapshot(index(i,1),index(i,2),sam)
         !$OMP END CRITICAL
         
         ! 
         call preprocess_snapshot(sam,sam_sel)
         
         if (size(sam)>0) then
            allocate(sam_replica(size(sam)))
            sam_replica = 0
            do itile = 1,size(tile)
               if ((snapshot(index(i,1))%dmax>=tile(itile)%dmin).and.(snapshot(index(i,1))%dmin<=tile(itile)%dmax)) then
               
                  call write_subvolume_into_tile(itile,index(i,1),index(i,2),sam,sam_sel,sam_replica,n_galaxies,n_groups)
                  
                  ! Only do the following by one thread at a time
                  !$OMP CRITICAL
                  n_galaxies_tot = n_galaxies_tot+n_galaxies
                  n_groups_tot = n_groups_tot+n_groups
                  strt = 'Thread 1/1'
                  !$ write(strt,'(A,I0,A,I0)') 'Thread ',OMP_GET_THREAD_NUM()+1,'/',OMP_GET_NUM_THREADS()
                  write(str,'(A,A,I0,A,I0,A,I0,A,I0,A,I0,A)') trim(strt),': Write snapshot ',index(i,1),', subvolume ', &
                     & index(i,2),' into tile ',itile,' (',n_galaxies,'/',n_galaxies_tot,' galaxies)'
                  call out(trim(str))
                  if ((n_galaxies_tot>=1e8_8).and.(n_galaxies_tot-n_galaxies<1e8_8)) call error(&
                  &'The cone has reached 1e8 galaxies. Consider using a more restrictive selection function.')
                  !$OMP END CRITICAL
                  
               end if
            end do
            n_distinct_galaxies = n_distinct_galaxies+count(sam_replica>0)
            n_replica_max = max(n_replica_max,maxval(sam_replica))
            deallocate(sam_replica)
         end if
         
      end if
   end do
   !$OMP END DO NOWAIT
   !$OMP END PARALLEL       
   
   ! deallocate arrays
   if (allocated(tile)) deallocate(tile)
   if (allocated(snapshot)) deallocate(snapshot)
   if (allocated(sam)) deallocate(sam)
   if (allocated(sam_sel)) deallocate(sam_sel)
   
   ! check number of galaxies
   if (n_galaxies_tot==0) then
      call warning('No galaxy in sky. Consider widening the sky geometry or relaxing the selection criteria.')
   end if
   
   ! add number of objects to beginning of file & close files
   n_replica_mean = real(real(n_galaxies_tot,8)/n_distinct_galaxies,4)
   write(1,pos=1) n_galaxies_tot,n_replica_mean,n_replica_max
   close(1)
   if (para%make_groups==1) then
      write(2,pos=1) n_groups_tot
      close(2)
   end if
   
   ! user output
   call out('Number of galaxies in mock sky:',n_galaxies_tot)
   if (para%make_groups==1) call out('Number of groups in mock sky:',n_groups_tot)
   call toc
   
end subroutine make_sky

subroutine convert_position_sam_to_sky(position,itile,dc,ra,dec,x,xbox)

   implicit none
   real*4,intent(in)    :: position(3) ! [box side-length] position in box
   integer*4,intent(in) :: itile
   real*4,intent(out)   :: dc,ra,dec   ! [box side-length,rad,rad] position on sky
   real*4,intent(out)   :: x(3)        ! [box side-length] position in box after symmetry operations
   real*4,intent(out)   :: xbox(3)     ! [box side-length] position in box after symmetry operations, before projection
   
   xbox = modulo(matmul(rot(:,:,tile(itile)%rotation),position)+tile(itile)%translation,1.0) ! apply symmetries
   x = xbox+tile(itile)%ix-0.5 ! translate coordinates to tile position
   x = matmul(para%sky_rotation,x) ! convert SAM-coordinates to Sky-coordinates
   call car2sph(x,dc,ra,dec)
   
end subroutine convert_position_sam_to_sky

subroutine preprocess_snapshot(sam,sam_sel)

   ! 1) determine for each galaxy, if it has been selected according to it sam_selection
   ! 2) sort galaxies in 'sam' by groups, making the central galaxy the first one in its group
   ! 3) remove all groups, where no galaxy was selected
   
   implicit none
   type(type_sam),allocatable,intent(inout)  :: sam(:)
   logical,allocatable,intent(inout)         :: sam_sel(:)
   integer*8,allocatable                     :: id(:,:)
   integer*4                                 :: group_n ! number of objects in current group
   integer*4                                 :: group_central_n ! number of central objects in current group
   integer*4                                 :: i,n_objects,n
   logical                                   :: group_selected
   
   n = size(sam)
   
   if (para%make_groups==1) then
   
      ! determine if the objects pass the SAM selection
      if (allocated(sam_sel)) deallocate(sam_sel)
      allocate(sam_sel(n))
      do i = 1,n
         sam_sel(i) = sam_selection(sam(i))
      end do
   
      ! order galaxies by group id
      allocate(id(n,2))
      do i = 1,n
         id(i,1) = sam(i)%get_groupid()*2+1
         if (sam(i)%is_group_center()) id(i,1)=id(i,1)-1 ! to make sure that this galaxy gets listed first
         id(i,2) = i
      end do
      call merge_sort_list(id)
   
      ! only retain groups, where at least one object is selected based on its SAM properties
      n_objects = 0
      group_selected = .false.
      group_n = 0
      group_central_n = 0
      do i = 1,n
   
         group_n = group_n+1
         if (mod(id(i,1),2)==0) group_central_n = group_central_n+1
         group_selected = group_selected.or.sam_sel(id(i,2))
      
         if ((id(i,1)/2.ne.id(min(i+1,n),1)/2).or.(i==n)) then ! check if group id differs from next group id
      
            ! close group
            if (group_central_n.ne.1) then
               write(*,*) group_central_n
               call error('Each group must have exactly one central member.')
            end if
            if (group_selected) then
               id(n_objects+1:n_objects+group_n,:) = id(i-group_n+1:i,:)
               n_objects = n_objects+group_n 
            end if
         
            ! start new group
            group_selected = .false.
            group_n = 0
            group_central_n = 0
         
         end if
      
      end do
   
      ! apply reordering and selection to arrays sam(:) and sam_sel(:)
      sam = sam(id(1:n_objects,2))
      sam_sel = sam_sel(id(1:n_objects,2))
      deallocate(id)
      
   else
   
      n_objects = 0
      do i = 1,n
         if (sam_selection(sam(i))) then
            n_objects = n_objects+1
            sam(n_objects) = sam(i)
         end if
      end do
      sam = sam(1:n_objects)
      
      if (allocated(sam_sel)) deallocate(sam_sel)
      allocate(sam_sel(n_objects))
      sam_sel = .true.
   
   end if
   
end subroutine preprocess_snapshot

subroutine write_subvolume_into_tile(itile,isnapshot,isubvolume,sam,sam_sel,sam_replica,n_galaxies,n_groups)

   ! NB: most of this routine deals with groups
   ! if no make_groups==0, this routine essentially places all the galaxies in 'sam' into the sky,
   ! making sure to only retain the objects that pass pos_selection, sam_selection, pre_selection, sky_selection

   implicit none
   integer*4,intent(in)                      :: itile
   integer*4,intent(in)                      :: isnapshot
   integer*4,intent(in)                      :: isubvolume
   type(type_sam),allocatable,intent(in)     :: sam(:)
   logical,allocatable,intent(in)            :: sam_sel(:)
   integer*4,allocatable,intent(inout)       :: sam_replica(:)
   integer*8,intent(out)                     :: n_galaxies
   integer*8,intent(out)                     :: n_groups
   type(type_base)                           :: base
   integer*4                                 :: i,j,n,jmin,k,d,group_nselected
   real*4                                    :: xbox(3)
   real*4                                    :: xmin(3),xmax(3),dx(3)
   integer*4                                 :: n_in_snapshot
   integer*4                                 :: n_in_survey_volume
   real*4,allocatable                        :: dc(:),ra(:),dec(:),x(:,:)
   logical,allocatable                       :: ok(:)
   logical                                   :: group_preselected
   logical                                   :: last_galaxy_in_group
   logical                                   :: wrapped
   integer*8                                 :: galaxyid,groupid,prefixid
   type(type_sky_galaxy),allocatable         :: sky_galaxy(:)
   type(type_sky_group)                      :: sky_group
   
   ! set tile
   base%tile = itile
   
   ! allocate arrays
   n = size(sam)
   allocate(dc(n),ra(n),dec(n),x(n,3),ok(n))
   ok = .true.
   
   ! initialise group flagging
   call reset_group_variables
   n_galaxies = 0
   n_groups = 0
   if (n>9999999) call error('There are more than 9,999,999 galaxies in a single subvolume.')
   if (itile>9999) call error('There are more than 9999 tiles, which is not allowed.')
   if (isubvolume>999) call error('There are more than 999 subvolumes, which is not allowed.')
   if (isubvolume>999) call error('There are more than 999 subvolumes, which is not allowed.')
   prefixid = 100000000000000_8*isnapshot+100000000000_8*isubvolume+10000000_8*itile
   
   ! iterate over all mock galaxies in the subvolume, ordered by group, and check how groups have been truncated
   do i = 1,n
      
      ! count number of galaxies in the group
      base%group_ntot = base%group_ntot+1
      
      ! check if this galaxy is the last in its group
      last_galaxy_in_group = (sam(i)%get_groupid().ne.sam(min(i+1,n))%get_groupid()).or.(i==n)
      
      ! compute sky position and determine range
      call convert_position_sam_to_sky(sam(i)%get_position()/para%box_side,itile,dc(i),ra(i),dec(i),x(i,:),xbox)
      do j = 1,3
         xmin(j) = min(xmin(j),xbox(j))
         xmax(j) = max(xmax(j),xbox(j))
      end do
   
      ! check if distance is in the range covered by snapshot isnapshot
      if ((dc(i)>=snapshot(isnapshot)%dmin).and.(dc(i)<snapshot(isnapshot)%dmax)) then
         n_in_snapshot = n_in_snapshot+1
      else
         ok(i) = .false. ! => this object will be rejected
      end if
      
      ! check full position
      if (is_in_fov(dc(i)*para%box_side,ra(i),dec(i)).and.pos_selection(dc(i)*para%box_side,ra(i)/degree,dec(i)/degree)) then
         n_in_survey_volume = n_in_survey_volume+1
      else
         ok(i) = .false. ! => this object will be rejected
      end if
      
      ! check if any galaxy satisfies the snapshot+position selection
      group_preselected = group_preselected.or.ok(i) ! true if at least one galaxy in the group has passed the previous tests
      
      ! close this group
      if (last_galaxy_in_group) then ! check if group id differs from next group id
      
         if (group_preselected) then
         
            ! allocate sky array
            jmin = i-base%group_ntot+1
            if (allocated(sky_galaxy)) deallocate(sky_galaxy)
            allocate(sky_galaxy(jmin:i))
            
            ! make group id
            if (base%group_ntot==1) then
               groupid = -1_8
            else
               groupid = prefixid+n_groups+1
            end if
            
            ! make group flag
            if (base%group_ntot==1) then
               base%group_flag = 0
            else
               wrapped = maxval(xmax-xmin)>0.5
               base%group_flag = 0
               if (n_in_survey_volume.ne.base%group_ntot) base%group_flag = base%group_flag+1   ! group truncated by survey volume
               if (n_in_snapshot.ne.base%group_ntot) base%group_flag = base%group_flag+2        ! group truncated by snapshot limit
               if (wrapped) base%group_flag = base%group_flag+4                                 ! group wrapped around
            end if
               
            ! re-iterate over all accepted group members to check if they also pass the sam_selection, pre_selection and sky_selection
            group_nselected = 0
            do j = jmin,i
               if (ok(j)) then
                  ok(j) = .false.
                  if (sam_sel(j)) then
                     galaxyid = prefixid+n_galaxies+1
                     ok(j) = pre_selection(sam(j),dc(j)*para%box_side,ra(j)/degree,dec(j)/degree)
                     if (ok(j)) then
                        base%ra = ra(j)
                        base%dec = dec(j)
                        base%dc = dc(j)
                        call sky_galaxy(j)%make_from_sam(sam(j),base,groupid,galaxyid)
                        ok(j) = sky_selection(sky_galaxy(j),sam(j))
                        if (ok(j)) then
                           n_galaxies = n_galaxies+1
                           group_nselected = group_nselected+1
                           sam_replica(j) = sam_replica(j)+1
                        end if
                     end if
                  end if
               end if
            end do
            
            ! if at least one group member has passed all the selections, write relevant galaxies
            if (group_nselected>0) then
               
               ! write group, if the "group" has intrinsically more than one member (and if user wants group)
               if ((para%make_groups==1).and.(base%group_ntot>1)) then
               
                  n_groups = n_groups+1
            
                  if (.not.ok(jmin)) then
                     ! jmin is a first group member (group centre), which has not been selected, but is saved
                     ! for the group catalogue. In case this galaxy has been mapped far away from *all* the
                     ! selected group members, due to wrapping at the edge of the tile, it is here moved back
                     ! next to the selected galaxies in the group, accepting a violation of the tile boundary
                     if (base%group_flag>=4) then ! group wrapped around
                        xmin = 2.0
                        xmax = -1.0
                        do k = jmin+1,i
                           if (ok(k)) then
                              do d = 1,3
                                 xmin(d) = min(xmin(d),x(k,d))
                                 xmax(d) = max(xmax(d),x(k,d))
                              end do
                           end if
                        end do
                        if (maxval(xmax-xmin)<0.5) then
                           dx = (xmax+xmin)/2.0-x(jmin,:)
                           if (sum(dx**2)>0.25) then
                              dx = nint(matmul(transpose(para%sky_rotation),dx))
                              x(jmin,:) = x(jmin,:)+matmul(para%sky_rotation,dx)
                              call car2sph(x(jmin,:),dc(jmin),ra(jmin),dec(jmin))
                           end if
                        end if
                     end if
                  end if
               
                  ! make group
                  base%dc = dc(jmin)
                  base%ra = ra(jmin)
                  base%dec = dec(jmin)
                  call sky_group%make_from_sam(sam(jmin:i),sky_galaxy(jmin:i),ok(jmin:i),base,groupid,group_nselected)
                  
               end if
               
               ! Only do the following by one thread at a time
               !$OMP CRITICAL
               do j = jmin,i
                  if (ok(j)) call sky_galaxy(j)%write_to_file(1)
               end do
               if ((para%make_groups==1).and.(base%group_ntot>1)) call sky_group%write_to_file(2)
               !$OMP END CRITICAL
               
            end if
         
         end if
         
         call reset_group_variables
            
      end if
      
   end do
   
contains
   
   subroutine reset_group_variables
      implicit none
      xmin = 2.0
      xmax = -1.0
      base%group_ntot = 0
      n_in_snapshot = 0
      n_in_survey_volume = 0
      group_preselected = .false.
   end subroutine reset_group_variables

end subroutine write_subvolume_into_tile

subroutine make_distance_ranges

   implicit none
   integer*4            :: i
   real*4,allocatable   :: d(:)
   
   allocate(d(para%snapshot_min:para%snapshot_max))
   do i = para%snapshot_min,para%snapshot_max
      d(i) = redshift_to_dc(snapshot(i)%redshift)*(Mpc/para%length_unit)/para%box_side ! [box side-length] comoving distance to redshift of the box
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

subroutine count_tiles_of_snapshot

   implicit none
   integer*4   :: isnapshot,itile

   snapshot%n_tiles = 0
   do isnapshot = para%snapshot_min,para%snapshot_max
      do itile = 1,size(tile)
         if ((snapshot(isnapshot)%dmax>=tile(itile)%dmin).and.(snapshot(isnapshot)%dmin<=tile(itile)%dmax)) then
            snapshot(isnapshot)%n_tiles = snapshot(isnapshot)%n_tiles+1
         end if
      end do
   end do
   
end subroutine count_tiles_of_snapshot
   
end module module_sky