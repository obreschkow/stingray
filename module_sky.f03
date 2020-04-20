module module_sky

   use shared_module_core
   use shared_module_constants
   use shared_module_parameters
   use shared_module_cosmology
   use shared_module_maths
   use shared_module_sort
   use module_global
   use module_user_routines
   use module_user_selection
   use module_tiling
   
   public   :: make_sky             ! positions galaxies into mock sky and produces their properties
   public   :: write_sky_to_hdf5    ! call user routine make_hdf5() to save mock sky in HDF5 format
   
   private
   
   integer*4,parameter :: fid = 1 ! default file unit used in this module to write mock sky data
   
contains

subroutine make_sky

   use omp_lib

   implicit none
   
   type type_task
      integer*4   :: isnapshot
      integer*4   :: isubvolume
   end type type_task
   
   integer*4                           :: isnapshot,isubvolume,itile
   type(type_sam),allocatable          :: sam(:)
   logical*4,allocatable               :: sam_sel(:)
   integer*4,allocatable               :: sam_replica(:),sam_replica_tile(:)
   type(type_sky_galaxy),allocatable   :: sky_galaxy(:)
   type(type_sky_group),allocatable    :: sky_group(:)
   character(255)                      :: strt = 'Thread 1/1'  ! default thread number, openmp not used
   character(255)                      :: str
   character(255)                      :: filename
   type(type_task),allocatable         :: task(:)
   type(type_skystats),allocatable     :: substats(:) ! statistics of mock sky for particular subvolumes
   type(type_skystats)                 :: totstats    ! statistics of total mock sky
   integer*4                           :: itask,n_tasks
   
   call tic('POSITION OBJECTS INTO SKY AND MAKE APPARENT PROPERTIES')
   
   ! make snapshot properties
   allocate(snapshot(para%snapshot_min:para%snapshot_max))
   call make_redshifts
   call make_distance_ranges
   call count_tiles_of_snapshot
   
   ! make task list
   n_tasks = (para%snapshot_max-para%snapshot_min+1)*(para%subvolume_max-para%subvolume_min+1)
   allocate(task(n_tasks))
   n_tasks = 0
   do isnapshot = para%snapshot_min,para%snapshot_max
      if (snapshot(isnapshot)%n_tiles>0) then ! snapshot is used for at least some tiles
         do isubvolume = para%subvolume_min,para%subvolume_max
            n_tasks = n_tasks+1
            task(n_tasks)%isnapshot = isnapshot
            task(n_tasks)%isubvolume = isubvolume
         end do
      end if
   end do
   task = task(1:n_tasks)
   
   ! initialze binary files
   do isubvolume = para%subvolume_min,para%subvolume_max
      filename = filename_sky_galaxies(isubvolume)
      open(fid,file=trim(filename),action='write',form='unformatted',status='replace',access='stream')
      write(fid) 0_4,0_4,0_4 ! place holder for number of galaxies and replication values
      close(fid)
      if (para%make_groups==1) then
         filename = filename_sky_groups(isubvolume)
         open(fid,file=trim(filename),action='write',form='unformatted',status='replace',access='stream')
         write(fid) 0_4 ! place holder for number of groups
         close(fid)
      end if
   end do
   
   ! fill galaxies with intrinsic properties into sky
   allocate(substats(para%subvolume_min:para%subvolume_max))
   
   !$ call OMP_set_nested(.true.) ! enables nested parallelism
   
   !$OMP PARALLEL PRIVATE(isnapshot,isubvolume,itile,sam,sam_sel,sam_replica)
   !$OMP DO
   do itask = 1,n_tasks
   
      isnapshot = task(itask)%isnapshot
      isubvolume = task(itask)%isubvolume
   
      !$OMP CRITICAL ! load snapshot, one thread at a time
      call load_sam_snapshot(isnapshot,isubvolume,sam)
      !$OMP END CRITICAL
         
      call preprocess_snapshot(sam,sam_sel)
         
      if (size(sam)>0) then
      
         allocate(sam_replica(size(sam)))
         sam_replica = 0
         
         !$OMP PARALLEL PRIVATE(sam_replica_tile,sky_galaxy,sky_group,filename,str,strt)
         !$OMP DO
         do itile = 1,size(tile)
            if ((snapshot(isnapshot)%dmax>=tile(itile)%dmin).and.(snapshot(isnapshot)%dmin<=tile(itile)%dmax)) then
         
               allocate(sam_replica_tile(size(sam)))
               sam_replica_tile = 0
               call place_subvolume_into_tile(itile,isnapshot,isubvolume,sam,sam_sel,sam_replica_tile,sky_galaxy,sky_group)
               
               !$OMP CRITICAL 
               ! count objects, one thread at a time
               sam_replica = sam_replica+sam_replica_tile
               deallocate(sam_replica_tile)
               substats(isubvolume)%n_galaxies = substats(isubvolume)%n_galaxies+size(sky_galaxy)
               substats(isubvolume)%n_groups = substats(isubvolume)%n_groups+size(sky_group)
               totstats%n_galaxies = totstats%n_galaxies+size(sky_galaxy)
               totstats%n_groups = totstats%n_groups+size(sky_group)
               if ((totstats%n_galaxies>=limit%n_galaxies_sky_max).and.&
               &(totstats%n_galaxies-size(sky_galaxy)<limit%n_galaxies_sky_max)) &
               &call error('Number of galaxies in the sky exceeds ',limit%n_galaxies_sky_max)
               
               ! progress update, one thread at a time
               !$ write(strt,'(A,I0,A,I0)') 'Thread ',OMP_GET_THREAD_NUM()+1,'/',OMP_GET_NUM_THREADS()
               write(str,'(A,A,I0,A,I0,A,I0,A,I0,A,I0,A)') trim(strt),': Write snapshot ',isnapshot,', subvolume ', &
                  & isubvolume,' into tile ',itile,' (',size(sky_galaxy),'/',totstats%n_galaxies,' galaxies)'
               call out(trim(str))
               
               ! save_sky_part, one thread at a time
               filename = filename_sky_galaxies(isubvolume)
               open(fid,file=trim(filename),action='write',status='old',position='append',form='unformatted',access='stream')
               write(fid) sky_galaxy
               close(fid)
               if (para%make_groups==1) then
                  filename = filename_sky_groups(isubvolume)
                  open(fid,file=trim(filename),action='write',status='old',position='append',form='unformatted',access='stream')
                  write(fid) sky_group
                  close(fid)
               end if
               !$OMP END CRITICAL
              
            end if
         end do
         !$OMP END DO
         !$OMP END PARALLEL  
               
         substats(isubvolume)%n_distinct = substats(isubvolume)%n_distinct+count(sam_replica>0)
         substats(isubvolume)%n_replica_max = max(substats(isubvolume)%n_replica_max,maxval(sam_replica))
         deallocate(sam_replica)
         
      end if
      
   end do ! task
   !$OMP END DO
   !$OMP END PARALLEL
      
   ! finalize binary files
   do isubvolume = para%subvolume_min,para%subvolume_max
      filename = filename_sky_galaxies(isubvolume)
      open(fid,file=trim(filename),action='write',form='unformatted',status='old',access='stream')
      write(fid,pos=1) substats(isubvolume)%n_galaxies,substats(isubvolume)%n_distinct,substats(isubvolume)%n_replica_max
      close(fid)
      if (para%make_groups==1) then
         filename = filename_sky_groups(isubvolume)
         open(fid,file=trim(filename),action='write',form='unformatted',status='old',access='stream')
         write(fid) substats(isubvolume)%n_groups
         close(fid)
      end if
   end do
   
   ! check number of galaxies
   if (totstats%n_galaxies==0) then
      call warning('No galaxy in sky. Consider widening the sky geometry or relaxing the selection criteria.')
   end if
   
   ! user output
   call out('Number of galaxies in mock sky: ',totstats%n_galaxies)
   if (para%make_groups==1) call out('Number of groups in mock sky: ',totstats%n_groups)
   call toc
   
end subroutine make_sky

function filename_sky_galaxies(isubvolume) result(fn)

   implicit none
   integer*4,intent(in) :: isubvolume
   character(255)       :: fn
   
   fn = dir(para%path_output,fn_galaxies//'.'//val2str(isubvolume))

end function filename_sky_galaxies

function filename_sky_groups(isubvolume) result(fn)

   implicit none
   integer*4,intent(in) :: isubvolume
   character(255)       :: fn
   
   fn = dir(para%path_output,fn_groups//'.'//val2str(isubvolume))

end function filename_sky_groups

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
   call car2sph(x,dc,ra,dec,astro=.true.)
   
end subroutine convert_position_sam_to_sky

subroutine preprocess_snapshot(sam,sam_sel)

   ! 1) determine for each galaxy, if it has been selected according to it sam_selection
   ! 2) sort galaxies in 'sam' by groups, making the central galaxy the first one in its group
   ! 3) remove all groups, where no galaxy was selected
   
   implicit none
   type(type_sam),allocatable,intent(inout)  :: sam(:)
   logical*4,allocatable,intent(inout)       :: sam_sel(:)
   integer*8,allocatable                     :: id(:)
   integer*4,allocatable                     :: index(:)
   integer*4                                 :: i,n
   integer*4                                 :: n_groups
   integer*4                                 :: i_group
   integer*8                                 :: groupid
   integer*4                                 :: n_galaxies_in_group
   logical*4                                 :: group_selected
   integer*4                                 :: n_keep
   
   n = size(sam)
   
   if (para%make_groups==1) then
   
      ! determine if the objects pass the SAM selection
      if (allocated(sam_sel)) deallocate(sam_sel)
      allocate(sam_sel(n))
      do i = 1,n
         sam_sel(i) = sam_selection(sam(i))
      end do
   
      ! order galaxies by group id, such that central galaxies come first
      allocate(id(n),index(n))
      n_groups = 0
      do i = 1,n
         id(i) = sam(i)%get_groupid()*2+1
         if (sam(i)%is_group_center()) then
            id(i)=id(i)-1 ! to make sure that this galaxy gets listed first
            n_groups = n_groups+1
         end if
      end do
      call sort(id,index)
      sam = sam(index)
      sam_sel = sam_sel(index)
      
      ! keep all galaxies of groups, where at least one object is sam-selected; reject all other galaxies
      n_keep = 0
      i = 1
      do i_group = 1,n_groups
      
         ! handle central galaxy
         if (.not.sam(i)%is_group_center()) call error('each group must have exactly one central member.')
         groupid = sam(i)%get_groupid() ! identical to sam(i)%get_groupid(), but faster
         
         ! count number of galaxies in group and number of selected galaxies in group
         n_galaxies_in_group = 1
         group_selected = sam_sel(i)
         index(n_keep+n_galaxies_in_group) = i
         i = i+1
         do while (i<=n.and.sam(min(i,n))%get_groupid()==groupid)
            if (sam(i)%is_group_center()) call error('each group must have exactly one central member.')
            n_galaxies_in_group = n_galaxies_in_group+1
            group_selected = group_selected.or.sam_sel(i)
            index(n_keep+n_galaxies_in_group) = i
            i = i+1
         end do
         
         ! keep group, if at least one galaxy is selected
         if (group_selected) n_keep = n_keep+n_galaxies_in_group
         
      end do
      
      if (i.ne.n+1) call error('something wrong with group indexing.')
      
      ! apply reordering and selection to arrays sam(:) and sam_sel(:)
      sam = sam(index(1:n_keep))
      sam_sel = sam_sel(index(1:n_keep))
      
   else
   
      n_keep = 0
      do i = 1,n
         if (sam_selection(sam(i))) then
            n_keep = n_keep+1
            sam(n_keep) = sam(i)
         end if
      end do
      sam = sam(1:n_keep)
      
      if (allocated(sam_sel)) deallocate(sam_sel)
      allocate(sam_sel(n_keep))
      sam_sel = .true.
   
   end if
   
end subroutine preprocess_snapshot

subroutine place_subvolume_into_tile(itile,isnapshot,isubvolume,sam,sam_sel,sam_replica,sky_galaxy_list,sky_group_list)

   ! NB: most of this routine deals with groups
   ! if no make_groups==0, this routine essentially places all the galaxies in 'sam' into the sky,
   ! making sure to only retain the objects that pass pos_selection, sam_selection, pre_selection, sky_selection

   implicit none
   integer*4,intent(in)                            :: itile
   integer*4,intent(in)                            :: isnapshot
   integer*4,intent(in)                            :: isubvolume
   type(type_sam),allocatable,intent(in)           :: sam(:)
   logical,allocatable,intent(in)                  :: sam_sel(:)
   integer*4,allocatable,intent(inout)             :: sam_replica(:)
   type(type_sky_galaxy),allocatable,intent(out)   :: sky_galaxy_list(:)
   type(type_sky_group),allocatable,intent(out)    :: sky_group_list(:)
   integer*8                                       :: n_galaxies,n_galaxies_saved
   integer*8                                       :: n_groups,n_groups_saved
   type(type_base)                                 :: base
   integer*4                                       :: i,j,n,jmin,k,d,group_nselected
   real*4                                          :: xbox(3)
   real*4                                          :: xmin(3),xmax(3),dx(3)
   integer*4                                       :: n_in_snapshot
   integer*4                                       :: n_in_survey_volume
   real*4,allocatable                              :: dc(:),ra(:),dec(:),x(:,:)
   logical,allocatable                             :: ok(:)
   logical                                         :: group_preselected
   logical                                         :: last_galaxy_in_group
   logical                                         :: wrapped
   integer*8                                       :: galaxyid,groupid,prefixid
   type(type_sky_galaxy),allocatable               :: sky_galaxy(:)
   type(type_sky_group)                            :: sky_group
   
   ! set tile
   base%tile = itile
   
   ! allocate arrays
   n = size(sam)
   allocate(dc(n),ra(n),dec(n),x(n,3),ok(n),sky_galaxy_list(n),sky_group_list(n))
   ok = .true.
   
   ! initialise group flagging
   call reset_group_variables
   n_galaxies = 0
   n_groups = 0
   n_galaxies_saved = 0
   n_groups_saved = 0
   if (n>limit%n_galaxies_sky_max) call error('number of galaxies in mock sky exceeds',limit%n_galaxies_sky_max)
   if (itile>limit%n_tiles_max) call error('number of tiles exceeds',limit%n_tiles_max)
   if (isubvolume>=limit%n_subvolumes_max) call error('number of subvolumes exceeds',limit%n_subvolumes_max)
   if (isnapshot>=limit%n_snapshots_max) call error('number of snapshots exceeds',limit%n_snapshots_max)
   prefixid = int(limit%n_galaxies_per_tile_max,8)*(int(limit%n_subvolumes_max,8)*int(isnapshot,8)+int(isubvolume,8))
   
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
      if (is_in_fov(dc(i)*para%box_side,ra(i),dec(i)).and. &
         & pos_selection(dc(i)*para%box_side,ra(i)/unit%degree,dec(i)/unit%degree)) then
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
                     ok(j) = pre_selection(sam(j),dc(j)*para%box_side,ra(j)/unit%degree,dec(j)/unit%degree)
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
                              call car2sph(x(jmin,:),dc(jmin),ra(jmin),dec(jmin),astro=.true.)
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
               
               ! save galaxies and groups
               do j = jmin,i
                  if (ok(j)) then
                     n_galaxies_saved = n_galaxies_saved+1
                     sky_galaxy_list(n_galaxies_saved) = sky_galaxy(j)
                  end if
               end do
               if ((para%make_groups==1).and.(base%group_ntot>1)) then
                  n_groups_saved = n_groups_saved+1
                  sky_group_list(n_groups_saved) = sky_group
               end if
               
            end if
         
         end if
         
         call reset_group_variables
            
      end if
      
   end do
   
   sky_galaxy_list = sky_galaxy_list(1:n_galaxies_saved)
   sky_group_list = sky_group_list(1:n_groups_saved)
   
   if (n_galaxies.ne.n_galaxies_saved) call deverror('n_galaxies.ne.n_galaxies_saved')
   if (n_groups.ne.n_groups_saved) call deverror('n_groups.ne.n_groups_saved')
   
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

end subroutine place_subvolume_into_tile

subroutine make_distance_ranges

   implicit none
   integer*4            :: i
   real*4,allocatable   :: d(:)
   
   allocate(d(para%snapshot_min:para%snapshot_max))
   do i = para%snapshot_min,para%snapshot_max
      d(i) = redshift_to_dc(snapshot(i)%redshift)*(unit%Mpc/para%length_unit)/para%box_side ! [box side-length] comoving distance to redshift of the box
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

subroutine write_sky_to_hdf5
   
      implicit none
      integer*4                           :: isubvolume
      type(type_sky_galaxy),allocatable   :: sky_galaxy(:)
      type(type_sky_group),allocatable    :: sky_group(:)
      character(255)                      :: filename
      character(255)                      :: str
      type(type_skystats)                 :: substats
      type(type_skystats)                 :: totstats
      
      call tic('WRITE SKY TO HDF5 FILE')
      
      if (para%merge_output==1) then
      
         ! count galaxies and groups
         do isubvolume = para%subvolume_min,para%subvolume_max
            open(fid,file=filename_sky_galaxies(isubvolume),action='read',form='unformatted',access='stream')
            read(fid) substats%n_galaxies
            close(fid)
            totstats%n_galaxies = totstats%n_galaxies+substats%n_galaxies
            if (para%make_groups==1) then
               open(fid,file=filename_sky_groups(isubvolume),action='read',form='unformatted',access='stream')
               read(fid) substats%n_groups
               close(fid)
               totstats%n_groups = totstats%n_groups+substats%n_groups
            end if
         end do
         
         ! load galaxies and groups
         allocate(sky_galaxy(totstats%n_galaxies))
         allocate(sky_group(totstats%n_groups))
         totstats%n_galaxies = 0
         totstats%n_groups = 0
         do isubvolume = para%subvolume_min,para%subvolume_max
         
            call out('Load sky subvolume ',isubvolume)
            
            ! galaxies
            open(fid,file=filename_sky_galaxies(isubvolume),action='read',form='unformatted',access='stream')
            read(fid) substats%n_galaxies,substats%n_distinct,substats%n_replica_max
            read(fid) sky_galaxy(totstats%n_galaxies+1:totstats%n_galaxies+substats%n_galaxies)
            close(fid)
            totstats%n_galaxies = totstats%n_galaxies+substats%n_galaxies
            totstats%n_distinct = totstats%n_distinct+substats%n_distinct
            totstats%n_replica_max = max(totstats%n_replica_max,substats%n_replica_max)
            
            ! groups
            if (para%make_groups==1) then
               open(fid,file=filename_sky_groups(isubvolume),action='read',form='unformatted',access='stream')
               read(fid) substats%n_groups
               read(fid) sky_group(totstats%n_groups+1:totstats%n_groups+substats%n_groups)
               close(fid)
               totstats%n_groups = totstats%n_groups+substats%n_groups
            end if
            
         end do
         
         ! write single new file
         filename = trim(para%path_output)//trim(para%filename_sky)//'.hdf5'
         call out('Write file ',trim(filename))
         call make_hdf5(trim(filename),sky_galaxy,sky_group,totstats)
         
         ! delete binary files
         do isubvolume = para%subvolume_min,para%subvolume_max
            if (para%keep_binaries==0) then
               call system('rm '//filename_sky_galaxies(isubvolume))
               if (para%make_groups==1) call system('rm '//filename_sky_groups(isubvolume))
            end if
         end do
      
      else
      
         do isubvolume = para%subvolume_min,para%subvolume_max
      
            ! user output
            write(filename,'(A,A,A,I0,A)') trim(para%path_output),trim(para%filename_sky),'.',isubvolume,'.hdf5'
            write(str,'(A,I0,A,A)') 'Write subvolume ',isubvolume,' to ',trim(filename)
            call out(str)
         
            ! load galaxies
            open(fid,file=filename_sky_galaxies(isubvolume),action='read',form='unformatted',access='stream')
            read(fid) substats%n_galaxies,substats%n_distinct,substats%n_replica_max
            allocate(sky_galaxy(totstats%n_galaxies+1:totstats%n_galaxies+substats%n_galaxies))
            read(fid) sky_galaxy
            close(fid)
            totstats%n_galaxies = totstats%n_galaxies+substats%n_galaxies
         
            ! load groups
            if (para%make_groups==1) then
               open(fid,file=filename_sky_groups(isubvolume),action='read',form='unformatted',access='stream')
               read(fid) substats%n_groups
               allocate(sky_group(totstats%n_groups+1:totstats%n_groups+substats%n_groups))
               read(fid) sky_group
               close(fid)
               totstats%n_groups = totstats%n_groups+substats%n_groups
            end if
         
            ! write new file
            call make_hdf5(trim(filename),sky_galaxy,sky_group,totstats,substats,isubvolume)
         
            ! free memory
            deallocate(sky_galaxy)
            deallocate(sky_group)
            
            ! delete binary files
            if (para%keep_binaries==0) then
               call system('rm '//filename_sky_galaxies(isubvolume))
               if (para%make_groups==1) call system('rm '//filename_sky_groups(isubvolume))
            end if
         
         end do
         
      end if
      
      call toc
   
   end subroutine write_sky_to_hdf5
   
end module module_sky