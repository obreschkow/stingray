module module_sky

   use shared_module_core
   use shared_module_constants
   use shared_module_parameters
   use shared_module_cosmology
   use shared_module_maths
   use shared_module_vectors
   use shared_module_hdf5
   use shared_module_sort
   use module_global
   use module_parameters
   use module_user_routines
   use module_selection_tools
   use module_user_selection
   use module_tiling
   
   public   :: make_sky             ! positions galaxies into mock sky and produces their properties
   public   :: write_sky_to_hdf5    ! call user routine write_hdf5() to save mock sky in HDF5 format
   
   private
   
   integer*4,parameter :: fid = 1 ! default file unit used in this module to write mock sky data
   
   type type_skystats
   
      integer*4   :: n_galaxies = 0    ! number of galaxies in mock sky (or some part thereof)
      integer*4   :: n_distinct = 0    ! number of distinct SAM galaxies in mock sky (or some part thereof)
      integer*4   :: n_replica_max = 0 ! maximum number of replications of the same SAM galaxy
      integer*4   :: n_groups = 0      ! number of groups
   
   end type type_skystats
   
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
   character(50)                       :: fmt
   character(255)                      :: str
   character(255)                      :: filename
   character(10)                       :: progress
   type(type_task),allocatable         :: task(:)
   type(type_skystats),allocatable     :: substats(:) ! statistics of mock sky for particular subvolumes
   type(type_skystats)                 :: totstats    ! statistics of total mock sky
   integer*4                           :: itask,n_tasks,n_stamps,i_stamps,stamp_id
   
   call tic('POSITION OBJECTS INTO SKY AND MAKE APPARENT PROPERTIES')
   
   ! make task list
   n_tasks = (para%snapshot_max-para%snapshot_min+1)*(para%subvolume_max-para%subvolume_min+1)
   allocate(task(n_tasks))
   n_tasks = 0
   n_stamps = 0
   do isnapshot = para%snapshot_min,para%snapshot_max
      if (snapshot(isnapshot)%n_tiles>0) then ! snapshot is used for at least some tiles
         do isubvolume = para%subvolume_min,para%subvolume_max
            n_tasks = n_tasks+1
            if (n_stamps==huge(n_stamps)) call error('Number of stamps larger than largest integer.')
            n_stamps = n_stamps+snapshot(isnapshot)%n_tiles
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
      write(fid) 0_4,0_4,0_4,0_4,0_4,0_4 ! place holder for number of galaxies and replication values
      close(fid)
      if (para%make_groups) then
         filename = filename_sky_groups(isubvolume)
         open(fid,file=trim(filename),action='write',form='unformatted',status='replace',access='stream')
         write(fid) 0_4,0_4 ! place holder for number of groups
         close(fid)
      end if
   end do
   
   ! fill galaxies with intrinsic properties into sky
   allocate(substats(para%subvolume_min:para%subvolume_max))
   substats%n_galaxies = 0
   substats%n_groups = 0
   substats%n_distinct = 0
   substats%n_replica_max = 0
   totstats%n_galaxies = 0
   totstats%n_groups = 0
   totstats%n_distinct = 0
   totstats%n_replica_max = 0
   
   fmt = '(A,A,I0.'//val2str(ceiling(log10(para%snapshot_max+1.0)))//',A,I0.'//&
       & val2str(max(1,ceiling(log10(para%subvolume_max+1.0))))// &
       & ',A,I0.'//val2str(ceiling(log10(size(tile)+1.0)))//',A,A,A,I0,A,I0,A)'
       
   i_stamps = 0
   
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
         
         !$OMP PARALLEL PRIVATE(sam_replica_tile,sky_galaxy,sky_group,filename,str,strt,stamp_id)
         !$OMP DO
         do itile = 1,size(tile)
            if ((snapshot(isnapshot)%dmax>=tile(itile)%dmin).and.(snapshot(isnapshot)%dmin<=tile(itile)%dmax)) then
         
               allocate(sam_replica_tile(size(sam)))
               sam_replica_tile = 0
               stamp_id = size(tile)*itask+itile ! unique, independent of threading
               call place_subvolume_into_tile(stamp_id,itile,isnapshot,isubvolume,sam,sam_sel,sam_replica_tile,sky_galaxy,sky_group)
               
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
               i_stamps = i_stamps+1
               write(progress,'(F7.3)') 100.0*real(i_stamps)/n_stamps
               write(str,fmt) trim(strt),': Place snapshot ',isnapshot,', subvolume ', &
                  & isubvolume,' into tile ',itile,' (',trim(adjustl(progress)), &
                  & '%, ',size(sky_galaxy),'/',totstats%n_galaxies,' galaxies)'
               call out(trim(str))
               
               ! save_sky_part, one thread at a time
               filename = filename_sky_galaxies(isubvolume)
               open(fid,file=trim(filename),action='write',status='old',position='append',form='unformatted',access='stream')
               write(fid) sky_galaxy
               close(fid)
               if (para%make_groups) then
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
         
         !$OMP CRITICAL 
         substats(isubvolume)%n_distinct = substats(isubvolume)%n_distinct+count(sam_replica>0)
         substats(isubvolume)%n_replica_max = max(substats(isubvolume)%n_replica_max,maxval(sam_replica))
         !$OMP END CRITICAL
         
         deallocate(sam_replica)
         
      end if
      
   end do ! task
   !$OMP END DO
   !$OMP END PARALLEL
   
   totstats%n_distinct = sum(substats%n_distinct)
   totstats%n_replica_max = maxval(substats%n_replica_max)
      
   ! finalize binary files
   do isubvolume = para%subvolume_min,para%subvolume_max
      filename = filename_sky_galaxies(isubvolume)
      open(fid,file=trim(filename),action='write',form='unformatted',status='old',access='stream')
      write(fid,pos=1) totstats%n_galaxies,totstats%n_distinct,totstats%n_replica_max, &
                       & substats(isubvolume)%n_galaxies,substats(isubvolume)%n_distinct,substats(isubvolume)%n_replica_max
      close(fid)
      if (para%make_groups) then
         filename = filename_sky_groups(isubvolume)
         open(fid,file=trim(filename),action='write',form='unformatted',status='old',access='stream')
         write(fid,pos=1) totstats%n_groups,substats(isubvolume)%n_groups
         close(fid)
      end if
   end do
   
   ! check number of galaxies
   if (totstats%n_galaxies==0) then
      call warning('No galaxy in sky. Consider widening the sky geometry or relaxing the selection criteria.')
   end if
   
   ! user output
   call out('Number of galaxies in mock sky: ',totstats%n_galaxies)
   if (para%make_groups) call out('Number of groups in mock sky: ',totstats%n_groups)
   call toc
   
end subroutine make_sky

function filename_sky_galaxies(isubvolume) result(fn)

   implicit none
   integer*4,intent(in) :: isubvolume
   character(255)       :: fn
   
   fn = dir(path_tmp,fn_galaxies//'.'//val2str(isubvolume))

end function filename_sky_galaxies

function filename_sky_groups(isubvolume) result(fn)

   implicit none
   integer*4,intent(in) :: isubvolume
   character(255)       :: fn
   
   fn = dir(path_tmp,fn_groups//'.'//val2str(isubvolume))

end function filename_sky_groups

subroutine preprocess_snapshot(sam,sam_sel)

   ! 1) determine for each galaxy, if it has been selected according to its selection_function(sam)
   ! 2) sort galaxies in 'sam' by groups, making the central galaxy the first one in its group
   ! 3) remove all groups, where no galaxy was selected
   
   implicit none
   type(type_sam),allocatable,intent(inout)  :: sam(:)
   logical*4,allocatable,intent(inout)       :: sam_sel(:)
   integer*8,allocatable                     :: id(:)
   integer*4,allocatable                     :: index(:)
   integer*4                                 :: i,n
   integer*4                                 :: n_groups
   integer*8                                 :: groupid
   integer*4                                 :: n_galaxies_in_group
   logical*4                                 :: group_selected
   integer*4                                 :: n_keep
   logical*4                                 :: selected
   
   n = size(sam)
   
   if (para%make_groups) then
   
      ! order galaxies by group id, such that central galaxies come first (if it exists)
      allocate(id(n),index(n))
      do i = 1,n
         id(i) = sam(i)%get_groupid()*2
         if (sam(i)%is_group_center()) then
            id(i)=id(i)-1 ! to make sure that this galaxy gets listed first
         end if
      end do
      call sort(id,index)
      sam = sam(index)
      deallocate(id,index)
      
      ! determine if the objects pass the SAM selection
      if (allocated(sam_sel)) deallocate(sam_sel)
      allocate(sam_sel(n))
      do i = 1,n
         sam_sel(i) = .true.
         call selection_function(sam=sam(i),selected=sam_sel(i))
      end do
            
      ! keep all galaxies of groups, where at least one object is sam-selected; reject all other galaxies
      allocate(index(n))
      n_groups = 0
      n_keep = 0
      i = 1
      do while (i<=n)
      
         ! start new group
         n_groups = n_groups+1
         groupid = sam(i)%get_groupid()
         
         ! count number of galaxies in group and number of selected galaxies in group
         group_selected = .false.
         n_galaxies_in_group = 0
         do while (i<=n.and.sam(min(i,n))%get_groupid()==groupid)
            if (n_galaxies_in_group>0) then
               if (sam(i)%is_group_center()) call error('each group can have at most one central object.')
            end if
            n_galaxies_in_group = n_galaxies_in_group+1
            group_selected = group_selected.or.sam_sel(i)
            index(n_keep+n_galaxies_in_group) = i
            i = i+1
         end do
         
         ! keep group, if at least one galaxy is selected
         if (group_selected) n_keep = n_keep+n_galaxies_in_group
         
      end do
      
      if (i/=n+1) call error('something wrong with group indexing.')
      
      ! apply reordering and selection to arrays sam(:) and sam_sel(:)
      sam = sam(index(1:n_keep))
      sam_sel = sam_sel(index(1:n_keep))
      deallocate(index)
      
   else
   
      n_keep = 0
      do i = 1,n
         selected = .true.
         call selection_function(sam=sam(i),selected=selected)
         if (selected) then
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

subroutine place_subvolume_into_tile(stamp_index,itile,isnapshot,isubvolume,sam,sam_sel,sam_replica,sky_galaxy_list,sky_group_list)

   ! NB: most of this routine deals with groups
   ! if make_groups==0, this routine essentially places all the galaxies in 'sam' into the sky,
   ! making sure to only retain the objects that pass the selection_function of the user module

   implicit none
   integer*4,intent(in)                            :: stamp_index
   integer*4,intent(in)                            :: itile
   integer*4,intent(in)                            :: isnapshot
   integer*4,intent(in)                            :: isubvolume
   integer*4                                       :: ishell
   type(type_sam),allocatable,intent(in)           :: sam(:)
   logical,allocatable,intent(in)                  :: sam_sel(:)
   integer*4,allocatable,intent(inout)             :: sam_replica(:)
   type(type_sky_galaxy),allocatable,intent(out)   :: sky_galaxy_list(:)
   type(type_sky_group),allocatable,intent(out)    :: sky_group_list(:)
   integer*8                                       :: n_galaxies
   integer*8                                       :: n_groups
   type(type_base)                                 :: base
   integer*4                                       :: nsam,d,isam,k
   logical                                         :: selected
   integer*8                                       :: prefixid
   type(type_sky_galaxy),allocatable               :: sky_galaxy(:)
   logical,allocatable                             :: sky_galaxy_selected(:)
   type(type_sky_group)                            :: sky_group
   logical                                         :: position_selected
   real*4                                          :: dx
   integer*4                                       :: group_flag
   real*4                                          :: xsam(3),xtile_min,xtile_max
   real*4,allocatable                              :: xtile(:,:)
   logical                                         :: tiling_mode
   logical                                         :: devoption_yzgrid
   
   type type_lim
      logical  :: survey = .false.
      logical  :: snapshot = .false.
      logical  :: tile = .false.
      logical  :: shell = .false.
   end type type_lim
   
   type type_group
      integer*8            :: id
      integer*4            :: n
      integer*4            :: n_selected = 0
      integer*4            :: first_index
      integer*4            :: last_index
      type(type_lim)       :: has_galaxy_outside = type_lim()
      type(type_lim)       :: is_at_edge_of = type_lim()
      type(type_base)      :: central_galaxy_base = type_base()
   end type type_group
   
   type(type_group)  :: group
   
   ! input checks
   nsam = size(sam)
   if (nsam>limit%n_galaxies_sky_max) call error('number of galaxies in mock sky potentially exceeds ',limit%n_galaxies_sky_max)
   if (itile>limit%n_tiles_max) call error('number of tiles exceeds ',limit%n_tiles_max)
   if (isubvolume>=limit%n_subvolumes_max) call error('number of subvolumes exceeds ',limit%n_subvolumes_max)
   if (isnapshot>=limit%n_snapshots_max) call error('number of snapshots exceeds ',limit%n_snapshots_max)
   
   ! initialise variables
   ishell = tile(itile)%shell
   n_galaxies = 0
   n_groups = 0
   !prefixid = int(limit%n_galaxies_per_tile_max,8)*(int(limit%n_subvolumes_max,8)*int(isnapshot,8)+int(isubvolume,8)) ! up to v0.37
   prefixid = int(limit%n_galaxies_per_tile_max,8)*int(stamp_index,8)
   
   tiling_mode = trim(para%randomisation)=='tiles'
   devoption_yzgrid = devoption('yzgrid')
   
   ! make fixed base properties
   base%index%shell = ishell
   base%index%tile = itile
   base%index%snapshot = isnapshot
   base%index%subvolume = isubvolume
   base%snapshot_redshift = snapshot(isnapshot)%redshift
   base%transformation%inverted = xor(shell(ishell)%transformation%inverted,tile(itile)%transformation%inverted)
   base%transformation%rotation = matmul(shell(ishell)%transformation%rotation,tile(itile)%transformation%rotation)
   
   ! allocate galaxy-arrays
   allocate(sky_galaxy_list(nsam),sky_group_list(nsam))

   ! iterate over all mock galaxies in the subvolume, group-by-group
   isam = 0
   do while (isam<nsam)
      
      ! initialise new group
      call reset_group(group,isam+1)
      if (allocated(sky_galaxy)) deallocate(sky_galaxy)
      allocate(sky_galaxy(group%n))
      if (allocated(sky_galaxy_selected)) deallocate(sky_galaxy_selected)
      allocate(sky_galaxy_selected(group%n))
      if (allocated(xtile)) deallocate(xtile)
      allocate(xtile(group%n,3))
      base%index%group = group%id
      
      ! make galaxy positions in tile
      do k = 1,group%n
      
         ! evaluate normalised position of object in sam-output
         xsam = sam(isam+k)%get_position()/para%box_side
         if (devoption_yzgrid) then
            if (mod(sam(isam+k)%get_groupid(),11)==0_8) then
               xsam(2) = (xsam(2)-0.5)*0.4+0.65
               xsam(3) = 0.5
            else if (mod(sam(isam+k)%get_groupid(),11)==1_8) then
               xsam(3) = (xsam(3)-0.5)*0.4+0.65
               xsam(2) = 0.5
            end if
            xsam(1) = 0.5
         end if
      
         ! apply tile-symmetry to normalised coordinates
         xtile(k,:) = apply_tile_symmetry(xsam,itile)
         
      end do
      
      ! handle wrapped groups, unwrap if appropriate
      group%is_at_edge_of%tile = .false.
      do d = 1,3
         xtile_min = minval(xtile(:,d))
         xtile_max = maxval(xtile(:,d))
         dx = xtile_max-xtile_min
         if ((dx>limit%group_diameter_max).and.(dx<1.0-limit%group_diameter_max)) then
            if (.not.devoption_yzgrid) call error('group wider than box side times ',limit%group_diameter_max)
         end if
         if (dx>0.5) then ! group is wrapped
            ! translate all members to side of the first galaxy in the group (which is the central if it exists)
            xtile(:,d) = modulo(xtile(:,d)+0.5,1.0)+xtile(1,d)-modulo(xtile(1,d)+0.5,1.0)
            group%is_at_edge_of%tile = tiling_mode
         end if
      end do
      
      ! map galaxies onto sky and check full selection
      do k = 1,group%n
      
         selected = .true.
         
         ! remove galaxies outside tile (after unwrapping) in tiling-based randomisation mode
         if ((group%is_at_edge_of%tile).and.tiling_mode) then
            if ((minval(xtile(k,:))<0.0).or.(maxval(xtile(k,:))<1.0)) then
               selected = .false.
            end if
         end if
         
         ! make sky coordinates
         base%cartesian = map_tile_onto_sky(xtile(k,:),itile)
         base%transformation%translation = components(base%cartesian)-xtile(k,:)
         call car2sph(base%cartesian,base%spherical)
         
         ! check if distance is in the range covered by snapshot isnapshot
         if ((base%spherical%dc<snapshot(isnapshot)%dmin).or.(base%spherical%dc>=snapshot(isnapshot)%dmax).or.&
            & (base%spherical%dc<=epsilon(0.0))) then ! exclude objects that exactly coincide with the observer
            group%has_galaxy_outside%snapshot = .true.
            selected = .false.
         end if
      
         ! check if distance is in the range covered by shell
         if ((base%spherical%dc<shell(ishell)%dmin).or.(base%spherical%dc>=shell(ishell)%dmax)) then
            group%has_galaxy_outside%shell = .true.
            selected = .false.
         end if
   
         ! check full position
         position_selected = is_in_fov(base%spherical)
         if (position_selected) then
            call selection_function(pos=sph_deg_l(base%spherical),selected=position_selected)
         end if
         if (.not.position_selected) then
            group%has_galaxy_outside%survey = .true.
            selected = .false.
         end if
      
         ! check sam(+pos(+sky))
         selected = selected.and.sam_sel(isam+k)
         if (selected) then
            call selection_function(pos=sph_deg_l(base%spherical),sam=sam(isam+k),selected=selected)
            if (selected) then
               base%index%galaxy = prefixid+n_galaxies+1
               call make_sky_object(sky_galaxy(k),sam(isam+k),convert_units(base))
               call make_sky_galaxy(sky_galaxy(k),sam(isam+k),convert_units(base))
               call selection_function(pos=sph_deg_l(base%spherical),sam=sam(isam+k),sky=sky_galaxy(k),selected=selected)
            end if
         end if
   
         ! save / count selected galaxy
         if (selected) then
            n_galaxies = n_galaxies+1
            group%n_selected = group%n_selected+1
            sam_replica(isam+k) = sam_replica(isam+k)+1
            sky_galaxy_list(n_galaxies) = sky_galaxy(k)
         else
            base%index%galaxy = -1 ! just to make sure that the id of the central group galaxy is -1, if this galaxy is not selected 
         end if
   
         ! save some properties for group
         sky_galaxy_selected(k) = selected
         if (k==1) group%central_galaxy_base = base ! base-data of first galaxy, which is the group's central (if it exists)
         
      end do
   
      ! make group
      if ((para%make_groups).and.(group%n>1).and.(group%n_selected>0)) then ! by choice, we only consider a group as such, if it has at least 2 members
      
         n_groups = n_groups+1

         ! make truncation flag
         group%is_at_edge_of%survey = group%has_galaxy_outside%survey
         group%is_at_edge_of%snapshot = group%has_galaxy_outside%snapshot
         group%is_at_edge_of%shell = group%has_galaxy_outside%shell
         group_flag = log2int(group%is_at_edge_of%survey)+2*log2int(group%is_at_edge_of%snapshot)+ &
         & 4*log2int(group%is_at_edge_of%tile)+8*log2int(group%is_at_edge_of%shell)
      
         ! make group
         call make_sky_object(sky_group,sam(group%first_index),convert_units(group%central_galaxy_base))
         call make_sky_group(sky_group,sam(group%first_index:group%last_index),sky_galaxy, &
         & sky_galaxy_selected,convert_units(group%central_galaxy_base),group_flag)

         ! save group
         sky_group_list(n_groups) = sky_group
   
      end if
      
      isam = isam+group%n
   
   end do
   
   sky_galaxy_list = sky_galaxy_list(1:n_galaxies)
   sky_group_list = sky_group_list(1:n_groups)
   
contains

   subroutine reset_group(group,first_index)
   
      implicit none
      type(type_group),intent(inout)   :: group
      integer*4,intent(in)             :: first_index
      integer*4                        :: i
      logical                          :: last_galaxy_in_group
      integer*8                        :: id
      
      do i = first_index,nsam
         if (i==nsam) then
            last_galaxy_in_group = .true.
         else
            last_galaxy_in_group = sam(i)%get_groupid()/=sam(i+1)%get_groupid()
         end if
         if (last_galaxy_in_group) exit
      end do
      
      if ((para%make_groups).and.(i>first_index)) then ! group has at least 2 members
         id = prefixid+n_groups+1
      else
         id = -1_8
      end if
      
      group = type_group(first_index=first_index,last_index=i,n=i-first_index+1,id=id)
      
   end subroutine reset_group
   
   function convert_units(bin) result(bout)

      ! converts rad to deg
   
      implicit none
      type(type_base),intent(in) :: bin
      type(type_base) :: bout
   
      bout = bin
      bout%spherical%dc = bin%spherical%dc*para%box_side
      bout%cartesian = bout%cartesian*para%box_side
      
   end function convert_units

end subroutine place_subvolume_into_tile

subroutine write_sky_to_hdf5
   
   implicit none
   integer*4                           :: isubvolume,i_galaxy,i_group
   type(type_sky_galaxy),allocatable   :: sky_galaxy(:)
   type(type_sky_group),allocatable    :: sky_group(:)
   character(255)                      :: filename
   character(255)                      :: str
   type(type_skystats)                 :: substats
   type(type_skystats)                 :: totstats
   
   call tic('WRITE SKY TO HDF5 FILE')
   
   i_galaxy = 0
   i_group = 0
   
   if (para%merge_output) then
   
      ! count galaxies and groups
      isubvolume = para%subvolume_min
      open(fid,file=filename_sky_galaxies(isubvolume),action='read',form='unformatted',access='stream')
      read(fid) totstats%n_galaxies
      close(fid)
      if (para%make_groups) then
         open(fid,file=filename_sky_groups(isubvolume),action='read',form='unformatted',access='stream')
         read(fid) totstats%n_groups
         close(fid)
      end if
      
      ! load galaxies and groups
      allocate(sky_galaxy(totstats%n_galaxies))
      allocate(sky_group(totstats%n_groups))
      do isubvolume = para%subvolume_min,para%subvolume_max
      
         call out('Process subvolume ',isubvolume)
         
         ! galaxies
         open(fid,file=filename_sky_galaxies(isubvolume),action='read',form='unformatted',access='stream')
         read(fid) totstats%n_galaxies,totstats%n_distinct,totstats%n_replica_max
         read(fid) substats%n_galaxies,substats%n_distinct,substats%n_replica_max
         read(fid) sky_galaxy(i_galaxy+1:i_galaxy+substats%n_galaxies)
         close(fid)
         i_galaxy = i_galaxy+substats%n_galaxies
         
         ! groups
         if (para%make_groups) then
            open(fid,file=filename_sky_groups(isubvolume),action='read',form='unformatted',access='stream')
            read(fid) totstats%n_groups,substats%n_groups
            read(fid) sky_group(i_group+1:i_group+substats%n_groups)
            close(fid)
            i_group = i_group+substats%n_groups
         end if
         
      end do
      
      ! write single new file
      filename = trim(para%path_output)//trim(para%filename_sky)//'.hdf5'
      call out('Write file ',trim(filename))
      call initialize_hdf5(trim(filename),totstats)
      call write_hdf5_parameters(trim(filename))
      call write_hdf5_mapping(trim(filename))
      call write_hdf5(trim(filename),sky_galaxy,sky_group)
   
   else
   
      do isubvolume = para%subvolume_min,para%subvolume_max
   
         ! user output
         write(filename,'(A,A,A,I0,A)') trim(para%path_output),trim(para%filename_sky),'.',isubvolume,'.hdf5'
         write(str,'(A,I0,A,A)') 'Write subvolume ',isubvolume,' to ',trim(filename)
         call out(str)
      
         ! load galaxies
         open(fid,file=filename_sky_galaxies(isubvolume),action='read',form='unformatted',access='stream')
         read(fid) totstats%n_galaxies,totstats%n_distinct,totstats%n_replica_max
         read(fid) substats%n_galaxies,substats%n_distinct,substats%n_replica_max
         allocate(sky_galaxy(i_galaxy+1:i_galaxy+substats%n_galaxies))
         read(fid) sky_galaxy
         close(fid)
         i_galaxy = i_galaxy+substats%n_galaxies
      
         ! load groups
         if (para%make_groups) then
            open(fid,file=filename_sky_groups(isubvolume),action='read',form='unformatted',access='stream')
            read(fid) totstats%n_groups,substats%n_groups
            allocate(sky_group(i_group+1:i_group+substats%n_groups))
            read(fid) sky_group
            close(fid)
            i_group = i_group+substats%n_groups
         end if
      
         ! write new file
         call initialize_hdf5(trim(filename),totstats,substats,isubvolume)
         call write_hdf5_parameters(trim(filename))
         call write_hdf5_mapping(trim(filename))
         call write_hdf5(trim(filename),sky_galaxy,sky_group)
      
         ! free memory
         deallocate(sky_galaxy)
         deallocate(sky_group)
      
      end do
      
   end if
   
   if (i_galaxy/=totstats%n_galaxies) call deverror('i_galaxy/=totstats%n_galaxies')
   if ((para%make_groups).and.(i_group/=totstats%n_groups)) call deverror('i_group/=totstats%n_groups')
   
   ! user output
   call out('Number of galaxies in mock sky: ',totstats%n_galaxies)
   if (para%make_groups) call out('Number of groups in mock sky: ',totstats%n_groups)
   call toc

end subroutine write_sky_to_hdf5

subroutine initialize_hdf5(filename_hdf5,totstats,substats,isubvolume)

   implicit none
   character(*),intent(in)                      :: filename_hdf5  ! output filename
   type(type_skystats),intent(in)               :: totstats
   type(type_skystats),intent(in),optional      :: substats ! only provided if the sky corresponds to a subvolume
   integer*4,intent(in),optional                :: isubvolume
   character(:),allocatable                     :: name
   
   allocate(character(1)::name) ! empty allocation to avoid compiler flags

  ! create and open HDF5 file
   call hdf5_create(filename_hdf5)
   call hdf5_open(filename_hdf5,.true.)
   
   ! write group "run_info"
   name = 'run_info/'
   call hdf5_add_group(name)
   call hdf5_write_data(name//'stingray_version',version,'Version of Stingray used to produce this mock sky')
   call hdf5_write_data(name//'stingray_timestamp',timestamp(),'Time at which this mockfile was written')
   
   ! write group "statistics"
   name = 'statistics/'
   call hdf5_add_group(name)
   call hdf5_write_data(name//'n_galaxies',totstats%n_galaxies,'Number of galaxies')
   call hdf5_write_data(name//'n_groups',totstats%n_groups,'Number of groups')
   call hdf5_write_data(name//'n_replica_mean',real(totstats%n_galaxies,4).safedivide.real(totstats%n_distinct,4),&
   & 'Mean number a galaxy was replicated)')
   call hdf5_write_data(name//'n_replica_max',totstats%n_replica_max,'Maximum number a galaxy was replicated')
   if (present(substats).and.present(isubvolume)) then
      call hdf5_write_data(name//'subvolume_index',isubvolume,'Index of the SAM subvolume used in this file')
      call hdf5_write_data(name//'subvolume_n_galaxies',substats%n_galaxies,'Number of galaxies in this subvolume')
      call hdf5_write_data(name//'subvolume_n_groups',substats%n_groups,'Number of groups in this subvolume')
      call hdf5_write_data(name//'subvolume_n_replica_mean',real(substats%n_galaxies,4).safedivide.real(substats%n_distinct,4),&
      & 'Mean number a galaxy was replicated in this subvolume)')
      call hdf5_write_data(name//'subvolume_n_replica_max',substats%n_replica_max,&
      &'Maximum number a galaxy was replicated in this subvolume')   
   end if

   ! close HDF5 file
   call hdf5_close()

end subroutine initialize_hdf5
   
end module module_sky