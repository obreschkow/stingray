module module_tiling

   use shared_module_core
   use shared_module_parameters
   use shared_module_maths
   use shared_module_vectors
   use shared_module_constants
   use shared_module_cosmology
   use shared_module_hdf5
   use module_global
   use module_parameters
   use module_user_routines
   use module_user_selection
   use module_selection_tools
   
   public   :: make_tiling
   public   :: is_in_fov
   public   :: apply_tile_symmetry
   public   :: map_tile_onto_sky
   public   :: sph_deg_l
   public   :: write_hdf5_mapping
   public   :: snapshot,tile,shell ! read-only
   
   private

   type type_snapshot

      real*4      :: redshift
      real*4      :: dmin           ! [box side length] minimum comoving distance at which galaxies are drawn from this snapshot
      real*4      :: dmax           ! [box side length] maximum comoving distance at which galaxies are drawn from this snapshot
      integer*4   :: n_tiles        ! Number of tiles this snapshot has been considered for, irrespective of whether a galaxy was selected

   end type type_snapshot
   
   type type_shell
   
      real*4      :: dmin        ! [box side length] minimum comoving distance to observer
      real*4      :: dmax        ! [box side length] maximum comoving distance to observer
      type(type_transformation)  :: transformation
      
   end type type_shell
   
   type type_tile
   
      integer*4   :: shell          ! shell index
      integer*4   :: ix(3)          ! integer position, where ix=(0,0,0) is the central box with the observer in the middle
      real*4      :: dmin           ! [box side length] minimum comoving distance to observer
      real*4      :: dmax           ! [box side length] maximum comoving distance to observer
      type(type_transformation)  :: transformation
      
   end type type_tile
   
   type(type_tile),allocatable,protected         :: tile(:)
   type(type_shell),allocatable,protected        :: shell(:)
   type(type_snapshot),allocatable,protected     :: snapshot(:)
   
   type(type_fov),protected   :: fov   ! survey range in [box side length], [rad]
   
   integer*4,allocatable      :: intersection(:,:,:,:) ! ==  0 : not checked
                                                      ! == -1 : does not intersect
                                                      ! >= +1 : id of intersecting tile
   integer*4                  :: imax
   integer*4                  :: ntiles
   integer*4                  :: nchecks
   
contains

subroutine make_tiling

   implicit none
   integer*4   :: starting_point(3)
   integer*4   :: ishell,isnapshot
   
   call tic
   call out('MAKE 3D TILING')
   
   call set_seed(para%seed)
   call get_position_range
   
   imax = ceiling(fov%dc(2))
   if (abs(imax)>limit%n_tiles_max) call error('maximum comoving distance too large for this box side length')
   
   ! make snapshot properties
   allocate(snapshot(para%snapshot_min:para%snapshot_max))
   do isnapshot = para%snapshot_min,para%snapshot_max
      snapshot(isnapshot)%redshift = get_redshift(isnapshot)
   end do
   call make_distance_ranges
   
   ! make shells
   call make_shell_list
   call out('Number of shells = ',size(shell))
   
   ! make tiles
   ntiles = 0
   nchecks = 0
   allocate(intersection(size(shell),-imax:imax,-imax:imax,-imax:imax))
   intersection = 0
   do ishell = 1,size(shell)
      call make_starting_point(ishell,starting_point)
      call check_tile(ishell,starting_point,0)
   end do
   call make_tile_list
   call count_tiles_of_snapshot
   call out('Number of tiles = ',size(tile))
   call out('Number of points checked = ',nchecks)
   
   call toc
   
end subroutine make_tiling

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
         snapshot(i)%dmax = huge(0.0_4)
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

logical function is_in_fov(sph)

   implicit none
   
   type(type_spherical),intent(in)  :: sph   ! spherical coordinates in rad and box side lengths
   
   if (fov%ra(1)<fov%ra(2)) then
   
      is_in_fov = (sph%dc>=fov%dc(1)).and.(sph%dc<=fov%dc(2)).and. &
                & (sph%ra>=fov%ra(1)).and.(sph%ra<=fov%ra(2)).and. &
                & (sph%dec>=fov%dec(1)).and.(sph%dec<=fov%dec(2))
                
   else
   
      is_in_fov = (sph%dc>=fov%dc(1)).and.(sph%dc<=fov%dc(2)).and. &
                & (sph%ra>=fov%ra(1)).or.(sph%ra<=fov%ra(2)).and. &
                & (sph%dec>=fov%dec(1)).and.(sph%dec<=fov%dec(2))
      
   end if
             
end function is_in_fov

subroutine get_position_range

   implicit none
   real*4,parameter  :: x = huge(0.0_4)
   
   ! set default values
   fov%dc = x
   fov%ra = x
   fov%dec = x
   
   ! extract ranges from selection function
   call selection_function(range=fov)
   
   ! check initialization
   if (any(fov%dc==x)) call error('range of comoving distances (range%dc) not provided in selection function')
   if (any(fov%ra==x)) call error('range of right ascension (range%ra) not provided in selection function')
   if (any(fov%dec==x)) call error('range of declination (range%dec) not provided in selection function')
   
   ! check values
   if (fov%dc(1)<0.0) call error('dc_min must be >=0')
   if (fov%dc(2)<=0.0) call error('dc_max must be >0')
   if (fov%dc(2)<=fov%dc(1)) call error('dc_min must be smaller than dc_max')
   
   if (fov%ra(1)<0.0) call error('ra_min must be >=0')
   if (fov%ra(2)>360.0) call error('ra_max must be <=360')
   if (fov%ra(2)==fov%ra(1)) call error('ra_min cannot be equal to ra_max')
   
   if (fov%dec(1)<-90.0) call error('dec_min must be >=-90')
   if (fov%dec(2)>+90.0) call error('dec_max must be <=90')
   if (fov%dec(2)<=fov%dec(1)) call error('dec_min must be smaller than dec_max')
   
   ! convert degrees to radian
   fov%ra = fov%ra*unit%degree
   fov%dec = fov%dec*unit%degree
   
   ! convert simulation units to box side lengths
   fov%dc = fov%dc/para%box_side

end subroutine get_position_range

logical function is_tile_in_shell(ishell,ix)

   implicit none
   integer*4,intent(in) :: ishell
   integer*4,intent(in) :: ix(3)
   real*4               :: xmin(3),rmin
   real*4               :: xmax(3),rmax
   integer*4            :: d
   
   do d = 1,3
      xmin(d) = max(0.0,real(abs(ix(d)),4)-0.5)
      xmax(d) = real(abs(ix(d)),4)+0.5
   end do
   rmin = norm(xmin)
   rmax = norm(xmax)
   is_tile_in_shell = (rmin<shell(ishell)%dmax).and.(rmax>shell(ishell)%dmin)
            
end function is_tile_in_shell

subroutine make_starting_point(ishell,ix)

   implicit none
   integer*4,intent(in)    :: ishell
   integer*4,intent(out)   :: ix(3)
   integer*4               :: i,j,k
   
   ix = (/0,0,0/)
   if (is_tile_in_shell(ishell,ix).and.is_tile_in_survey(ishell,ix,.false.).and.is_tile_in_survey(ishell,ix,.true.)) return
   
   do i = -imax,imax !xxx inside out might be faster
      do j = -imax,imax
         do k = -imax,imax
            ix = (/i,j,k/)
            if (is_tile_in_shell(ishell,ix).and.is_tile_in_survey(ishell,ix,.false.).and.is_tile_in_survey(ishell,ix,.true.)) return
         end do
      end do
   end do
      
   call error('no object passes position selection criterion')
   
end subroutine make_starting_point

recursive subroutine check_tile(ishell,ix,nrecursion)

   implicit none
   integer*4,intent(in) :: ishell
   integer*4,intent(in) :: ix(3)
   integer*4,intent(in) :: nrecursion
   
   if (nrecursion>1e4) call error('too many recursions in tiling; please restrict the ranges of dc, ra, dec in selection_function')
   
   if (maxval(abs(ix))>imax) return
   if (intersection(ishell,ix(1),ix(2),ix(3)).ne.0) return
   
   if (is_tile_in_survey(ishell,ix,.false.).and.is_tile_in_shell(ishell,ix)) then
      
      ! tile intersects with global survey volume and shell
   
      if (is_tile_in_survey(ishell,ix,.true.)) then
         ! tile also intersects with detailed survey volume
         ntiles = ntiles+1
         if (ntiles>limit%n_tiles_max) call error ('number of required tiles exceeds ',limit%n_tiles_max)
         intersection(ishell,ix(1),ix(2),ix(3)) = ntiles
      else
         ! tile does not intersect with survey volume specified by user file
         intersection(ishell,ix(1),ix(2),ix(3)) = -1
      end if
      
      ! check neighboring tiles
      call check_tile(ishell,ix+(/+1,0,0/),nrecursion+1)
      call check_tile(ishell,ix+(/-1,0,0/),nrecursion+1)
      call check_tile(ishell,ix+(/0,+1,0/),nrecursion+1)
      call check_tile(ishell,ix+(/0,-1,0/),nrecursion+1)
      call check_tile(ishell,ix+(/0,0,+1/),nrecursion+1)
      call check_tile(ishell,ix+(/0,0,-1/),nrecursion+1)
      
   else
   
      ! tile does not intersect with survey volume specified by parameter file
      intersection(ishell,ix(1),ix(2),ix(3)) = -1
      
   end if

end subroutine check_tile

subroutine make_tile_list

   implicit none
   real*4                     :: dmin,dmax
   integer*4                  :: i,j,k,itile,ishell
   
   if (ntiles/=count(intersection>0)) call deverror('tile number mismatch 1')
   if (ntiles/=maxval(intersection)) call deverror('tile number mismatch 2')
   if (ntiles==0) call error('no galaxy passes the selection criterion')
   
   allocate(tile(ntiles))
   
   do ishell = 1,size(shell)
   
      do i = -imax,imax
         do j = -imax,imax
            do k = -imax,imax
               if (intersection(ishell,i,j,k)>0) then
      
                  itile = intersection(ishell,i,j,k)
         
                  ! position
                  tile(itile)%ix = (/i,j,k/)
      
                  ! distance range of tile
                  call get_distance_range_of_cube(tile(itile)%ix,dmin,dmax)
                  tile(itile)%dmin = max(fov%dc(1),dmin)
                  tile(itile)%dmax = min(fov%dc(2),dmax)
               
                  ! shell index
                  tile(itile)%shell = ishell
               
                  ! random transformation
                  if (trim(para%randomisation)=='tiles') then
                     call assign_random_transformation(tile(itile)%transformation,.true.)
                     if ((para%fix_observer_position).and.(i==0.and.j==0.and.k==0)) then
                        tile(itile)%transformation%translation = 0.0
                        tile(itile)%transformation%inverted = .false.
                     end if
                     if ((para%fix_observer_rotation).and.(i==0.and.j==0.and.k==0)) tile(itile)%transformation%rotation = 0.0
                  end if
                  
                  ! NB: for the other randomisation schemes ("single" and "shell"), the symmetry operations (translations,
                  !     rotations and inversion) are specified at the level of the shell(s). See "make_shell_list" below.
                  !     In this case the tile transformations will just be identities as specified in the definition
                  !     of type_transformation (see module_global).
               
               end if
            end do
         end do
      end do
   end do

end subroutine make_tile_list

subroutine make_shell_list

   implicit none
   integer*4         :: ishell,nshell,isnapshot
   real*4            :: dr_first          ! thickness/radius of inner most shell
   real*4            :: r
   real*4,parameter  :: dtweak = 0.25
   
   if (trim(para%randomisation)=='shells') then
      
      if (fov%dc(1)<0.5) then
         dr_first = 0.5
      else
         dr_first = 1.0
      end if
   
      nshell = ceiling(fov%dc(2)-fov%dc(1)-dr_first)+1
      
      if (nshell<1) call error('distance range of selection function is too narrow')
   
      allocate(shell(nshell))
      
      r = fov%dc(1)
   
      do ishell = 1,nshell
   
         ! distance range of shell, interspaced by one box side length
         shell(ishell)%dmin = r
         if (ishell==1) then
            shell(ishell)%dmax = shell(ishell)%dmin+dr_first
         else
            shell(ishell)%dmax = shell(ishell)%dmin+1.0
         end if
         r = shell(ishell)%dmax
         
         ! truncate last shell to maximum distance
         if (shell(ishell)%dmax>fov%dc(2)) then
            if (ishell/=nshell) call deverror('shell too large')
            shell(ishell)%dmax = fov%dc(2)
         end if
      
         ! random transformation
         call assign_random_transformation(shell(ishell)%transformation,.false.)
      
      end do
      
      if (para%shell_tweaking) then
      
         ! if radius of closest snapshot boundary is less than dtweak away, change shell radius to snapshot boundary
         
         do ishell = 1,nshell-1
         
            if (minval(abs(snapshot%dmax-shell(ishell)%dmax))<dtweak) then
               isnapshot = minloc(abs(snapshot%dmax-shell(ishell)%dmax),1)+para%snapshot_min-1
               shell(ishell)%dmax = snapshot(isnapshot)%dmax
               shell(ishell+1)%dmin = snapshot(isnapshot)%dmax
            end if
         
         end do
         
      end if
      
   else ! para%randomisation is either 'single' or 'tiles'
   
      allocate(shell(1))
      
      ! distance range of shell
      shell(1)%dmin = fov%dc(1)
      shell(1)%dmax = fov%dc(2)
      
      call assign_random_transformation(shell(1)%transformation,.false.)
   
   end if
   
   if (shell(1)%dmin<0.5) then
   
      if (para%fix_observer_position) then
         shell(1)%transformation%translation = para%observer_translation
         shell(1)%transformation%inverted = .false.
      end if
      if (para%fix_observer_rotation) shell(1)%transformation%rotation = para%observer_rotation
      
   end if
      
end subroutine make_shell_list

subroutine assign_random_transformation(tr,discrete_rotation)

   ! NB: random_number() uses a better prng than rand(), but depends on the system and fortran version
   !     we therefore chose to use rand(), accessed with the option modern = .false.

   implicit none
   type(type_transformation),intent(inout)  :: tr
   logical,intent(in) :: discrete_rotation
   integer*4 :: k,d
   real*4   :: axis(3)
   real*4   :: angle
   real*4,parameter :: rxx(24) = real((/+1,+1,+1,+1,-1,-1,-1,-1,+0,+0,+0,+0,+0,+0,+0,+0,+0,+0,+0,+0,+0,+0,+0,+0/),4)
   real*4,parameter :: rxy(24) = real((/+0,+0,+0,+0,+0,+0,+0,+0,+1,-1,+0,+0,+1,-1,+0,+0,+1,-1,+0,+0,+1,-1,+0,+0/),4)
   real*4,parameter :: rxz(24) = real((/+0,+0,+0,+0,+0,+0,+0,+0,+0,+0,+1,-1,+0,+0,-1,+1,+0,+0,-1,+1,+0,+0,+1,-1/),4)
   real*4,parameter :: ryx(24) = real((/+0,+0,+0,+0,+0,+0,+0,+0,+1,+1,+1,+1,-1,-1,-1,-1,+0,+0,+0,+0,+0,+0,+0,+0/),4)
   real*4,parameter :: ryy(24) = real((/+1,-1,+0,+0,+1,-1,+0,+0,+0,+0,+0,+0,+0,+0,+0,+0,+0,+0,+1,-1,+0,+0,+1,-1/),4)
   real*4,parameter :: ryz(24) = real((/+0,+0,-1,+1,+0,+0,+1,-1,+0,+0,+0,+0,+0,+0,+0,+0,+1,-1,+0,+0,-1,+1,+0,+0/),4)
   real*4,parameter :: rzx(24) = real((/+0,+0,+0,+0,+0,+0,+0,+0,+0,+0,+0,+0,+0,+0,+0,+0,+1,+1,+1,+1,-1,-1,-1,-1/),4)
   real*4,parameter :: rzy(24) = real((/+0,+0,+1,-1,+0,+0,+1,-1,+0,+0,+1,-1,+0,+0,+1,-1,+0,+0,+0,+0,+0,+0,+0,+0/),4)
   real*4,parameter :: rzz(24) = real((/+1,-1,+0,+0,-1,+1,+0,+0,-1,+1,+0,+0,+1,-1,+0,+0,+0,+0,+0,+0,+0,+0,+0,+0/),4)
   
   ! choose random proper rotation
   if (para%rotate) then
      if (discrete_rotation) then
         if (devoption('xrotation')) then
            k = get_random_integer(1,4,modern=para%modern_prng)
         else
            k = get_random_integer(1,24,modern=para%modern_prng)
         end if
         tr%rotation = reshape((/rxx(k),rxy(k),rxz(k),ryx(k),ryy(k),ryz(k),rzx(k),rzy(k),rzz(k)/),(/3,3/))
      else
         if (devoption('xrotation')) then
            axis = (/1.0,0.0,0.0/)
         else
            axis = get_random_unit_vector(modern=para%modern_prng)
         end if
         angle = get_random_uniform_number(0.0,2*pi,modern=para%modern_prng)
         tr%rotation = rotation_matrix(axis,angle)
      end if
   else
      tr%rotation = const%identity3
   end if

   ! choose random inversion
   if (para%invert) then
      tr%inverted = get_random_uniform_number(0.0,1.0,modern=para%modern_prng)>0.5
   else
      tr%inverted = .false.
   end if

   ! choose random translation
   if (para%translate) then
      do d = 1,3
         tr%translation(d) = get_random_uniform_number(0.0,1.0,modern=para%modern_prng)
      end do
      if (devoption('yztranslation')) tr%translation(1) = 0.0
   else
      tr%translation = 0
   end if
   
end subroutine assign_random_transformation

subroutine get_distance_range_of_cube(ix,dmin,dmax)

   implicit none
   integer*4,intent(in) :: ix(3)
   real*4,intent(out)   :: dmin,dmax
   real*4               :: dmin_sq(6)
   integer*4            :: i,j,k
   
   if (maxval(abs(ix))==0) then
      dmin = 0
   else
      dmin_sq(1) = get_min_distance_to_square(ix+(/-0.5,0.0,0.0/),(/0.0,1.0,0.0/),(/0.0,0.0,1.0/))
      dmin_sq(2) = get_min_distance_to_square(ix+(/+0.5,0.0,0.0/),(/0.0,1.0,0.0/),(/0.0,0.0,1.0/))
      dmin_sq(3) = get_min_distance_to_square(ix+(/0.5,-0.5,0.0/),(/1.0,0.0,0.0/),(/0.0,0.0,1.0/))
      dmin_sq(4) = get_min_distance_to_square(ix+(/0.5,+0.5,0.0/),(/1.0,0.0,0.0/),(/0.0,0.0,1.0/))
      dmin_sq(5) = get_min_distance_to_square(ix+(/0.0,0.0,-0.5/),(/1.0,0.0,0.0/),(/0.0,1.0,0.0/))
      dmin_sq(6) = get_min_distance_to_square(ix+(/0.0,0.0,+0.5/),(/1.0,0.0,0.0/),(/0.0,1.0,0.0/))
      dmin = minval(dmin_sq)
   end if
   
   dmax = 0
   do i = -1,1,2
      do j = -1,1,2
         do k = -1,1,2
            dmax = max(dmax,norm(real(ix)+(/i*0.5,j*0.5,k*0.5/)))
         end do
      end do
   end do
               
end subroutine get_distance_range_of_cube

function get_min_distance_to_square(p,e1,e2) result(dmin)
   
   implicit none
   real*4,intent(in)    :: p(3)        ! center of square
   real*4,intent(in)    :: e1(3),e2(3) ! orthonormal basis in the plane of the square, parallel to the edges
   real*4               :: dmin        ! minimal distance from the origin (0,0,0) to the square
   real*4               :: p1,p2       ! projections of p onto square
   real*4               :: s1,s2
   
   ! project p onto square
   p1 = sum(p*e1)
   p2 = sum(p*e2)
   s1 = sign(1.0,p1) ! sign of p1
   s2 = sign(1.0,p2) ! sign of p2
   
   ! handle the situation where projection lands inside the square => closest point lies inside square
   if (max(abs(p1),abs(p2))<=0.5) then
      dmin = sqrt(sum(p**2)-sum(p*e1)**2-sum(p*e2)**2) ! equation for the minimum distance to a plane
      return
   end if
   
   ! handle the situation where the closes point lies along an edge
   if ((min(abs(p1),abs(p2))<=0.5).and.(max(abs(p1),abs(p2))>0.5)) then
      if (abs(p1)>0.5) then
         dmin = get_min_distance_to_line(p-s1*e1/2,e2)
         return
      else if (abs(p2)>0.5) then
         dmin = get_min_distance_to_line(p-s2*e2/2,e1)
         return
      end if
   end if
   
   ! handle the situation where the closest point is a corner
   if (min(abs(p1),abs(p2))>0.5) then
      dmin = norm(p-s1*e1/2-s2*e2/2)
      return
   else
      call deverror('unknown error in get_min_distance_to_square')
   end if
   
   ! code never gets here, this is just to avoid compiler warning of uninitialised variable
   dmin = 0.0 
   return
   
end function get_min_distance_to_square

function get_min_distance_to_line(p,e) result(dmin)

   implicit none
   real*4,intent(in)    :: p(3)  ! position on the line
   real*4,intent(in)    :: e(3)  ! unit vector parallel to the line
   real*4               :: dmin  ! minimal distance from the origin (0,0,0) to the line
   
   dmin = sqrt(sum(p**2)-sum(p*e)**2)
   
end function get_min_distance_to_line

logical function is_tile_in_survey(ishell,ix,user)

   implicit none
   integer*4,intent(in)    :: ishell
   integer*4,intent(in)    :: ix(3) ! [box side-length] tile center in tiling coordinates
   logical,intent(in)      :: user
   integer*4               :: n2D,n3D
   real*4                  :: x(3),dx(3),dy(3),dz(3),sx(3),sy(3),sz(3),fi,fj,d
   integer*4               :: i,j,k
   integer*4,allocatable   :: index(:)
   
   sx = matmul(shell(ishell)%transformation%rotation,(/1.0,0.0,0.0/))
   sy = matmul(shell(ishell)%transformation%rotation,(/0.0,1.0,0.0/))
   sz = matmul(shell(ishell)%transformation%rotation,(/0.0,0.0,1.0/))
   x  = matmul(shell(ishell)%transformation%rotation,real(ix,4))

   d = max(0.5,sqrt(real(sum(ix**2),4))) ! [box side-length] approximate distance from origin to nearest tile face
   n2D = max(10,2*nint(0.5/(d*para%search_angle)))
   allocate(index(0:n2D))
   do i = 0,n2D/2-1
      index(i*2) = i
      index(i*2+1) = n2D-i
   end do
   index(n2D) = n2D/2

   do i = 0,n2D
      fi = real(index(i),4)/n2D-0.5
      do j = 0,n2D
         fj = real(index(j),4)/n2D-0.5
         if (is_point_in_survey(x+fi*sx+fj*sy+0.5*sz,user)) then; is_tile_in_survey = .true.; return; end if
         if (is_point_in_survey(x+fi*sx+fj*sy-0.5*sz,user)) then; is_tile_in_survey = .true.; return; end if
         if (is_point_in_survey(x+fi*sx+0.5*sy+fj*sz,user)) then; is_tile_in_survey = .true.; return; end if
         if (is_point_in_survey(x+fi*sx-0.5*sy+fj*sz,user)) then; is_tile_in_survey = .true.; return; end if
         if (is_point_in_survey(x+0.5*sx+fi*sy+fj*sz,user)) then; is_tile_in_survey = .true.; return; end if
         if (is_point_in_survey(x-0.5*sx+fi*sy+fj*sz,user)) then; is_tile_in_survey = .true.; return; end if
      end do
   end do
   
   n3D = int(2**para%volume_search_level-1,4)
   do i = 0,n3D
      dx = ((real(i,4)+0.5)/(n3D+1)-0.5)*sx
      do j = 0,n3D
         dy = ((real(j,4)+0.5)/(n3D+1)-0.5)*sy
         do k = 0,n3D
            dz = ((real(k,4)+0.5)/(n3D+1)-0.5)*sz
            if (is_point_in_survey(x+dx+dy+dz,user)) then; is_tile_in_survey = .true.; return; end if
         end do
      end do
   end do
   
   is_tile_in_survey = .false. 

end function is_tile_in_survey

logical function is_point_in_survey(x,user)

   implicit none
   real*4,intent(in)    :: x(3)     ! [box side length] position of point in tiling coordinates
   logical,intent(in)   :: user
   type(vector4)        :: v
   type(type_spherical) :: sph
   
   nchecks = nchecks+1
   
   v = x
   call car2sph(v,sph)
   
   if (user) then
   
      is_point_in_survey = .true.
      call selection_function(pos=sph_deg_l(sph),selected=is_point_in_survey)
      return
      
   else
   
      if (sph%dc<=epsilon(sph%dc)) then
         is_point_in_survey = (fov%dc(1)<=0)
         return
      end if
      is_point_in_survey = is_in_fov(sph)
      
   end if
   
end function is_point_in_survey

function apply_tile_symmetry(x,itile) result(y)

   ! converts the internal tile symmetries to the normalised coordinates [0...1]^3

   real*4,intent(in)    :: x(3)     ! [box side length] position in N-body box
   integer*4,intent(in) :: itile    ! tile index
   real*4               :: y(3)     ! [box side length] position in cartesian sky coordinates
   integer*4            :: ishell   ! shell index
   real*4               :: sgn      ! sign of axes transformations
   
   ishell = tile(itile)%shell
   
   ! apply translations of shell and tile (imposing periodic boundary conditions)
   y = modulo(x+shell(ishell)%transformation%translation+tile(itile)%transformation%translation,1.0)
   
   ! apply inversion of shell and tile
   sgn = (1-2*log2int(shell(ishell)%transformation%inverted))*(1-2*log2int(tile(itile)%transformation%inverted))
   y = (y-0.5)*sgn+0.5
   
   ! apply *discrete* rotation of tile about the tile centre
   y = matmul(tile(itile)%transformation%rotation,y-0.5)+0.5

end function apply_tile_symmetry

function map_tile_onto_sky(x,itile) result(y)

   ! converts the normalised tile coordinates [0...1]^3 into sky coordinates in R^3

   implicit none
   real*4,intent(in)    :: x(3)     ! [box side length] position in N-body box
   integer*4,intent(in) :: itile    ! tile index
   real*4               :: y(3)     ! [box side length] position in cartesian sky coordinates
   integer*4            :: ishell   ! shell index
   
   ishell = tile(itile)%shell
   
   ! translate nbody box to tile position
   y = x+tile(itile)%ix-0.5
   
   ! shell rotation
   y = matmul(shell(ishell)%transformation%rotation,y)
   
end function map_tile_onto_sky

function sph_deg_l(sph) result(s)

   ! converts rad to deg
   ! and box side lengths to user length units
   
   implicit none
   type(type_spherical),intent(in) :: sph
   type(type_spherical) :: s
   
   s = type_spherical(dc=sph%dc*para%box_side,ra=sph%ra/unit%degree,dec=sph%dec/unit%degree)
   
end function sph_deg_l

subroutine write_hdf5_mapping(filename_hdf5)

   implicit none
   character(*),intent(in)                      :: filename_hdf5  ! output filename
   character(:),allocatable                     :: name
   integer*4                                    :: i

   allocate(character(1)::name) ! empty allocation to avoid compiler flags

   ! open HDF5 file
   call hdf5_open(filename_hdf5,.true.)
   call hdf5_add_group('mapping/')
   
   ! write group "snapshots"
   name = 'mapping/snapshots/'
   call hdf5_add_group(name)
   call hdf5_write_data(name//'id',(/(i,i=para%snapshot_min,para%snapshot_max)/), &
   & 'snapshot number')
   call hdf5_write_data(name//'z',snapshot%redshift, &
   & 'redshift corresponding to the cosmic time of this snapshot')
   call hdf5_write_data(name//'dc_min',snapshot%dmin, &
   & '[box side length] minimal comoving distance at which this snapshot is used')
   call hdf5_write_data(name//'dc_max',snapshot%dmax, &
   & '[box side length] maximal comoving distance at which this snapshot is used')
   call hdf5_write_data(name//'n_tiles',snapshot%n_tiles, &
   & 'Number of tiles this snapshot has been considered for, irrespective of whether a galaxy was selected')

   ! write group "tiles"
   name = 'mapping/tiles/'
   call hdf5_add_group(name)
   call hdf5_write_data(name//'tile_id',(/(i,i=1,size(tile),1)/),'index of cubic tile')
   call hdf5_write_data(name//'shell_id',tile%shell,'index of the shell containing the cubic tile')
   call hdf5_write_data(name//'center_x',tile%ix(1),'[box side length] x-coordinate of tile centre')
   call hdf5_write_data(name//'center_y',tile%ix(2),'[box side length] y-coordinate of tile centre')
   call hdf5_write_data(name//'center_z',tile%ix(3),'[box side length] z-coordinate of tile centre')
   call hdf5_write_data(name//'dc_min',tile%dmin,'[box side length] minimum comoving distance of tile')
   call hdf5_write_data(name//'dc_max',tile%dmax,'[box side length] maximum comoving distance of tile')
   call hdf5_add_group(name//'transformation')
   call hdf5_write_data(name//'transformation/rotation_xx',tile%transformation%rotation(1,1),'xx-element of rotation matrix')
   call hdf5_write_data(name//'transformation/rotation_xy',tile%transformation%rotation(1,2),'xy-element of rotation matrix')
   call hdf5_write_data(name//'transformation/rotation_xz',tile%transformation%rotation(1,3),'xz-element of rotation matrix')
   call hdf5_write_data(name//'transformation/rotation_yx',tile%transformation%rotation(2,1),'yx-element of rotation matrix')
   call hdf5_write_data(name//'transformation/rotation_yy',tile%transformation%rotation(2,2),'yy-element of rotation matrix')
   call hdf5_write_data(name//'transformation/rotation_yz',tile%transformation%rotation(2,3),'yz-element of rotation matrix')
   call hdf5_write_data(name//'transformation/rotation_zx',tile%transformation%rotation(3,1),'zx-element of rotation matrix')
   call hdf5_write_data(name//'transformation/rotation_zy',tile%transformation%rotation(3,2),'zy-element of rotation matrix')
   call hdf5_write_data(name//'transformation/rotation_zz',tile%transformation%rotation(3,3),'zz-element of rotation matrix')
   call hdf5_write_data(name//'transformation/translation_x',tile%transformation%translation(1), & 
   & '[box side length] x-component of translation vector')
   call hdf5_write_data(name//'transformation/translation_y',tile%transformation%translation(2), & 
   & '[box side length] y-component of translation vector')
   call hdf5_write_data(name//'transformation/translation_z',tile%transformation%translation(3), & 
   & '[box side length] z-component of translation vector')
   call hdf5_write_data(name//'transformation/inverted',tile%transformation%inverted, &
   & 'logical flag for axis inversion (0 = no inversion, 1 = all three axes inverted)')

   ! write group "shells"
   name = 'mapping/shells/'
   call hdf5_add_group(name)
   call hdf5_write_data(name//'shell_id',(/(i,i=1,size(shell),1)/),'index of spherical shell')
   call hdf5_write_data(name//'dc_min',shell%dmin,'[box side length] minimum comoving distance of shell')
   call hdf5_write_data(name//'dc_max',shell%dmax,'[box side length] maximum comoving distance of shell')
   call hdf5_add_group(name//'transformation')
   call hdf5_write_data(name//'transformation/rotation_xx',shell%transformation%rotation(1,1),'xx-element of rotation matrix')
   call hdf5_write_data(name//'transformation/rotation_xy',shell%transformation%rotation(1,2),'xy-element of rotation matrix')
   call hdf5_write_data(name//'transformation/rotation_xz',shell%transformation%rotation(1,3),'xz-element of rotation matrix')
   call hdf5_write_data(name//'transformation/rotation_yx',shell%transformation%rotation(2,1),'yx-element of rotation matrix')
   call hdf5_write_data(name//'transformation/rotation_yy',shell%transformation%rotation(2,2),'yy-element of rotation matrix')
   call hdf5_write_data(name//'transformation/rotation_yz',shell%transformation%rotation(2,3),'yz-element of rotation matrix')
   call hdf5_write_data(name//'transformation/rotation_zx',shell%transformation%rotation(3,1),'zx-element of rotation matrix')
   call hdf5_write_data(name//'transformation/rotation_zy',shell%transformation%rotation(3,2),'zy-element of rotation matrix')
   call hdf5_write_data(name//'transformation/rotation_zz',shell%transformation%rotation(3,3),'zz-element of rotation matrix')
   call hdf5_write_data(name//'transformation/translation_x',shell%transformation%translation(1), & 
   & '[box side length] x-component of translation vector')
   call hdf5_write_data(name//'transformation/translation_y',shell%transformation%translation(2), & 
   & '[box side length] y-component of translation vector')
   call hdf5_write_data(name//'transformation/translation_z',shell%transformation%translation(3), & 
   & '[box side length] z-component of translation vector')
   call hdf5_write_data(name//'transformation/inverted',shell%transformation%inverted, &
   & 'logical flag for axis inversion (0 = no inversion, 1 = all three axes inverted)')
   
   ! write group "fov"
   name = 'mapping/fov/'
   call hdf5_add_group(name)
   call hdf5_write_data(name//'dc_min',fov%dc(1)*para%box_side, &
   & '[simulation length units] minimum comoving distance as defined in selection function')
   call hdf5_write_data(name//'ra_min',fov%ra(1)/unit%degree, &
   & '[deg] minimum right ascension as defined in selection function')
   call hdf5_write_data(name//'dec_min',fov%dec(1)/unit%degree, &
   & '[deg] minimum declination as defined in selection function')
   call hdf5_write_data(name//'dc_max',fov%dc(2)*para%box_side, &
   & '[simulation length units] maximum comoving distance as defined in selection function')
   call hdf5_write_data(name//'ra_max',fov%ra(2)/unit%degree, &
   & '[deg] maximum right ascension as defined in selection function')
   call hdf5_write_data(name//'dec_max',fov%dec(2)/unit%degree, &
   & '[deg] maximum declination as defined in selection function')

   ! close HDF5 file
   call hdf5_close()

end subroutine write_hdf5_mapping
             
end module module_tiling