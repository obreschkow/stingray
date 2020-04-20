! This module defines global derived types, used by other modules

module module_global

   public

   ! File names
   character(*),parameter  :: fn_parameters_default = 'parameters.txt'
   character(*),parameter  :: fn_parameters = 'parameters'
   character(*),parameter  :: fn_snapshots = 'snapshots'
   character(*),parameter  :: fn_tiles = 'tiles'
   character(*),parameter  :: fn_galaxies = 'galaxies'
   character(*),parameter  :: fn_groups = 'groups'
   
   ! Limits
   type type_limit
      integer*4   :: n_tiles_max = 10000
      integer*4   :: n_tiles_on_a_line_max = 100
      integer*4   :: n_snapshots_max = 10000                 ! use factors of 10 for readable galaxy indices
      integer*4   :: n_subvolumes_max = 10000                ! use factors of 10 for readable galaxy indices
      integer*4   :: n_galaxies_per_tile_max = int(1e8,4)    ! use factors of 10 for readable galaxy indices
      integer*4   :: n_galaxies_sky_max = int(1e9,4)
   end type type_limit
   
   type(type_limit),protected   :: limit
   
   ! 90-degree rotation matrices
   real*4,protected  :: rot(3,3,-6:6)
   
   ! default values for type_para
   integer*4,private,parameter      :: i4 = -huge(0_4)
   real*4,private,parameter         :: r4 = -huge(0.0_4)
   character(1),private,parameter   :: c1 = ' '

   type type_para

      ! name of mock survey
      character(len=255)   :: survey = c1
      
      ! paths & file names
      character(len=255)   :: path_output = c1
      character(len=255)   :: path_input = c1
      character(len=255)   :: filename_sky = c1
   
      ! simulation box
      real*4               :: box_side = r4 ! [length unit of simulation] box side length
      real*4               :: length_unit = r4 ! [m] simulation length unit expressed in SI unit
      integer*4            :: snapshot_min = i4
      integer*4            :: snapshot_max = i4
      integer*4            :: subvolume_min = i4
      integer*4            :: subvolume_max = i4
   
      ! cosmology
      real*4               :: h = r4
      real*4               :: omega_l = r4
      real*4               :: omega_m = r4
      real*4               :: omega_b = r4
      
      ! distance range
      real*4               :: dc_min = r4  ! [length unit of simulation]
      real*4               :: dc_max = r4  ! [length unit of simulation]
      
      ! fov
      real*4               :: ra_min = r4  ! [rad]
      real*4               :: ra_max = r4  ! [rad]
      real*4               :: dec_min = r4 ! [rad]
      real*4               :: dec_max = r4 ! [rad]
      
      ! mapping of SAM coordinates onto survey coordinates
      real*4               :: zaxis_ra = r4   ! [rad]
      real*4               :: zaxis_dec = r4   ! [rad]
      real*4               :: xy_angle = r4    ! [rad]

      ! sky parameters
      integer*4            :: seed  ! seed of random number generator (integer >=1)
      integer*4            :: translate
      integer*4            :: rotate
      integer*4            :: invert
      
      ! observer velocity relative to CMB
      real*4               :: velocity_ra = r4    ! [rad]
      real*4               :: velocity_dec = r4   ! [rad]
      real*4               :: velocity_norm = r4  ! [km/s] peculiar velocity of observer with respect to Hubble flow
      
      ! advanced options
      real*4               :: search_angle = r4         ! [rad] minimal angular separation of points on faces
      real*4               :: volume_search_level = r4
      
      ! galaxy options
      integer*4            :: make_groups
      integer*4            :: line_parameters
      
      ! I/O options
      integer*4            :: merge_output = i4
      integer*4            :: keep_binaries = i4
      integer*4            :: keep_log = i4
      
      ! derived parameters, not directly specified by the user (no default needed)
      real*4               :: velocity_car(3)   ! [km/s] velocity of observer cartesian survey-coordinates
      real*4               :: sky_rotation(3,3) ! rotation matrix to move the (x,y,z)-sky axis onto the central (RA,dec)-sky
      
   end type type_para
   
   type type_base

      real*4      :: dc,ra,dec      ! [box lengths,rad,rad] position in spherical Sky-coords
      integer*4   :: tile           ! unique identifier of box in mock sky
      integer*4   :: group_ntot     ! total number of members in group
      integer*4   :: group_flag     ! group flag (0 if group unclipped, >0 if clipped by survey edge (+1), snapshot limit (+2), box limit (+4))
      
   end type type_base

   type type_tile
   
      integer*4   :: ix(3)          ! integer position, where ix=(0,0,0) is the central box with the observer in the middle
      real*4      :: dmin           ! minimum distance to observer in units of box side-length
      real*4      :: dmax           ! maximum distance ...
      integer*4   :: rotation       ! 1...6, describing the type of proper 90-rotation, 1 being the identity; if negative with inversion
      real*4      :: Rvector(3,3)   ! matrix of full rotation = 90 tiling rotation, followed by sky rotation
      real*4      :: Rpseudo(3,3)   ! matrix of full rotation = 90 tiling rotation without inversion, followed by sky rotation
      real*4      :: translation(3) ! [units of side-length] translation vector [0...1]
   
   end type type_tile

   type type_snapshot

      real*4      :: redshift
      real*4      :: dmin           ! [units of side-length] minimum comoving distance at which galaxies are drawn from this redshift
      real*4      :: dmax           ! [units of side-length] maximum ...
      integer*4   :: n_tiles        ! Number of tiles this snapshot has been considered for, irrespective of whether a galaxy was selected

   end type type_snapshot
   
   type type_skystats
   
      integer*4   :: n_galaxies = 0    ! number of galaxies in mock sky (or some part thereof)
      integer*4   :: n_distinct = 0    ! number of distinct SAM galaxies in mock sky (or some part thereof)
      integer*4   :: n_replica_max = 0 ! maximum number of replications of the same SAM galaxy
      integer*4   :: n_groups = 0      ! number of groups
   
   end type type_skystats
   
   type(type_para)                     :: para
   type(type_tile),allocatable         :: tile(:)
   type(type_snapshot),allocatable     :: snapshot(:)
   
contains

   subroutine initialize_global_variables
   
      implicit none
      integer*4 :: i
      
      rot(:,:,0) = 0
      
      ! proper rotations
      rot(:,:,1) = reshape((/+1,+0,+0,+0,+1,+0,+0,+0,+1/),(/3,3/))   ! identity
      rot(:,:,2) = reshape((/-1,+0,+0,+0,+0,+1,+0,+1,+0/),(/3,3/))   ! invert x-axis, while permuting y and z
      rot(:,:,3) = reshape((/+0,+0,+1,+0,-1,+0,+1,+0,+0/),(/3,3/))   ! invert y-axis, while permuting z and x
      rot(:,:,4) = reshape((/+0,+1,+0,+1,+0,+0,+0,+0,-1/),(/3,3/))   ! invert z-axis, while permuting x and y
      rot(:,:,5) = reshape((/+0,+1,+0,+0,+0,+1,+1,+0,+0/),(/3,3/))   ! permute (x,y,z) -> (y,z,x)
      rot(:,:,6) = reshape((/+0,+0,+1,+1,+0,+0,+0,+1,+0/),(/3,3/))   ! permute (x,y,z) -> (z,x,y)
      
      ! improper rotations (with axis flip)
      do i = 1,6
         rot(:,:,-i) = -rot(:,:,i)
      end do
      
   end subroutine initialize_global_variables

end module module_global