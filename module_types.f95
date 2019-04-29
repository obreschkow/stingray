! This module defines global derived types, used by other modules, as well as the routines
! to save and load them

module module_types

   use module_constants
   use module_system

   public

   type type_para

      ! NOTE: if changing this type, also update the following routines in module_parameters:
      !       subroutine reset_parameters
      !       subroutine check_parameters
      !       subroutine adjust_parameters
      !       subroutine load_user_parameters
      !       subroutine save_parameters
      !       [automatic parameter initialization in module_user_###]
      
      ! name of mock survey
      character(len=255)   :: survey
      
      ! paths
      character(len=255)   :: path_output
      character(len=255)   :: path_input
   
      ! simulation box
      real*4               :: L ! box side length in simulation units
      real*4               :: length_unit ! [m]
      integer*4            :: snapshot_min
      integer*4            :: snapshot_max
      integer*4            :: subvolume_min
      integer*4            :: subvolume_max
   
      ! cosmology
      real*4               :: h
      real*4               :: omega_l
      real*4               :: omega_m
      real*4               :: omega_b
      
      ! distance range
      real*4               :: dc_min   ! [length unit of simulation]
      real*4               :: dc_max   ! [length unit of simulation]
      
      ! fov
      real*4               :: ra_min   ! [rad]
      real*4               :: ra_max   ! [rad]
      real*4               :: dec_min  ! [rad]
      real*4               :: dec_max  ! [rad]
      
      ! mapping of SAM coordinates onto survey coordinates
      real*4               :: zaxis_ra    ! [rad]
      real*4               :: zaxis_dec   ! [rad]
      real*4               :: xy_angle    ! [rad]

      ! sky parameters
      integer*4            :: seed  ! seed of random number generator (integer >=1)
      integer*4            :: translate
      integer*4            :: rotate
      integer*4            :: invert
      
      ! observer velocity relative to CMB
      real*4               :: velocity_ra    ! [rad]
      real*4               :: velocity_dec   ! [rad]
      real*4               :: velocity_norm  ! [km/s] peculiar velocity of observer with respect to Hubble flow
      
      ! advanced options
      real*4               :: search_angle         ! [rad] minimal angular separation of points on faces
      real*4               :: volume_search_level
      
      ! derived parameters, not directly specified by the user
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
   
   type(type_para)                     :: para
   type(type_tile),allocatable         :: tile(:)
   type(type_snapshot),allocatable     :: snapshot(:)

end module module_types