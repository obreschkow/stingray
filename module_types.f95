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
      character(len=255)   :: name
      
      ! paths
      character(len=255)   :: path_output
      character(len=255)   :: path_input
   
      ! simulation box
      real*4               :: L ! box side length in simulation units
      real*4               :: length_unit ! [m]
      integer*4            :: snapshot_min
      integer*4            :: snapshot_max
      integer*4            :: subsnapshot_min
      integer*4            :: subsnapshot_max
   
      ! cosmology
      real*4               :: h
      real*4               :: OmegaL
      real*4               :: OmegaM
      real*4               :: OmegaB
      
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

      real*4      :: dc,ra,dec      ! [simulation length unit,rad,rad] position in spherical Sky-coords
      integer*4   :: tile           ! unique identifier of box in mock sky
      integer*4   :: group_ntot     ! total number of members in group
      integer*4   :: group_nsel     ! selected number of members in group
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

      real*4   :: redshift
      real*4   :: dmin        ! [units of side-length] minimum comoving distance at which galaxies are drawn from this redshift
      real*4   :: dmax        ! [units of side-length] maximum ...

   end type type_snapshot
   
   type(type_para)                     :: para
   type(type_tile),allocatable         :: tile(:)
   type(type_snapshot),allocatable     :: snapshot(:)
   
contains

   subroutine nil(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20)
      class(*),optional :: x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20
      if (.false.) then
         select type(x1);  end select;    select type(x2);  end select
         select type(x3);  end select;    select type(x4);  end select
         select type(x5);  end select;    select type(x6);  end select
         select type(x7);  end select;    select type(x8);  end select
         select type(x9);  end select;    select type(x10); end select
         select type(x11); end select;    select type(x12); end select
         select type(x13); end select;    select type(x14); end select
         select type(x15); end select;    select type(x16); end select
         select type(x17); end select;    select type(x18); end select
         select type(x19); end select;    select type(x20); end select
      end if
   end subroutine nil

end module module_types