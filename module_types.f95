! This module defines global derived types, used by other modules

module module_types

   public

   type type_para

      ! NOTE: if changing this type, also update the following routines in module_parameters:
      !       subroutine reset_parameters
      !       subroutine check_parameters
      !       subroutine adjust_parameters
      !       subroutine load_user_parameters
      !       subroutine save_parameters
      !       [automatic parameter initialization in module_user_###]

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
      
      ! footprint on the sky
      real*4               :: ra       ! [rad]
      real*4               :: dec      ! [rad]
      real*4               :: angle    ! [rad] half opening angle
      real*4               :: ra_min   ! [rad]
      real*4               :: ra_max   ! [rad]
      real*4               :: dec_min  ! [rad]
      real*4               :: dec_max  ! [rad]
      
      ! distance range
      real*4               :: dc_min   ! [length unit of simulation]
      real*4               :: dc_max   ! [length unit of simulation]
      
      ! direction and orientation in the tiling grid
      real*4               :: axis(3)
      real*4               :: turn     ! [rad] 

      ! cone parameters
      integer*4            :: seed  ! seed of random number generator (integer >=1)
      integer*4            :: translate
      integer*4            :: rotate
      integer*4            :: invert
      integer*4            :: preserve_groups
   
      ! observer
      real*4               :: velocity(3)    ! [km/s] peculiar velocity of observer with respect to Hubble flow
      
      ! derived parameters, not directly specified by the user
      real*4               :: sky_rotation(3,3) ! rotation matrix to move the (x,y,z)-cone axis onto the central (RA,dec)-cone
   
   end type type_para

   type type_galaxy_base

      integer*8            :: groupid        ! unique identifier of group
      real*4               :: xbox(3)        ! [simulation length unit] position in SAM snapshot
      integer*4            :: snapshot       ! snapshot index
      integer*4            :: subsnapshot    ! sub-snapshot index
      integer*4            :: box            ! unique identifier of box in mock cone
      real*4               :: xcone(3)       ! [units of side-length L] position in mock cone after tilying
                                             ! (translations + 90deg-rotations), but before sky-rotation
      real*4               :: dc             ! [simulation length unit] comoving distance in cone
      real*4               :: ra,dec         ! [rad] position on sky
   
   end type type_galaxy_base

   type type_box
   
      integer*4            :: ix(3)          ! integer position, where ix=(0,0,0) is the central box with the observer in the middle
      real*4               :: dmin           ! minimum distance to observer in units of box side-length
      real*4               :: dmax           ! maximum distance ...
      integer*4            :: rotation       ! 1...6, describing the type of proper 90-rotation, 1 being the identity; if negative with inversion
      real*4               :: Rvector(3,3)   ! matrix of full rotation = 90 tiling rotation, followed by sky rotation
      real*4               :: Rpseudo(3,3)   ! matrix of full rotation = 90 tiling rotation without inversion, followed by sky rotation
      real*4               :: translation(3) ! [units of side-length] translation vector [0...1]
   
   end type type_box

   type type_snapshot

      real*4               :: redshift
      real*4               :: dmin        ! [units of side-length] minimum comoving distance at which galaxies are drawn from this redshift
      real*4               :: dmax        ! [units of side-length] maximum ...

   end type type_snapshot

end module module_types