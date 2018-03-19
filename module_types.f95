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
   
      ! cone geometry
      real*4               :: angle
      integer*4            :: square_base ! integer (0/1) representation of logical number
      real*4               :: dc_min
      real*4               :: dc_max
      real*4               :: axis(3)

      ! cone parameters
      integer*4            :: translate
      integer*4            :: rotate
      integer*4            :: invert
      integer*4            :: preserve_groups
   
      ! observer
      real*4               :: velocity(3)    ! [km/s] peculiar velocity of observer with respect to Hubble flow
   
   end type type_para

   type type_galaxy_base

      integer*8            :: id             ! unique identifier 
      integer*8            :: groupid        ! unique identifier of group
      real*4               :: xbox(3)        ! [simulation length unit] position in snapshot
      integer*4            :: snapshot       ! snapshot index
      integer*4            :: subsnapshot    ! sub-snapshot index
      integer*4            :: box            ! unique identifier of box in mock cone
      integer*4            :: rotation       ! 1...6, describing the type of proper 90-rotation, 1 being the identity; if negative with inversion
      real*4               :: xcone(3)       ! [units of side-length] position in mock cone
   
   end type type_galaxy_base

   type type_box
   
      integer*4            :: ix(3)          ! integer position, where ix=(0,0,0) is the central box with the observer in the middle
      real*4               :: dmin           ! minimum distance to observer in units of box side-length
      real*4               :: dmax           ! maximum distance ...
      integer*4            :: rotation       ! 1...6, describing the type of proper 90-rotation, 1 being the identity; if negative with inversion
      real*4               :: translation(3) ! [units of side-length] translation vector [0...1]
   
   end type type_box

   type type_snapshot

      real*4               :: redshift
      real*4               :: dmin        ! [units of side-length] minimum comoving distance at which galaxies are drawn from this redshift
      real*4               :: dmax        ! [units of side-length] maximum ...

   end type type_snapshot

end module module_types