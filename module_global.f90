! This module defines global derived types, used by other modules

module module_global

   use shared_module_core
   use shared_module_constants
   use shared_module_vectors
   use shared_module_maths

   public

   ! File names
   character(*),parameter  :: fn_parameters_default = 'parameters.txt'
   character(*),parameter  :: fn_parameters = 'parameters'
   character(*),parameter  :: fn_snapshots = 'snapshots'
   character(*),parameter  :: fn_tiles = 'tiles'
   character(*),parameter  :: fn_galaxies = 'galaxies'
   character(*),parameter  :: fn_groups = 'groups'
   character(255)          :: path_tmp
   
   ! Limits
   type type_limit
      integer*4   :: n_tiles_max = 100000
      integer*4   :: n_tiles_on_a_line_max = 100
      integer*4   :: n_snapshots_max = 10000                 ! use factors of 10 for readable galaxy indices
      integer*4   :: n_subvolumes_max = 10000                ! use factors of 10 for readable galaxy indices
      integer*4   :: n_galaxies_per_tile_max = int(1e8,4)    ! use factors of 10 for readable galaxy indices
      integer*4   :: n_galaxies_sky_max = int(1e9,4)
      real*4      :: group_diameter_max = 0.25               ! [box side length] maximum allow group diameter
   end type type_limit
   
   ! Spherical coordinates
   type type_spherical
      real*4   :: dc = 0.0    ! [simulation units/box lengths] comoving distance from observer
      real*4   :: ra = 0.0    ! [rad/deg] right ascension
      real*4   :: dec = 0.0   ! [rad/deg] declination
   end type type_spherical
   
   ! Field of view ranges
   type type_fov
      real*4   :: dc(2)    ! [simulation units/box lengths] range of comoving distance from observer
      real*4   :: ra(2)    ! [deg/rad] range of right ascension
      real*4   :: dec(2)   ! [deg/rad] range of declination 
   end type type_fov
   
   ! Standard indices for all sky objects
   type type_index
      integer*8   :: galaxy = 0_8   ! unique galaxy index in mock sky (for groups: central galaxy or -1 if central not selected)
      integer*8   :: group = 0_8    ! unique group id in the mock sky (-1 for isolated galaxies)
      integer*4   :: tile = 0       ! tile index
      integer*4   :: shell = 0      ! shell index
      integer*4   :: snapshot = 0   ! snapshot index
      integer*4   :: subvolume = 0  ! subvolume index
   end type type_index
   
   ! Transformations of tiles and shells
   type type_transformation
      real*4      :: translation(3) = (/0.0,0.0,0.0/) ! [box side length] translation vector [0...1]
      real*4      :: rotation(3,3)  = const%identity3 ! proper rotation matrix of tile/shell
      logical     :: inverted       = .false.         ! logical flag for axis inversion (0 = no inversion, 1 = all three axes inverted)
   end type type_transformation

   type type_base
      type(type_spherical)       :: spherical = type_spherical()           ! spherical coordinates (dc,ra,dec) [simulation length units, rad, rad]
      type(vector4)              :: cartesian = vector4()                  ! cartesian coordinates (x,y,z) [simulation length units]
      type(type_index)           :: index = type_index()                   ! various indices
      type(type_transformation)  :: transformation = type_transformation() ! comoving coordinate transformation from N-body box to sky
      real*4                     :: snapshot_redshift = 0.0                ! redshift of this snapshot, z=1/a-1
   end type type_base
   
   type(type_limit),parameter          :: limit = type_limit()
   
   interface car2sph
      procedure car2sph_compact
   end interface car2sph
   
   interface sph2car
      procedure sph2car_compact
   end interface sph2car

   
contains

   subroutine car2sph_compact(car,sph)
   
      implicit none
      type(vector4),intent(in)         :: car
      type(type_spherical),intent(out) :: sph
      real*4                           :: x(3),radius,longitude,lattitude
      x = car
      call car2sph(x,radius,longitude,lattitude)
      sph = type_spherical(dc=radius,ra=longitude,dec=lattitude)
   
   end subroutine car2sph_compact

   subroutine sph2car_compact(sph,car)
   
      implicit none
      type(type_spherical),intent(in)  :: sph
      type(vector4),intent(out)        :: car
      real*4                           :: x(3)
      call sph2car(sph%dc,sph%ra,sph%dec,x)
      car = x
   
   end subroutine sph2car_compact

end module module_global