module module_global

character(*),parameter  :: version = '0.1'

type type_para

   ! file names
   character(len=255)   :: file_snapshot_base
   character(len=5)     :: file_snapshot_extension
   character(len=255)   :: file_redshifts
   character(len=255)   :: path_output
   
   ! simulation box
   real*4               :: L ! box side length in simulation units
   real*4               :: length_unit ! [m]
   integer*4            :: snapshot_min
   integer*4            :: snapshot_max
   
   ! cosmology
   real*4               :: h
   real*4               :: OmegaL
   real*4               :: OmegaM
   
   ! cone geometry
   real*4               :: angle
   real*4               :: dcmax
   real*4               :: dcmin
   real*4               :: axis(3)

   ! cone parameters
   logical              :: translate
   logical              :: rotate
   logical              :: preserve_groups
   
end type type_para

type type_galaxy

   integer*8            :: id          ! unique identifier 
   integer*8            :: groupid     ! unique identifier of group
   real*4               :: x(3)        ! position in snapshot
   
end type type_galaxy

type type_mockgalaxy

   integer*8            :: id          ! unique identifier
   integer*8            :: groupid     ! unique identifier of group
   integer*8            :: box         ! unique identifier of box in mock cone
   real*4               :: x(3)        ! position in mock cone
   real*4               :: dc          ! [Mpc/h] comoving distance = sqrt(sum(x**3))
   real*4               :: R(6)        ! rotation matrix of box: R11, R12, R13, R22, R23, R33

end type type_mockgalaxy

type type_box
   
   integer*4            :: ix(3)          ! integer position, where ix=(0,0,0) is the central box with the observer in the middle
   real*4               :: dmin           ! minimum distance to observer in units of box side-length
   real*4               :: dmax           ! maximum distance ...
   integer*4            :: rotation       ! 1...6, describing the type of proper 90-rotation, 1 being the identity
   real*4               :: translation(3) ! [units of side-length] translation vector [0...1]
   
end type type_box

type(type_para)                     :: para
type(type_galaxy),allocatable       :: galaxy(:)
type(type_mockgalaxy),allocatable   :: mockgalaxy(:)
type(type_box),allocatable          :: box(:)

real*4,allocatable                  :: redshift(:)

end module module_global