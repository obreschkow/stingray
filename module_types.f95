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

      real*4      :: dc,ra,dec         ! [simulation length unit,rad,rad] position in spherical Sky-coords
      integer*4   :: tile              ! unique identifier of box in mock sky
      integer*4   :: sam_selection     ! selection by sam-properties (selected if >0)
      integer*4   :: pos_selection     ! selection by sky-position (selected if >0)
      integer*4   :: sky_selection     ! selection by sky-properties other than position (selected if >0)
   
   end type type_base

   type type_tile
   
      integer*4            :: ix(3)          ! integer position, where ix=(0,0,0) is the central box with the observer in the middle
      real*4               :: dmin           ! minimum distance to observer in units of box side-length
      real*4               :: dmax           ! maximum distance ...
      integer*4            :: rotation       ! 1...6, describing the type of proper 90-rotation, 1 being the identity; if negative with inversion
      real*4               :: Rvector(3,3)   ! matrix of full rotation = 90 tiling rotation, followed by sky rotation
      real*4               :: Rpseudo(3,3)   ! matrix of full rotation = 90 tiling rotation without inversion, followed by sky rotation
      real*4               :: translation(3) ! [units of side-length] translation vector [0...1]
   
   end type type_tile

   type type_snapshot

      real*4               :: redshift
      real*4               :: dmin        ! [units of side-length] minimum comoving distance at which galaxies are drawn from this redshift
      real*4               :: dmax        ! [units of side-length] maximum ...

   end type type_snapshot
   
   type(type_para)                  :: para
   type(type_tile),allocatable      :: tile(:)
   type(type_snapshot),allocatable  :: snapshot(:)
   
   interface line
      module procedure line_string
      module procedure line_int4
      module procedure line_real4
   end interface line
   
contains

   subroutine save_parameters

      implicit none
   
      character(len=255)   :: filename
   
      ! Binary file
      filename = trim(para%path_output)//'parameters.bin'
      open(1,file=trim(filename),action='write',form="unformatted",status='replace')
      write(1) para
      close(1)
   
      ! ASCII file
      ! open ascii file (only for non-derived parameters)
      filename = trim(para%path_output)//'parameters.txt'
      open(1,file=trim(filename),action='write',form="formatted",status='replace')
      
      ! write data
      call line('path_output',para%path_output)
      call line('path_input',para%path_input)
      call line('L',para%L)
      call line('length_unit',para%length_unit)
      call line('snapshot_min',para%snapshot_min)
      call line('snapshot_max',para%snapshot_max)
      call line('subsnapshot_min',para%subsnapshot_min)
      call line('subsnapshot_max',para%subsnapshot_max)
      call line('h',para%h)
      call line('OmegaL',para%OmegaL)
      call line('OmegaM',para%OmegaM)
      call line('OmegaB',para%OmegaB)
      call line('dc_min',para%dc_min)
      call line('dc_max',para%dc_max)
      call line('ra_min',para%ra_min/degree)
      call line('ra_max',para%ra_max/degree)
      call line('dec_min',para%dec_min/degree)
      call line('dec_max',para%dec_max/degree)
      call line('zaxis_ra',para%zaxis_ra)
      call line('zaxis_dec',para%zaxis_dec)
      call line('xy_angle',para%xy_angle)
      call line('seed',para%seed)
      call line('translate',para%translate)
      call line('rotate',para%rotate)
      call line('invert',para%invert)
      call line('velocity_ra',para%velocity_ra)
      call line('velocity_dec',para%velocity_dec)
      call line('velocity_norm',para%velocity_norm)
      call line('search_angle',para%search_angle)
      call line('volume_search_level',para%volume_search_level)
      
      ! close files
      close(1)

   end subroutine save_parameters
   
   subroutine load_parameters

      implicit none
      character(len=255)   :: filename
      filename = trim(para%path_output)//'parameters.bin'
      call check_exists(filename)
      open(1,file=trim(filename),action='read',form="unformatted")
      read(1) para
      close(1)
   
   end subroutine load_parameters

   subroutine save_snapshot_list
   
      implicit none
      character(len=255)   :: filename
      integer*4            :: i
      
      ! write to ascii file
      filename = trim(para%path_output)//'snapshots.txt'
      open(1,file=trim(filename),action='write',form="formatted",status='replace')
      write(1,'(A)') 'Stingray snapshot list'
      write(1,'(A)') '------------------------------------------------------------------------------------------'
      write(1,'(A)') 'Col  1:  Snapshot ID'
      write(1,'(A)') 'Col  2:  Redshift corresponding to the cosmic time of the snapshot'
      write(1,'(A)') 'Col  3:  [sim units] min comoving distance at which galaxies are drawn from this snapshot.'
      write(1,'(A)') 'Col  4:  [sim units] max comoving distance at which galaxies are drawn from this snapshot.'
      write(1,'(A)') '------------------------------------------------------------------------------------------'
      do i = lbound(snapshot,1),ubound(snapshot,1)
         write(1,'(I6,3F14.7,I3,3F9.5)') i,snapshot(i)%redshift,snapshot(i)%dmin,snapshot(i)%dmax
      end do
      close(1)
      
      ! write to binary file
      filename = trim(para%path_output)//'snapshots.bin'
      open(1,file=trim(filename),action='write',form='unformatted',status='replace')
      write(1) lbound(snapshot,1),ubound(snapshot,1)
      do i = lbound(snapshot,1),ubound(snapshot,1)
         write(1) snapshot(i)
      end do
      close(1)
      
   end subroutine save_snapshot_list
   
   subroutine load_snapshot_list
   
      implicit none
      character(len=255)   :: filename
      integer*4            :: i,sn_min,sn_max
   
         ! write to binary file
      filename = trim(para%path_output)//'snapshots.bin'
      open(1,file=trim(filename),action='read',form='unformatted')
      read(1) sn_min,sn_max
      if (allocated(snapshot)) deallocate(snapshot)
      allocate(snapshot(sn_min:sn_max))
      do i = lbound(snapshot,1),ubound(snapshot,1)
         read(1) snapshot(i)
      end do
      close(1)
      
   end subroutine load_snapshot_list

   subroutine save_box_list

      implicit none
      character(len=255)   :: filename
      integer*4            :: i
   
      ! write to ascii file
      filename = trim(para%path_output)//'tiling.txt'
      open(1,file=trim(filename),action='write',form="formatted",status='replace')
      write(1,'(A)') 'Stingray 3D tiling geometry'
      write(1,'(A)') '--------------------------------------------------------------------------------'
      write(1,'(A)') 'Col  1:  Box index, starting at 1'
      write(1,'(A)') 'Col  2:  x-position of box-centre in units of box side lengths (L)'
      write(1,'(A)') 'Col  3:  y-position of box-centre in units of L'
      write(1,'(A)') 'Col  4:  z-position of box-centre in units of L'
      write(1,'(A)') 'Col  5:  min comoving distance to be considered to fill this box in units of L'
      write(1,'(A)') 'Col  6:  max comoving distance to be considered to fill this box in units of L'
      write(1,'(A)') 'Col  7:  index [1,...,6] of rotation, where 1 is the identity (negative if inversion)'
      write(1,'(A)') 'Col  8:  x-component of translation vector in units of L'
      write(1,'(A)') 'Col  9:  y-component of translation vector in units of L'
      write(1,'(A)') 'Col 10:  z-component of translation vector in units of L'
      write(1,'(A)') '--------------------------------------------------------------------------------'
      do i = 1,size(tile)
         write(1,'(I6,3I6,2F13.6,I3,3F9.5)') i,tile(i)%ix,tile(i)%dmin,tile(i)%dmax,tile(i)%rotation,tile(i)%translation
      end do
      close(1)
   
      ! write to binary file
      filename = trim(para%path_output)//'tiling.bin'
      open(1,file=trim(filename),action='write',form='unformatted',status='replace')
      write(1) size(tile)
      do i = 1,size(tile)
         write(1) tile(i)
      end do
      close(1)

   end subroutine save_box_list

   subroutine load_box_list

      implicit none
      character(len=255)                     :: filename
      integer*4                              :: i,ntile
   
      filename = trim(para%path_output)//'tiling.bin'
      call check_exists(filename)
      open(1,file=trim(filename),action='read',form='unformatted')
      read(1) ntile
      if (allocated(tile)) deallocate(tile)
      allocate(tile(ntile))
      do i = 1,ntile
         read(1) tile(i)
      end do
      close(1)

   end subroutine load_box_list

   ! subroutines used by interface line
   
   subroutine line_string(name,value)

      implicit none
      character(*),intent(in) :: name
      character(*),intent(in) :: value
      character(len=30)       :: txtname
   
      write(txtname,'(A30)') trim(name)
      write(1,'(A30,A)') adjustl(txtname),trim(adjustl(value))

   end subroutine line_string

   subroutine line_int4(name,value)

      implicit none
      character(*),intent(in) :: name
      integer*4,intent(in)    :: value
      character(len=30)       :: txtname
   
      write(txtname,'(A30)') trim(name)
      write(1,'(A30,I6)') adjustl(txtname),value
   
   end subroutine line_int4

   subroutine line_real4(name,value)

      implicit none
      character(*),intent(in) :: name
      real*4,intent(in)       :: value
      character(len=30)       :: txtname
   
      write(txtname,'(A30)') trim(name)
      write(1,'(A30,E14.7)') adjustl(txtname),value
   
   end subroutine line_real4

end module module_types