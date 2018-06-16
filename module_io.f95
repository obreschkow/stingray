module module_io

   use module_types
   
   public
   private :: line
   
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

end module module_io