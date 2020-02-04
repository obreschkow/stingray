module module_io

   use module_constants
   use module_types
   use module_system
   
   public
   
contains

   subroutine save_parameters

      implicit none
      
      open(1,file=trim(para%path_output)//fn_parameters,action='write',form="unformatted",status='replace')
      write(1) para
      close(1)

   end subroutine save_parameters
   
   subroutine load_parameters

      implicit none
      character(len=255)   :: filename
      filename = trim(para%path_output)//fn_parameters
      call check_exists(filename)
      open(1,file=trim(filename),action='read',form="unformatted")
      read(1) para
      close(1)
   
   end subroutine load_parameters

   subroutine save_snapshot_list
   
      implicit none
      character(len=255)   :: filename
      integer*4            :: i
      
      filename = trim(para%path_output)//fn_snapshots
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
      filename = trim(para%path_output)//fn_snapshots
      open(1,file=trim(filename),action='read',form='unformatted')
      read(1) sn_min,sn_max
      if (allocated(snapshot)) deallocate(snapshot)
      allocate(snapshot(sn_min:sn_max))
      do i = lbound(snapshot,1),ubound(snapshot,1)
         read(1) snapshot(i)
      end do
      close(1)
      
   end subroutine load_snapshot_list

   subroutine save_tile_list

      implicit none
      character(len=255)   :: filename
      integer*4            :: i
   
      filename = trim(para%path_output)//fn_tiles
      open(1,file=trim(filename),action='write',form='unformatted',status='replace')
      write(1) size(tile)
      do i = 1,size(tile)
         write(1) tile(i)
      end do
      close(1)

   end subroutine save_tile_list

   subroutine load_tile_list

      implicit none
      character(len=255)                     :: filename
      integer*4                              :: i,ntile
   
      filename = trim(para%path_output)//fn_tiles
      call check_exists(filename)
      open(1,file=trim(filename),action='read',form='unformatted')
      read(1) ntile
      if (allocated(tile)) deallocate(tile)
      allocate(tile(ntile))
      do i = 1,ntile
         read(1) tile(i)
      end do
      close(1)

   end subroutine load_tile_list
   
   subroutine clean_up_files
   
      implicit none
      
      if (para%keep_binaries==0) then 
         call system('rm '//trim(para%path_output)//fn_parameters)
         call system('rm '//trim(para%path_output)//fn_snapshots)
         call system('rm '//trim(para%path_output)//fn_tiles)
         call system('rm '//trim(para%path_output)//fn_galaxies)
         if (para%make_groups==1) call system('rm '//trim(para%path_output)//fn_groups)
      end if
      
      if (para%keep_log==0) then 
         call system('rm '//trim(para%path_output)//fn_log)
      end if
   
   end subroutine clean_up_files

end module module_io