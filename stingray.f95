program stingray

   use module_constants
   use module_types
   use module_system
   use module_io
   use module_user_routines
   use module_parameters
   use module_tiling
   use module_sky

   implicit none

   character(len=255)   :: parameter_filename_custom
   character(len=255)   :: arg_task
   character(len=255)   :: arg_option
   character(len=255)   :: arg_value
   character(len=255)   :: custom_option
   integer*4            :: i,narg
   
   narg = iargc() ! number of arguments
   
   if (narg==1) then
      call getarg(1,arg_task)
      if ((trim(arg_task)=='version').or.(trim(arg_task)=='-version').or.(trim(arg_task)=='-v')) then
         call task_version
         stop
      else if ((trim(arg_task)=='help').or.(trim(arg_task)=='-help').or.(trim(arg_task)=='-h')) then
         call task_help
         stop
      else
         write(*,'(A)') 'Unknown argument or wrong use of argument.'
         call task_help
         stop
      end if
   end if
   
   ! Change default options
   parameter_filename_custom = ''
   opt_logscreen = .true.
   opt_logfile = .true.
   custom_option = ''
   if (narg>1) then
      if (mod(narg,2)==1) then
         call out('ERROR: wrong input format. Each option needs to be given an argument.')
         stop
      end if
      do i = 1,narg-1,2
         call getarg(i,arg_option)
         call getarg(i+1,arg_value)
         select case (trim(arg_option))
         case ('-parameterfile')
            parameter_filename_custom = trim(arg_value)
         case ('-verbose')
            opt_logscreen = (arg_value=='1')
         case default
            call out('ERROR: '//trim(arg_option)//' is an unknown option.')
            stop
         end select
      end do
   end if
   
   ! Load/create paths
   call load_paths(parameter_filename_custom)
   
   ! Initialize verbose
   call out_open(version)
   
   ! Initialize constants
   call initialize_constants
   
   ! Load parameter file
   call make_parameters(parameter_filename_custom)
   
   ! Execute tasks
   call make_tiling
   call make_sky
   call make_hdf5
   
   ! Finalize verbose
   call out_close
   
   ! Clean up files
   call clean_up_files
   
   contains

   subroutine task_version
      implicit none
      write(*,'(A,A,A)') 'This is Stingray Version ',version,'.'
   end subroutine task_version
   
   subroutine task_help
      implicit none
      write(*,'(A)') 'The general use is stingray [-option ...]'
      write(*,'(A)') 'Consult the README file for additional information.'
   end subroutine
    
end program stingray