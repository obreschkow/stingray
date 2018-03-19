program surfsuite

   use module_system
   use module_parameters
   use module_user
   use module_geometry
   use module_cone_intrinsic
   use module_cone_apparent

   implicit none
   
   character(*),parameter  :: version = '0.1'

   character(len=255)      :: parameter_filename_custom
   character(len=255)      :: arg_task
   character(len=255)      :: arg_option
   character(len=255)      :: arg_value
   character(len=255)      :: custom_option
   integer*4               :: i,narg
   integer*4               :: seed
   logical                 :: success
      
   narg = iargc() ! number of arguments
   
   if (narg==0) then
      write(*,'(A)') 'The general use is stingray task [-option ...]'
      write(*,'(A)') 'Consult the README file for additional information.'
      stop
   end if
   
   if (narg==1) then
      call getarg(1,arg_task)
      if ((trim(arg_task)=='version').or.(trim(arg_task)=='-version').or.(trim(arg_task)=='-v')) then
         call task_version
         stop
      end if
   end if
   
   ! Change default options
   parameter_filename_custom = ''
   opt_logscreen = .true.
   opt_logfile = .true.
   seed = 1
   custom_option = ''
   if (narg>1) then
      if (mod(narg,2)==0) then
         call out('ERROR: wrong input format. Each option needs to be given an argument.')
         stop
      end if
      do i = 2,narg-1,2
         call getarg(i,arg_option)
         call getarg(i+1,arg_value)
         select case (trim(arg_option))
         case ('-parameterfile')
            parameter_filename_custom = trim(arg_value)
         case ('-verbose')
            opt_logscreen = (arg_value=='0')
         case ('-seed')
            read(arg_value,*) seed
         case ('-option')
            custom_option = trim(arg_value)
         case default
            call out('ERROR: '//trim(arg_option)//' is an unknown option.')
            stop
         end select
      end do
   end if
   
   ! Load paths
   call load_paths(parameter_filename_custom)
   logfilename = trim(para%path_output)//'log.txt'
   
   ! Initialize verbose
   call out_open(version)
   
   ! Execute tasks
   call getarg(1,arg_task)
   select case (trim(arg_task))
   case ('make.all')
      call make_parameters(parameter_filename_custom)
      call make_geometry(seed)
      call make_cone_intrinsic
      call make_cone_apparent
   case ('make.parameters')
      call make_parameters(parameter_filename_custom)
   case ('make.geometry')
      call make_parameters(parameter_filename_custom)
      call make_geometry(seed)
   case ('make.intrinsic.cone')
      call make_cone_intrinsic
   case ('make.apparent.cone')
      call make_cone_apparent
   case default
      call set_seed(seed)
      call handle_custom_arguments(trim(arg_task),trim(custom_option),success)
      if (.not.success) then
         call out('ERROR: '//trim(arg_task)//' is an unknown task.')
         stop
      end if
   end select
   
   ! Finalize verbose
   call out_close
   
   contains

   subroutine  task_version
      implicit none
      write(*,'(A,A,A)') 'This is SurfSuite Version ',version,'.'
      write(*,'(A)') 'Consult the README file for additional information.'
      stop
   end subroutine task_version
    
end program surfsuite