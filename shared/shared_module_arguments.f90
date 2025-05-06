! **********************************************************************************************************************************
! Shared Fortran module to interpret user arguments supplied with the program call
! Developed by Danail Obreschkow
!
! The module interprets program calls in the form
! ./program [task_name [task_value]] [-option_name option_value] [-option_name option_value] ...
! **********************************************************************************************************************************

module shared_module_arguments

   use shared_module_core

   private
   
   ! argument handling routines
   public   :: handle_arguments ! analyses arguments and processes -version, -logfile, -verbose
   public   :: istask ! logical function checking if the task_name is equal to the argument; can also check existence of task_value and options
   public   :: get_task_value ! subroutine returning the task value (as output argument, not as a function)
   public   :: get_option_value ! function returning option value (see function description for further details)
   public   :: require_no_options_left ! check if get_option_value has been used on all available options, if not produce error
   public   :: unknown_task ! display error message for unknown task_name
   
   ! less frequently used
   public   :: task_exists ! read-only logical flag specifying if a task name is given; use require_task argument in handle_arguments
   public   :: task_has_value ! read-only logical flag specifying if the task has a value; normally this is checked with istask()
   public   :: noptions ! read-only integer giving the number of optional arguments provided by the user
   
   ! handling arguments
   character(len=255),protected  :: task_name
   character(len=255),protected  :: task_value
   logical*4,protected           :: task_has_value
   logical*4,protected           :: task_exists
   integer*4,parameter           :: n_options_max = 10
   integer*4,protected           :: n_options
   character(len=255),protected  :: option_name(n_options_max)
   character(len=255),protected  :: option_value(n_options_max)
   logical*4,protected           :: option_used(n_options_max)
   
   interface get_option_value
      module procedure get_option_value_string
      module procedure get_option_value_int4
      module procedure get_option_value_int8
      module procedure get_option_value_real4
      module procedure get_option_value_real8
      module procedure get_option_value_logical
   end interface get_option_value
   
contains
   
! master routines ***************************************************************************************************************

subroutine handle_arguments(require_task,require_options)

   implicit none
   logical*4,intent(in),optional    :: require_task
   logical*4,intent(in),optional    :: require_options
   integer*4                        :: narg
   character(255)                   :: argument
   logical*4                        :: value_verbose
   character(255)                   :: value_logfile_name
   
   if (time_start.ne.-1) call deverror('call handle_arguments before calling start_output')
   
   ! return version text, if the first argument is version/-version/-v/help/-help
   narg = iargc() ! get number of arguments
   if (narg>0) then
      call getarg(1,argument)
      if ((argument=='version').or.(argument=='-version').or.(argument=='-v').or.(argument=='-help').or.(argument=='help')) then
         call hline
         call out('Running '//execname()//' version '//trim(version)//'.')
         if (.not.isempty(help)) call out(help)
         if (.not.isempty(copyright)) call out(copyright)
         call hline
         stop
      end if
   end if
   
   ! extract arguments and split into tasks and options
   call get_arguments
   
   ! checks if task existence
   if (present(require_task)) then
      if (require_task.and.(.not.task_exists)) call error('Task missing', showhelp=.true.)
      if ((.not.require_task).and.task_exists) call error('Unknown argument "'//trim(task_name)//'".', showhelp=.true.)
   end if
   
   ! check option existence
   if (present(require_options)) then
      if (require_options.and.(n_options==0)) call error('Options missing', showhelp=.true.)
      if ((.not.require_options).and.(n_options>0)) call error('Unknown argument "'//trim(option_name(1))//'".', showhelp=.true.)
   end if
   
   ! start logfile if requested with option -logfile
   call get_option_value(value_logfile_name,'-logfile','')
   call set_logfile_name(value_logfile_name)
   
   ! start output on screen, as specified by options -verbose (default: false if logfile set, otherwise true)
   call get_option_value(value_verbose,'-verbose',isempty(value_logfile_name))
   call set_verbose(value_verbose)

end subroutine handle_arguments


! handle input arguments ********************************************************************************************************

subroutine get_arguments

   ! extract all the user arguments and group them into the global variables
   ! task_name, task_value, option_name(:), option_value(:)

   implicit none
   
   integer*4         :: narg
   character(255)    :: argument
   integer*4         :: i,j
   
   ! default values
   task_name = ''
   task_value = ''
   task_has_value = .false.
   task_exists = .false.
   n_options = 0
   option_name = ''
   option_value = ''
   option_used = .false.
   
   ! get number of arguments
   narg = iargc()

   if (narg>0) then
   
      ! handle task
      call getarg(1,argument)
      if (argument(1:1)=='-') then
         i = 1
      else
         task_exists = .true.
         task_name = argument
         if (narg>1) then
            call getarg(2,argument)
            if (argument(1:1)=='-') then
               i = 2
            else
               task_has_value = .true.
               task_value = argument
               i = 3
            end if
         else
            i = 2
         end if
      end if
      
      ! handle options
      if ((i<=narg).and.(mod(narg-i,2).ne.1)) call error('each option must have exactly one argument.')
      do while (i<=narg)
         n_options = n_options+1
         if (n_options>n_options_max) call error('too many options.')
         call getarg(i,argument)
         if (argument(1:1).ne.'-') call error('option starting with "-" expected instead of "'//trim(argument)//'".')
         option_name(n_options) = argument
         call getarg(i+1,argument)
         if (argument(1:1)=='-') call error('option "'//trim(option_name(n_options))//'" must have a value not starting with "-".')
         option_value(n_options) = argument
         i = i+2
      end do
      if (i.ne.(narg+1)) call deverror('unknown error in task handler.')
      
      ! check if any options appear multiple times
      do i = 1,n_options-1
         do j = i+1,n_options
            if (trim(option_name(i))==trim(option_name(j))) then
               call error('option '//trim(option_name(i))//' appears multiple times.')
            end if
         end do
      end do
      
   end if

end subroutine get_arguments

logical function istask(name,require_value,require_options)

   ! check if task_name is equal to name and optionally check for the existence of values and options
   ! if require_value is not given, the task can have a value or not
   ! ir require_options is not given *or set to true*, the task can have options or not -> use get_option_values to require options
   
   implicit none
   character(*),intent(in)       :: name
   logical*4,intent(in),optional :: require_value ! if set, the task must have a value (true) or *must not* have a value (false)
   logical*4,intent(in),optional :: require_options ! if false, the task *must not* have options; true is ignored -> use get_option_values
   integer*4                     :: i
   
   if (task_exists) then
   
      istask = trim(name)==trim(task_name)
   
      if (istask) then
         if (present(require_value)) then
            if (require_value.and.(.not.task_has_value)) call error('task "'//trim(task_name)//'" requires an argument.')
            if ((.not.require_value).and.task_has_value) call error('task "'//trim(task_name)//'" requires no argument.')
         end if
         if (present(require_options)) then
            if (.not.require_options) then
               do i = 1,n_options
                  if (.not.option_used(i)) call error('invalid argument "'//trim(option_name(i))//'" for task "'//&
                  & trim(task_name)//'"')
               end do
            end if
         end if
      end if
      
   else
   
      istask = .false.
   
   end if

end function istask

subroutine unknown_task

   implicit none
   if (.not.task_exists) call deverror('do not call unknown_task if task does not exist')
   call error('"'//trim(task_name)//'" is an unknown task.')
   
end subroutine unknown_task

subroutine get_task_value(value)

   use,intrinsic :: iso_fortran_env, only : int8, int16, int32, int64, real32, real64, real128
   
   implicit none
   class(*),intent(out) :: value
   integer*4            :: status = 0
   
   if (.not.task_has_value) call deverror('do not call get_task_value if no task_value exists; '//&
   & 'use istask(...,require_value=.false.)')
   
   select type (value)
   type is (integer(kind=int8));    read(task_value,*,iostat=status) value
   type is (integer(kind=int16));   read(task_value,*,iostat=status) value
   type is (integer(kind=int32));   read(task_value,*,iostat=status) value
   type is (integer(kind=int64));   read(task_value,*,iostat=status) value
   type is (real(kind=real32));     read(task_value,*,iostat=status) value
   type is (real(kind=real64));     read(task_value,*,iostat=status) value
   !type is (real(kind=real128));    read(task_value,*,iostat=status) value ! not supported by all compilers
   type is (logical);               read(task_value,*,iostat=status) value
   type is (character(*));          read(task_value,*,iostat=status) value
   type is (complex(kind=real32));  read(task_value,*,iostat=status) value
   type is (complex(kind=real64));  read(task_value,*,iostat=status) value
   !type is (complex(kind=real128)); read(task_value,*,iostat=status) value
   class default
      call deverror('unknown variable type')
   end select
   
   if (status.ne.0) call error('invalid argument "'//trim(task_value)//'" for task "'//trim(task_name)//'"')
      
end subroutine get_task_value

subroutine require_no_options_left

   ! check if all options have been used and produce an error if not

   implicit none
   integer*4   :: i
   
   do i = 1,n_options
      if (.not.option_used(i)) call error('option "'//trim(option_name(i))//'" unknown or not used in this task.')
   end do
   
end subroutine require_no_options_left

subroutine get_option_value_string(value,name,preset)
   implicit none
   character(*),intent(out)         :: value
   character(*),intent(in)          :: name     ! name of option
   character(*),intent(in),optional :: preset   ! default value
   integer*4                        :: i
   i = option_index(name,required=.not.present(preset))
   if (i==0) then
      if (present(preset)) value = preset
   else
      value = trim(option_value(i))
   end if
end subroutine get_option_value_string

subroutine get_option_value_int4(value,name,preset,min,max)
   implicit none
   integer*4,intent(out)         :: value
   character(*),intent(in)       :: name     ! name of option
   integer*4,intent(in),optional :: preset   ! default value
   integer*4,intent(in),optional :: min,max  ! optional range required for value
   integer*4                     :: i
   integer*4                     :: status
   i = option_index(name,required=.not.present(preset))
   if (i==0) then
      if (present(preset)) value = preset
   else
      read(option_value(i),*,iostat=status) value
      if (status.ne.0) call error('non-integer value given for option "'//trim(name)//'"')
   end if
   if (present(min)) then
      if (value<min) call error('option "'//trim(name)//'" cannot be smaller than '//val2str(min))
   end if
   if (present(max)) then
      if (value>max) call error('option "'//trim(name)//'" cannot be larger than '//val2str(max))
   end if
end subroutine get_option_value_int4

subroutine get_option_value_int8(value,name,preset,min,max)
   implicit none
   integer*8,intent(out)         :: value
   character(*),intent(in)       :: name     ! name of option
   integer*8,intent(in),optional :: preset   ! default value
   integer*8,intent(in),optional :: min,max  ! optional range required for value
   integer*4                     :: i
   integer*4                     :: status
   i = option_index(name,required=.not.present(preset))
   if (i==0) then
      if (present(preset)) value = preset
   else
      read(option_value(i),*,iostat=status) value
      if (status.ne.0) call error('non-integer value given for option "'//trim(name)//'"')
   end if
   if (present(min)) then
      if (value<min) call error('option "'//trim(name)//'" cannot be smaller than '//val2str(min))
   end if
   if (present(max)) then
      if (value>max) call error('option "'//trim(name)//'" cannot be larger than '//val2str(max))
   end if
end subroutine get_option_value_int8

subroutine get_option_value_real4(value,name,preset,min,max)
   implicit none
   real*4,intent(out)            :: value
   character(*),intent(in)       :: name     ! name of option
   real*4,intent(in),optional    :: preset   ! default value
   real*4,intent(in),optional    :: min,max  ! optional range required for value
   integer*4                     :: i
   integer*4                     :: status
   i = option_index(name,required=.not.present(preset))
   if (i==0) then
      if (present(preset)) value = preset
   else
      read(option_value(i),*,iostat=status) value
      if (status.ne.0) call error('non-numeric value given for option "'//trim(name)//'"')
   end if
   if (present(min)) then
      if (value<min) call error('option "'//trim(name)//'" cannot be smaller than '//val2str(min))
   end if
   if (present(max)) then
      if (value>max) call error('option "'//trim(name)//'" cannot be larger than '//val2str(max))
   end if
end subroutine get_option_value_real4

subroutine get_option_value_real8(value,name,preset,min,max)
   implicit none
   real*8,intent(out)            :: value
   character(*),intent(in)       :: name     ! name of option
   real*8,intent(in),optional    :: preset   ! default value
   real*8,intent(in),optional    :: min,max  ! optional range required for value
   integer*4                     :: i
   integer*4                     :: status
   i = option_index(name,required=.not.present(preset))
   if (i==0) then
      if (present(preset)) value = preset
   else
      read(option_value(i),*,iostat=status) value
      if (status.ne.0) call error('non-numeric value given for option "'//trim(name)//'"')
   end if
   if (present(min)) then
      if (value<min) call error('option "'//trim(name)//'" cannot be smaller than '//val2str(min))
   end if
   if (present(max)) then
      if (value>max) call error('option "'//trim(name)//'" cannot be larger than '//val2str(max))
   end if
end subroutine get_option_value_real8

subroutine get_option_value_logical(value,name,preset)
   implicit none
   logical,intent(out)           :: value
   character(*),intent(in)       :: name     ! name of option
   logical*4,intent(in),optional :: preset   ! default value
   integer*4                     :: i
   i = option_index(name,required=.not.present(preset))
   if (i==0) then
      if (present(preset)) value = preset
   else
      if (any(trim(option_value(i))==(/'0','F','f','N','n'/))) then
         value = .false.
      else if (any(trim(option_value(i))==(/'1','T','t','Y','y'/))) then
         value = .true.
      else
         call error('option "'//trim(name)//'" only takes logical arguments (0/f/F/n/N and 1/t/T/y/Y).')
      end if
   end if
end subroutine get_option_value_logical

integer*4 function option_index(name,required)

   ! check of option exists, and if so, write it's value into opt_val

   implicit none
   character(*),intent(in)       :: name ! name of option
   logical*4,intent(in),optional :: required
   integer*4                     :: i
   
   option_index = 0
   do i = 1,n_options
      if (trim(option_name(i))==trim(name)) then
         option_index = i
         exit
      end if
   end do
   
   if (option_index>0) option_used(option_index)=.true.
   
   if (present(required)) then
      if (required.and.(option_index==0)) call error('argument "'//trim(name)//'" must be specified.')
   end if
   
end function option_index

end module shared_module_arguments