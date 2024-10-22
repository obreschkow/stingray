! **********************************************************************************************************************************
! Shared Fortran module to read parameter files
! Developed by Danail Obreschkow
!
! Basic structure
! The parameter file has to be an ascii-file, with one parameter per line, structured as:
! ----------------------------------------------------------------------
! # Example file
! parameter1_name    parameter1_value  optional description
! parameter2_name    parameter2_value  optional description
! ...
! ----------------------------------------------------------------------
! Empty lines and lines starting with "#" are ignored, as well as all the text that follows the parameter values in the same line
! (e.g. "optional description" in the example above). 
!
! Advanced structure
! It is possible do define different parameter-sets in a single parameter file, using the following syntax:
! ----------------------------------------------------------------------
! parameterset abc
! parameter1_name    parameter1_value
! parameter2_name    parameter2_value
! end
!
! parameterset xyz
! parameter1_name    parameter1_value
! parameter3_name    parameter3_value
! end
!
! # default parameters
! parameter2_name    parameter2_value
! parameter3_name    parameter3_value
! parameter4_name    parameter4_value
! ----------------------------------------------------------------------
! If the protected variable parameterset is non-empty, e.g. parameterset='abc', the parameters listed in this parameterset, e.g.
! between the lines "parameterset abc" and "end" overwrite the default parameters specified outside the parameterset. The user can
! mark at most one parameterset in the parameterfile as the default parameterset using an asterix, e.g. "parameterset* abc". This
! parameterset is taken as the default, if the variable parameterset is empty. If no parameterset is marked as default and if the
! variable parameterset is empty, all parameters in parametersets are ignored. Multiple parametersets of the same name are allowed,
! as long as they do not repeat parameters.
!
! Example code:
! character(255)  :: parameter_filename
! character(255)  :: parameter_set
! real            :: x
! integer         :: y
! character(20)   :: s
! logical         :: l
! call get_option_value(parameter_filename,'-parameterfile','parameters.txt')
! call set_parameterfile(trim(parameter_filename))
! call get_option_value(parameter_set,'-parameterset','')
! call set_parameterset(parameter_set)
! call read_parameters
! call set_auto_string('auto')
! call set_auto_parameter('parameter1_name',0.5) ! this value will be assigned if parameter1 is set to 'auto' in the parameter file
! call get_parameter_value(x,'parameter1_name',0.3,min=0.0,max=1.0) ! 0.3 is the default if the parameter is not present in the file
! call get_parameter_value(y,'parameter2_name',max=10) ! no default value given => parameter is required
! call get_parameter_value(s,'parameter3_name','na')
! call get_parameter_value(l,'parameter4_name',.true.)
! call require_no_parameters_left ! produces an error, if there are additional parameters
! **********************************************************************************************************************************

module shared_module_parameters

   use shared_module_core

   private
   
   public   :: set_parameterfile
   public   :: set_parameterset
   public   :: read_parameters
   public   :: get_parameter_value
   public   :: require_no_parameters_left
   public   :: parameterfile ! read-only
   public   :: parameterset ! read-only, optional name of the user-selected parameterset, set via set_parameterset() or "*"
   public   :: n_parameters ! read-only, number of read parameters, accessible once read_parameters has been called
   public   :: set_auto_string ! routine to set a default parameter value, e.g. "auto", for parameters specified elsewhere
   public   :: set_auto_parameter
   
   ! handling arguments
   character(len=255),protected  :: parameterfile = ''
   character(len=255),protected  :: parameterset = ''
   integer*4,parameter           :: n_parameters_max = 100
   integer*4,protected           :: n_parameters = 0
   character(len=255),protected  :: parameter_name(n_parameters_max)
   character(len=255),protected  :: parameter_value(n_parameters_max)
   logical*4,protected           :: parameter_used(n_parameters_max)
   logical*4,protected           :: parameter_has_auto_value(n_parameters_max)
   integer*4,protected           :: n_used = 0
   character(len=255),protected  :: used_name(n_parameters_max)
   character(len=255),protected  :: auto_string = ''
   logical,protected             :: parameters_handled = .false.
   
   interface get_parameter_value
      module procedure get_parameter_value_string
      module procedure get_parameter_value_int4
      module procedure get_parameter_value_int8
      module procedure get_parameter_value_real4
      module procedure get_parameter_value_real8
      module procedure get_parameter_value_logical
   end interface get_parameter_value
   
contains

subroutine set_parameterfile(txt) ! sets parameter file name

   implicit none
   character(*),intent(in) :: txt
   if (parameters_handled) call deverror('call set_parameterset before calling read_parameters')
   parameterfile = txt

end subroutine set_parameterfile

subroutine set_parameterset(txt) ! sets parameter set name

   implicit none
   character(*),intent(in) :: txt
   if (parameters_handled) call deverror('call set_parameterset before calling read_parameters')
   parameterset = txt

end subroutine set_parameterset

subroutine set_auto_string(txt)

   implicit none
   character(*),intent(in) :: txt
   auto_string = trim(adjustl(txt))

end subroutine set_auto_string

subroutine set_auto_parameter(name,value)

   implicit none
   character(*),intent(in) :: name
   class(*),intent(in)     :: value
   integer*4               :: i
   
   if (.not.parameters_handled) call deverror('read_parameters must be called before set_auto_parameter')
   if (isempty(auto_string)) call deverror('auto_string must be set before calling set_auto_parameter')
   
   i = parameter_index(name,.false.)
   
   if (i==0) then
      call error('attempting to assign an automatic value to the parameter "'//trim(name)//&
      &'", which does not exist in parameterfile')
   else
      if (parameter_has_auto_value(i)) then
         call deverror('the parameter "'//trim(name)//'" can only be set once using set_auto_parameter')
      else if (parameter_used(i)) then
         call deverror('the parameter "'//trim(name)//'" has already been used before its call of set_auto_parameter')
      else
         if (trim(parameter_value(i))==trim(auto_string)) then
            parameter_value(i) = val2str(value)
            parameter_has_auto_value(i) = .true.
         end if
      end if
   end if
   
end subroutine set_auto_parameter

subroutine read_parameters

   implicit none
   character(255)    :: line
   character(255)    :: name
   character(255)    :: value
   integer*4         :: io
   integer*4         :: i
   logical           :: reading = .true.
   logical           :: inside_set = .false.
   character(255)    :: current_set = ''
   logical           :: set_found = .false.
   logical           :: default_parameterset_found = .false.
   logical           :: parameter_in_default(n_parameters_max)
   logical           :: parameter_in_set(n_parameters_max)
   
   if (parameters_handled) call deverror('only call read_parameters once')

   ! check if parameter file exists and if the user has read-access
   call check_file(parameterfile,'r')
   
   ! reset
   parameter_used = .false.
   parameter_has_auto_value = .false.
   parameter_in_default = .false.
   parameter_in_set = .false.
   
   ! read parameter file line-by-line
   open(1,file=trim(parameterfile),action='read',form='formatted')
   do
      
      read(1,'(A)',IOSTAT=io) line
      if (io/=0) exit
      if (.not.(isempty(line).or.(line(1:1)=='#'))) then
      
         line = tabs2spaces(line)
      
         read(line,*) name
         
         if (trim(name)=='end') then
         
            if (.not.inside_set) call error('parameterfile: "end" must be preceded by "parameterset"')
            inside_set = .false.
            reading = .true.
         
         else
         
            value = trim(adjustl(line(len(trim(name))+1:len(line))))
            do i = 2,len(trim(value))
               if (any(value(i:i)==(/' ','#'/))) then
                  value = value(1:i-1)
                  exit
               end if
            end do
            
            if ((trim(name)=='parameterset').or.(trim(name)=='parameterset*')) then
            
               if (inside_set) call error('parameterfile: a new parameterset appears before the parameterset '// &
               & trim(current_set)//' has been terminated with "end".')
               if (isempty(value)) call error('parameterfile: each parameterset must have a name')
               if (trim(name)=='parameterset*') then
                  if (default_parameterset_found) call error('parameterfile: several parametersets have been '//&
                  &'specified as default using "*"')
                  default_parameterset_found = .true.
                  if (isempty(parameterset)) parameterset = trim(value)
               end if
               inside_set = .true.
               current_set = trim(value)
               if (trim(value)==trim(parameterset)) then
                  reading = .true.
                  set_found = .true.
               else
                  reading = .false.
               end if
               
            else
            
               if (reading) then
               
                  i = parameter_index(trim(name))
                  if (i==0) then
                     n_parameters = n_parameters+1
                     if (n_parameters>n_parameters_max) call error('parameterfile: too many parameters')
                     i = n_parameters
                  end if
                  
                  if (inside_set) then
                     if (parameter_in_set(i)) call error('parameterfile: the parameter "'//trim(name)//'" appears more '//&
                     &'than once in the set '//trim(current_set))
                     parameter_in_set(i) = .true.
                  else
                     if (parameter_in_default(i)) call error('parameterfile: the parameter "'//trim(name)//'" appears more '//&
                     &'than once (outside parametersets)')
                     parameter_in_default(i) = .true.
                  end if
               
                  parameter_name(i) = name
                  if (inside_set) then
                     parameter_value(i) = value
                  else
                     if (.not.parameter_in_set(i)) parameter_value(i) = value
                  end if
               
               end if
            
            end if
         
         end if
         
      end if
   
   end do
   
   close(1)
   
   if (inside_set) call error('parameterfile: the parameterset "'//trim(current_set)//'" must be terminated with "end".')
   if ((.not.set_found).and.(.not.isempty(parameterset))) call error('parameterset "'//trim(parameterset)//'" not found in '//&
   &'parameterfile')
   
   parameters_handled = .true.
   
end subroutine read_parameters

subroutine require_no_parameters_left

   ! check if all parameters have been used

   implicit none
   integer*4   :: i
   
   do i = 1,n_parameters
      if (.not.parameter_used(i)) call error('parameter "'//trim(parameter_name(i))//'" not used')
   end do
   
end subroutine require_no_parameters_left

function parameter_index(name,required) result(index)

   implicit none

   character(*),intent(in)       :: name
   logical*4,intent(in),optional :: required
   integer*4                     :: index
   integer*4                     :: i
   
   index = 0
   
   do i = 1,n_parameters
      if (trim(parameter_name(i))==trim(name)) then
         if (index>0) call deverror('parameterfile: non-captured repetition of identical parameters')
         index = i
      end if
   end do
   
   if (present(required)) then
      if (required.and.(index==0)) call error('parameter "'//trim(name)//'" is required in parameterfile')
   end if

end function parameter_index

subroutine get_parameter_value_string(value,name,preset,allowmultiuse,min,max)
   implicit none
   character(*),intent(out)            :: value    ! value of parameter
   character(*),intent(in)             :: name     ! name of parameter
   character(*),intent(in),optional    :: preset   ! default value, if parameter name not found; if not given, parameter is required
   logical*4,intent(in),optional       :: allowmultiuse ! if set to true, the same name can be queried multiple times
   integer*4,intent(in),optional       :: min,max  ! optional values defining the min and max length of the trimmed value
   integer*4                           :: i
   character(255)                      :: val
   
   call check_multiuse(name,allowmultiuse)
   i = parameter_index(name,required=.not.present(preset))
   if (i==0) then
      val = preset
   else
      val = parameter_value(i)
      parameter_used(i) = .true.
   end if
   if (present(min)) then
      if (len(trim(val))<min) call error('parameter "'//trim(name)//'" must have at least '//val2str(min)//' characters')
   end if
   if (present(max)) then
      if (len(trim(val))>max) call error('parameter "'//trim(name)//'" must have at most '//val2str(max)//' characters')
   end if
   value = trim(val)
   
end subroutine get_parameter_value_string

subroutine get_parameter_value_int4(value,name,preset,allowmultiuse,min,max)
   implicit none
   integer*4,intent(out)         :: value    ! value of parameter
   character(*),intent(in)       :: name     ! name of parameter
   integer*4,intent(in),optional :: preset   ! default value, if parameter name not found; if not given, the parameter is required
   logical*4,intent(in),optional :: allowmultiuse ! if set to true, the same name can be queried multiple times
   integer*4,intent(in),optional :: min,max  ! optional range required for value
   integer*4                     :: i,status
   call check_multiuse(name,allowmultiuse)
   i = parameter_index(name,required=.not.present(preset))
   if (i==0) then
      value = preset
   else
      read(parameter_value(i),*,iostat=status) value
      if (status/=0) call error('non-integer value found for parameter "'//trim(name)//'"')
      parameter_used(i) = .true.
   end if
   if (present(min)) then
      if (value<min) call error('parameter "'//trim(name)//'" cannot be smaller than '//val2str(min))
   end if
   if (present(max)) then
      if (value>max) call error('parameter "'//trim(name)//'" cannot be larger than '//val2str(max))
   end if
end subroutine get_parameter_value_int4

subroutine get_parameter_value_int8(value,name,preset,allowmultiuse,min,max)
   implicit none
   integer*8,intent(out)         :: value    ! value of parameter
   character(*),intent(in)       :: name     ! name of parameter
   integer*8,intent(in),optional :: preset   ! default value, if parameter name not found; if not given, the parameter is required
   logical*4,intent(in),optional :: allowmultiuse ! if set to true, the same name can be queried multiple times
   integer*8,intent(in),optional :: min,max  ! optional range required for value
   integer*4                     :: i,status
   call check_multiuse(name,allowmultiuse)
   i = parameter_index(name,required=.not.present(preset))
   if (i==0) then
      value = preset
   else
      read(parameter_value(i),*,iostat=status) value
      if (status/=0) call error('non-integer value found for parameter "'//trim(name)//'"')
      parameter_used(i) = .true.
   end if
   if (present(min)) then
      if (value<min) call error('parameter "'//trim(name)//'" cannot be smaller than '//val2str(min))
   end if
   if (present(max)) then
      if (value>max) call error('parameter "'//trim(name)//'" cannot be larger than '//val2str(max))
   end if
end subroutine get_parameter_value_int8

subroutine get_parameter_value_real4(value,name,preset,allowmultiuse,min,max)
   implicit none
   real*4,intent(out)            :: value    ! value of parameter
   character(*),intent(in)       :: name     ! name of parameter
   real*4,intent(in),optional    :: preset   ! default value, if parameter name not found; if not given, the parameter is required
   logical*4,intent(in),optional :: allowmultiuse ! if set to true, the same name can be queried multiple times
   real*4,intent(in),optional    :: min,max  ! optional range required for value
   integer*4                     :: i,status
   call check_multiuse(name,allowmultiuse)
   i = parameter_index(name,required=.not.present(preset))
   if (i==0) then
      value = preset
   else
      read(parameter_value(i),*,iostat=status) value
      if (status/=0) call error('non-numeric value found for parameter "'//trim(name)//'"')
      parameter_used(i) = .true.
   end if
   if (present(min)) then
      if (value<min) call error('parameter "'//trim(name)//'" cannot be smaller than '//val2str(min))
   end if
   if (present(max)) then
      if (value>max) call error('parameter "'//trim(name)//'" cannot be larger than '//val2str(max))
   end if
end subroutine get_parameter_value_real4

subroutine get_parameter_value_real8(value,name,preset,allowmultiuse,min,max)
   implicit none
   real*8,intent(out)            :: value    ! value of parameter
   character(*),intent(in)       :: name     ! name of parameter
   real*8,intent(in),optional    :: preset   ! default value, if parameter name not found; if not given, the parameter is required
   logical*4,intent(in),optional :: allowmultiuse ! if set to true, the same name can be queried multiple times
   real*8,intent(in),optional    :: min,max  ! optional range required for value
   integer*4                     :: i,status
   call check_multiuse(name,allowmultiuse)
   i = parameter_index(name,required=.not.present(preset))
   if (i==0) then
      value = preset
   else
      read(parameter_value(i),*,iostat=status) value
      if (status/=0) call error('non-numeric value found for parameter "'//trim(name)//'"')
      parameter_used(i) = .true.
   end if
   if (present(min)) then
      if (value<min) call error('parameter "'//trim(name)//'" cannot be smaller than '//val2str(min))
   end if
   if (present(max)) then
      if (value>max) call error('parameter "'//trim(name)//'" cannot be larger than '//val2str(max))
   end if
end subroutine get_parameter_value_real8

subroutine get_parameter_value_logical(value,name,preset,allowmultiuse)
   implicit none
   logical*4,intent(out)         :: value    ! value of parameter
   character(*),intent(in)       :: name     ! name of parameter
   logical*4,intent(in),optional :: preset   ! default value, if parameter name not found; if not given, the parameter is required
   logical*4,intent(in),optional :: allowmultiuse ! if set to true, the same name can be queried multiple times
   integer*4                     :: i
   call check_multiuse(name,allowmultiuse)
   i = parameter_index(name,required=.not.present(preset))
   if (i==0) then
      value = preset
   else
      if (any(trim(parameter_value(i))==(/'0','F','f','N','n'/))) then
         value = .false.
      else if (any(trim(parameter_value(i))==(/'1','T','t','Y','y'/))) then
         value = .true.
      else
         call error('parameter "'//trim(name)//'" only takes logical arguments (0/f/F/n/N and 1/t/T/y/Y).')
      end if
      parameter_used(i) = .true.
   end if
end subroutine get_parameter_value_logical

subroutine check_multiuse(name,allowmultiuse)
   implicit none
   character(*),intent(in)       :: name  ! name of parameter
   logical*4,intent(in),optional :: allowmultiuse ! if set true, the same name can be queried multiple times
   integer*4                     :: i,j
   if (.not.parameters_handled) call deverror('must call read_parameters before using get_parameter_value')
   j = 0
   do i = 1,n_used
      if (trim(used_name(i))==trim(adjustl(name))) then
         j = i
         exit
      end if
   end do
   if (j==0) then
      n_used = n_used+1
      if (n_used>n_parameters_max) call deverror('unknown error in check_multiused')
      used_name(n_used) = trim(adjustl(name))
   else
      if (present(allowmultiuse)) then
         if (.not.allowmultiuse) call deverror('Parameter '//trim(adjustl(name))//' used multiple times.')
      else
         call deverror('Parameter '//trim(adjustl(name))//' used multiple times.')
      end if
   end if
end subroutine check_multiuse

end module shared_module_parameters