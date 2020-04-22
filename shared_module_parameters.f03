! **********************************************************************************************************************************
! Shared Fortran module to read parameter files
! Developed by Danail Obreschkow
!
! Basic structure
! The parameter file has to be an ascii-file, with one parameter per line, structured as:
! ----------------------------------------------------------------------
! # Example file
! parameter1_name    parameter1_value  # number of marbles
! parameter2_name    parameter2_value  # name of player
! ...
! ----------------------------------------------------------------------
! Empty lines are ignored, as well as all text following the #-symbol in non-empty lines. Use this symbol for comments.
!
! Advanced structure
! It is possible do define different parameter-sets in a single parameter file, using the syntax:
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
! If the protected variable parameterset is specified, e.g. parameterset='abc', only the parameters listed in this parameterset,
! e.g. the parameters between the lines 'parameterset abc' and 'end', as well as the parameters outside any parameterset, e.g.
! those listed at the end under '# default parameters' are read. If a parameter appears in the parameterset *and* outside, the value
! in the parameterset is used. Multiple parametersets of the same name are allowed, as long as they do not repeat parameters.
! If the variable parameterset is not specified, i.e. parameterset='', the effect is as if this variable is set to the first set,
! e.g. parameterset='abc'.
! **********************************************************************************************************************************

module shared_module_parameters

   use shared_module_core

   private
   
   public   :: set_parameterfile
   public   :: set_parameterset
   public   :: handle_parameters
   public   :: get_parameter_value
   public   :: require_no_parameters_left
   public   :: parameterfile ! read-only
   public   :: parameterset ! read-only
   public   :: n_parameters ! read-only
   public   :: set_autostring
   
   ! handling arguments
   character(len=255),protected  :: parameterfile = ''
   character(len=255),protected  :: parameterset = ''
   integer*4,parameter           :: n_parameters_max = 100
   integer*4,protected           :: n_parameters = 0
   character(len=255),protected  :: parameter_name(n_parameters_max)
   character(len=255),protected  :: parameter_value(n_parameters_max)
   logical*4,protected           :: parameter_used(n_parameters_max)
   integer*4,protected           :: n_used = 0
   character(len=255),protected  :: used_name(n_parameters_max)
   character(len=255),protected  :: autostring = ''
   
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
   parameterfile = txt

end subroutine set_parameterfile

subroutine set_parameterset(txt) ! sets parameter set name

   implicit none
   character(*),intent(in) :: txt
   parameterset = txt

end subroutine set_parameterset

subroutine set_autostring(txt)

   implicit none
   character(*),intent(in) :: txt
   autostring = trim(adjustl(txt))

end subroutine set_autostring

subroutine handle_parameters

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
   logical           :: parameter_in_default(n_parameters_max)
   logical           :: parameter_in_set(n_parameters_max)
   
   if (n_parameters.ne.0) call deverror('only call handle_parameters once')

   ! check if parameter file exists and if the user has read-access
   call check_file(parameterfile,'r')
   
   ! reset
   parameter_used = .false.
   parameter_in_default = .false.
   parameter_in_set = .false.
   
   ! read parameter file line-by-line
   open(1,file=trim(parameterfile),action='read',form='formatted')
   do
      
      read(1,'(A)',IOSTAT=io) line
      if (io.ne.0) exit
      if (.not.(isempty(line).or.(line(1:1)=='#'))) then
      
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
            
            if (isempty(value)) call error('parameterfile: the parameter "'//trim(name)//'" is empty')
            
            if (trim(name)=='parameterset') then
            
               if (inside_set) call error('parameterfile: a new parameterset appears before the parameterset '// &
               & trim(current_set)//' has been terminated with "end".')
               inside_set = .true.
               current_set = trim(value)
               if (isempty(parameterset)) parameterset = trim(current_set)
               reading = (trim(current_set)==trim(parameterset))
               if (reading) set_found = .true.
               
            else
            
               if (reading) then
               
                  i = parameter_index(trim(name))
                  if (i==0) then
                     n_parameters = n_parameters+1
                     if (n_parameters>n_parameters_max) call error('parameterfile: too many parameters')
                     i = n_parameters
                  end if
                  
                  if (inside_set) then
                     if (parameter_in_set(i)) call error('parameterfile: the parameter '//trim(name)//' appears more '//&
                     &'than once in the set '//trim(current_set))
                     parameter_in_set(i) = .true.
                  else
                     if (parameter_in_default(i)) call error('parameterfile: the parameter '//trim(name)//' appears more '//&
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
         
end subroutine handle_parameters

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
      if (required.and.(index==0)) call error('parameter '//trim(name)//' is required in parameterfile')
   end if

end function parameter_index

subroutine get_parameter_value_string(value,name,preset,auto,allowmultiuse)
   implicit none
   character(*),intent(inout)          :: value    ! value of parameter, must be initialised to an impossible value
   character(*),intent(in)             :: name     ! name of parameter
   character(*),intent(in),optional    :: preset   ! default value, if parameter name not found; if not given, parameter is required
   character(*),intent(in),optional    :: auto     ! value assigned, if parameter is identical to autostring, set via set_autostring
   logical*4,intent(in),optional       :: allowmultiuse ! if set to true, the same name can be queried multiple times
   integer*4                           :: i
   call check_multiuse(name,allowmultiuse)
   i = parameter_index(name,required=.not.present(preset))
   if (i==0) then
      value = preset
   else
      if ((.not.isempty(autostring)).and.(trim(parameter_value(i))==trim(autostring))) then
         if (present(auto)) then
            if (value==auto) call error('parameter '//trim(name)//' has not been set automatically')
            value = auto
         else
            call error('parameter '//trim(name)//' cannot be set automatically')
         end if
      else
         value = parameter_value(i)
      end if
      parameter_used(i) = .true.
   end if
end subroutine get_parameter_value_string

subroutine get_parameter_value_int4(value,name,preset,auto,allowmultiuse,min,max)
   implicit none
   integer*4,intent(inout)       :: value    ! value of parameter, must be initialised to an impossible value
   character(*),intent(in)       :: name     ! name of parameter
   integer*4,intent(in),optional :: preset   ! default value, if parameter name not found; if not given, the parameter is required
   integer*4,intent(in),optional :: auto     ! value assigned, if the parameter is identical to autostring, set via set_autostring
   logical*4,intent(in),optional :: allowmultiuse ! if set to true, the same name can be queried multiple times
   integer*4,intent(in),optional :: min,max  ! optional range required for value
   integer*4                     :: i,status
   call check_multiuse(name,allowmultiuse)
   i = parameter_index(name,required=.not.present(preset))
   if (i==0) then
      value = preset
   else
      if ((.not.isempty(autostring)).and.(trim(parameter_value(i))==trim(autostring))) then
         if (present(auto)) then
            if (value==auto) call error('parameter '//trim(name)//' has not been set automatically')
            value = auto
         else
            call error('parameter '//trim(name)//' cannot be set automatically')
         end if
      else
         read(parameter_value(i),*,iostat=status) value
         if (status.ne.0) call error('non-integer value found for parameter "'//trim(name)//'"')
      end if
      parameter_used(i) = .true.
   end if
   if (present(min)) then
      if (value<min) call error('parameter '//trim(name)//' cannot be smaller than '//val2str(min))
   end if
   if (present(max)) then
      if (value>max) call error('parameter '//trim(name)//' cannot be larger than '//val2str(max))
   end if
end subroutine get_parameter_value_int4

subroutine get_parameter_value_int8(value,name,preset,auto,allowmultiuse,min,max)
   implicit none
   integer*8,intent(inout)       :: value    ! value of parameter, must be initialised to an impossible value
   character(*),intent(in)       :: name     ! name of parameter
   integer*8,intent(in),optional :: preset   ! default value, if parameter name not found; if not given, the parameter is required
   integer*8,intent(in),optional :: auto     ! value assigned, if the parameter is identical to autostring, set via set_autostring
   logical*4,intent(in),optional :: allowmultiuse ! if set to true, the same name can be queried multiple times
   integer*8,intent(in),optional :: min,max  ! optional range required for value
   integer*4                     :: i,status
   call check_multiuse(name,allowmultiuse)
   i = parameter_index(name,required=.not.present(preset))
   if (i==0) then
      value = preset
   else
      if ((.not.isempty(autostring)).and.(trim(parameter_value(i))==trim(autostring))) then
         if (present(auto)) then
            if (value==auto) call error('parameter '//trim(name)//' has not been set automatically')
            value = auto
         else
            call error('parameter '//trim(name)//' cannot be set automatically')
         end if
      else
         read(parameter_value(i),*,iostat=status) value
         if (status.ne.0) call error('non-integer value found for parameter "'//trim(name)//'"')
      end if
      parameter_used(i) = .true.
   end if
   if (present(min)) then
      if (value<min) call error('parameter '//trim(name)//' cannot be smaller than '//val2str(min))
   end if
   if (present(max)) then
      if (value>max) call error('parameter '//trim(name)//' cannot be larger than '//val2str(max))
   end if
end subroutine get_parameter_value_int8

subroutine get_parameter_value_real4(value,name,preset,auto,allowmultiuse,min,max)
   implicit none
   real*4,intent(inout)          :: value    ! value of parameter, must be initialised to an impossible value
   character(*),intent(in)       :: name     ! name of parameter
   real*4,intent(in),optional    :: preset   ! default value, if parameter name not found; if not given, the parameter is required
   real*4,intent(in),optional    :: auto     ! value assigned, if the parameter is identical to autostring, set via set_autostring
   logical*4,intent(in),optional :: allowmultiuse ! if set to true, the same name can be queried multiple times
   real*4,intent(in),optional    :: min,max  ! optional range required for value
   integer*4                     :: i,status
   call check_multiuse(name,allowmultiuse)
   i = parameter_index(name,required=.not.present(preset))
   if (i==0) then
      value = preset
   else
      if ((.not.isempty(autostring)).and.(trim(parameter_value(i))==trim(autostring))) then
         if (present(auto)) then
            if (value==auto) call error('parameter '//trim(name)//' has not been set automatically')
            value = auto
         else
            call error('parameter '//trim(name)//' cannot be set automatically')
         end if
      else
         read(parameter_value(i),*,iostat=status) value
         if (status.ne.0) call error('non-numeric value found for parameter "'//trim(name)//'"')
      end if
      parameter_used(i) = .true.
   end if
   if (present(min)) then
      if (value<min) call error('parameter '//trim(name)//' cannot be smaller than '//val2str(min))
   end if
   if (present(max)) then
      if (value>max) call error('parameter '//trim(name)//' cannot be larger than '//val2str(max))
   end if
end subroutine get_parameter_value_real4

subroutine get_parameter_value_real8(value,name,preset,auto,allowmultiuse,min,max)
   implicit none
   real*8,intent(inout)          :: value    ! value of parameter, must be initialised to an impossible value
   character(*),intent(in)       :: name     ! name of parameter
   real*8,intent(in),optional    :: preset   ! default value, if parameter name not found; if not given, the parameter is required
   real*8,intent(in),optional    :: auto     ! value assigned, if the parameter is identical to autostring, set via set_autostring
   logical*4,intent(in),optional :: allowmultiuse ! if set to true, the same name can be queried multiple times
   real*8,intent(in),optional    :: min,max  ! optional range required for value
   integer*4                     :: i,status
   call check_multiuse(name,allowmultiuse)
   i = parameter_index(name,required=.not.present(preset))
   if (i==0) then
      value = preset
   else
      if ((.not.isempty(autostring)).and.(trim(parameter_value(i))==trim(autostring))) then
         if (present(auto)) then
            if (value==auto) call error('parameter '//trim(name)//' has not been set automatically')
            value = auto
         else
            call error('parameter '//trim(name)//' cannot be set automatically')
         end if
      else
         read(parameter_value(i),*,iostat=status) value
         if (status.ne.0) call error('non-numeric value found for parameter "'//trim(name)//'"')
      end if
      parameter_used(i) = .true.
   end if
   if (present(min)) then
      if (value<min) call error('parameter '//trim(name)//' cannot be smaller than '//val2str(min))
   end if
   if (present(max)) then
      if (value>max) call error('parameter '//trim(name)//' cannot be larger than '//val2str(max))
   end if
end subroutine get_parameter_value_real8

subroutine get_parameter_value_logical(value,name,preset,auto,allowmultiuse)
   implicit none
   logical*4,intent(inout)       :: value    ! value of parameter, must be initialised to an impossible value
   character(*),intent(in)       :: name     ! name of parameter
   logical*4,intent(in),optional :: preset   ! default value, if parameter name not found; if not given, the parameter is required
   logical*4,intent(in),optional :: auto     ! value assigned, if the parameter is identical to autostring, set via set_autostring
   logical*4,intent(in),optional :: allowmultiuse ! if set to true, the same name can be queried multiple times
   integer*4                     :: i
   call check_multiuse(name,allowmultiuse)
   i = parameter_index(name,required=.not.present(preset))
   if (i==0) then
      value = preset
   else
      if ((.not.isempty(autostring)).and.(trim(parameter_value(i))==trim(autostring))) then
         if (present(auto)) then
            if (value.eqv.auto) call error('parameter '//trim(name)//' has not been set automatically')
            value = auto
         else
            call error('parameter '//trim(name)//' cannot be set automatically')
         end if
      else
         if (any(trim(parameter_value(i))==(/'0','F','f','N','n'/))) then
            value = .false.
         else if (any(trim(parameter_value(i))==(/'1','T','t','Y','y'/))) then
            value = .true.
         else
            call error('parameter '//trim(name)//' only takes logical arguments (0/f/F/n/N and 1/t/T/y/Y).')
         end if
      end if
      parameter_used(i) = .true.
   end if
end subroutine get_parameter_value_logical

subroutine check_multiuse(name,allowmultiuse)
   implicit none
   character(*),intent(in)       :: name  ! name of parameter
   logical*4,intent(in),optional :: allowmultiuse ! if set true, the same name can be queried multiple times
   integer*4                     :: i,j
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