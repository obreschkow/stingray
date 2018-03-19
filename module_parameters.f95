module module_parameters

   use module_constants
   use module_types
   use module_system
   use module_user
   
   private
   public   :: make_parameters, load_parameters, load_paths
   
contains

subroutine make_parameters(parameter_filename_custom)

   implicit none
   character(255),intent(in)  :: parameter_filename_custom
   character(255)             :: parameter_filename
   
   if (len(trim(parameter_filename_custom))>0) then
      parameter_filename = parameter_filename_custom
   else
      parameter_filename = parameter_filename_default
   end if
   call check_exists(parameter_filename)
   
   call reset_parameters
   call make_automatic_parameters
   call load_user_parameters(parameter_filename)
   call check_parameters
   call adjust_parameters
   call save_parameters

end subroutine make_parameters

subroutine reset_parameters

   implicit none
   para%L = huge(para%L)
   para%length_unit = huge(para%length_unit)
   para%snapshot_min = huge(para%snapshot_min)
   para%snapshot_max = huge(para%snapshot_max)
   para%subsnapshot_min = huge(para%subsnapshot_min)
   para%subsnapshot_max = huge(para%subsnapshot_max)
   para%h = huge(para%h)
   para%OmegaL = huge(para%OmegaL)
   para%OmegaM = huge(para%OmegaM)
   para%angle = huge(para%angle)
   para%square_base = huge(para%square_base)
   para%dc_min = huge(para%dc_min)
   para%dc_max = huge(para%dc_max)
   para%axis = huge(para%axis)
   para%translate = huge(para%translate)
   para%rotate = huge(para%rotate)
   para%invert = huge(para%invert)
   para%preserve_groups = huge(para%preserve_groups)
   para%velocity = huge(para%velocity)

end subroutine reset_parameters

subroutine check_parameters

   implicit none
   
   ! check if all parameters have been initialized
   if (para%L == huge(para%L)) call wrong('L')
   if (para%length_unit == huge(para%length_unit)) call wrong('length_unit')
   if (para%snapshot_min == huge(para%snapshot_min)) call wrong('snapshot_min')
   if (para%snapshot_max == huge(para%snapshot_max)) call wrong('snapshot_max')
   if (para%subsnapshot_min == huge(para%subsnapshot_min)) call wrong('subsnapshot_min')
   if (para%subsnapshot_max == huge(para%subsnapshot_max)) call wrong('subsnapshot_max')
   if (para%h == huge(para%h)) call wrong('h')
   if (para%OmegaL == huge(para%OmegaL)) call wrong('OmegaL')
   if (para%OmegaM == huge(para%OmegaM)) call wrong('OmegaM')
   if (para%angle == huge(para%angle)) call wrong('angle')
   if (para%square_base == huge(para%square_base)) call wrong('square_base')
   if (para%dc_min == huge(para%dc_min)) call wrong('dc_min')
   if (para%dc_max == huge(para%dc_max)) call wrong('dc_max')
   if (para%axis(1) == huge(para%axis)) call wrong('axis_x')
   if (para%axis(2) == huge(para%axis)) call wrong('axis_y')
   if (para%axis(3) == huge(para%axis)) call wrong('axis_z')
   if (para%translate == huge(para%translate)) call wrong('translate')
   if (para%rotate == huge(para%rotate)) call wrong('rotate')
   if (para%invert == huge(para%invert)) call wrong('invert')
   if (para%preserve_groups == huge(para%preserve_groups)) call wrong('preserve_groups')
   if (para%velocity(1) == huge(para%velocity)) call wrong('velocity_x')
   if (para%velocity(2) == huge(para%velocity)) call wrong('velocity_y')
   if (para%velocity(3) == huge(para%velocity)) call wrong('velocity_x')
   
   ! check parameter ranges
   if (para%L<=0) call error('L must be larger than 0.')
   if (para%length_unit<=0) call error('length_unit must be larger than 0.')
   if (para%snapshot_min>para%snapshot_max) call error('snapshot_min must be smaller than snapshot_max.')
   if (para%subsnapshot_min>para%subsnapshot_max) call error('subsnapshot_min must be smaller than subsnapshot_max.')
   if (para%h<=0) call error('h must be larger than 0.')
   if (para%OmegaL<0) call error('OmegaL must be >=0.')
   if (para%OmegaM<0) call error('OmegaM must be >=0.')
   if (para%OmegaL>1) call error('OmegaL must be <=1.')
   if (para%OmegaM>1) call error('OmegaM must be <=1.')
   if (para%angle<=0) call error('angle must be larger than 0.')
   if (.not.islogical(para%square_base)) call error('square_base can only be 0 or 1.')
   if (para%dc_min<0) call error('dc_min must be >=0.')
   if (para%dc_max<=0) call error('dc_min must be >0.')
   if (para%dc_max<=para%dc_min) call error('dc_min must smaller than dc_max.')
   if (sqrt(sum(para%axis**2))<epsilon(para%axis)) call error('The cone axis is ill-defined.')
   if (.not.islogical(para%translate)) call error('translate can only be 0 or 1.')
   if (.not.islogical(para%rotate)) call error('rotate can only be 0 or 1.')
   if (.not.islogical(para%invert)) call error('invert can only be 0 or 1.')
   if (.not.islogical(para%preserve_groups)) call error('preserve_groups can only be 0 or 1.')
   if (sqrt(sum(para%velocity**2))>0.1*c) call error('The observer velocity cannot exceed 0.1*c.')
   
   contains
   
   subroutine wrong(str)
      implicit none
      character(*),intent(in) :: str
      call error('Parameter '//str//' not initialized')
   end subroutine wrong
   
   logical function islogical(i)
      implicit none
      integer*4,intent(in) :: i
      islogical = ((i==0).or.(i==1))
   end function islogical

end subroutine check_parameters

subroutine adjust_parameters

   implicit none
   para%angle = min(2*pi,para%angle)
   para%dc_min = max(0.0,para%dc_min)
   para%axis = para%axis/sqrt(sum(para%axis**2))

end subroutine adjust_parameters

subroutine load_user_parameters(parameter_filename)

   ! loads all parameters from text file, except for paths

   implicit none
   character(255),intent(in)  :: parameter_filename
   character(255)             :: line
   character(50)              :: var_name
   character(205)             :: var_value,var_value_first
   integer                    :: io
   logical                    :: manual
   
   call check_exists(parameter_filename)
   
   open(1,file=trim(parameter_filename),action='read',form='formatted')
   do
      read(1,'(A)',IOSTAT=io) line
      if (io.ne.0) exit
      if (.not.((len(trim(line))==0).or.(line(1:1)=='#'))) then
         read(line,*) var_name
         var_value = adjustl(line(len(trim(var_name))+1:len(line)))
         read(var_value,*) var_value_first
         manual = (trim(var_value_first).ne.'auto')
         select case (trim(var_name))
            case ('L')
               if (manual) read(var_value,*) para%L
            case ('length_unit')
               if (manual) read(var_value,*) para%length_unit
            case ('snapshot_min')
               if (manual) read(var_value,*) para%snapshot_min
            case ('snapshot_max')
               if (manual) read(var_value,*) para%snapshot_max
            case ('subsnapshot_min')
               if (manual) read(var_value,*) para%subsnapshot_min
            case ('subsnapshot_max')
               if (manual) read(var_value,*) para%subsnapshot_max
            case ('h')
               if (manual) read(var_value,*) para%h
            case ('OmegaL')
               if (manual) read(var_value,*) para%OmegaL
            case ('OmegaM')
               if (manual) read(var_value,*) para%OmegaM
            case ('angle')
               if (manual) read(var_value,*) para%angle
            case ('square_base')
               read(var_value,*) para%square_base
            case ('dc_max')
               if (manual) read(var_value,*) para%dc_max
            case ('dc_min')
               if (manual) read(var_value,*) para%dc_min
            case ('axis.x')
               if (manual) read(var_value,*) para%axis(1)
            case ('axis.y')
               if (manual) read(var_value,*) para%axis(2)
            case ('axis.z')
               if (manual) read(var_value,*) para%axis(3)
            case ('translate')
               if (manual) read(var_value,*) para%translate
            case ('rotate')
               if (manual) read(var_value,*) para%rotate
            case ('invert')
               if (manual) read(var_value,*) para%invert
            case ('preserve_groups')
               if (manual) read(var_value,*) para%preserve_groups
            case ('velocity.x')
               if (manual) read(var_value,*) para%velocity(1)
            case ('velocity.y')
               if (manual) read(var_value,*) para%velocity(2)
            case ('velocity.z')
               if (manual) read(var_value,*) para%velocity(3)
            case default
               if ((trim(var_name).ne.'path_input').and.(trim(var_name).ne.'path_output')) then
                  call out('ERROR: '//trim(var_name)//' is an unknown parameter.')
                  stop
               end if
         end select
      end if
   end do
   close(1)
   
end subroutine load_user_parameters

subroutine load_paths(parameter_filename_custom)

   ! loads all parameters from text file, except for paths

   implicit none
   character(255),intent(in)  :: parameter_filename_custom
   character(255)             :: parameter_filename
   character(255)             :: line
   character(50)              :: var_name
   character(205)             :: var_value
   integer                    :: io
   
   if (len(trim(parameter_filename_custom))>0) then
      parameter_filename = parameter_filename_custom
   else
      parameter_filename = parameter_filename_default
   end if
   call check_exists(parameter_filename)
   
   para%path_output = ''
   para%path_input = ''
   
   open(1,file=trim(parameter_filename),action='read',form='formatted')
   do
      read(1,'(A)',IOSTAT=io) line
      if (io.ne.0) exit
      if (.not.((len(trim(line))==0).or.(line(1:1)=='#'))) then
         read(line,*) var_name
         var_value = adjustl(line(len(trim(var_name))+1:len(line)))
         select case (trim(var_name))
            case ('path_output')
               para%path_output = trim(noslash(var_value))//'/'
               if (len(trim(para%path_output))==0) then
                  call out('ERROR: Output path cannot be found in parameter file.')
                  stop
               end if
               call system('mkdir -p '//trim(para%path_output))
            case ('path_input')
               para%path_input = trim(noslash(var_value))//'/'
               if (len(trim(para%path_input))==0) then
                  call out('ERROR: Input path cannot be found in parameter file.')
                  stop
               end if
               call system('mkdir -p '//trim(para%path_input))
         end select
      end if
   end do
   close(1)
   
end subroutine load_paths

subroutine load_parameters

   implicit none
   character(len=255)   :: filename
   filename = trim(para%path_output)//'parameters.bin'
   call check_exists(filename)
   open(1,file=trim(filename),action='read',form="unformatted")
   read(1) para
   close(1)
   
end subroutine load_parameters

subroutine save_parameters

   implicit none
   
   character(len=255)   :: filename
   character(len=255)   :: txt
   
   ! open binary file
   filename = trim(para%path_output)//'parameters.bin'
   open(1,file=trim(filename),action='write',form="unformatted",status='replace')
   write(1) para
   close(1)
   
   ! open ascii file
   filename = trim(para%path_output)//'parameters.txt'
   open(1,file=trim(filename),action='write',form="formatted",status='replace')
   
   ! file names
   write(txt,*)  trim(para%path_output);        call line('path_output',txt)
   write(txt,*)  trim(para%path_input);         call line('path_input',txt)
   
   ! simulation box
   write(txt,'(E14.7)') para%L;                 call line('L',txt)
   write(txt,'(E14.7)') para%length_unit;       call line('length_unit',txt)
   write(txt,'(I6)')    para%snapshot_min;      call line('snapshot_min',txt)
   write(txt,'(I6)')    para%snapshot_max;      call line('snapshot_max',txt)
   write(txt,'(I6)')    para%subsnapshot_min;   call line('subsnapshot_min',txt)
   write(txt,'(I6)')    para%subsnapshot_max;   call line('subsnapshot_max',txt)
   
   ! cosmology
   write(txt,'(E14.7)') para%h;                 call line('h',txt)
   write(txt,'(E14.7)') para%OmegaL;            call line('OmegaL',txt)
   write(txt,'(E14.7)') para%OmegaM;            call line('OmegaM',txt)
   
   ! cone geometry
   write(txt,'(E14.7)') para%angle;             call line('angle',txt)
   write(txt,'(I1)')    para%square_base;       call line('square_base',txt)
   write(txt,'(E14.7)') para%dc_min;            call line('dc_min',txt)
   write(txt,'(E14.7)') para%dc_max;            call line('dc_max',txt)
   write(txt,'(E14.7)') para%axis(1);           call line('axis.x',txt)
   write(txt,'(E14.7)') para%axis(2);           call line('axis.y',txt)
   write(txt,'(E14.7)') para%axis(3);           call line('axis.z',txt)

   ! cone parameters
   write(txt,'(I1)')    para%translate;         call line('translate',txt)
   write(txt,'(I1)')    para%rotate;            call line('rotate',txt)
   write(txt,'(I1)')    para%invert;            call line('invert',txt)
   write(txt,'(I1)')    para%preserve_groups;   call line('preserve_groups',txt)
   
   ! observer
   write(txt,'(E14.7)') para%velocity(1);       call line('velocity.x',txt)
   write(txt,'(E14.7)') para%velocity(2);       call line('velocity.y',txt)
   write(txt,'(E14.7)') para%velocity(3);       call line('velocity.z',txt)
   
   ! close file
   close(1)
   
   contains
   
   subroutine line(name,value)
   
      implicit none
      character(*),intent(in) :: name,value
      character(len=30)       :: txtname
      write(txtname,'(A30)') trim(name)
      write(1,'(A30,A)') adjustl(txtname),trim(adjustl(value))
   
   end subroutine line
   
   integer*4 function log2int(l)
      implicit none
      logical,intent(in) :: l
      if (l) then
         log2int = 1
      else
         log2int = 0
      end if
   end function log2int

end subroutine save_parameters

end module module_parameters