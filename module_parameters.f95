module module_parameters

   use module_constants
   use module_types
   use module_system
   use module_user
   
   private
   public   :: make_parameters, load_parameters, load_paths
   
contains

! Note: the routine make_automatic_parameters needs to be specified in the user module

subroutine make_parameters(parameter_filename_custom)

   implicit none
   character(255),intent(in)  :: parameter_filename_custom
   character(255)             :: parameter_filename
   
   call tic
   call out('MAKE PARAMETERS')
   
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
   call make_derived_parameters
   call save_parameters
   
   call toc

end subroutine make_parameters

subroutine reset_parameters

   ! resets all parameters other than path_output and path_input

   implicit none
   para%name = ''
   para%L = huge(para%L)
   para%length_unit = huge(para%length_unit)
   para%snapshot_min = huge(para%snapshot_min)
   para%snapshot_max = huge(para%snapshot_max)
   para%subsnapshot_min = huge(para%subsnapshot_min)
   para%subsnapshot_max = huge(para%subsnapshot_max)
   para%h = huge(para%h)
   para%OmegaL = huge(para%OmegaL)
   para%OmegaM = huge(para%OmegaM)
   para%OmegaB = huge(para%OmegaB)
   para%ra = huge(para%ra)
   para%dec = huge(para%dec)
   para%angle = huge(para%angle)
   para%dc_min = huge(para%dc_min)
   para%dc_max = huge(para%dc_max)
   para%axis = huge(para%axis)
   para%turn = huge(para%turn)
   para%seed = huge(para%seed)
   para%translate = huge(para%translate)
   para%rotate = huge(para%rotate)
   para%invert = huge(para%invert)
   para%velocity = huge(para%velocity)

end subroutine reset_parameters

subroutine check_parameters

   implicit none
   
   ! check if all parameters have been initialized
   if (trim(para%name)=='') call wrong('name')
   if (trim(para%path_output)=='') call wrong('path_output')
   if (trim(para%path_input)=='') call wrong('path_input')
   if (para%L == huge(para%L)) call wrong('L')
   if (para%length_unit == huge(para%length_unit)) call wrong('length_unit')
   if (para%snapshot_min == huge(para%snapshot_min)) call wrong('snapshot_min')
   if (para%snapshot_max == huge(para%snapshot_max)) call wrong('snapshot_max')
   if (para%subsnapshot_min == huge(para%subsnapshot_min)) call wrong('subsnapshot_min')
   if (para%subsnapshot_max == huge(para%subsnapshot_max)) call wrong('subsnapshot_max')
   if (para%h == huge(para%h)) call wrong('h')
   if (para%OmegaL == huge(para%OmegaL)) call wrong('OmegaL')
   if (para%OmegaM == huge(para%OmegaM)) call wrong('OmegaM')
   if (para%OmegaB == huge(para%OmegaB)) call wrong('OmegaB')
   if (para%dc_min == huge(para%dc_min)) call wrong('dc_min')
   if (para%dc_max == huge(para%dc_max)) call wrong('dc_max')
   if (para%axis(1) == huge(para%axis)) call wrong('axis_x')
   if (para%axis(2) == huge(para%axis)) call wrong('axis_y')
   if (para%axis(3) == huge(para%axis)) call wrong('axis_z')
   if (para%turn == huge(para%turn)) call wrong('turn')
   if (para%seed == huge(para%seed)) call wrong('seed')
   if (para%translate == huge(para%translate)) call wrong('translate')
   if (para%rotate == huge(para%rotate)) call wrong('rotate')
   if (para%invert == huge(para%invert)) call wrong('invert')
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
   if (para%OmegaB<0) call error('OmegaB must be >=0.')
   if (para%OmegaL>1) call error('OmegaL must be <=1.')
   if (para%OmegaM>1) call error('OmegaM must be <=1.')
   if (para%OmegaB>para%OmegaM) call error('OmegaB must be <= OmegaM.')
   if (para%dc_min<0) call error('dc_min must be >=0.')
   if (para%dc_max<=0) call error('dc_min must be >0.')
   if ((para%ra<0).or.(para%ra>=360.0)) call error('ra must lie between 0 and 360')
   if ((para%dec<-180.0).or.(para%dec>180.0)) call error('dec must lie between -180 and 180')
   if (para%angle<=0) call error('angle must be larger than 0.')
   if (para%dc_max<=para%dc_min) call error('dc_min must smaller than dc_max.')
   if (norm(para%axis)<epsilon(para%axis)) call error('The cone axis is ill-defined.')
   if (para%seed<=0) call error('seed must be a positive integer.')
   if (.not.islogical(para%translate)) call error('translate can only be 0 or 1.')
   if (.not.islogical(para%rotate)) call error('rotate can only be 0 or 1.')
   if (.not.islogical(para%invert)) call error('invert can only be 0 or 1.')
   if (norm(para%velocity)>0.1*c) call error('The observer velocity cannot exceed 0.1*c.')
   
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
   
   para%dc_min = max(0.0,para%dc_min)
   para%axis = para%axis/norm(para%axis)
   para%turn = para%turn*degree
   para%ra = para%ra*degree
   para%dec = para%dec*degree
   para%angle = para%angle*degree
   para%angle = min(pi,para%angle)

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
            case ('name')
               if (manual) read(var_value,*) para%name
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
            case ('OmegaB')
               if (manual) read(var_value,*) para%OmegaB
            case ('ra')
               if (manual) read(var_value,*) para%ra
            case ('dec')
               if (manual) read(var_value,*) para%dec
            case ('angle')
               if (manual) read(var_value,*) para%angle
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
            case ('turn')
               if (manual) read(var_value,*) para%turn
            case ('seed')
               if (manual) read(var_value,*) para%seed
            case ('translate')
               if (manual) read(var_value,*) para%translate
            case ('rotate')
               if (manual) read(var_value,*) para%rotate
            case ('invert')
               if (manual) read(var_value,*) para%invert
            case ('velocity.x')
               if (manual) read(var_value,*) para%velocity(1)
            case ('velocity.y')
               if (manual) read(var_value,*) para%velocity(2)
            case ('velocity.z')
               if (manual) read(var_value,*) para%velocity(3)
            case default
               if ((trim(var_name).ne.'path_input').and.(trim(var_name).ne.'path_output')) then
                  call error(trim(var_name)//' is an unknown parameter.')
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
               if (len(trim(var_value))==0) then
                  write(*,*) 'ERROR: parameter path_output empty.'
                  stop
               end if
               para%path_output = trim(noslash(var_value))//'/'
               call system('mkdir -p '//trim(para%path_output))
            case ('path_input')
               if (len(trim(var_value))==0) then
                  write(*,*) 'ERROR: parameter path_input empty.'
                  stop
               end if
               para%path_input = trim(noslash(var_value))//'/'
               if (.not.exists(trim(para%path_input),.true.)) then
                  write(*,*) 'ERROR: Input path cannot be found: '//trim(para%path_input)
                  stop
               end if
         end select
      end if
   end do
   close(1)
   
   ! check if paths set   
   if (trim(para%path_output)=='') then
      write(*,*) 'ERROR: parameter path_output missing in parameter file.'
      stop
   end if
   if (trim(para%path_input)=='') then
      write(*,*) 'ERROR: parameter path_input missing in parameter file.'
      stop
   end if
   
end subroutine load_paths

subroutine make_derived_parameters

   call make_sky_rotation

contains

   subroutine make_sky_rotation

      implicit none
      real*4      :: rotationvector(3)
      real*4      :: coneaxis(3)
      real*4      :: angle
      real*4      :: nrot
      
      ! turn around cone axis
      para%sky_rotation = rotation_matrix(para%axis,para%turn)
      
      ! rotate cone axis onto central RA and DEC
      coneaxis = (/cos(para%dec)*sin(para%ra),sin(para%dec),cos(para%dec)*cos(para%ra)/)
      rotationvector = cross_product(para%axis,coneaxis)
      nrot = norm(rotationvector)
      
      if ((nrot>epsilon(nrot)).and.(sum(coneaxis*para%axis)<1.0)) then
         rotationvector = rotationvector/nrot
         angle = acos(sum(coneaxis*para%axis))
         para%sky_rotation = matmul(rotation_matrix(rotationvector,angle),para%sky_rotation)
      end if
   
   end subroutine make_sky_rotation

end subroutine make_derived_parameters

end module module_parameters