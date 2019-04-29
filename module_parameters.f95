module module_parameters

   use module_constants
   use module_system
   use module_types
   use module_io
   use module_linalg
   use module_user_routines
   
   private
   public   :: make_parameters, load_paths
   
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
      parameter_filename = get_parameter_filename_default()
   end if
   call out('File: '//trim(parameter_filename))
   call check_exists(parameter_filename)
   call reset_parameters
   call make_automatic_parameters
   call load_user_parameters(parameter_filename)
   call check_parameters
   call adjust_parameters
   call make_derived_parameters
   call save_parameters(parameter_filename)
   
   call toc

end subroutine make_parameters

subroutine reset_parameters

   ! resets all parameters other than path_output and path_input

   implicit none
   para%survey = ''
   para%L = huge(para%L)
   para%length_unit = huge(para%length_unit)
   para%snapshot_min = huge(para%snapshot_min)
   para%snapshot_max = huge(para%snapshot_max)
   para%subvolume_min = huge(para%subvolume_min)
   para%subvolume_max = huge(para%subvolume_max)
   para%h = huge(para%h)
   para%omega_l = huge(para%omega_l)
   para%omega_m = huge(para%omega_m)
   para%omega_b = huge(para%omega_b)
   para%dc_min = huge(para%dc_min)
   para%dc_max = huge(para%dc_max)
   para%ra_min = huge(para%ra_min)
   para%ra_max = huge(para%ra_max)
   para%dec_min = huge(para%dec_min)
   para%dec_max = huge(para%dec_max)
   para%zaxis_ra = huge(para%zaxis_ra)
   para%zaxis_dec = huge(para%zaxis_dec)
   para%xy_angle = huge(para%xy_angle)
   para%seed = huge(para%seed)
   para%translate = huge(para%translate)
   para%rotate = huge(para%rotate)
   para%invert = huge(para%invert)
   para%velocity_ra = huge(para%velocity_ra)
   para%velocity_dec = huge(para%velocity_dec)
   para%velocity_norm = huge(para%velocity_norm)
   para%search_angle = huge(para%search_angle)
   para%volume_search_level= huge(para%volume_search_level)

end subroutine reset_parameters

subroutine check_parameters

   implicit none
   
   ! check if all parameters have been initialized
   if (trim(para%survey)=='') call wrong('survey')
   if (trim(para%path_output)=='') call wrong('path_output')
   if (trim(para%path_input)=='') call wrong('path_input')
   if (para%L == huge(para%L)) call wrong('L')
   if (para%length_unit == huge(para%length_unit)) call wrong('length_unit')
   if (para%snapshot_min == huge(para%snapshot_min)) call wrong('snapshot_min')
   if (para%snapshot_max == huge(para%snapshot_max)) call wrong('snapshot_max')
   if (para%subvolume_min == huge(para%subvolume_min)) call wrong('subvolume_min')
   if (para%subvolume_max == huge(para%subvolume_max)) call wrong('subvolume_max')
   if (para%h == huge(para%h)) call wrong('h')
   if (para%omega_l == huge(para%omega_l)) call wrong('omega_l')
   if (para%omega_m == huge(para%omega_m)) call wrong('omega_m')
   if (para%omega_b == huge(para%omega_b)) call wrong('omega_b')
   if (para%dc_min == huge(para%dc_min)) call wrong('dc_min')
   if (para%dc_max == huge(para%dc_max)) call wrong('dc_max')
   if (para%ra_min == huge(para%dc_min)) call wrong('ra_min')
   if (para%ra_max == huge(para%dc_max)) call wrong('ra_max')
   if (para%dec_min == huge(para%dc_min)) call wrong('dec_min')
   if (para%dec_max == huge(para%dc_max)) call wrong('dec_max')
   if (para%zaxis_ra == huge(para%zaxis_ra)) call wrong('zaxis_ra')
   if (para%zaxis_dec == huge(para%zaxis_dec)) call wrong('zaxis_dec')
   if (para%xy_angle == huge(para%xy_angle)) call wrong('xy_angle')
   if (para%seed == huge(para%seed)) call wrong('seed')
   if (para%translate == huge(para%translate)) call wrong('translate')
   if (para%rotate == huge(para%rotate)) call wrong('rotate')
   if (para%invert == huge(para%invert)) call wrong('invert')
   if (para%velocity_ra == huge(para%velocity_ra)) call wrong('velocity_ra')
   if (para%velocity_dec == huge(para%velocity_dec)) call wrong('velocity_dec')
   if (para%velocity_norm == huge(para%velocity_norm)) call wrong('velocity_norm')
   if (para%search_angle == huge(para%search_angle)) call wrong('search_angle')
   if (para%volume_search_level == huge(para%volume_search_level)) call wrong('volume_search_level')
   
   ! check parameter ranges
   if (para%L<=0) call error('L must be larger than 0')
   if (para%length_unit<=0) call error('length_unit must be larger than 0')
   if (para%snapshot_min>para%snapshot_max) call error('snapshot_min must be smaller than snapshot_max')
   if (para%subvolume_min>para%subvolume_max) call error('subvolume_min must be smaller than subvolume_max')
   if (para%h<=0) call error('h must be larger than 0')
   if (para%omega_l<0) call error('omega_l must be >=0')
   if (para%omega_m<0) call error('omega_m must be >=0')
   if (para%omega_b<0) call error('omega_b must be >=0')
   if (para%omega_l>1) call error('omega_l must be <=1')
   if (para%omega_m>1) call error('omega_m must be <=1')
   if (para%omega_b>para%omega_m) call error('omega_b must be <= omega_m')
   if (para%dc_min<0) call error('dc_min must be >=0')
   if (para%dc_max<=0) call error('dc_min must be >0')
   if (para%dc_max<=para%dc_min) call error('dc_min must smaller than dc_max')
   if ((para%ra_min<0).or.(para%ra_min>360.0)) call error('ra_min must lie between 0 and 360')
   if ((para%ra_max<0).or.(para%ra_max>360.0)) call error('ra_max must lie between 0 and 360')
   if (para%ra_min>=para%ra_max) call error('ra_max must be larger than ra_min')
   if ((para%dec_min<-180.0).or.(para%dec_min>180.0)) call error('dec_min must lie between -180 and 180')
   if ((para%dec_max<-180.0).or.(para%dec_max>180.0)) call error('dec_max must lie between -180 and 180')
   if (para%dec_min>=para%dec_max) call error('dec_max must be larger than dec_min')
   if ((para%zaxis_ra<0).or.(para%zaxis_ra>360.0)) call error('zaxis_ra must lie between 0 and 360')
   if ((para%zaxis_dec<-180.0).or.(para%zaxis_dec>180.0)) call error('zaxis_dec must lie between -180 and 180')
   if ((para%xy_angle<0).or.(para%xy_angle>360.0)) call error('xy_angle must lie between 0 and 360')
   if (para%seed<=0) call error('seed must be a positive integer.')
   if (.not.islogical(para%translate)) call error('translate can only be 0 or 1')
   if (.not.islogical(para%rotate)) call error('rotate can only be 0 or 1')
   if (.not.islogical(para%invert)) call error('invert can only be 0 or 1')
   if ((para%velocity_ra<0).or.(para%velocity_ra>360.0)) call error('velocity_ra must lie between 0 and 360')
   if ((para%velocity_dec<-180.0).or.(para%velocity_dec>180.0)) &
   & call error('velocity_dec must lie between -180 and 180')
   if (para%velocity_norm>0.1*c) call error('The observer velocity cannot exceed 0.1*c')
   if (para%search_angle>=10.0) call error('search_angle must be smaller than 10')
   if (para%search_angle<=0.0) call error('search_angle must be larger than 0')
   if (para%volume_search_level<0) call error('volume_search_level must be equal to or larger than 0')
   if (para%volume_search_level>10) call error('volume_search_level cannot be larger than 10')
   
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
   
   ! convert degrees to radian
   para%ra_min = para%ra_min*degree
   para%ra_max = para%ra_max*degree
   para%dec_min = para%dec_min*degree
   para%dec_max = para%dec_max*degree
   para%zaxis_ra = para%zaxis_ra*degree
   para%zaxis_dec = para%zaxis_dec*degree
   para%xy_angle = para%xy_angle*degree
   para%velocity_ra = para%velocity_ra*degree
   para%velocity_dec = para%velocity_dec*degree
   para%search_angle = para%search_angle*degree

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
            case ('survey')
               if (manual) read(var_value,*) para%survey
            case ('L')
               if (manual) read(var_value,*) para%L
            case ('length_unit')
               if (manual) read(var_value,*) para%length_unit
            case ('snapshot_min')
               if (manual) read(var_value,*) para%snapshot_min
            case ('snapshot_max')
               if (manual) read(var_value,*) para%snapshot_max
            case ('subvolume_min')
               if (manual) read(var_value,*) para%subvolume_min
            case ('subvolume_max')
               if (manual) read(var_value,*) para%subvolume_max
            case ('h')
               if (manual) read(var_value,*) para%h
            case ('omega_l')
               if (manual) read(var_value,*) para%omega_l
            case ('omega_m')
               if (manual) read(var_value,*) para%omega_m
            case ('omega_b')
               if (manual) read(var_value,*) para%omega_b
            case ('dc_min')
               if (manual) read(var_value,*) para%dc_min
            case ('dc_max')
               if (manual) read(var_value,*) para%dc_max
            case ('ra_min')
               if (manual) read(var_value,*) para%ra_min
            case ('ra_max')
               if (manual) read(var_value,*) para%ra_max
            case ('dec_min')
               if (manual) read(var_value,*) para%dec_min
            case ('dec_max')
               if (manual) read(var_value,*) para%dec_max
            case ('zaxis_ra')
               if (manual) read(var_value,*) para%zaxis_ra
            case ('zaxis_dec')
               if (manual) read(var_value,*) para%zaxis_dec
            case ('xy_angle')
               if (manual) read(var_value,*) para%xy_angle
            case ('seed')
               if (manual) read(var_value,*) para%seed
            case ('translate')
               if (manual) read(var_value,*) para%translate
            case ('rotate')
               if (manual) read(var_value,*) para%rotate
            case ('invert')
               if (manual) read(var_value,*) para%invert
            case ('velocity_ra')
               if (manual) read(var_value,*) para%velocity_ra
            case ('velocity_dec')
               if (manual) read(var_value,*) para%velocity_dec
            case ('velocity_norm')
               if (manual) read(var_value,*) para%velocity_norm
            case ('search_angle')
               if (manual) read(var_value,*) para%search_angle
            case ('volume_search_level')
               if (manual) read(var_value,*) para%volume_search_level
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
      parameter_filename = get_parameter_filename_default()
   end if
   
   if (.not.exists(parameter_filename,.true.)) then
      write(*,*) 'ERROR: parameter file does not exist. Please update the variable parameter_filename_default'
      write(*,*) '       in the user module or specify a correct filename using the option -parameterfile.'
      stop
   end if
   
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

   call sph2car(para%velocity_norm,para%velocity_ra,para%velocity_dec,para%velocity_car)
   call make_sky_rotation

contains

   subroutine make_sky_rotation

      implicit none
      real*4            :: rotationvector(3)
      real*4            :: angle
      real*4            :: nrot
      real*4,parameter  :: ezsam(3) = (/0.0,0.0,1.0/) ! unit vector of the SAM z-axis inSAM coordinates
      real*4            :: ezsky(3) ! unit vector of the SAM z-axis in Survey coordinates
      
      ! rotate SAM coordinates in the SAM xy-plane
      para%sky_rotation = rotation_matrix(ezsam,para%xy_angle)
      
      ! rotate SAM coordinates onto Survey coordinates
      call sph2car(1.0,para%zaxis_ra,para%zaxis_dec,ezsky)
      rotationvector = cross_product(ezsam,ezsky)
      nrot = norm(rotationvector)
      
      if ((nrot>epsilon(nrot)).and.(sum(ezsam*ezsky)<1.0)) then
         rotationvector = rotationvector/nrot
         angle = acos(sum(ezsam*ezsky))
         para%sky_rotation = matmul(rotation_matrix(rotationvector,angle),para%sky_rotation)
      end if
   
   end subroutine make_sky_rotation

end subroutine make_derived_parameters

end module module_parameters