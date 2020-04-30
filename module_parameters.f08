module module_parameters

   use shared_module_core
   use shared_module_arguments
   use shared_module_parameters
   use shared_module_maths
   use shared_module_constants
   use module_global
   use module_user_routines
   
   private
   public   :: make_parameters
   
contains

! Note: the routine make_automatic_parameters needs to be specified in the user module

subroutine make_parameters

   implicit none
   character(255)    :: parameter_filename
   character(255)    :: parameter_set
   type(type_para)   :: auto,empty

   
   call tic('MAKE PARAMETERS')
   
   ! process user arguments
   call get_option_value(parameter_filename,'-parameterfile',fn_parameters_default)
   call get_option_value(parameter_set,'-parameterset','')
   call set_parameterset(parameter_set)
   
   ! read parameter file
   call out('File: '//trim(parameter_filename))
   call set_parameterfile(trim(parameter_filename))
   call set_autostring('auto')
   call handle_parameters
   
   ! get input path
   call get_parameter_value(para%path_input,'path_input')
   para%path_input = dir(para%path_input,ispath=.true.)
   call check_file(para%path_input,permission='r')
   
   ! make automatic parameters, as specified in the module "module_user_routines_..."
   empty = para
   call make_automatic_parameters
   auto = para
   para = empty
   
   ! interpret parameters of parameter file
   call get_parameter_value(para%survey,'survey',auto=auto%survey,min=1)
   call get_parameter_value(para%path_output,'path_output',auto=auto%path_output,min=1)
   call get_parameter_value(para%filename_sky,'filename_sky',auto=auto%filename_sky,min=1)
   call get_parameter_value(para%box_side,'box_side',auto=auto%box_side,min=0.0)
   call get_parameter_value(para%length_unit,'length_unit',auto=auto%length_unit,min=0.0)
   call get_parameter_value(para%snapshot_min,'snapshot_min',auto=auto%snapshot_min,min=0)
   call get_parameter_value(para%snapshot_max,'snapshot_max',auto=auto%snapshot_max,min=0)
   call get_parameter_value(para%subvolume_min,'subvolume_min',auto=auto%subvolume_min,min=0)
   call get_parameter_value(para%subvolume_max,'subvolume_max',auto=auto%subvolume_max,min=0)
   call get_parameter_value(para%h,'h',auto=auto%h,min=0.1)
   call get_parameter_value(para%omega_l,'omega_l',auto=auto%omega_l,min=0.0,max=1.0)
   call get_parameter_value(para%omega_m,'omega_m',auto=auto%omega_m,min=0.0,max=1.0)
   call get_parameter_value(para%omega_b,'omega_b',auto=auto%omega_b,min=0.0,max=1.0)
   call get_parameter_value(para%zaxis_ra,'zaxis_ra',auto=auto%zaxis_ra,min=0.0,max=360.0)
   call get_parameter_value(para%zaxis_dec,'zaxis_dec',auto=auto%zaxis_dec,min=-90.0,max=90.0)
   call get_parameter_value(para%xy_angle,'xy_angle',auto=auto%xy_angle,min=0.0,max=360.0)
   call get_parameter_value(para%seed,'seed',auto=auto%seed,min=1)
   call get_parameter_value(para%translate,'translate',auto=auto%translate,min=0,max=1)
   call get_parameter_value(para%rotate,'rotate',auto=auto%rotate,min=0,max=1)
   call get_parameter_value(para%invert,'invert',auto=auto%invert,min=0,max=1)
   call get_parameter_value(para%velocity_ra,'velocity_ra',auto=auto%velocity_ra,min=0.0,max=360.0)
   call get_parameter_value(para%velocity_dec,'velocity_dec',auto=auto%velocity_dec,min=-90.0,max=90.0)
   call get_parameter_value(para%velocity_norm,'velocity_norm',auto=auto%velocity_norm,min=0.0)
   call get_parameter_value(para%search_angle,'search_angle',auto=auto%search_angle,min=0.0,max=10.0)
   call get_parameter_value(para%volume_search_level,'volume_search_level',auto=auto%volume_search_level,min=0.0,max=10.0)
   call get_parameter_value(para%make_groups,'make_groups',auto=auto%make_groups,min=0,max=1)
   call get_parameter_value(para%options,'options',auto=auto%options)
   call get_parameter_value(para%merge_output,'merge_output',auto=auto%merge_output,min=0,max=1)
   call get_parameter_value(para%keep_binaries,'keep_binaries',auto=auto%keep_binaries,min=0,max=1)
   call get_parameter_value(para%keep_log,'keep_log',auto=auto%keep_log,min=0,max=1)
   call require_no_parameters_left
   
   ! checks other than simple ranges
   if (para%snapshot_min>para%snapshot_max) call error('snapshot_min must be <= snapshot_max')
   if (para%subvolume_min>para%subvolume_max) call error('subvolume_min must be <= subvolume_max')
   if (para%omega_b>para%omega_m) call error('omega_b must be <= omega_m')
   if (para%velocity_norm>0.1*const%c) call error('observer velocity cannot exceed 0.1*c')
   
   ! make output path
   para%path_output = dir(para%path_output,ispath=.true.)
   call make_path(para%path_output)
   call check_file(para%path_output,permission='rw')
   
   ! convert user units to internal units
   call convert_parameter_units
   
   ! compute derived parameters, such as the sky rotation matrix
   call make_derived_parameters
   
   call toc

end subroutine make_parameters

subroutine convert_parameter_units

   implicit none
   
   ! convert degrees to radian
   para%zaxis_ra = para%zaxis_ra*unit%degree
   para%zaxis_dec = para%zaxis_dec*unit%degree
   para%xy_angle = para%xy_angle*unit%degree
   para%velocity_ra = para%velocity_ra*unit%degree
   para%velocity_dec = para%velocity_dec*unit%degree
   para%search_angle = para%search_angle*unit%degree

end subroutine convert_parameter_units

subroutine make_derived_parameters

   call sph2car(para%velocity_norm,para%velocity_ra,para%velocity_dec,para%velocity_car,astro=.true.)
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
      call sph2car(1.0,para%zaxis_ra,para%zaxis_dec,ezsky,astro=.true.)
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