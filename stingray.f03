program stingray

   use shared_module_core
   use shared_module_arguments
   use shared_module_parameters
   use shared_module_cosmology
   use module_global
   use module_user_routines
   use module_user_selection
   use module_parameters
   use module_tiling
   use module_sky

   implicit none

   character(255) :: parameter_filename
   
   ! start user interface
   call set_version('0.22')
   call handle_arguments(require_task=.false.)
   call start_output
   
   ! read parameters
   call get_option_value(parameter_filename,'-parameterfile',fn_parameters_default)
   call make_parameters(parameter_filename)
   
   ! initialise global variables
   call initialize_global_variables
   call set_cosmology('stingray',para%h,para%omega_m,para%omega_l)
   call assign_selection_function
   
   ! main tasks
   call make_tiling
   call make_sky
   call write_sky_to_hdf5
   
   ! finalize output on screen/logfile
   call stop_output(delete_logfile=.not.int2log(para%keep_log))
    
end program stingray