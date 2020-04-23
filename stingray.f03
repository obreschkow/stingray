program stingray

   use shared_module_core
   use shared_module_arguments
   use shared_module_parameters
   use shared_module_cosmology
   use module_global
   use module_interface
   use module_user_routines
   use module_user_selection
   use module_parameters
   use module_tiling
   use module_sky

   implicit none
   
   ! start user interface
   call set_version('0.24')
   call handle_arguments(require_task=.false.)
   call start_output
   
   ! read user arguments and parameters
   call make_parameters
   call require_no_options_left
   
   ! initialise variables
   call initialize_global_variables
   call get_keywords
   call set_cosmology('stingray',para%h,para%omega_m,para%omega_l)
   call assign_selection_function
   
   ! main tasks
   call make_tiling
   call make_sky
   call write_sky_to_hdf5
   
   ! finalize output on screen/logfile
   call stop_output(delete_logfile=.not.int2log(para%keep_log))
    
end program stingray