program stingray

   use shared_module_core
   use shared_module_arguments
   use shared_module_parameters
   use shared_module_cosmology
   use shared_module_maths
   use module_global
   use module_parameters
   use module_user_routines
   use module_user_selection
   use module_tiling
   use module_sky

   implicit none
   
   ! start user interface
   call set_version('0.43')
   call handle_arguments(require_task=.false.)
   call start_output
   
   ! handle parameters
   call read_parameterfile
   call require_no_options_left
   call set_auto_string('auto')
   call make_auto_parameters
   
   ! initialise variables
   call initialize_parameters
   call set_cosmology('stingray',para%h,para%omega_m,para%omega_l)
   call assign_selection_function
   call make_tmp_path
   
   ! main tasks
   call make_tiling
   call make_sky
   call write_sky_to_hdf5
   
   ! finalize
   call remove_path(path_tmp)
   call stop_output(delete_logfile=.not.para%keep_log)
    
end program stingray