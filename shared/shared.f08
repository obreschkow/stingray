program shared

   use shared_module_core
   use shared_module_constants
   use shared_module_arguments
   use shared_module_assignable_types
   use shared_module_parameters
   use shared_module_hdf5
   use shared_module_vectors
   use shared_module_maths
   use shared_module_lapack
   use shared_module_graphics
   use shared_module_cosmology

   implicit none
   
   ! start user interface
   call handle_arguments(require_task=.true.) ! also deals with -version, -verbose, -logfile
   call start_output
   
   ! run task
   if (istask('core',require_value=.false.,require_options=.false.)) then
   
   else if  (istask('constants',.false.,.false.)) then
      call example_constants
   else if (istask('vectors',.false.,.false.)) then
      call example_vectors
   else if (istask('maths',.false.,.false.)) then
      call out('4.safedivide.5                   = ',4.safedivide.5.0_8)
      call out('4.safedivide.0 (+Inf)            = ',4.safedivide.0)
      call out('-4.safedivide.0 (-Inf)           = ',-4.safedivide.0)
      call out('0.safedivide.0 (undefined)       = ',0.safedivide.0)
      call out('1e20.safedivide.1e-20 (overflow) = ',1e20.safedivide.1e-20)
      call out('1e20_8.safedivide.1e-20          = ',1e20_8.safedivide.1e-20_8)
   else
      call unknown_task
   end if
   
   ! finalize output on screen/logfile
   call stop_output
   
end program shared