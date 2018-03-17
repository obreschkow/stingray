module module_parameters

   use module_constants
   use module_types
   use module_system
   
   type(type_para)   :: para
   
contains

subroutine initialize_default_parameters

   implicit none
   para%x_min = -1e20
   para%x_max = +1e20
   para%y_min = -1e20
   para%y_max = +1e20
   para%angle = pi
   para%axis = (/0,0,1/)
   para%translate = .true.
   para%rotate = .true.
   para%preserve_groups = .true.
   para%velocity = (/0,0,0/)
   para%subsnapshot_min = 0
   para%subsnapshot_max = 0

end subroutine initialize_default_parameters

subroutine check_and_adjust_parameters

   implicit none
   para%axis = para%axis/sqrt(sum(para%axis**2))
   para%dc_min = max(0.0,para%dc_min)

end subroutine check_and_adjust_parameters

subroutine load_parameters(parameterfile)

   implicit none
   character(255),intent(in)  :: parameterfile
   character(255)             :: line
   character(50)              :: var_name
   character(205)             :: var_value
   real*4                     :: nb
   integer                    :: io
   logical                    :: file_exists
   
   inquire(file=trim(parameterfile),exist=file_exists)
   if (.not.file_exists) then
      call out('Error: The parameter file does not exist:')
      call out(trim(parameterfile))
      stop
   end if
   
   open(1,file=trim(parameterfile),action='read',form='formatted')
   do
      read(1,'(A)',IOSTAT=io) line
      if (io.ne.0) exit
      if (.not.((len(trim(line))==0).or.(line(1:1)=='#'))) then
         read(line,*) var_name
         var_value = adjustl(line(len(trim(var_name))+1:len(line)))
         select case (trim(var_name))
            case ('path_output')
               para%path_output = trim(noslash(var_value))//'/'
               call system('mkdir -p '//trim(para%path_output))
            case ('path_input')
               para%path_input = trim(noslash(var_value))//'/'
               call system('mkdir -p '//trim(para%path_input))
            case ('L')
               read(var_value,*) para%L
            case ('length_unit')
               read(var_value,*) para%length_unit
            case ('snapshot_min')
               read(var_value,*) para%snapshot_min
            case ('snapshot_max')
               read(var_value,*) para%snapshot_max
            case ('subsnapshot_min')
               read(var_value,*) para%subsnapshot_min
            case ('subsnapshot_max')
               read(var_value,*) para%subsnapshot_max
            case ('h')
               read(var_value,*) para%h
            case ('OmegaL')
               read(var_value,*) para%OmegaL
            case ('OmegaM')
               read(var_value,*) para%OmegaM
            case ('angle')
               read(var_value,*) para%angle
            case ('x_min')
               read(var_value,*) para%x_min
            case ('x_max')
               read(var_value,*) para%x_max
            case ('y_min')
               read(var_value,*) para%y_min
            case ('y_max')
               read(var_value,*) para%y_max
            case ('dc_max')
               read(var_value,*) para%dc_max
            case ('dc_min')
               read(var_value,*) para%dc_min
            case ('axis.x')
               read(var_value,*) para%axis(1)
            case ('axis.y')
               read(var_value,*) para%axis(2)
            case ('axis.z')
               read(var_value,*) para%axis(3)
            case ('translate')
               read(var_value,*) nb
               para%translate = (nb==1)
            case ('rotate')
               read(var_value,*) nb
               para%rotate = (nb==1)
            case ('invert')
               read(var_value,*) nb
               para%invert = (nb==1)
            case ('preserve_groups')
               read(var_value,*) nb
               para%preserve_groups = (nb==1)
            case ('velocity.x')
               read(var_value,*) para%velocity(1)
            case ('velocity.y')
               read(var_value,*) para%velocity(2)
            case ('velocity.z')
               read(var_value,*) para%velocity(3)
            case default
               call out('ERROR: '//trim(var_name)//' is an unknown parameter.')
               stop
         end select
      end if
   end do
   close(1)
   
end subroutine load_parameters

subroutine save_parameters

   implicit none
   
   character(len=255)   :: filename
   character(len=255)   :: txt
   
   ! open file
   filename = trim(para%path_output)//'parameters.txt'
   open(1,file=trim(filename),action='write',form="formatted",status='replace')
   
   ! file names
   write(txt,*)  trim(para%path_output);              call line('path_output',txt)
   write(txt,*)  trim(para%path_input);               call line('path_input',txt)
   
   ! simulation box
   write(txt,'(E14.7)')  para%L;                      call line('L',txt)
   write(txt,'(E14.7)')  para%length_unit;            call line('length_unit',txt)
   write(txt,'(I6)')  para%snapshot_min;              call line('snapshot_min',txt)
   write(txt,'(I6)')  para%snapshot_max;              call line('snapshot_max',txt)
   write(txt,'(I6)')  para%subsnapshot_min;           call line('subsnapshot_min',txt)
   write(txt,'(I6)')  para%subsnapshot_max;           call line('subsnapshot_max',txt)
   
   ! cosmology
   write(txt,'(E14.7)')  para%h;                      call line('h',txt)
   write(txt,'(E14.7)')  para%OmegaL;                 call line('OmegaL',txt)
   write(txt,'(E14.7)')  para%OmegaM;                 call line('OmegaM',txt)
   
   ! cone geometry
   write(txt,'(E14.7)')  para%angle;                  call line('angle',txt)
   write(txt,'(E14.7)')  para%dc_min;                 call line('dc_min',txt)
   write(txt,'(E14.7)')  para%dc_max;                 call line('dc_max',txt)
   write(txt,'(E14.7)')  para%axis(1);                call line('axis.x',txt)
   write(txt,'(E14.7)')  para%axis(2);                call line('axis.y',txt)
   write(txt,'(E14.7)')  para%axis(3);                call line('axis.z',txt)

   ! cone parameters
   write(txt,'(I1)')  log2int(para%translate);        call line('translate',txt)
   write(txt,'(I1)')  log2int(para%rotate);           call line('rotate',txt)
   write(txt,'(I1)')  log2int(para%invert);           call line('invert',txt)
   write(txt,'(I1)')  log2int(para%preserve_groups);  call line('preserve_groups',txt)
   
   ! observer
   write(txt,'(E14.7)')  para%velocity(1);            call line('velocity.x',txt)
   write(txt,'(E14.7)')  para%velocity(2);            call line('velocity.y',txt)
   write(txt,'(E14.7)')  para%velocity(3);            call line('velocity.z',txt)
   
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