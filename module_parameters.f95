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
   para%ra_min = huge(para%ra_min)
   para%ra_max = huge(para%ra_max)
   para%dec_min = huge(para%dec_min)
   para%dec_max = huge(para%dec_max)
   para%dc_min = huge(para%dc_min)
   para%dc_max = huge(para%dc_max)
   para%axis = huge(para%axis)
   para%turn = huge(para%turn)
   para%seed = huge(para%seed)
   para%translate = huge(para%translate)
   para%rotate = huge(para%rotate)
   para%invert = huge(para%invert)
   para%preserve_groups = huge(para%preserve_groups)
   para%velocity = huge(para%velocity)

end subroutine reset_parameters

subroutine check_parameters

   implicit none
   integer*4   :: mode1,mode2
   
   ! check if all parameters have been initialized
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
   if (para%OmegaB<0) call error('OmegaB must be >=0.')
   if (para%OmegaL>1) call error('OmegaL must be <=1.')
   if (para%OmegaM>1) call error('OmegaM must be <=1.')
   if (para%OmegaB>para%OmegaM) call error('OmegaB must be <= OmegaM.')
   if (para%dc_min<0) call error('dc_min must be >=0.')
   if (para%dc_max<=0) call error('dc_min must be >0.')
   if (para%dc_max<=para%dc_min) call error('dc_min must smaller than dc_max.')
   if (norm(para%axis)<epsilon(para%axis)) call error('The cone axis is ill-defined.')
   if (para%seed<=0) call error('seed must be a positive integer.')
   if (.not.islogical(para%translate)) call error('translate can only be 0 or 1.')
   if (.not.islogical(para%rotate)) call error('rotate can only be 0 or 1.')
   if (.not.islogical(para%invert)) call error('invert can only be 0 or 1.')
   if (.not.islogical(para%preserve_groups)) call error('preserve_groups can only be 0 or 1.')
   if (norm(para%velocity)>0.1*c) call error('The observer velocity cannot exceed 0.1*c.')
   
   ! handle foot print on the sky
   mode1 = 0
   if (para%ra .ne. huge(para%ra)) mode1 = mode1+1
   if (para%dec .ne. huge(para%dec)) mode1 = mode1+1
   if (para%angle .ne. huge(para%angle)) mode1 = mode1+1
   
   mode2 = 0
   if (para%ra_min .ne. huge(para%ra_min)) mode2 = mode2+1
   if (para%ra_max .ne. huge(para%ra_max)) mode2 = mode2+1
   if (para%dec_min .ne. huge(para%dec_min)) mode2 = mode2+1
   if (para%dec_max .ne. huge(para%dec_max)) mode2 = mode2+1
   
   if ((mode1==3).and.(mode2==0)) then
      if ((para%ra<0).or.(para%ra>=360.0)) call error('ra must lie between 0 and 360')
      if ((para%dec<-180.0).or.(para%dec>180.0)) call error('dec must lie between -180 and 180')
      if (para%angle<=0) call error('angle must be larger than 0.')
   else if ((mode1==0).and.(mode2==4)) then
      if ((para%ra_min<0).or.(para%ra_min>360.0)) call error('ra_min must lie between 0 and 360')
      if ((para%ra_max<0).or.(para%ra_max>360.0)) call error('ra_max must lie between 0 and 360')
      if (para%ra_min>=para%ra_max) call error('ra_min must smaller than ra_max')
      if ((para%dec_min<-180.0).or.(para%dec_min>180.0)) call error('dec_min must lie between -180 and 180')
      if ((para%dec_max<-180.0).or.(para%dec_max>180.0)) call error('dec_max must lie between -180 and 180')
      if (para%dec_min>=para%dec_max) call error('dec_min must smaller than dec_max')
   else
      call error('Specify all of ra/dec/angle or all of ra_min,ra_max,dec_min,dec_max, but not both.')
   end if
   
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
   real*4   :: angle(4)
   real*4   :: xcenter(3)
   real*4   :: xcorner(3)
   
   para%dc_min = max(0.0,para%dc_min)
   para%axis = para%axis/norm(para%axis)
   para%turn = para%turn*degree
   
   ! make ra,dec,angle if ra/dec range are provide as input
   if (para%dec_min<=180.0) then
      ! convert deg to rad
      para%ra_min = para%ra_min*degree
      para%ra_max = para%ra_max*degree
      para%dec_min = para%dec_min*degree
      para%dec_max = para%dec_max*degree
      
      ! make cone parameters
      para%ra = (para%ra_min+para%ra_max)/2
      para%dec = (para%dec_min+para%dec_max)/2
      xcenter = (/cos(para%dec)*sin(para%ra),sin(para%dec),cos(para%dec)*cos(para%ra)/)
      xcorner = (/cos(para%dec_min)*sin(para%ra_min),sin(para%dec_min),cos(para%dec_min)*cos(para%ra_min)/)
      angle(1) = acos(min(1.0,sum(xcenter*xcorner)))
      xcorner = (/cos(para%dec_min)*sin(para%ra_max),sin(para%dec_min),cos(para%dec_min)*cos(para%ra_max)/)
      angle(2) = acos(min(1.0,sum(xcenter*xcorner)))
      xcorner = (/cos(para%dec_max)*sin(para%ra_min),sin(para%dec_max),cos(para%dec_max)*cos(para%ra_min)/)
      angle(3) = acos(min(1.0,sum(xcenter*xcorner)))
      xcorner = (/cos(para%dec_max)*sin(para%ra_max),sin(para%dec_max),cos(para%dec_max)*cos(para%ra_max)/)
      angle(4) = acos(min(1.0,sum(xcenter*xcorner)))
      para%angle = maxval(angle)
   else
      ! convert deg to rad
      para%ra = para%ra*degree
      para%dec = para%dec*degree
      para%angle = para%angle*degree
      para%angle = min(pi,para%angle)
   end if

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
         manual = (trim(var_value_first).ne.'auto').and.(trim(var_value_first).ne.'Auto') &
         & .and.(trim(var_value_first).ne.'AUTO').and.(trim(var_value_first).ne.'a') &
         & .and.(trim(var_value_first).ne.'A')
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
            case ('OmegaB')
               if (manual) read(var_value,*) para%OmegaB
            case ('ra')
               if (manual) read(var_value,*) para%ra
            case ('dec')
               if (manual) read(var_value,*) para%dec
            case ('angle')
               if (manual) read(var_value,*) para%angle
            case ('ra_min')
               if (manual) read(var_value,*) para%ra_min
            case ('ra_max')
               if (manual) read(var_value,*) para%ra_max
            case ('dec_min')
               if (manual) read(var_value,*) para%dec_min
            case ('dec_max')
               if (manual) read(var_value,*) para%dec_max
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
   
   ! open ascii file (only for non-derived parameters)
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
   write(txt,'(E14.7)') para%OmegaB;            call line('OmegaB',txt)
   
   ! cone geometry
   if (para%dec_max<=180.0) then
      write(txt,'(E14.7)') para%ra_min/degree;  call line('ra_min',txt)
      write(txt,'(E14.7)') para%ra_max/degree;  call line('ra_max',txt)
      write(txt,'(E14.7)') para%dec_min/degree; call line('dec_min',txt)
      write(txt,'(E14.7)') para%dec_max/degree; call line('dec_max',txt)
   else
      write(txt,'(E14.7)') para%ra/degree;      call line('ra',txt)
      write(txt,'(E14.7)') para%dec/degree;     call line('dec',txt)
      write(txt,'(E14.7)') para%angle/degree;   call line('angle',txt)
   end if
   write(txt,'(E14.7)') para%dc_min;            call line('dc_min',txt)
   write(txt,'(E14.7)') para%dc_max;            call line('dc_max',txt)
   write(txt,'(E14.7)') para%axis(1);           call line('axis.x',txt)
   write(txt,'(E14.7)') para%axis(2);           call line('axis.y',txt)
   write(txt,'(E14.7)') para%axis(3);           call line('axis.z',txt)
   write(txt,'(E14.7)') para%turn/degree;       call line('turn',txt)

   ! cone parameters
   write(txt,'(I1)')    para%seed;              call line('seed',txt)
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