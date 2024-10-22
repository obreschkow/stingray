module module_parameters

   use shared_module_core
   use shared_module_arguments
   use shared_module_parameters
   use shared_module_maths
   use shared_module_vectors
   use shared_module_constants
   use shared_module_hdf5
   use module_global
   
   public   :: read_parameterfile
   public   :: initialize_parameters
   public   :: write_hdf5_parameters
   public   :: option
   public   :: devoption
   public   :: para
   public   :: make_tmp_path
   
   private
   
   integer*4,parameter           :: n_keywords_max = 20
   integer*4,protected           :: n_keywords
   character(len=255),protected  :: keyword_name(n_keywords_max)
   integer*4,protected           :: n_devwords
   character(len=255),protected  :: devword_name(n_keywords_max)

   type type_para

      ! name of mock survey
      character(len=255)   :: survey
      
      ! paths & file names
      character(len=255)   :: path_input
      character(len=255)   :: path_output
      character(len=255)   :: filename_sky
   
      ! simulation box
      real*4               :: box_side     ! [length unit of simulation] box side length
      real*4               :: length_unit  ! [m] simulation length unit expressed in SI unit
      integer*4            :: snapshot_min
      integer*4            :: snapshot_max
      integer*4            :: subvolume_min
      integer*4            :: subvolume_max
   
      ! cosmology
      real*4               :: h
      real*4               :: omega_l
      real*4               :: omega_m
      real*4               :: omega_b

      ! large-scale structure randomisation
      character(7)         :: randomisation ! randomisation type
      logical*4            :: shell_tweaking
      character(4)         :: prng
      integer*4            :: seed  ! seed of random number generator (integer >=1)
      logical*4            :: translate
      logical*4            :: rotate
      logical*4            :: invert
      
      ! observer position
      logical*4            :: fix_observer_position
      real*4               :: observer_x  ! [length unit]
      real*4               :: observer_y  ! [length unit]
      real*4               :: observer_z  ! [length unit]
      
      ! observer rotation
      logical*4            :: fix_observer_rotation
      real*4               :: xaxis_ra  ! [rad]
      real*4               :: xaxis_dec ! [rad]
      real*4               :: yz_angle  ! [rad]
      
      ! observer velocity relative to CMB
      real*4               :: velocity_ra     ! [rad]
      real*4               :: velocity_dec    ! [rad]
      real*4               :: velocity_norm   ! [km/s] peculiar velocity of observer with respect to Hubble flow
      
      ! advanced options
      real*4               :: search_angle         ! [rad] minimal angular separation of points on faces
      real*4               :: volume_search_level
      
      ! galaxy options
      logical*4            :: make_groups
      character(len=255)   :: options
      
      ! I/O options
      logical*4            :: merge_output
      logical*4            :: keep_tmp
      logical*4            :: keep_log
      
      ! derived parameters, not directly specified by the user (no default needed)
      real*4               :: velocity_car(3)   ! [km/s] velocity of observer cartesian survey-coordinates
      real*4               :: observer_translation(3) ! [box side length] fixed observer translation
      real*4               :: observer_rotation(3,3) ! rotation matrix to move the (x,y,z)-sky axis onto the central (RA,dec)-sky
      logical*4            :: modern_prng ! logical flag, true iff prng==F95
      
      ! hidden developer parameters
      character(len=255)   :: devoptions
      
   end type type_para
   
   type(type_para),protected     :: para
   
contains

subroutine read_parameterfile

   implicit none
   character(255)    :: parameter_filename
   character(255)    :: parameter_set
   
   call tic('MAKE PARAMETERS')
   
   ! set parameter file name
   call get_option_value(parameter_filename,'-parameterfile',fn_parameters_default)
   call set_parameterfile(trim(parameter_filename))
   
   ! set parameter set name
   call get_option_value(parameter_set,'-parameterset','')
   call set_parameterset(parameter_set)
   
   ! read parameter file
   call out('File: '//trim(parameter_filename))
   call read_parameters
   
   ! get input path
   call get_parameter_value(para%path_input,'path_input')
   para%path_input = dir(para%path_input,ispath=.true.)
   call check_file(para%path_input,permission='r')
   
end subroutine read_parameterfile

subroutine initialize_parameters

   implicit none
   
   ! interpret parameters of parameter file
   call get_parameter_value(para%survey,'survey',min=1)
   call get_parameter_value(para%path_output,'path_output',min=1)
   call get_parameter_value(para%filename_sky,'filename_sky',min=1)
   call get_parameter_value(para%box_side,'box_side',min=0.0)
   call get_parameter_value(para%length_unit,'length_unit',min=0.0)
   call get_parameter_value(para%snapshot_min,'snapshot_min',min=0)
   call get_parameter_value(para%snapshot_max,'snapshot_max',min=0)
   call get_parameter_value(para%subvolume_min,'subvolume_min',min=0)
   call get_parameter_value(para%subvolume_max,'subvolume_max',min=0)
   call get_parameter_value(para%h,'h',min=0.1)
   call get_parameter_value(para%omega_l,'omega_l',min=0.0,max=1.0)
   call get_parameter_value(para%omega_m,'omega_m',min=0.0,max=1.0)
   call get_parameter_value(para%omega_b,'omega_b',min=0.0,max=1.0)
   call get_parameter_value(para%randomisation,'randomisation')
   call get_parameter_value(para%shell_tweaking,'shell_tweaking')
   call get_parameter_value(para%prng,'prng')
   call get_parameter_value(para%seed,'seed',min=1)
   call get_parameter_value(para%translate,'translate')
   call get_parameter_value(para%rotate,'rotate')
   call get_parameter_value(para%invert,'invert')
   call get_parameter_value(para%fix_observer_position,'fix_observer_position')
   call get_parameter_value(para%observer_x,'observer_x',min=0.0)
   call get_parameter_value(para%observer_y,'observer_y',min=0.0)
   call get_parameter_value(para%observer_z,'observer_z',min=0.0)
   call get_parameter_value(para%fix_observer_rotation,'fix_observer_rotation')
   call get_parameter_value(para%xaxis_ra,'xaxis_ra',min=0.0,max=360.0)
   call get_parameter_value(para%xaxis_dec,'xaxis_dec',min=-90.0,max=90.0)
   call get_parameter_value(para%yz_angle,'yz_angle',min=0.0,max=360.0)
   call get_parameter_value(para%velocity_ra,'velocity_ra',min=0.0,max=360.0)
   call get_parameter_value(para%velocity_dec,'velocity_dec',min=-90.0,max=90.0)
   call get_parameter_value(para%velocity_norm,'velocity_norm',min=0.0)
   call get_parameter_value(para%search_angle,'search_angle',min=0.0,max=10.0)
   call get_parameter_value(para%volume_search_level,'volume_search_level',min=0.0,max=10.0)
   call get_parameter_value(para%make_groups,'make_groups')
   call get_parameter_value(para%options,'options')
   call get_parameter_value(para%merge_output,'merge_output')
   call get_parameter_value(para%keep_tmp,'keep_tmp')
   call get_parameter_value(para%keep_log,'keep_log')
   call get_parameter_value(para%devoptions,'devoptions','')
   call require_no_parameters_left
   
   ! checks other than simple ranges
   if (para%snapshot_min>para%snapshot_max) call error('snapshot_min must be <= snapshot_max')
   if (para%subvolume_min>para%subvolume_max) call error('subvolume_min must be <= subvolume_max')
   if (para%omega_b>para%omega_m) call error('omega_b must be <= omega_m')
   if (para%velocity_norm>0.1*const%c) call error('observer velocity cannot exceed 0.1*c')
   if (para%observer_x>para%box_side) call error('observer_x must not be larger than box_side')
   if (para%observer_y>para%box_side) call error('observer_y must not be larger than box_side')
   if (para%observer_z>para%box_side) call error('observer_z must not be larger than box_side')
   if (.not.any(para%randomisation==(/'single ','tiles  ','shells '/))) &
   & call error('parameter "randomisation" must be either "single", "tiles", "shells"')
   if (.not.any(para%prng==(/'F77 ','F95 '/))) call error('parameter "prng" must be either "F77" or "F95"')
   
   ! make output path
   para%path_output = dir(para%path_output,ispath=.true.)
   call make_path(para%path_output)
   call check_file(para%path_output,permission='rw')
   
   ! convert user units to internal units
   call convert_parameter_units
   
   ! compute derived parameters, such as the sky rotation matrix
   call make_derived_parameters
   call get_keywords
   call get_devwords
   
   call toc

end subroutine initialize_parameters

subroutine convert_parameter_units

   implicit none
   
   ! convert degrees to radian
   para%xaxis_ra = para%xaxis_ra*unit%degree
   para%xaxis_dec = para%xaxis_dec*unit%degree
   para%yz_angle = para%yz_angle*unit%degree
   para%velocity_ra = para%velocity_ra*unit%degree
   para%velocity_dec = para%velocity_dec*unit%degree
   para%search_angle = para%search_angle*unit%degree

end subroutine convert_parameter_units

subroutine make_derived_parameters

   call sph2car(para%velocity_norm,para%velocity_ra,para%velocity_dec,para%velocity_car)
   call fix_observer
   para%modern_prng = trim(para%prng)=="F95"

contains

   subroutine fix_observer

      implicit none
      real*4            :: rotationaxis(3)
      real*4            :: angle
      real*4,parameter  :: exsam(3) = (/1.0,0.0,0.0/) ! unit vector along the x-axis
      real*4            :: exsky(3) ! expression of exsam in sky coordinates
      
      if (para%fix_observer_position) then
      
         para%observer_translation = modulo(0.5-(/para%observer_x,para%observer_y,para%observer_z/)/para%box_side,1.0)
      
      else
      
         para%observer_translation = 0.0
      
      end if
      
      if (para%fix_observer_rotation) then
      
         ! rotate SAM coordinates in the SAM xy-plane
         para%observer_rotation = rotation_matrix(exsam,para%yz_angle)
      
         ! rotate SAM coordinates onto survey coordinates
         call sph2car(1.0,para%xaxis_ra,para%xaxis_dec,exsky)
         rotationaxis = exsam.cross.exsky
      
         if ((norm(rotationaxis)>epsilon(1.0)).and.(abs(exsam.dot.exsky)<1.0)) then
            angle = acos(exsam.dot.exsky)
            para%observer_rotation = matmul(rotation_matrix(rotationaxis,angle),para%observer_rotation)
         end if
         
      else
      
         para%observer_rotation = 0.0
         
      end if
   
   end subroutine fix_observer

end subroutine make_derived_parameters

subroutine get_keywords
   
   implicit none
   integer*4   :: i,imin,imax,n
   
   n_keywords = 0
   
   if (.not.isempty(para%options)) then
   
      n = len(trim(para%options))
      imin = 1
      do i = 1,n
         if ((para%options(i:i)==',').or.(i==n)) then
            if (para%options(i:i)==',') then
               imax = i-1
            else
               imax = i
            end if
            if (imax<imin) call error('cannot interpret parameter "options"')
            n_keywords = n_keywords+1
            keyword_name(n_keywords) = para%options(imin:imax)
            imin = imax+2
         end if
      end do
      
   end if
   
end subroutine get_keywords

subroutine get_devwords
   
   implicit none
   integer*4   :: i,imin,imax,n
   
   n_devwords = 0
   
   if (.not.isempty(para%devoptions)) then
   
      n = len(trim(para%devoptions))
      imin = 1
      do i = 1,n
         if ((para%devoptions(i:i)==',').or.(i==n)) then
            if (para%devoptions(i:i)==',') then
               imax = i-1
            else
               imax = i
            end if
            if (imax<imin) call error('cannot interpret parameter "devoptions"')
            n_devwords = n_devwords+1
            devword_name(n_devwords) = para%devoptions(imin:imax)
            imin = imax+2
         end if
      end do
      
   end if
   
end subroutine get_devwords

logical function option(string)

   implicit none
   character(*),intent(in) :: string
   integer*4               :: i
   
   option = .false.
   do i = 1,n_keywords
      if (trim(keyword_name(i))==trim(string)) then
         option = .true.
         exit
      end if
   end do
   
end function option

logical function devoption(string)

   implicit none
   character(*),intent(in) :: string
   integer*4               :: i
   
   devoption = .false.
   do i = 1,n_devwords
      if (trim(devword_name(i))==trim(string)) then
         devoption = .true.
         exit
      end if
   end do
   
end function devoption

subroutine make_tmp_path
   
   ! should be called before PRNG seed is initialised
   
   implicit none
   real  :: x
   call cpu_time(x)
   write(path_tmp,'(A,A,I0.15,A)') trim(para%path_output),'tmp_', &
   & int(mod(x,1.0)*1e10,8)+10000000000_8*get_random_integer(0,99999,.true.),'/'
   call make_path(path_tmp)
   
end subroutine make_tmp_path

subroutine write_hdf5_parameters(filename_hdf5)

   implicit none
   character(*),intent(in)    :: filename_hdf5  ! output filename
   character(:),allocatable   :: name
   
   allocate(character(1)::name) ! empty allocation to avoid compiler flags

   ! open HDF5 file
   call hdf5_open(filename_hdf5,.true.)
   
   ! write group "parameters"
   name = 'parameters/'
   call hdf5_add_group(name)
   call hdf5_write_data(name//'survey',trim(para%survey),'name of simulated survey')
   call hdf5_write_data(name//'path_output',trim(para%path_output))
   call hdf5_write_data(name//'path_input',trim(para%path_input))
   call hdf5_write_data(name//'para%filename_sky',trim(para%filename_sky), &
   & 'filename of mock sky file(s), without extension .hdf5 and without subvolume index')
   call hdf5_write_data(name//'box_side',para%box_side, &
   & '[length_unit] comoving side-length of simulation box in multiples of length_unit')
   call hdf5_write_data(name//'length_unit',para%length_unit,'[m] SI-value of comoving length unit')
   call hdf5_write_data(name//'snapshot_min',para%snapshot_min,'index of earliest snapshot used for the mock sky')
   call hdf5_write_data(name//'snapshot_max',para%snapshot_max,'index of latest snapshot used for the mock sky')
   call hdf5_write_data(name//'subvolume_min',para%subvolume_min,'index of first subvolume used for the mock sky')
   call hdf5_write_data(name//'subvolume_max',para%subvolume_max,'index of last subvolume used for the mock sky')
   call hdf5_write_data(name//'h',para%h,'[-] Hubble parameter H0=h*100 km/s/Mpc')
   call hdf5_write_data(name//'omega_l',para%omega_l,'energy density of dark energy relative to closure density')
   call hdf5_write_data(name//'omega_m',para%omega_m,'energy density of all matter relative to closure density')
   call hdf5_write_data(name//'omega_b',para%omega_b,'energy density of baryonic matter relative to closure density')
   call hdf5_write_data(name//'randomisation',trim(para%randomisation),'method used to randomize cosmic large-scale structure')
   call hdf5_write_data(name//'shell_tweaking',para%shell_tweaking, &
   & 'logical flag specifying if the shell radii are tweaked to match snapshot boundaries; used only if randomisation=shells')
   call hdf5_write_data(name//'prng',trim(para%prng),'type of pseudo-random number generator')
   call hdf5_write_data(name//'seed',para%seed,'seed for the random number generator of symmetry operations')
   call hdf5_write_data(name//'translate',para%translate, &
   & 'logical flag specifying if random translations are applied (0=false, 1=true)')
   call hdf5_write_data(name//'rotate',para%rotate, &
   & 'logical flag specifying if random rotations are applied (0=false, 1=true)')
   call hdf5_write_data(name//'invert',para%invert, &
   & 'logical flag specifying if random inversions are applied (0=false, 1=true)')
   call hdf5_write_data(name//'fixed_observer_position',para%fix_observer_rotation, &
   & 'logical flag specifying if the position of the observer in the nbody box is fixed (0=false, 1=true)')
   call hdf5_write_data(name//'observer_x',para%observer_x, &
   & '[length units] x-coordinate of fixed observere position in the N-body box')
   call hdf5_write_data(name//'observer_y',para%observer_y, &
   & '[length units] y-coordinate of fixed observere position in the N-body box')
   call hdf5_write_data(name//'observer_z',para%observer_z, &
   & '[length units] z-coordinate of fixed observere position in the N-body box')
   call hdf5_write_data(name//'fixed_observer_rotation',para%fix_observer_rotation, &
   & 'logical flag specifying if the rotation between the N-body and sky coordinates at the observer is fixed (0=false, 1=true)')
   call hdf5_write_data(name//'xaxis_ra',para%xaxis_ra/unit%degree, &
   & '[deg] RA coordinate of the SAM x-axis in spherical survey coordinates')
   call hdf5_write_data(name//'xaxis_dec',para%xaxis_dec/unit%degree, &
   & '[deg] Dec coordinate the SAM x-axis in spherical survey coordinates')
   call hdf5_write_data(name//'yz_angle',para%yz_angle/unit%degree,'[deg] Rotation of the SAM (y,z)-plane on the sky')
   call hdf5_write_data(name//'velocity_norm',para%velocity_norm, &
   & '[km/s] observer velocity relative to CMB rest-frame')
   call hdf5_write_data(name//'velocity_ra',para%velocity_ra, &
   & '[deg] RA coordinate to which the observer is moving relative to CMB rest-frame')
   call hdf5_write_data(name//'velocity_dec',para%velocity_dec, &
   & '[deg] Dec coordinate to which the observer is moving relative to CMB rest-frame')
   call hdf5_write_data(name//'search_angle',para%search_angle, &
   & '[deg] typical angle in which overlaps between survey volume and tiling grid are searched')
   call hdf5_write_data(name//'volume_search_level',para%volume_search_level, &
   & 'specifies the number of search points (2^#)^3 inside each tile')
   call hdf5_write_data(name//'options',trim(para%options), &
   & 'string of options specifying what properties are generated')
   call hdf5_write_data(name//'keep_tmp',para%keep_tmp, &
   & 'logical flag specifying if binary output files are kept in additino to this HDF5 (0=false, 1=true)')
   call hdf5_write_data(name//'keep_log',para%keep_tmp, &
   & 'logical flag specifying if the logfile is kept after successful runs (0=false, 1=true)')
   if (.not.isempty(para%devoptions)) call hdf5_write_data(name//'devoptions',trim(para%devoptions), &
   & 'string of developer options')

   ! close HDF5 file
   call hdf5_close()

end subroutine write_hdf5_parameters

end module module_parameters