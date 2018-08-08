module module_user

! ==============================================================================================================
! IMPORT MODULES
! ==============================================================================================================

! default modules, do not edit
use module_constants
use module_system
use module_types
use module_io
use module_linalg
use module_cosmology
use module_conversion

! custom modules
use module_hdf5

private

public :: type_sam   ! SAM class, requires procedures get_position, get_groupid, is_selected, rotate_vectors
public :: type_sky_galaxy
public :: type_sky_group
public :: pos_selection
public :: get_parameter_filename_default
public :: make_automatic_parameters
public :: make_redshifts
public :: load_sam_snapshot
public :: custom_routines

! ==============================================================================================================
! SET DEFAULT PARAMETER FILENAME
! ==============================================================================================================

! Here, specify the default parameter filename, used when stingray is called without the optional argument
! "parameterfile"

character(len=255),parameter  :: parameter_filename_default = &
   & '/Users/do/Dropbox/Code/Fortran/stingray/parameters.txt'


! ==============================================================================================================
! TYPE DECLARATIONS
! ==============================================================================================================

! Here, specify the class of SAM-properties to be loaded for making the mock sky. How these properties are
! read from files is specified in the subroutine load_sam_snapshot below
! The class requires four class functions:
! get_position: returns xyz-position of the galaxy
! get_groupid: returns the group id of the galaxy
! is_selected: logical function specifying if the galaxy is selected
! is_group_center: logical function specifying if the galaxy is the group_center
! rotate_vectors: specification of how objects are rotated

type type_sam

   integer*8   :: id_galaxy      ! unique galaxy ID
   integer*8   :: id_halo        ! unique ID of parent halo
   integer*4   :: snapshot       ! snapshot ID
   integer*4   :: subvolume      ! subvolume index
   integer*4   :: typ            ! galaxy type (0=central, 1=satellite in halo, 2=orphan)
   real*4      :: position(3)    ! [Mpc/h] position of galaxy centre in simulation box
   real*4      :: velocity(3)    ! [proper km/s] peculiar velocity
   real*4      :: J(3)           ! [proper Msun/h pMpc/h km/s] angular momentum
   real*4      :: mstars_disk    ! [Msun/h] stellar mass disk
   real*4      :: mstars_bulge   ! [Msun/h] stellar mass bulge
   real*4      :: mgas_disk      ! [Msun/h] gas mass disk
   real*4      :: mgas_bulge     ! [Msun/h] gas mass bulge
   real*4      :: matom_disk     ! [Msun/h] atomic gas mass disk
   real*4      :: matom_bulge    ! [Msun/h] atomic gas mass bulge
   real*4      :: mmol_disk      ! [Msun/h] molecular gas mass disk
   real*4      :: mmol_bulge     ! [Msun/h] molecular gas mass bulge
   real*4      :: rstar_disk     ! [cMpc/h] half-mass radius of stars in the disk
   real*4      :: rstar_bulge    ! [cMpc/h] half-mass radius of stars in the bulge
   real*4      :: rgas_disk      ! [cMpc/h] half-mass radius of gas in the disk
   real*4      :: rgas_bulge     ! [cMpc/h] half-mass radius of gas in the bulge
   real*4      :: mvir_hosthalo  ! [Msun/h]
   real*4      :: mvir_subhalo   ! [Msun/h]
   real*4      :: cnfw_subhalo   ! concentration of NFW fit to subhalo
   
contains

   procedure   :: get_position      => sam_get_position     ! required function
   procedure   :: get_groupid       => sam_get_groupid      ! required function               
   procedure   :: is_selected       => sam_is_selected      ! required function
   procedure   :: is_group_center   => sam_is_group_center  ! required function
   procedure   :: rotate_vectors    => sam_rotate_vectors   ! required subroutine

end type type_sam

! Here, specify the class or classes of mock object properties in the sky, such as apparent galaxy properties.
! See the existing routines below for clarifications. The class must contain the following class functions:
! make_from_sam
! write_to_file
! is_selected
   
type type_sky

   integer*4   :: snapshot       ! snapshot ID
   integer*4   :: subvolume      ! subvolume index
   integer*4   :: tile           ! tile ID
   real*4      :: zobs           ! redshift in observer-frame
   real*4      :: zcmb           ! redshift in CMB frame
   real*4      :: zcos           ! cosmological redshift without peculiar motions
   real*4      :: dc             ! [simulation units = Mpc/h] comoving distance
   real*4      :: ra             ! [rad] right ascension
   real*4      :: dec            ! [rad] declination
   
   contains

   procedure   :: make_from_sam  => sky_make_from_sam    ! required subroutine
   procedure   :: write_to_file  => sky_write_to_file    ! required subroutine
   procedure   :: is_selected    => sky_is_selected      ! required logical function

end type type_sky

type,extends(type_sky) :: type_sky_galaxy ! must exist

   integer*8   :: id_galaxy      ! unique ID in the mock sky
   integer*8   :: id_group       ! unique group ID in the mock sky
   integer*8   :: id_galaxy_sam  ! galaxy ID in the SAM
   integer*8   :: id_halo_sam    ! galaxy parent halo ID in the SAM
   integer*4   :: typ            ! galaxy type (0=central, 1=satellite, 2=orphan)
   real*4      :: inclination    ! [rad] inclination
   real*4      :: pa             ! [rad] position angle from North to East
   real*4      :: mag            ! apparent magnitude (generic: M/L ratio of 1, no k-correction)
   real*8      :: SHI            ! [W/m^2] integrated HI line flux
   real*4      :: vpecrad        ! [proper km/s] radial peculiar velocity
   real*4      :: mstars         ! [Msun/h] total stellar mass
   real*4      :: mvir_hosthalo  ! [Msun/h] mass of 1st generation halo (i.e. direct host of type 0 galaxies)
   real*4      :: mvir_subhalo   ! [Msun/h] mass of 1st generation halo (i.e. direct host of type 0 galaxies)
   real*4      :: rstar_disk     ! [arcsec] apparent half-mass radius of stars in the disk
   real*4      :: rstar_bulge    ! [arcsec] apparent half-mass radius of stars in the bulge
   real*4      :: rgas_disk      ! [arcsec] apparent half-mass radius of stars in the disk
   
end type type_sky_galaxy

type,extends(type_sky) :: type_sky_group ! must exist
   
   integer*8   :: id_halo_sam    ! galaxy parent halo ID in the SAM
   integer*8   :: id_group       ! unique group ID in the mock sky
   real*4      :: mvir           ! [Msun/h] virial mass of group
   integer*4   :: group_ntot     ! total number of galaxies in group
   integer*4   :: group_nsel     ! number of selected galaxies in group
   integer*4   :: group_flag     ! 0=group complete, >0 group truncated
   
end type type_sky_group

! In order to place the objects in the mock sky, the class type_sam must have the following two functions
! enabling stingray to extract the position and group id of each object.

contains

function sam_get_position(self) result(position)
   class(type_sam) :: self
   real*4 :: position(3)
   position = self%position ! [length_unit of simulation] position if the galaxy in the box
end function sam_get_position

integer*8 function sam_get_groupid(self)
   class(type_sam) :: self
   sam_get_groupid = self%id_halo ! unique group identifier
end function sam_get_groupid

logical function sam_is_group_center(self)
   class(type_sam) :: self
   sam_is_group_center = self%typ==0
end function sam_is_group_center

! For instances of the user-defined sky-types to be saved, add all sky types in the following subroutine

subroutine sky_write_to_file(self,fileID)
   class(type_sky) :: self
   integer*4,intent(in) :: fileID
   select type (self)
   type is (type_sky_galaxy); write(fileID) self
   type is (type_sky_group); write(fileID) self
   class default
   call error('Unknown class for I/O.')
   end select
end subroutine sky_write_to_file


! ==============================================================================================================
! SELECTION OF GALAXIES IN THE MOCK SKY
! ==============================================================================================================

! selection function acting on comoving position, applied when making the tiling and intrinsic sky

logical function pos_selection(dc,ra,dec)
   
   implicit none
   real*4,intent(in) :: dc    ! [simulation units] comoving distance
   real*4,intent(in) :: ra    ! [deg] right ascension
   real*4,intent(in) :: dec   ! [deg] right ascension
   
   call nil(dc,ra,dec) ! dummy to avoid compiler warnings for unused arguments
   
   select case (trim(para%name))
   case ('test')
      pos_selection = ((ra>=5.0).and.(ra<=10.0).and.(dec>=-1.0).and.(dec<=3.0))
   case ('devils')
      pos_selection = ((ra>= 34.000).and.(ra<= 37.050).and.(dec>= -5.200).and.(dec<= -4.200)).or. &
                    & ((ra>= 52.263).and.(ra<= 53.963).and.(dec>=-28.500).and.(dec<=-27.500)).or. &
                    & ((ra>=149.380).and.(ra<=150.700).and.(dec>= +1.650).and.(dec<= +2.790))
   case ('gama')
      pos_selection = ((ra>= 30.200).and.(ra<= 38.800).and.(dec>=-10.250).and.(dec<= -3.720)).or. &
                    & ((ra> 129.000).and.(ra<=141.000).and.(dec>= -2.000).and.(dec<= +3.000)).or. &
                    & ((ra>=174.000).and.(ra<=186.000).and.(dec>= -3.000).and.(dec<= +2.000)).or. &
                    & ((ra>=211.500).and.(ra<=223.500).and.(dec>= -2.000).and.(dec<= +3.000)).or. &
                    & ((ra>=339.000).and.(ra<=351.000).and.(dec>=-35.000).and.(dec<=-30.000))
   case default
      pos_selection = .false.
   end select

end function pos_selection

! pre-selection function applied to the sam properties when placing the galaxies. More selections can be applied
! later in the function 'selected'. The function 'sam_selection' is simply to avoid carrying too many galaxies
! through the whole analysis

logical function sam_is_selected(sam) result(selected)

   implicit none
   class(type_sam),intent(in) :: sam
   
   call nil(sam)  ! dummy to avoid compiler warnings for unused arguments
   
   select case (trim(para%name))
   case ('test')
      selected = (sam%mstars_disk>1e8)
   case ('devils')
      selected = (sam%mstars_disk>1e8).or.((sam%mvir_hosthalo>1e11).and.(sam%typ==0))
   case ('gama')
      selected = (sam%mstars_disk>1e8).or.((sam%mvir_hosthalo>1e11).and.(sam%typ==0))
   case default
      selected = .true.
   end select
   
end function sam_is_selected

! selection applied to sky classes

logical function sky_is_selected(sky,sam) result(selected)

   class(type_sky)            :: sky ! self
   type(type_sam),intent(in)  :: sam
   real*4,parameter           :: dmag = 2.0
   
   call nil(sky,sam) ! dummy to avoid compiler warnings for unused arguments
   
   select type (sky)
   type is (type_sky_galaxy)
   
      select case (trim(para%name))
      case ('devils')
         selected = sky%mag<=21.2+dmag
      case ('gama')
         selected = ((sky%mag<=19.8+dmag).and.(sky%ra<330.0*degree)).or.(sky%mag<=19.2+dmag)
      case default
         selected = .true.
      end select
      
   type is (type_sky_group)
   
      selected = (sam%typ==0)
   
   class default
      call error('Unknown class for selection.')
   end select
   
end function sky_is_selected


! ==============================================================================================================
! CONVERSION BETWEEN INTRINSIC AND APPARENT GALAXY PROPERTIES
! ==============================================================================================================

subroutine sam_rotate_vectors(sam)

   ! This function rotates all the vector-properties of the galaxy, specified in type_sam, except for the
   ! position, which has already been processed when producing the intrinsic sky. For each such vector-property
   ! of the SAM, call
   ! sam%vector = rotate(sam%vector,ps)
   ! where ps is a logical argument that specifies whether the vector transforms like a normal vector,
   ! or like a pseudo-vector (also known as axial vector)

   implicit none
   class(type_sam),intent(inout)    :: sam
   
   sam%velocity   = rotate(sam%velocity,.false.)
   sam%J          = rotate(sam%J,.true.)
   
end subroutine sam_rotate_vectors

subroutine sky_make_from_sam(sky,sam,base,galaxyid,groupid,group_nselected)

   implicit none
   class(type_sky),intent(out)   :: sky
   type(type_sam),intent(in)     :: sam
   type(type_base),intent(in)    :: base              ! basic properties of the position of this galaxy in the sky
   integer*8,intent(in),optional :: galaxyid          ! (only for galaxy types) unique galaxy in sky index
   integer*8,intent(in),optional :: groupid           ! unique group in sky index
   integer*4,intent(in),optional :: group_nselected   ! (only for group types) number of galaxies selected in group
   real*4                        :: pos(3)            ! [simulation length units] position vector of galaxy
   real*4                        :: dl                ! [simulation length units] luminosity distance to observer
   real*4                        :: elos(3)           ! unit vector pointing from the observer to the object in comoving space
   real*4                        :: mHI
   
   call nil(sky,sam,base,galaxyid,groupid,group_nselected) ! dummy statement to avoid compiler warnings
   
   ! sky coordinates
   sky%dc  = base%dc*para%L   ! [Mpc/h]
   sky%ra  = base%ra          ! [rad]
   sky%dec = base%dec         ! [rad]
   
   ! position vector
   call sph2car(sky%dc,sky%ra,sky%dec,pos)
   elos = pos/norm(pos)
   
   ! make IDs
   sky%tile          = base%tile
   sky%snapshot      = sam%snapshot
   sky%subvolume     = sam%subvolume
   
   ! make redshift, provided the galaxy position [simulation units] galaxy-velocity [km/s]
   call make_redshift(pos*(para%length_unit/Mpc),sam%velocity,zobs=sky%zobs,zcmb=sky%zcmb,zcos=sky%zcos)
   
   select type(sky)
   type is (type_sky_galaxy)
   
      ! make IDs
      sky%id_galaxy_sam    = sam%id_galaxy   ! copy other properties
      sky%id_galaxy        = galaxyid        ! unique galaxy id
      sky%id_halo_sam      = sam%id_halo     ! ...
      sky%id_group         = groupid         ! unique group id
      
      ! copy properties
      sky%typ              = sam%typ
      sky%mstars           = sam%mstars_disk+sam%mstars_bulge
      sky%mvir_hosthalo    = sam%mvir_hosthalo
      sky%mvir_subhalo     = sam%mvir_subhalo
      
      ! make inclination and position angle [rad]
      call make_inclination_and_pa(pos,sam%J,inclination=sky%inclination,pa=sky%pa)
      
      ! appaprent propoerties
      dl = sky%dc*(1+sky%zobs) ! [Mpc/h]
      sky%mag = convert_absmag2appmag(convert_stellarmass2absmag(sky%mstars/para%h,1.0),dl/para%h)
      mHI = (sam%matom_disk+sam%matom_bulge)/1.35/para%h ! [Msun] HI mass
      sky%SHI = convert_luminosity2flux(real(mHI,8)*real(L2MHI,8)*Lsun,dl/para%h)
      sky%vpecrad = sum(sam%velocity*elos)
      
      ! apparent radii (note division by comoving distance, because intrinsic radii given in comoving units)
      sky%rstar_disk = sam%rstar_disk/sky%dc/degree*3600.0 ! [arcsec]
      sky%rstar_bulge = sam%rstar_bulge/sky%dc/degree*3600.0 ! [arcsec]
      sky%rgas_disk = sam%rgas_disk/sky%dc/degree*3600.0 ! [arcsec]
   
   type is (type_sky_group)
   
      sky%id_halo_sam  = sam%id_halo
      sky%id_group     = groupid          ! unique group id
      sky%mvir         = sam%mvir_hosthalo
      sky%group_ntot   = base%group_ntot
      sky%group_nsel   = group_nselected
      sky%group_flag   = base%group_flag
   
   class default
   end select
   
end subroutine sky_make_from_sam


! ==============================================================================================================
! IO routines
! ==============================================================================================================

! Set parameters automatically (e.g. from information provided with the snapshot files).
! These automatic parameters are only adopted if the values in the parameter files are set to 'auto'.
! Otherwise the parameter file overwrites these values.

subroutine make_automatic_parameters
   
   implicit none
   
   character(len=255)   :: filename
   integer*4            :: isnapshot,nsub
   
   isnapshot = -1
   filename = ''
   do while (.not.exists(filename,.true.).and.(isnapshot<10000))
      isnapshot = isnapshot+1
      write(filename,'(A,I0,A)') trim(para%path_input),isnapshot,'/0/galaxies.hdf5'
   end do
   if (isnapshot==10000) call error('No snapshots found in '//trim(para%path_input))
   para%snapshot_min = isnapshot
   
   do while (exists(filename,.true.))
      isnapshot = isnapshot+1
      write(filename,'(A,I0,A)') trim(para%path_input),isnapshot,'/0/galaxies.hdf5'
   end do
   para%snapshot_max = isnapshot-1
   
   write(filename,'(A,I0,A)') trim(para%path_input),para%snapshot_min,'/0/galaxies.hdf5'
   call hdf5_open(filename)
   call hdf5_read_data('/run_info/tot_n_subvolumes',nsub)
   para%subvolume_min = 0
   para%subvolume_max = nsub-1
   call hdf5_read_data('/cosmology/h',para%h)
   call hdf5_read_data('/cosmology/omega_l',para%omega_l)
   call hdf5_read_data('/cosmology/omega_m',para%omega_m)
   call hdf5_read_data('/cosmology/omega_b',para%omega_b)
   call hdf5_read_data('/run_info/lbox',para%L)
   para%length_unit = Mpc/para%h
   call hdf5_close()
   
end subroutine make_automatic_parameters

! load redshifts
! this routine must allocate the array snapshot and fill in its real*4-valued property 'redshift'

subroutine make_redshifts

   implicit none
   character(len=255)                              :: filename
   integer*4                                       :: isnapshot
   real*8                                          :: z
   
   do isnapshot = para%snapshot_min,para%snapshot_max
      write(filename,'(A,I0,A)') trim(para%path_input),isnapshot,'/0/galaxies.hdf5'
      call hdf5_open(filename) ! NB: this routine also checks if the file exists
      call hdf5_read_data('/run_info/redshift',z)
      snapshot(isnapshot)%redshift = real(z,4)
      call hdf5_close()
   end do
   
end subroutine make_redshifts

! load SAM snapshot file

subroutine load_sam_snapshot(index,subindex,sam)

   ! variable declaration
   implicit none
   integer*4,intent(in)                            :: index             ! snapshot index
   integer*4,intent(in)                            :: subindex          ! subindex, if the snapshot is split into several files
   type(type_sam),allocatable,intent(out)          :: sam(:)            ! derived type of all the relevant SAM properties
   character(len=255)                              :: filename
   integer*8                                       :: n
   character(*),parameter                          :: g = '/galaxies/'  ! group name
   
   ! open file
   write(filename,'(A,I0,A,I0,A)') trim(para%path_input),index,'/',subindex,'/galaxies.hdf5'
   call hdf5_open(filename) ! NB: this routine also checks if the file exists
   
   ! determine number of galaxies in this (sub)snapshot
   n = hdf5_dataset_size(g//'id_galaxy')

   ! allocate array
   if (allocated(sam)) deallocate(sam)
   allocate(sam(n))
   
   ! read file
   call hdf5_read_data(g//'id_galaxy',sam%id_galaxy)
   call hdf5_read_data(g//'id_halo',sam%id_halo)
   call hdf5_read_data(g//'type',sam%typ)
   call hdf5_read_data(g//'position_x',sam%position(1))
   call hdf5_read_data(g//'position_y',sam%position(2))
   call hdf5_read_data(g//'position_z',sam%position(3))
   call hdf5_read_data(g//'velocity_x',sam%velocity(1))
   call hdf5_read_data(g//'velocity_y',sam%velocity(2))
   call hdf5_read_data(g//'velocity_z',sam%velocity(3))
   call hdf5_read_data(g//'l_x',sam%J(1))
   call hdf5_read_data(g//'l_y',sam%J(2))
   call hdf5_read_data(g//'l_z',sam%J(3))
   call hdf5_read_data(g//'mstars_disk',sam%mstars_disk)
   call hdf5_read_data(g//'mstars_bulge',sam%mstars_bulge)
   call hdf5_read_data(g//'mgas_disk',sam%mgas_disk)
   call hdf5_read_data(g//'mgas_bulge',sam%mgas_bulge)
   call hdf5_read_data(g//'matom_disk',sam%matom_disk)
   call hdf5_read_data(g//'matom_bulge',sam%matom_bulge)
   call hdf5_read_data(g//'mmol_disk',sam%mmol_disk)
   call hdf5_read_data(g//'mmol_bulge',sam%mmol_bulge)
   call hdf5_read_data(g//'rstar_disk',sam%rstar_disk)
   call hdf5_read_data(g//'rstar_bulge',sam%rstar_bulge)
   call hdf5_read_data(g//'rgas_disk',sam%rgas_disk)
   call hdf5_read_data(g//'rgas_bulge',sam%rgas_bulge)
   call hdf5_read_data(g//'mvir_hosthalo',sam%mvir_hosthalo)
   call hdf5_read_data(g//'mvir_subhalo',sam%mvir_subhalo)
   call hdf5_read_data(g//'cnfw_subhalo',sam%cnfw_subhalo)
   
   ! assign other properties
   sam%snapshot = index
   sam%subvolume = subindex
   
   ! close file
   call hdf5_close()
   
end subroutine load_sam_snapshot

! return default parameter file name
character(255) function get_parameter_filename_default()
   get_parameter_filename_default = trim(parameter_filename_default)
end function get_parameter_filename_default


! ==============================================================================================================
! HANDLE CUSTOM TASKS (= optional subroutines and arguments)
! ==============================================================================================================

subroutine custom_routines(task,custom_option,success)

   implicit none
   character(*),intent(in) :: task
   character(*),intent(in) :: custom_option
   logical,intent(out)     :: success
   
   ! custom task handler
   success = .true.
   select case (trim(task))
   case ('my_task') ! dymmy argument to illustrate the use
      call out('Here, specify what to do as "my_task"')
      if (len(custom_option)>0) call out('Using the custom option: '//custom_option)
   case ('make_hdf5')
      call make_hdf5
   case ('my_additions_to_make_all')
      ! here add optional commands to be executed automatically at the end of a make.all run
      call make_hdf5
   case default
      success = .false.
   end select
   
contains

subroutine make_hdf5
   
   implicit none
   logical,parameter                   :: show_test_sum = .false.
   character(len=255)                  :: filename_bin
   character(len=255)                  :: filename_hdf5
   type(type_sky_galaxy),allocatable   :: sky_galaxy(:)
   type(type_sky_group),allocatable    :: sky_group(:)
   integer*8                           :: n,i
   character(len=255)                  :: name,str
   character(len=255)                  :: filename
   character(len=255)                  :: shark_version,shark_git_revision,shark_timestamp
   real*8                              :: test(10),expected_sum
   integer*8                           :: empty_int8
   real*4                              :: n_replica_mean
   integer*4                           :: n_replica_max
   
   call tic
   call out('CONVERT MOCK SKY FROM BINARY TO HDF5')
   
   ! load auxilary data
   call load_parameters
   call load_box_list
   call load_snapshot_list
   
   ! load auxilary data from shark output
   write(filename,'(A,I0,A)') trim(para%path_input),para%snapshot_min,'/0/galaxies.hdf5'
   call hdf5_open(filename)
   call hdf5_read_data('/run_info/shark_version',shark_version)
   call hdf5_read_data('/run_info/shark_git_revision',shark_git_revision)
   call hdf5_read_data('/run_info/timestamp',shark_timestamp)
   call hdf5_close()
   
   ! load auxiliary data on replication
   filename = trim(para%path_output)//'mocksky_intrinsic.bin'
   open(1,file=trim(filename),action='read',form='unformatted',status='old',access='stream')
   read(1) empty_int8,n_replica_mean,n_replica_max ! number of galaxies
   close(1)
   
   ! create HDF5 file
   filename_hdf5 = trim(para%path_output)//'mocksky.hdf5'
   call hdf5_create(filename_hdf5)
   
   ! open HDF5 file
   call hdf5_open(filename_hdf5,.true.)
   
   ! Group "Parameters"
   call hdf5_add_group('parameters')
   call hdf5_write_data('parameters/path_output',para%path_output)
   call hdf5_write_data('parameters/path_input',para%path_input)
   call hdf5_write_data('parameters/box_length',para%L)
   call hdf5_write_data('parameters/length_unit',para%length_unit)
   call hdf5_write_data('parameters/snapshot_min',para%snapshot_min)
   call hdf5_write_data('parameters/snapshot_max',para%snapshot_max)
   call hdf5_write_data('parameters/subvolume_min',para%subvolume_min)
   call hdf5_write_data('parameters/subvolume_max',para%subvolume_max)
   call hdf5_write_data('parameters/h',para%h,'Normalization of Hubble parameter H0 = h * 100 km/s/Mpc')
   call hdf5_write_data('parameters/omega_l',para%omega_l)
   call hdf5_write_data('parameters/omega_m',para%omega_m)
   call hdf5_write_data('parameters/omega_b',para%omega_b)
   call hdf5_write_data('parameters/dc_min',para%dc_min)
   call hdf5_write_data('parameters/dc_max',para%dc_max)
   call hdf5_write_data('parameters/ra_min',para%ra_min/degree)
   call hdf5_write_data('parameters/ra_max',para%ra_max/degree)
   call hdf5_write_data('parameters/dec_min',para%dec_min/degree)
   call hdf5_write_data('parameters/dec_max',para%dec_max/degree)
   call hdf5_write_data('parameters/zaxis_ra',para%zaxis_ra/degree)
   call hdf5_write_data('parameters/zaxis_dec',para%zaxis_dec/degree)
   call hdf5_write_data('parameters/xy_angle',para%xy_angle/degree)
   call hdf5_write_data('parameters/seed',para%seed)
   call hdf5_write_data('parameters/translate',para%translate)
   call hdf5_write_data('parameters/rotate',para%rotate)
   call hdf5_write_data('parameters/invert',para%invert)
   call hdf5_write_data('parameters/velocity_ra',para%velocity_ra)
   call hdf5_write_data('parameters/velocity_dec',para%velocity_dec)
   call hdf5_write_data('parameters/velocity_norm',para%velocity_norm)
   call hdf5_write_data('parameters/search_angle',para%search_angle)
   call hdf5_write_data('parameters/volume_search_level',para%volume_search_level)
   call hdf5_write_data('parameters/skyrotation',para%sky_rotation, &
   & 'Rotation matrix to map xyz-coordinates of the tiling structure onto sky coordinates.')
   
   ! Group "galaxies"
   name = 'galaxies'
   filename_bin = trim(para%path_output)//'mocksky_galaxies.bin'
   open(1,file=trim(filename_bin),action='read',form='unformatted',access='stream')
   read(1) n
   allocate(sky_galaxy(n))
   read(1) sky_galaxy
   close(1)
   call hdf5_add_group(trim(name))
   call hdf5_write_data(trim(name)//'/snapshot',sky_galaxy%snapshot,'snapshot ID')
   call hdf5_write_data(trim(name)//'/subvolume',sky_galaxy%subvolume,'subvolume index')
   call hdf5_write_data(trim(name)//'/tile',sky_galaxy%tile,'tile ID in tiling array')
   call hdf5_write_data(trim(name)//'/zobs',sky_galaxy%zobs,'redshift in observer-frame')
   call hdf5_write_data(trim(name)//'/zcmb',sky_galaxy%zcmb,'redshift in CMB frame')
   call hdf5_write_data(trim(name)//'/zcos',sky_galaxy%zcos,'cosmological redshift without peculiar motions')
   call hdf5_write_data(trim(name)//'/dc',sky_galaxy%dc,'[Mpc/h] comoving distance')
   call hdf5_write_data(trim(name)//'/ra',sky_galaxy%ra/degree,'[deg] right ascension')
   call hdf5_write_data(trim(name)//'/dec',sky_galaxy%dec/degree,'[deg] declination')
   call hdf5_write_data(trim(name)//'/id_galaxy',sky_galaxy%id_galaxy,'unique galaxy ID in mock sky')
   call hdf5_write_data(trim(name)//'/id_group',sky_galaxy%id_group,'unique group id if galaxy is in a group, -1 otherwise')
   call hdf5_write_data(trim(name)//'/id_galaxy_sam',sky_galaxy%id_galaxy_sam,'galaxy ID in SAM')
   call hdf5_write_data(trim(name)//'/id_halo_sam',sky_galaxy%id_halo_sam,'host halo ID in SAM')
   call hdf5_write_data(trim(name)//'/type',sky_galaxy%typ,'galaxy type (0=central, 1=satellite in halo, 2=orphan)')
   call hdf5_write_data(trim(name)//'/inclination',sky_galaxy%inclination/degree, &
   & '[deg] inclination = angle between line-of-sight and spin axis')
   call hdf5_write_data(trim(name)//'/pa',sky_galaxy%pa/degree,'[deg] position angle from north to east')
   call hdf5_write_data(trim(name)//'/mag',sky_galaxy%mag, &
   & 'apparent magnitude (generic: M/L ratio of 1, no k-correction)')
   call hdf5_write_data(trim(name)//'/s_hi',sky_galaxy%SHI,'[W/m^2] integrated HI line flux')
   call hdf5_write_data(trim(name)//'/vpecrad',sky_galaxy%vpecrad,'[proper km/s] radial peculiar velocity')
   call hdf5_write_data(trim(name)//'/mstars',sky_galaxy%mstars,'[Msun/h] stellar mass')
   call hdf5_write_data(trim(name)//'/mvir_hosthalo',sky_galaxy%mstars,'[Msun/h] host halo mass')
   call hdf5_write_data(trim(name)//'/mvir_subhalo',sky_galaxy%mstars,'[Msun/h] subhalo mass')
   call hdf5_write_data(trim(name)//'/rstar_disk',sky_galaxy%rstar_disk,'[arcsec] half-mass radius of stellar disk')
   call hdf5_write_data(trim(name)//'/rstar_bulge',sky_galaxy%rstar_bulge,'[arcsec] half-mass radius of stellar bulge')
   call hdf5_write_data(trim(name)//'/rgas_disk',sky_galaxy%rgas_disk,'[arcsec] half-mass radius of gas disk')
   test(1) = n
   test(3) = sum(sky_galaxy%tile)
   test(4) = sum(sky_galaxy%inclination)
   test(5) = sum(sky_galaxy%zobs)
   deallocate(sky_galaxy)
   
   ! Group "groups"
   name = 'groups'
   filename_bin = trim(para%path_output)//'mocksky_groups.bin'
   open(1,file=trim(filename_bin),action='read',form='unformatted',access='stream')
   read(1) n
   allocate(sky_group(n))
   read(1) sky_group
   close(1)
   call hdf5_add_group(trim(name))
   call hdf5_write_data(trim(name)//'/snapshot',sky_group%snapshot,'snapshot ID')
   call hdf5_write_data(trim(name)//'/subvolume',sky_group%subvolume,'subvolume index')
   call hdf5_write_data(trim(name)//'/tile',sky_group%tile,'tile ID in tiling array')
   call hdf5_write_data(trim(name)//'/zobs',sky_group%zobs,'redshift in observer-frame')
   call hdf5_write_data(trim(name)//'/zcmb',sky_group%zcmb,'redshift in CMB frame')
   call hdf5_write_data(trim(name)//'/zcos',sky_group%zcos,'cosmological redshift without peculiar motions')
   call hdf5_write_data(trim(name)//'/dc',sky_group%dc,'[Mpc/h] comoving distance')
   call hdf5_write_data(trim(name)//'/ra',sky_group%ra/degree,'[deg] right ascension')
   call hdf5_write_data(trim(name)//'/dec',sky_group%dec/degree,'[deg] declination')
   call hdf5_write_data(trim(name)//'/id_group',sky_group%id_group,'unique parent halo ID in mock sky')
   call hdf5_write_data(trim(name)//'/id_halo_sam',sky_group%id_halo_sam,'parent halo ID in SAM')
   call hdf5_write_data(trim(name)//'/mvir',sky_group%mvir,'[Msun/h] virial mass')
   call hdf5_write_data(trim(name)//'/n_galaxies_total',sky_group%group_ntot, &
   & 'total number of galaxies that live in the same group (host halo)')
   call hdf5_write_data(trim(name)//'/n_galaxies_selected',sky_group%group_nsel, &
   & 'number of galaxies that live in the same group (host halo) and are present in the mock survey')
   call hdf5_write_data(trim(name)//'/flag',sky_group%group_flag, &
   & 'group flag (0 if group complete, >0 if truncated by survey edge (+1), snapshot limit (+2), tile edge (+4))')
   test(6) = sum(sky_group%tile)
   test(7) = sum(sky_group%zcmb)
   test(2) = sum(sky_group%group_flag)+sum(sky_group%group_flag**2)
   deallocate(sky_group)
   
   ! Group "Tiling"
   call hdf5_add_group('tiling')
   call hdf5_write_data('tiling/tile_id',(/(i,i=1,size(tile),1)/),'unique ID of cubic tile')
   call hdf5_write_data('tiling/center_x',tile%ix(1),'x-coordinate of box-centre in units of box side length')
   call hdf5_write_data('tiling/center_y',tile%ix(2),'y-coordinate of box-centre in units of box side length')
   call hdf5_write_data('tiling/center_z',tile%ix(3),'z-coordinate of box-centre in units of box side length')
   call hdf5_write_data('tiling/dc_min',tile%dmin,'minimum comoving distance in units of box side length')
   call hdf5_write_data('tiling/dc_max',tile%dmax,'maximum comoving distance in units of box side length')
   call hdf5_write_data('tiling/rotation',tile%rotation, & 
   & 'index [1,...,6] of rotation, where 1 is the identity (negative if inversion)')
   call hdf5_write_data('tiling/translation_x',tile%translation(1), & 
   & 'x-component of translation vector in units of box side length')
   call hdf5_write_data('tiling/translation_y',tile%translation(2), & 
   & 'y-component of translation vector in units of box side length')
   call hdf5_write_data('tiling/translation_z',tile%translation(3), & 
   & 'z-component of translation vector in units of box side length')
   
   ! Group "runInfo"
   call hdf5_add_group('run_info')
   call hdf5_write_data('run_info/survey_name',para%name)
   call hdf5_write_data('run_info/stingray_version',version,'Version of Stingray used to produce this mock survey')
   call hdf5_write_data('run_info/shark_version',trim(shark_version),'Version of Shark used to produce this mock survey')
   call hdf5_write_data('run_info/shark_git_revision',trim(shark_git_revision),'Git revision of shark used to produce this data')
   call hdf5_write_data('run_info/shark_timestamp',trim(shark_timestamp),'Time at which this shark execution started')
   call hdf5_write_data('run_info/n_replica_mean',n_replica_mean,'Mean number a galaxy was replicated in the intrinsic cone '&
   &//'(before selection by apparent properties)')
   call hdf5_write_data('run_info/n_replica_max',n_replica_max,'Maximum number a galaxy was replicated in the intrinsic cone '&
   &//'(before selection by apparent properties)')
   
   ! Group "Snapshots"
   call hdf5_add_group('snapshots')
   call hdf5_write_data('snapshots/id',(/(i,i=para%snapshot_min,para%snapshot_max,1)/), &
   & 'snapshot number')
   call hdf5_write_data('snapshots/z',snapshot%redshift, &
   & 'redshift corresponding to the cosmic time of this snapshot')
   call hdf5_write_data('snapshots/dc_min',snapshot%dmin*para%L, &
   & '[Mpc/h] minimal comoving distance at which this snapshot is used')
   call hdf5_write_data('snapshots/dc_max',snapshot%dmax*para%L, &
   & '[Mpc/h] maximal comoving distance at which this snapshot is used')
   call hdf5_write_data('snapshots/n_replication',snapshot%n_replication, &
   & 'Number of tiles this snapshot has been considered for, irrespective of whether a galaxy was selected')
   
   ! close HDF5 file
   call hdf5_close()
   
   ! check test sum
   if (show_test_sum) then
      expected_sum = 17520.368820190
      if (abs(sum(test)/expected_sum-1)>1e-5) then
         write(str,'(F16.9)') sum(test)
         call out('Test sum FAILED, expected_sum = '//trim(str))
      else
         call out('Test sum OK.')
      end if
   end if
   
   call toc

end subroutine make_hdf5
   
end subroutine custom_routines

end module module_user