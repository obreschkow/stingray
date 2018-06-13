! ==============================================================================================================
! PERSONAL NOTES
! ==============================================================================================================

! Download/upload data from/to hyades
! rsync -av --exclude 'star_formation_histories.hdf5' dobreschkow@hyades.icrar.org:/mnt/su3ctm/clagos/SHArk_Out/medi-SURFS/SHArk/19* ~/Data/SURFS/stingray/shark/input/
! rsync -av --include '*/0/g*' --exclude '*/*/*' dobreschkow@hyades.icrar.org:/mnt/su3ctm/clagos/SHArk_Out/medi-SURFS/SHArk/ ~/Data/SURFS/stingray/shark/input/
! scp -r ~/Data/SURFS/stingray/shark/outputDEVILS/mocksky_DEVILS.hdf5 dobreschkow@hyades.icrar.org:/mnt/su3ctm/dobreschkow/Stingray/

module module_user

! ==============================================================================================================
! IMPORT MODULES
! ==============================================================================================================

! default modules, do not edit
use module_constants
use module_system
use module_types
use module_linalg
use module_cosmology
use module_conversion

! custom modules
use module_hdf5

public
private :: parameter_filename_default ! must be accessed via get_parameter_filename_default()

! ==============================================================================================================
! SET DEFAULT PARAMETER FILENAME
! ==============================================================================================================

! Here, specify the default parameter filename, used when stingray is called without the optional argument
! "parameterfile"

character(len=255),parameter  :: parameter_filename_default = &
   & '/Users/do/Data/SURFS/stingray/shark/parameters_shark_test.txt'


! ==============================================================================================================
! TYPE DECLARATIONS
! ==============================================================================================================

! Here, specify the galaxy properties output by the semi-analytic model (SAM), which should be loaded
! to make the mock sky. How these galaxies will be loaded will be specified in the subroutine
! load_sam_snapshot below.

type type_sam

   integer*8   :: id_galaxy      ! unique galaxy ID
   integer*8   :: id_halo        ! unique ID of parent halo
   integer*4   :: snapshot       ! snapshot ID
   integer*4   :: subsnapshot    ! subsnapshot ID
   integer*4   :: typ            ! galaxy type (0=central, 1=satellite in halo, 2=orphan)
   real*4      :: position(3)    ! [Mpc/h] position of galaxy centre in simulation box
   real*4      :: velocity(3)    ! [proper km/s] peculiar velocity
   real*4      :: J(3)           ! [proper Msun/h Mpc/h km/s ?] angular momentum
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

   procedure   :: getPosition    ! required function
   procedure   :: getGroupID     ! required function

end type type_sam

! Here, specify the class or classes of mock object properties in the sky, such as apparent galaxy properties.
! Importantly, these classes must all be sub-classes of type_empty_sky to inherit some core functionality.
! All classes must have the three class-procedures convertSam, writeToFile and selected.
! See the existing routines below for clarifications.
   
type type_sky

   integer*4   :: snapshot       ! snapshot ID
   integer*4   :: subsnapshot    ! subsnapshot ID
   integer*4   :: tile           ! tile ID
   real*4      :: zobs           ! redshift in observer-frame
   real*4      :: zcmb           ! redshift in CMB frame
   real*4      :: zcos           ! cosmological redshift without peculiar motions
   real*4      :: dc             ! [simulation units = Mpc/h] comoving distance
   real*4      :: ra             ! [rad] right ascension
   real*4      :: dec            ! [rad] declination
   
   contains

   procedure   :: convertSam     ! required subroutine
   procedure   :: writeToFile    ! required subroutine
   procedure   :: selected       ! required logical function
   procedure   :: name           ! required character function, used for filenames

end type type_sky

type,extends(type_sky) :: type_sky_galaxy

   integer*8   :: id_halo_sky    ! unique parent halo ID in the mock sky
   integer*8   :: id_halo_sam    ! galaxy parent halo ID in the SAM
   integer*8   :: id_galaxy_sky  ! unique ID in the mock sky
   integer*8   :: id_galaxy_sam  ! galaxy ID in the SAM
   integer*4   :: typ            ! galaxy type (0=central, 1=satellite, 2=orphan)
   integer*4   :: group_ntot     ! total number of galaxies in group
   integer*4   :: group_nsel     ! number of selected galaxies in group
   integer*4   :: group_flag     ! 0=group complete, >0 group truncated
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

type,extends(type_sky) :: type_sky_group
   
   integer*8   :: id_halo_sky    ! unique parent halo ID in the mock sky
   integer*8   :: id_halo_sam    ! galaxy parent halo ID in the SAM
   real*4      :: mvir
   
end type type_sky_group

type,extends(type_sky) :: type_sky_lens

   integer*8   :: id_galaxy_sky  ! unique ID in the mock sky
   real*4      :: mhalo
   real*4      :: mdisk
   real*4      :: mbulge
   real*4      :: chalo

end type type_sky_lens

! Now link up the types of mock properties to the pointers 'skyclass()'
! It is important that each type is is given exactly one statement
! allocate([type name] :: skyclass(1)%ptr)
! in the subroutine set_class_pointers

type skyclass_pointer
   class(type_sky),allocatable :: ptr
end type skyclass_pointer

type(skyclass_pointer),allocatable  :: skyclass(:)

contains

subroutine set_class_pointers
   implicit none
   allocate(skyclass(3))
   allocate(type_sky_galaxy :: skyclass(1)%ptr)
   allocate(type_sky_group :: skyclass(2)%ptr)
   allocate(type_sky_lens :: skyclass(3)%ptr)
end subroutine set_class_pointers

! In order to place the objects in the mock sky, the class type_sam must have the following two functions
! enabling stingray to extract the position and group id of each object.

function getPosition(sam) result(position)
   class(type_sam) :: sam
   real*4 :: position(3)
   position = sam%position ! unique group identifier
end function getPosition

integer*8 function getGroupID(sam)
   class(type_sam) :: sam
   getGroupID = sam%id_halo ! [length_unit of simulation] position if the galaxy in the box
end function getGroupID

! For instances of the user-defined sky-types to be saved, add all sky types in the following subroutine

subroutine writeToFile(self,fileID)
   class(type_sky) :: self
   integer*4,intent(in) :: fileID
   select type (self)
   type is (type_sky_galaxy); write(fileID) self
   type is (type_sky_group); write(fileID) self
   type is (type_sky_lens); write(fileID) self
   class default
   call error('Unknown class for I/O.')
   end select
end subroutine writeToFile

! Give a name without spaces and special characters to each type

character(255) function name(self)
   class(type_sky) :: self
   select type (self)
   type is (type_sky_galaxy); name = 'Galaxies'
   type is (type_sky_group); name = 'Groups'
   type is (type_sky_lens); name = 'Lenses'
   class default
   call error('Unknown class for name.')
   end select
end function name


! ==============================================================================================================
! SELECTION OF GALAXIES IN THE MOCK OBSERVING sky FILE
! ==============================================================================================================

! selection function acting on comoving position, applied when making the tiling and intrinsic sky

logical function pos_selection(dc,ra,dec)
   
   implicit none
   real*4,intent(in) :: dc    ! [simulation units] comoving distance
   real*4,intent(in) :: ra    ! [deg] right ascension
   real*4,intent(in) :: dec   ! [deg] right ascension
   
   call nil(dc,ra,dec) ! dummy to avoid compiler warnings for unused arguments
   
   select case (trim(para%name))
   case ('DEVILS')
      pos_selection = ((ra>= 34.000).and.(ra<= 37.050).and.(dec>= -5.200).and.(dec<= -4.200)).or. &
                    & ((ra>= 52.263).and.(ra<= 53.963).and.(dec>=-28.500).and.(dec<=-27.500)).or. &
                    & ((ra>=149.380).and.(ra<=150.700).and.(dec>= +1.650).and.(dec<= +2.790))
   case ('GAMA')
      pos_selection = ((ra>= 30.200).and.(ra<= 38.800).and.(dec>=-10.250).and.(dec<= -3.720)).or. &
                    & ((ra> 129.000).and.(ra<=141.000).and.(dec>= -2.000).and.(dec<= +3.000)).or. &
                    & ((ra>=174.000).and.(ra<=186.000).and.(dec>= -2.000).and.(dec<= +3.000)).or. &
                    & ((ra>=211.500).and.(ra<=223.500).and.(dec>= -2.000).and.(dec<= +3.000)).or. &
                    & ((ra>=339.000).and.(ra<=351.000).and.(dec>=-35.000).and.(dec<=-30.000))
   case default
      pos_selection = .false.
   end select

end function pos_selection

! pre-selection function applied to the sam properties when placing the galaxies. More selections can be applied
! later in the function 'selected'. The function 'sam_selection' is simply to avoid carrying too many galaxies
! through the whole analysis

logical function sam_selection(sam)

   implicit none
   type(type_sam),intent(in)   :: sam
   
   call nil(sam)  ! dummy to avoid compiler warnings for unused arguments
   
   select case (trim(para%name))
   case ('DEVILS')
      sam_selection = (sam%mstars_disk>1e8).or.((sam%mvir_hosthalo>1e11).and.(sam%typ==0))
   case ('GAMA')
      sam_selection = (sam%mstars_disk>1e8).or.((sam%mvir_hosthalo>1e11).and.(sam%typ==0))
   case default
      sam_selection = .true.
   end select
   
end function sam_selection

! selection applied to sky classes

logical function selected(sky,sam)

   class(type_sky)            :: sky ! self
   type(type_sam),intent(in)  :: sam
   real*4,parameter           :: dmag = 0.0
   
   call nil(sky,sam) ! dummy to avoid compiler warnings for unused arguments
   
   select type (sky)
   type is (type_sky_galaxy)
   
      select case (trim(para%name))
      case ('DEVILS')
         selected = sky%mag<=21.2+dmag
      case ('GAMA')
         selected = ((sky%mag<=19.8+dmag).and.(sky%ra<330.0*degree)).or.(sky%mag<=19.2+dmag)
      case default
         selected = .true.
      end select
      
   type is (type_sky_group)
   
      selected = (sam%mvir_hosthalo>1e12).and.(sam%typ==0)
      
   type is (type_sky_lens)
   
      selected = sam%mvir_subhalo>1e13
   
   class default
      call error('Unknown class for selection.')
   end select
   
end function selected


! ==============================================================================================================
! CONVERSION BETWEEN INTRINSIC AND APPARENT GALAXY PROPERTIES
! ==============================================================================================================

subroutine rotate_vectors(sam)

   ! This function rotates all the vector-properties of the galaxy, specified in type_sam, except for the
   ! position, which has already been processed when producing the intrinsic sky. For each such vector-property
   ! of the SAM, call
   ! sam%vector = rotate(sam%vector,ps)
   ! where ps is a logical argument that specifies whether the vector transforms like a normal vector,
   ! or like a pseudo-vector (also known as axial vector)

   implicit none
   type(type_sam),intent(inout)    :: sam
   
   sam%velocity   = rotate(sam%velocity,.false.)
   sam%J          = rotate(sam%J,.true.)
   
end subroutine rotate_vectors

subroutine convertSam(sky,sam,base,nobj,nobjtot)

   implicit none
   class(type_sky)            :: sky
   type(type_sam),intent(in)  :: sam
   type(type_base),intent(in) :: base
   integer*8,intent(in)       :: nobj,nobjtot
   
   real*4                     :: pos(3)   ! [simulation length units] position vector of galaxy
   real*4                     :: dl       ! [simulation length units] luminosity distance to observer
   real*4                     :: elos(3)  ! unit vector pointing from the observer to the object in comoving space
   real*4                     :: mHI
   integer*8                  :: id_halo_sky,id_galaxy_sky
   
   call nil(sky,sam,base,nobj,nobjtot) ! dummy statement to avoid compiler warnings
   
   ! make unique IDs
   id_halo_sky   = sam%id_halo+base%tile*int(1e10,8)
   id_galaxy_sky = sam%id_galaxy+base%tile*int(1e10,8)
   
   
   ! sky coordinates
   sky%dc  = base%dc*para%L   ! [Mpc/h]
   sky%ra  = base%ra          ! [rad]
   sky%dec = base%dec         ! [rad]
   
   ! position vector
   call sph2car(sky%dc,sky%ra,sky%dec,pos)
   elos = pos/norm(pos)
   
   ! make IDs
   sky%tile            = base%tile
   sky%snapshot        = sam%snapshot
   sky%subsnapshot     = sam%subsnapshot
   
   ! make redshift, provided the galaxy position [simulation units] galaxy-velocity [km/s]
   call make_redshift(pos*(para%length_unit/Mpc),sam%velocity,zobs=sky%zobs,zcmb=sky%zcmb,zcos=sky%zcos)
   
   select type(sky)
   type is (type_sky_galaxy)
   
      ! make IDs
      sky%id_galaxy_sam    = sam%id_galaxy   ! copy other properties
      sky%id_galaxy_sky    = id_galaxy_sky
      sky%id_halo_sam      = sam%id_halo
      sky%id_halo_sky      = id_halo_sky   ! unique parent halo ID in the mock sky
      
      ! copy properties
      sky%typ              = sam%typ
      sky%mstars           = sam%mstars_disk+sam%mstars_bulge
      sky%mvir_hosthalo    = sam%mvir_hosthalo
      sky%mvir_subhalo     = sam%mvir_subhalo
      
      ! group flag
      sky%group_ntot       = base%group_ntot
      sky%group_nsel       = base%group_nsel
      sky%group_flag       = base%group_flag
      
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
      sky%id_halo_sky  = id_halo_sky   ! unique parent halo ID in the mock sky
      sky%mvir         = sam%mvir_hosthalo
   
   type is (type_sky_lens)
   
      sky%id_galaxy_sky = id_galaxy_sky
      sky%mhalo         = sam%mvir_subhalo
      sky%mdisk         = sam%mstars_disk+sam%mgas_disk
      sky%mbulge        = sam%mstars_bulge+sam%mgas_bulge
      sky%chalo         = sam%cnfw_subhalo
   
   class default
   end select
   
end subroutine convertSam


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
   call hdf5_read_data('/runInfo/tot_n_subvolumes',nsub)
   para%subsnapshot_min = 0
   para%subsnapshot_max = nsub-1
   call hdf5_read_data('/Cosmology/h',para%h)
   call hdf5_read_data('/Cosmology/OmegaL',para%OmegaL)
   call hdf5_read_data('/Cosmology/OmegaM',para%OmegaM)
   call hdf5_read_data('/Cosmology/OmegaB',para%OmegaB)
   call hdf5_read_data('/runInfo/lbox',para%L)
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
      call hdf5_read_data('/runInfo/redshift',z)
      snapshot(isnapshot)%redshift = real(z,4)
      call hdf5_close()
   end do
   
end subroutine make_redshifts

! load SAM snapshot file

subroutine load_sam_snapshot(index,subindex,sam,snapshotname)

   ! variable declaration
   implicit none
   integer*4,intent(in)                            :: index             ! snapshot index
   integer*4,intent(in)                            :: subindex          ! subindex, if the snapshot is split into several files
   type(type_sam),allocatable,intent(out)          :: sam(:)            ! derived type of all the relevant SAM properties
   character(len=100),intent(out)                  :: snapshotname      ! snapshot name to be returned for user display
   character(len=255)                              :: filename
   integer*8                                       :: n
   character(*),parameter                          :: g = '/Galaxies/'  ! group name
   
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
   call hdf5_read_data(g//'L_x',sam%J(1))
   call hdf5_read_data(g//'L_y',sam%J(2))
   call hdf5_read_data(g//'L_z',sam%J(3))
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
   sam%subsnapshot = subindex
   
   ! close file
   call hdf5_close()
   
   ! return snapshot name for screen output
   write(snapshotname,'(A,I0,A,I0,A,I0,A)') 'snapshot ',index,', subindex ',subindex,' (',n,' galaxies)'
   
end subroutine load_sam_snapshot

! return default parameter file name
character(255) function get_parameter_filename_default()
   get_parameter_filename_default = trim(parameter_filename_default)
end function get_parameter_filename_default


! ==============================================================================================================
! HANDLE CUSTOM TASKS (= optional subroutines and arguments)
! ==============================================================================================================

subroutine handle_custom_arguments(task,custom_option,success)

   implicit none
   character(*),intent(in) :: task
   character(*),intent(in) :: custom_option
   logical,intent(out)     :: success
   
   ! custom task handler
   success = .true.
   select case (trim(task))
   case ('my.task') ! dymmy argument to illustrate the use
      call out('Here, specify what to do as "my.task"')
      if (len(custom_option)>0) call out('Using the custom option: '//custom_option)
   case ('make.hdf5')
      call make_hdf5
   case ('my.additions.to.make.all')
      ! here add optional commands to be executed automatically at the end of a make.all run
      call make_hdf5
   case default
      success = .false.
   end select
   
contains

subroutine make_hdf5
   
   implicit none
   character(len=255)                  :: filename_bin
   character(len=255)                  :: filename_hdf5
   type(type_sky_galaxy),allocatable   :: sky_galaxy(:)
   type(type_sky_group),allocatable    :: sky_group(:)
   type(type_sky_lens),allocatable     :: sky_lens(:)
   integer*8                           :: n,i
   character(len=255)                  :: name
   
   ! load auxilary data
   call load_parameters
   call load_box_list
   call load_snapshot_list
   
   ! create HDF5 file
   filename_hdf5 = trim(para%path_output)//'mocksky_'//trim(para%name)//'.hdf5'
   call hdf5_create(filename_hdf5)
   
   ! open HDF5 file
   call hdf5_open(filename_hdf5)
   
   ! Group "Parameters"
   call hdf5_add_group('Parameters')
   call hdf5_write_data('Parameters/path_output',para%path_output)
   call hdf5_write_data('Parameters/path_input',para%path_input)
   call hdf5_write_data('Parameters/L',para%L)
   call hdf5_write_data('Parameters/length_unit',para%length_unit)
   call hdf5_write_data('Parameters/snapshot_min',para%snapshot_min)
   call hdf5_write_data('Parameters/snapshot_max',para%snapshot_max)
   call hdf5_write_data('Parameters/subsnapshot_min',para%subsnapshot_min)
   call hdf5_write_data('Parameters/subsnapshot_max',para%subsnapshot_max)
   call hdf5_write_data('Parameters/h',para%h,'Normalization of Hubble parameter H0 = h * 100 km/s/Mpc')
   call hdf5_write_data('Parameters/OmegaL',para%OmegaL)
   call hdf5_write_data('Parameters/OmegaM',para%OmegaM)
   call hdf5_write_data('Parameters/OmegaB',para%OmegaB)
   call hdf5_write_data('Parameters/dc_min',para%dc_min)
   call hdf5_write_data('Parameters/dc_max',para%dc_max)
   call hdf5_write_data('Parameters/ra_min',para%ra_min/degree)
   call hdf5_write_data('Parameters/ra_max',para%ra_max/degree)
   call hdf5_write_data('Parameters/dec_min',para%dec_min/degree)
   call hdf5_write_data('Parameters/dec_max',para%dec_max/degree)
   call hdf5_write_data('Parameters/zaxis_ra',para%zaxis_ra/degree)
   call hdf5_write_data('Parameters/zaxis_dec',para%zaxis_dec/degree)
   call hdf5_write_data('Parameters/xy_angle',para%xy_angle/degree)
   call hdf5_write_data('Parameters/seed',para%seed)
   call hdf5_write_data('Parameters/translate',para%translate)
   call hdf5_write_data('Parameters/rotate',para%rotate)
   call hdf5_write_data('Parameters/invert',para%invert)
   call hdf5_write_data('Parameters/velocity_ra',para%velocity_ra)
   call hdf5_write_data('Parameters/velocity_dec',para%velocity_dec)
   call hdf5_write_data('Parameters/velocity_norm',para%velocity_norm)
   call hdf5_write_data('Parameters/search_angle',para%search_angle)
   call hdf5_write_data('Parameters/volume_search_level',para%volume_search_level)
   call hdf5_write_data('Parameters/skyrotation',para%sky_rotation, &
   & 'Rotation matrix to map xyz-coordinates of the tiling structure onto sky coordinates.')
   
   ! Group "Galaxies"
   allocate(sky_galaxy(1)); name = sky_galaxy(1)%name(); deallocate(sky_galaxy)
   filename_bin = trim(para%path_output)//'mocksky_'//trim(name)//'.bin'
   open(1,file=trim(filename_bin),action='read',form='unformatted',access='stream')
   read(1) n
   allocate(sky_galaxy(n))
   name = sky_galaxy(1)%name()
   read(1) sky_galaxy
   call hdf5_add_group(trim(name))
   call hdf5_write_data(trim(name)//'/snapshot',sky_galaxy%snapshot,'snapshot ID')
   call hdf5_write_data(trim(name)//'/subsnapshot',sky_galaxy%subsnapshot,'subsnapshot ID')
   call hdf5_write_data(trim(name)//'/tile',sky_galaxy%tile,'tile ID in tiling array')
   call hdf5_write_data(trim(name)//'/zobs',sky_galaxy%zobs,'redshift in observer-frame')
   call hdf5_write_data(trim(name)//'/zcmb',sky_galaxy%zcmb,'redshift in CMB frame')
   call hdf5_write_data(trim(name)//'/zcos',sky_galaxy%zcos,'cosmological redshift without peculiar motions')
   call hdf5_write_data(trim(name)//'/dc',sky_galaxy%dc,'[Mpc/h] comoving distance')
   call hdf5_write_data(trim(name)//'/RA',sky_galaxy%ra/degree,'[deg] right ascension')
   call hdf5_write_data(trim(name)//'/Dec',sky_galaxy%dec/degree,'[deg] declination')
   call hdf5_write_data(trim(name)//'/id_halo_sky',sky_galaxy%id_halo_sky,'unique parent halo ID in mock sky')
   call hdf5_write_data(trim(name)//'/id_halo_sam',sky_galaxy%id_halo_sam,'parent halo ID in SAM')
   call hdf5_write_data(trim(name)//'/id_galaxy_sky',sky_galaxy%id_galaxy_sky,'unique galaxy ID in mock sky')
   call hdf5_write_data(trim(name)//'/id_galaxy_sam',sky_galaxy%id_galaxy_sam,'galaxy ID in SAM')
   call hdf5_write_data(trim(name)//'/type',sky_galaxy%typ,'galaxy type (0=central, 1=satellite in halo, 2=orphan)')
   call hdf5_write_data(trim(name)//'/group_ntot',sky_galaxy%group_ntot, &
   & 'total number of galaxies that live in the same host halo')
   call hdf5_write_data(trim(name)//'/group_nsel',sky_galaxy%group_nsel, &
   & 'number of galaxies that live in the same host halo and are present in the mock survey')
   call hdf5_write_data(trim(name)//'/group_flag',sky_galaxy%group_flag, &
   & 'group flag (0 if group complete, >0 if truncated by survey edge (+1), snapshot limit (+2), box limit (+4))')
   call hdf5_write_data(trim(name)//'/inclination',sky_galaxy%inclination/degree, &
   & '[deg] inclination = angle between line-of-sight and spin axis')
   call hdf5_write_data(trim(name)//'/pa',sky_galaxy%pa/degree,'[deg] position angle from north to east')
   call hdf5_write_data(trim(name)//'/mag',sky_galaxy%mag, &
   & 'apparent magnitude (generic: M/L ratio of 1, no k-correction)')
   call hdf5_write_data(trim(name)//'/SHI',sky_galaxy%SHI,'[W/m^2] integrated HI line flux')
   call hdf5_write_data(trim(name)//'/vpecrad',sky_galaxy%vpecrad,'[proper km/s] radial peculiar velocity')
   call hdf5_write_data(trim(name)//'/mstars',sky_galaxy%mstars,'[Msun/h] stellar mass')
   call hdf5_write_data(trim(name)//'/mvir_hosthalo',sky_galaxy%mstars,'[Msun/h] host halo mass')
   call hdf5_write_data(trim(name)//'/mvir_subhalo',sky_galaxy%mstars,'[Msun/h] subhalo mass')
   call hdf5_write_data(trim(name)//'/rstar_disk',sky_galaxy%rstar_disk,'[arcsec] half-mass radius of stellar disk')
   call hdf5_write_data(trim(name)//'/rstar_bulge',sky_galaxy%rstar_bulge,'[arcsec] half-mass radius of stellar bulge')
   call hdf5_write_data(trim(name)//'/rgas_disk',sky_galaxy%rgas_disk,'[arcsec] half-mass radius of gas disk')
   write(*,*) 'Test sums: ',sum(sky_galaxy%group_flag),sum(sky_galaxy%group_flag**2),sum(sky_galaxy%group_flag*sky_galaxy%zobs)
   deallocate(sky_galaxy)
   close(1)
   
   ! Group "Groups"
   allocate(sky_group(1)); name = sky_group(1)%name(); deallocate(sky_group)
   filename_bin = trim(para%path_output)//'mocksky_'//trim(name)//'.bin'
   open(1,file=trim(filename_bin),action='read',form='unformatted',access='stream')
   read(1) n
   allocate(sky_group(n))
   read(1) sky_group
   call hdf5_add_group(trim(name))
   call hdf5_write_data(trim(name)//'/snapshot',sky_group%snapshot,'snapshot ID')
   call hdf5_write_data(trim(name)//'/subsnapshot',sky_group%subsnapshot,'subsnapshot ID')
   call hdf5_write_data(trim(name)//'/tile',sky_group%tile,'tile ID in tiling array')
   call hdf5_write_data(trim(name)//'/zobs',sky_group%zobs,'redshift in observer-frame')
   call hdf5_write_data(trim(name)//'/zcmb',sky_group%zcmb,'redshift in CMB frame')
   call hdf5_write_data(trim(name)//'/zcos',sky_group%zcos,'cosmological redshift without peculiar motions')
   call hdf5_write_data(trim(name)//'/dc',sky_group%dc,'[Mpc/h] comoving distance')
   call hdf5_write_data(trim(name)//'/RA',sky_group%ra/degree,'[deg] right ascension')
   call hdf5_write_data(trim(name)//'/Dec',sky_group%dec/degree,'[deg] declination')
   call hdf5_write_data(trim(name)//'/id_halo_sky',sky_group%id_halo_sky,'unique parent halo ID in mock sky')
   call hdf5_write_data(trim(name)//'/id_halo_sam',sky_group%id_halo_sam,'parent halo ID in SAM')
   call hdf5_write_data(trim(name)//'/mvir',sky_group%mvir,'[Msun/h] virial mass')
   deallocate(sky_group)
   close(1)
   
   ! Group "Lenses"
   allocate(sky_lens(1)); name = sky_lens(1)%name(); deallocate(sky_lens)
   filename_bin = trim(para%path_output)//'mocksky_'//trim(name)//'.bin'
   open(1,file=trim(filename_bin),action='read',form='unformatted',access='stream')
   read(1) n
   allocate(sky_lens(n))
   read(1) sky_lens
   call hdf5_add_group(trim(name))
   call hdf5_write_data(trim(name)//'/snapshot',sky_lens%snapshot,'snapshot ID')
   call hdf5_write_data(trim(name)//'/subsnapshot',sky_lens%subsnapshot,'subsnapshot ID')
   call hdf5_write_data(trim(name)//'/tile',sky_lens%tile,'tile ID in tiling array')
   call hdf5_write_data(trim(name)//'/zobs',sky_lens%zobs,'redshift in observer-frame')
   call hdf5_write_data(trim(name)//'/zcmb',sky_lens%zcmb,'redshift in CMB frame')
   call hdf5_write_data(trim(name)//'/zcos',sky_lens%zcos,'cosmological redshift without peculiar motions')
   call hdf5_write_data(trim(name)//'/dc',sky_lens%dc,'[Mpc/h] comoving distance')
   call hdf5_write_data(trim(name)//'/RA',sky_lens%ra/degree,'[deg] right ascension')
   call hdf5_write_data(trim(name)//'/Dec',sky_lens%dec/degree,'[deg] declination')
   call hdf5_write_data(trim(name)//'/id_galaxy_sky',sky_lens%id_galaxy_sky,'unique galaxy ID in mock sky')
   call hdf5_write_data(trim(name)//'/mhalo',sky_lens%mhalo,'[Msun/h] halo mass')
   call hdf5_write_data(trim(name)//'/mdisk',sky_lens%mdisk,'[Msun/h] disk mass')
   call hdf5_write_data(trim(name)//'/mbluge',sky_lens%mbulge,'[Msun/h] bulge mass')
   call hdf5_write_data(trim(name)//'/chalo',sky_lens%mhalo,'halo concentration (NFW fit)')
   deallocate(sky_lens)
   close(1)
   
   ! Group "Tiling"
   call hdf5_add_group('Tiling')
   call hdf5_write_data('Tiling/BoxID',(/(i,i=1,size(tile),1)/),'unique ID of cubic box')
   call hdf5_write_data('Tiling/center.x',tile%ix(1),'x-coordinate of box-centre in units of box side length')
   call hdf5_write_data('Tiling/center.y',tile%ix(2),'y-coordinate of box-centre in units of box side length')
   call hdf5_write_data('Tiling/center.z',tile%ix(3),'z-coordinate of box-centre in units of box side length')
   call hdf5_write_data('Tiling/dc_min',tile%dmin,'minimum comoving distance in units of box side length')
   call hdf5_write_data('Tiling/dc_max',tile%dmax,'maximum comoving distance in units of box side length')
   call hdf5_write_data('Tiling/rotation',tile%rotation, & 
   & 'index [1,...,6] of rotation, where 1 is the identity (negative if inversion)')
   call hdf5_write_data('Tiling/translation.x',tile%translation(1), & 
   & 'x-component of translation vector in units of box side length')
   call hdf5_write_data('Tiling/translation.y',tile%translation(2), & 
   & 'y-component of translation vector in units of box side length')
   call hdf5_write_data('Tiling/translation.z',tile%translation(3), & 
   & 'z-component of translation vector in units of box side length')
   
   ! Group "runInfo"
   call hdf5_add_group('RunInfo')
   call hdf5_write_data('RunInfo/version',version,'Version of Stingray used to produce this mock survey')
   
   ! Group "Snapshots"
   call hdf5_add_group('Snapshots')
   call hdf5_write_data('Snapshots/id',(/(i,i=para%snapshot_min,para%snapshot_max,1)/), &
   & 'snapshot number')
   call hdf5_write_data('Snapshots/z',snapshot%redshift, &
   & 'redshift corresponding to the cosmic time of this snapshot')
   call hdf5_write_data('Snapshots/dc_min',snapshot%dmin*para%L, &
   & '[Mpc/h] minimal comoving distance at which this snapshot is used')
   call hdf5_write_data('Snapshots/dc_max',snapshot%dmax*para%L, &
   & '[Mpc/h] maximal comoving distance at which this snapshot is used')
   
   ! close HDF5 file
   call hdf5_close()

end subroutine make_hdf5
   
end subroutine handle_custom_arguments

end module module_user