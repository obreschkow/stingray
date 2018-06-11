! Notes:
! Download from hyades
! rsync -av --exclude 'star_formation_histories.hdf5' dobreschkow@hyades.icrar.org:/mnt/su3ctm/clagos/SHArk_Out/medi-SURFS/SHArk/19* ~/Data/SURFS/stingray/shark/input/
! rsync -av --include '*/0/g*' --exclude '*/*/*' dobreschkow@hyades.icrar.org:/mnt/su3ctm/clagos/SHArk_Out/medi-SURFS/SHArk/ ~/Data/SURFS/stingray/shark/input/

module module_user

! ==============================================================================================================
! IMPORT MODULES
! ==============================================================================================================

! default modules, do not edit
use module_constants
use module_types
use module_linalg
use module_system
use module_cosmology
use module_conversion

! custom modules
use module_hdf5

public

! ==============================================================================================================
! SET DEFAULT PARAMETER FILENAME
! ==============================================================================================================

! Here, specify the default parameter filename, used when stingray is called without the optional argument
! "parameterfile"

character(len=255),parameter  :: parameter_filename_default = &
   & '/Users/do/Data/SURFS/stingray/shark/parameters_shark_devils.txt'


! ==============================================================================================================
! VARIABLE TYPES
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
   real*4      :: matom_disk     ! [Msun/h] atomic gas mass disk
   real*4      :: matom_bulge    ! [Msun/h] atomic gas mass bulge
   real*4      :: rstar_disk     ! [cMpc/h] half-mass of stars in the disk
   real*4      :: rstar_bulge    ! [cMpc/h] half-mass of stars in the bulge
   real*4      :: rgas_disk      ! [cMpc/h] half-mass of gas in the disk
   real*4      :: rgas_bulge     ! [cMpc/h] half-mass of gas in the bulge
   real*4      :: mvir_hosthalo  ! [Msun/h]
   real*4      :: mvir_subhalo   ! [Msun/h]
   
contains

   procedure   :: getPosition
   procedure   :: getGroupID

end type type_sam

! Here, specify the class or classes of mock object properties in the sky, such as apparent galaxy properties.
! Importantly, these classes must all be sub-classes of type_empty_sky to inherit some core functionality.
! All classes must have the three class-procedures convertSAM, writeToFile and selected.
! See the existing routines below for clarifications.

type :: type_sky

   integer*8   :: id_galaxy_sky  ! unique ID in the mock sky
   integer*8   :: id_galaxy_sam  ! galaxy ID in the SAM
   integer*8   :: id_halo_sky    ! unique parent halo ID in the mock sky
   integer*8   :: id_halo_sam    ! galaxy parent halo ID in the SAM
   integer*4   :: snapshot       ! snapshot ID
   integer*4   :: subsnapshot    ! subsnapshot ID
   integer*4   :: tile           ! tile ID
   real*4      :: z              ! apparent redshift
   real*4      :: dc             ! [simulation units = Mpc/h] comoving distance
   real*4      :: ra             ! [rad] right ascension
   real*4      :: dec            ! [rad] declination
   
contains

   procedure   :: convertSAM
   procedure   :: writeToFile
   procedure   :: selected

end type type_sky

type,extends(type_sky) :: type_sky_galaxy

   integer*4   :: typ            ! galaxy type (0=central, 1=satellite, 2=orphan)
   real*4      :: inclination    ! [rad] inclination
   real*4      :: pa             ! [rad] position angle from North to East
   real*4      :: mag            ! apparent magnitude (generic: M/L ratio of 1, no k-correction)
   real*8      :: SHI            ! [W/m^2] integrated HI line flux
   real*4      :: vrad           ! [proper km/s] radial peculiar velocity
   real*4      :: mstars         ! [Msun/h] total stellar mass
   real*4      :: mvir           ! [Msun/h] mass of 1st generation halo (i.e. direct host of type 0 galaxies)
   
contains
   
   procedure   :: selected => sky_galaxy_selected
   
end type type_sky_galaxy

type,extends(type_sky) :: type_sky_group
   
   real*4      :: mvir
   
contains
   
   procedure   :: selected => sky_group_selected
   
end type type_sky_group

type,extends(type_sky) :: type_sky_lens

   real*4      :: mvir
   real*4      :: mdisk
   real*4      :: mbulge
   real*4      :: chalo

contains
   
   procedure   :: selected => sky_lens_selected

end type type_sky_lens

type(type_sky_galaxy),target   :: sky_galaxy
type(type_sky_group),target    :: sky_group
type(type_sky_lens),target     :: sky_lens
   
type sky_pointer
   class(type_sky),pointer :: sky
end type sky_pointer

type(sky_pointer),allocatable    :: ptr(:)

contains

subroutine set_class_pointers
   implicit none
   allocate(ptr(3))
   ptr(1)%sky => sky_galaxy
   ptr(2)%sky => sky_group
   ptr(3)%sky => sky_lens
end subroutine set_class_pointers

subroutine convertSAM(sky,sam,nobj,nobjtot,dc,ra,dec,tile)
   class(type_sky)            :: sky
   type(type_sam),intent(in)  :: sam
   integer*8,intent(in)       :: nobj,nobjtot
   real*4,intent(in)          :: dc,ra,dec
   integer*4,intent(in)       :: tile
end subroutine convertSAM

subroutine writeToFile(sky,fileID)
   class(type_sky)      :: sky
   integer*4,intent(in) :: fileID
   select type (sky)
   type is (type_sky_galaxy); write(fileID) sky
   type is (type_sky_group);  write(fileID) sky
   type is (type_sky_lens);   write(fileID) sky
   class default
   call error('Unknown class for I/O.')
   end select
end subroutine writeToFile

logical function selected(sky,sam)
   class(type_sky)            :: sky
   type(type_sam),intent(in)  :: sam
   selected = .true.
end function selected

logical function sky_galaxy_selected(sky,sam)
   class(type_sky_galaxy)     :: sky
   type(type_sam),intent(in)  :: sam
   sky_galaxy_selected = .true.
end function sky_galaxy_selected

logical function sky_group_selected(sky,sam)
   class(type_sky_group)     :: sky
   type(type_sam),intent(in)  :: sam
   sky_group_selected = .true.
end function sky_group_selected

logical function sky_lens_selected(sky,sam)
   class(type_sky_lens)     :: sky
   type(type_sam),intent(in)  :: sam
   sky_lens_selected = .true.
end function sky_lens_selected

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


! ==============================================================================================================
! SELECTION OF GALAXIES IN THE MOCK OBSERVING sky FILE
! ==============================================================================================================

! selection function acting on comoving position, applied when making the tiling and intrinsic sky

logical function pos_selection(dc,ra,dec)
   
   implicit none
   real*4,intent(in) :: dc    ! [simulation units] comoving distance
   real*4,intent(in) :: ra    ! [deg] right ascension
   real*4,intent(in) :: dec   ! [deg] right ascension
   
   if (.false.) then; write(*) dc,ra,dec; end if ! dummy statement to avoid compiler warnings for unused arguments
   
   select case (trim(para%name))
   case ('DEVILS')
      pos_selection = ((ra>= 34.000).and.(ra<= 37.050).and.(dec>= -5.200).and.(dec<= -4.200)).or. &
                    & ((ra>= 52.263).and.(ra<= 53.963).and.(dec>=-28.500).and.(dec<=-27.500)).or. &
                    & ((ra>=149.380).and.(ra<=150.700).and.(dec>= +1.650).and.(dec<= +2.790))
   case default
      pos_selection = .false.
   end select

end function pos_selection

! selection function applied to the sa properties (typically intrinsic properties in the sam)

logical function sam_selection(sam)

   implicit none
   type(type_sam),intent(in)   :: sam
   
   if (.false.) then; write(*) sam; end if ! dummy statement to avoid compiler warnings for unused arguments
   
   select case (trim(para%name))
   case ('DEVILS')
      sam_selection = (sam%mstars_disk>1e8).or.((sam%mvir_hosthalo>1e12).and.(sam%typ==0))
   case default
      sam_selection = .true.
   end select
   
end function sam_selection

! selection function applied to the sky properties (typically apparent properties in the mock sky)

!integer*4 function sky_selection(sky)
!
!   implicit none
!   type(type_sky),intent(in) :: sky
!   
!   if (.false.) then; write(*) sky; end if ! dummy statement to avoid compiler warnings for unused arguments
!   
!   select case (trim(para%name))
!   case ('DEVILS')
!      sky_selection = 0
!      if (sky%mag<=21.2) sky_selection = sky_selection+1 ! galaxy selection
!      if ((sky%mvir>1e12).and.(sky%typ==0)) sky_selection = sky_selection+2 ! halo selection
!   case default
!      sky_selection = 1
!   end select
!
!end function sky_selection


! ==============================================================================================================
! IO routines
! ==============================================================================================================

! Set parameters automatically (e.g. from information provided with the snapshot files).
! These automatic parameters are only adopted if the values in the parameter files are set to 'auto'.
! Otherwise the parameter file overwrites these values.

subroutine make_automatic_parameters
   
   implicit none
   
   character(len=255)   :: filename
   integer*4            :: isnapshot,isubsnapshot
   
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
   
   isubsnapshot = -1
   do while (.not.exists(filename,.true.))
      isubsnapshot = isubsnapshot+1
      write(filename,'(A,I0,A,I0,A)') trim(para%path_input),para%snapshot_max,'/',isubsnapshot,'/galaxies.hdf5'
   end do
   para%subsnapshot_min = isubsnapshot
   do while (exists(filename,.true.))
      isubsnapshot = isubsnapshot+1
      write(filename,'(A,I0,A,I0,A)') trim(para%path_input),para%snapshot_max,'/',isubsnapshot,'/galaxies.hdf5'
   end do
   para%subsnapshot_max = isubsnapshot-1
   
   write(filename,'(A,I0,A)') trim(para%path_input),para%snapshot_min,'/0/galaxies.hdf5'
   call hdf5_open(filename)
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
   call hdf5_read_data(g//'matom_disk',sam%matom_disk)
   call hdf5_read_data(g//'matom_bulge',sam%matom_bulge)
   call hdf5_read_data(g//'rstar_disk',sam%rstar_disk)
   call hdf5_read_data(g//'rstar_bulge',sam%rstar_bulge)
   call hdf5_read_data(g//'rgas_disk',sam%rgas_disk)
   call hdf5_read_data(g//'rgas_bulge',sam%rgas_bulge)
   call hdf5_read_data(g//'mvir_hosthalo',sam%mvir_hosthalo)
   call hdf5_read_data(g//'mvir_subhalo',sam%mvir_subhalo)
   
   ! assign other properties
   sam%snapshot = index
   sam%subsnapshot = subindex
   
   ! close file
   call hdf5_close()
   
   ! return snapshot name for screen output
   write(snapshotname,'(A,I0,A,I0,A,I0,A)') 'snapshot ',index,', subindex ',subindex,' (',n,' galaxies)'
   
end subroutine load_sam_snapshot


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
   integer*8                           :: n,i
   
   ! load auxilary data
   call load_parameters
   call load_box_list
   call load_snapshot_list
   
   ! create HDF5 file
   filename_hdf5 = trim(para%path_output)//'mocksky.hdf5'
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
   filename_bin = trim(para%path_output)//'mocksky_class01.bin'
   open(1,file=trim(filename_bin),action='read',form='unformatted',access='stream')
   read(1) n
   allocate(sky_galaxy(n))
   read(1) sky_galaxy
   call hdf5_add_group('Galaxies')
   call hdf5_write_data('Galaxies/id_galaxy_sky',sky_galaxy%id_galaxy_sky,'unique galaxy ID in mock sky')
   call hdf5_write_data('Galaxies/id_galaxy_sam',sky_galaxy%id_galaxy_sam,'galaxy ID in SAM')
   call hdf5_write_data('Galaxies/id_halo_sky',sky_galaxy%id_galaxy_sky,'unique parent halo ID in mock sky')
   call hdf5_write_data('Galaxies/id_halo_sam',sky_galaxy%id_galaxy_sam,'parent halo ID in SAM')
   call hdf5_write_data('Galaxies/snapshot',sky_galaxy%snapshot,'snapshot ID')
   call hdf5_write_data('Galaxies/subsnapshot',sky_galaxy%subsnapshot,'subsnapshot ID')
   call hdf5_write_data('Galaxies/box',sky_galaxy%tile,'tile ID in tiling array')
   call hdf5_write_data('Galaxies/type',sky_galaxy%typ,'galaxy type (0=central, 1=satellite in halo, 2=orphan)')
   call hdf5_write_data('Galaxies/z',sky_galaxy%z, &
   & 'apparent redshift (Hubble flow + peculiar motion of galaxy and observer)')
   call hdf5_write_data('Galaxies/dc',sky_galaxy%dc,'[Mpc/h] comoving distance')
   call hdf5_write_data('Galaxies/RA',sky_galaxy%ra/degree,'[deg] right ascension')
   call hdf5_write_data('Galaxies/Dec',sky_galaxy%dec/degree,'[deg] declination')
   call hdf5_write_data('Galaxies/inclination',sky_galaxy%inclination/degree, &
   & '[deg] inclination = angle between line-of-sight and spin axis')
   call hdf5_write_data('Galaxies/pa',sky_galaxy%pa/degree,'[deg] position angle from north to east')
   call hdf5_write_data('Galaxies/mag',sky_galaxy%mag, &
   & 'apparent magnitude (generic: M/L ratio of 1, no k-correction)')
   call hdf5_write_data('Galaxies/SHI',sky_galaxy%SHI,'[W/m^2] integrated HI line flux')
   call hdf5_write_data('Galaxies/vrad',sky_galaxy%vrad,'[proper km/s] radial peculiar velocity')
   call hdf5_write_data('Galaxies/mstars',sky_galaxy%mstars,'[Msun/h] stellar mass')
   deallocate(sky_galaxy)
   close(1)
   
   ! Group "Halos"
   filename_bin = trim(para%path_output)//'mocksky_class02.bin'
   open(1,file=trim(filename_bin),action='read',form='unformatted',access='stream')
   read(1) n
   allocate(sky_group(n))
   read(1) sky_group
   call hdf5_add_group('Groups')
   call hdf5_write_data('Groups/id_halo_sky',sky_group%id_galaxy_sky,'unique parent halo ID in mock sky')
   call hdf5_write_data('Groups/id_halo_sam',sky_group%id_galaxy_sam,'parent halo ID in SAM')
   call hdf5_write_data('Groups/dc',sky_group%dc,'[Mpc/h] comoving distance')
   call hdf5_write_data('Groups/RA',sky_group%ra/degree,'[deg] right ascension')
   call hdf5_write_data('Groups/Dec',sky_group%dec/degree,'[deg] declination')
   call hdf5_write_data('Groups/mvir',sky_group%mvir,'[Msun/h] virial mass')
   deallocate(sky_group)
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