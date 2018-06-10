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

type type_sam_galaxy

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

end type type_sam_galaxy

type type_sam_halo
   
   integer*8   :: id_halo        ! unique ID of parent halo
   real*4      :: position(3)    ! [Mpc/h] position of galaxy centre in simulation box
   real*4      :: spinparameter  ! [-] Peebles spin parameter
   real*4      :: rvir
   
end type type_sam_halo

type type_sam

   type(type_sam_galaxy)   :: galaxy
   type(type_sam_halo)     :: halo     ! parent halo properties, e.g. for group finders and lensing
   
contains

   procedure   :: getPosition
   procedure   :: getGroupID

end type type_sam

! Here, specify the galaxy properties in the mock sky. These are mostly apparent galaxy properties.

type type_sky

   integer*8   :: id_galaxy_sky  ! unique ID in the mock sky
   integer*8   :: id_galaxy_sam  ! galaxy ID in the SAM
   integer*8   :: id_halo_sky    ! unique parent halo ID in the mock sky
   integer*8   :: id_halo_sam    ! galaxy parent halo ID in the SAM
   integer*4   :: snapshot       ! snapshot ID
   integer*4   :: subsnapshot    ! subsnapshot ID
   integer*4   :: tile           ! tile ID
   integer*4   :: typ            ! galaxy type (0=central, 1=satellite, 2=orphan)
   real*4      :: z              ! apparent redshift
   real*4      :: dc             ! [simulation units = Mpc/h] comoving distance
   real*4      :: ra             ! [rad] right ascension
   real*4      :: dec            ! [rad] declination
   real*4      :: inclination    ! [rad] inclination
   real*4      :: pa             ! [rad] position angle from North to East
   real*4      :: mag            ! apparent magnitude (generic: M/L ratio of 1, no k-correction)
   real*8      :: SHI            ! [W/m^2] integrated HI line flux
   real*4      :: vrad           ! [proper km/s] radial peculiar velocity
   real*4      :: mstars         ! [Msun/h] total stellar mass
   
end type type_sky

contains

! In order to place the objects in the mock sky, the class type_sam must have the following two functions
! enabling stingray to extract the position and group id of each object.

function getPosition(sam) result(position)
   class(type_sam) :: sam
   real*4 :: position(3)
   position = sam%galaxy%position ! unique group identifier
end function getPosition

integer*8 function getGroupID(sam)
   class(type_sam) :: sam
   getGroupID = sam%galaxy%id_halo ! [length_unit of simulation] position if the galaxy in the box
end function getGroupID


! ==============================================================================================================
! SELECTION OF GALAXIES IN THE MOCK OBSERVING sky FILE
! ==============================================================================================================

! selection function acting on comoving position, applied when making the tiling and intrinsic sky

integer*4 function pos_selection(dc,ra,dec)
   
   implicit none
   real*4,intent(in) :: dc    ! [simulation units] comoving distance
   real*4,intent(in) :: ra    ! [deg] right ascension
   real*4,intent(in) :: dec   ! [deg] right ascension
   
   if (.false.) then; write(*) dc,ra,dec; end if ! dummy statement to avoid compiler warnings for unused arguments
   
   select case (trim(para%name))
   case ('DEVILS')
      if ((ra>= 34.000).and.(ra<= 37.050).and.(dec>= -5.200).and.(dec<= -4.200)) then
         pos_selection = 1
      else if ((ra>= 52.263).and.(ra<= 53.963).and.(dec>=-28.500).and.(dec<=-27.500)) then
         pos_selection = 2
      else if ((ra>=149.380).and.(ra<=150.700).and.(dec>= +1.650).and.(dec<= +2.790)) then
         pos_selection = 3
      else
         pos_selection = 0
      end if
   case default
      pos_selection = 1
   end select

end function pos_selection

! selection function applied to the sa properties (typically intrinsic properties in the sam)

integer*4 function sam_selection(sam)

   implicit none
   type(type_sam),intent(in)   :: sam
   
   if (.false.) then; write(*) sam; end if ! dummy statement to avoid compiler warnings for unused arguments
   
   select case (trim(para%name))
   case ('DEVILS')
      if (sam%galaxy%mstars_disk>1e8) then
         sam_selection = 1
      else
         sam_selection = 0
      end if
   case default
      sam_selection = 1
   end select
   
end function sam_selection

! selection function applied to the sky properties (typically apparent properties in the mock sky)

integer*4 function sky_selection(sky)

   implicit none
   type(type_sky),intent(in) :: sky
   
   if (.false.) then; write(*) sky; end if ! dummy statement to avoid compiler warnings for unused arguments
   
   select case (trim(para%name))
   case ('DEVILS')
      if (sky%mag<=21.2) then
         sky_selection = 1
      else
         sky_selection = 0
      end if
   case default
      sky_selection = 1
   end select

end function sky_selection


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
   call hdf5_read_data(g//'id_galaxy',sam%galaxy%id_galaxy)
   call hdf5_read_data(g//'id_halo',sam%galaxy%id_halo)
   call hdf5_read_data(g//'type',sam%galaxy%typ)
   call hdf5_read_data(g//'position_x',sam%galaxy%position(1))
   call hdf5_read_data(g//'position_y',sam%galaxy%position(2))
   call hdf5_read_data(g//'position_z',sam%galaxy%position(3))
   call hdf5_read_data(g//'velocity_x',sam%galaxy%velocity(1))
   call hdf5_read_data(g//'velocity_y',sam%galaxy%velocity(2))
   call hdf5_read_data(g//'velocity_z',sam%galaxy%velocity(3))
   call hdf5_read_data(g//'L_x',sam%galaxy%J(1))
   call hdf5_read_data(g//'L_y',sam%galaxy%J(2))
   call hdf5_read_data(g//'L_z',sam%galaxy%J(3))
   call hdf5_read_data(g//'mstars_disk',sam%galaxy%mstars_disk)
   call hdf5_read_data(g//'mstars_bulge',sam%galaxy%mstars_bulge)
   call hdf5_read_data(g//'matom_disk',sam%galaxy%matom_disk)
   call hdf5_read_data(g//'matom_bulge',sam%galaxy%matom_bulge)
   call hdf5_read_data(g//'rstar_disk',sam%galaxy%rstar_disk)
   call hdf5_read_data(g//'rstar_bulge',sam%galaxy%rstar_bulge)
   call hdf5_read_data(g//'rgas_disk',sam%galaxy%rgas_disk)
   call hdf5_read_data(g//'rgas_bulge',sam%galaxy%rgas_bulge)
   
   ! assign other properties
   sam%galaxy%snapshot = index
   sam%galaxy%subsnapshot = subindex
   
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
   ! sam%galaxy%vector = rotate(sam%galaxy%vector,ps)
   ! where ps is a logical argument that specifies whether the vector transforms like a normal vector,
   ! or like a pseudo-vector (also known as axial vector)

   implicit none
   type(type_sam),intent(inout)    :: sam
   
   sam%galaxy%velocity   = rotate(sam%galaxy%velocity,.false.)
   sam%galaxy%J          = rotate(sam%galaxy%J,.true.)
   
end subroutine rotate_vectors

function convert_properties(sam,id,dc,ra,dec,tile) result(sky)

   ! This is the central function of the user module. It makes the apparent properties of the galaxies
   ! based on the intrinsic properties and the basic positional properties stored in base.
   !
   ! NB: All vectors in sam have been rotated via rotate_vectors, before this function is called
   !     except for the galaxy position. The correct position (rotated+translated) is provided in spherical
   !     sky-coordinates: base%dc [comoving simulation length unit] base%ra [rad], base%dec [rad]. Only use this
   !     position, not the position vector of the  SAM to compute apparent properties.

   implicit none
   
   type(type_sam),intent(in)              :: sam      ! intrinsic object properties from SAM
   integer*8,intent(in)                   :: id       ! unique ID of object in mock survey
   real*4,intent(in)                      :: dc       ! [simulation length units] comoving distance
   real*4,intent(in)                      :: ra,dec   ! [rad] sky position
   integer*4,intent(in)                   :: tile     ! unique tile ID
   type(type_sky)                  :: sky      ! apparent galaxy properties
   
   real*4                                 :: pos(3)   ! [simulation length units] position vector of galaxy
   real*4                                 :: dl       ! [simulation length units] luminosity distance to observer
   real*4                                 :: elos(3)  ! unit vector pointing from the observer to the object in comoving space
   real*4                                 :: mHI
   
   if (.false.) then; write(*) sam,id,dc,ra,dec,tile; end if ! dummy statement to avoid compiler warnings
   
   ! position vector
   call sph2car(dc,ra,dec,pos)
   elos = pos/norm(pos)
   
   ! sky coordinates
   sky%dc  = dc         ! [Mpc/h]
   sky%ra  = ra         ! [rad]
   sky%dec = dec        ! [rad]
   
   ! make redshift, provided the galaxy position [simulation units] galaxy-velocity [km/s]
   call make_redshift(pos*(para%length_unit/Mpc),sam%galaxy%velocity,z=sky%z)
   
   ! make inclination and position angle [rad]
   call make_inclination_and_pa(pos,sam%galaxy%J,inclination=sky%inclination,pa=sky%pa)
   
   ! make IDs
   sky%id_galaxy_sky   = id
   sky%id_galaxy_sam   = sam%galaxy%id_galaxy
   sky%id_halo_sky     = sam%galaxy%id_halo+tile*int(1e10,8)
   sky%id_halo_sam     = sam%galaxy%id_halo
   sky%tile            = tile
   sky%snapshot        = sam%galaxy%snapshot
   sky%subsnapshot     = sam%galaxy%subsnapshot
   
   ! copy basic constants
   sky%typ             = sam%galaxy%typ
   
   ! convert intrinsic to apparent properties
   dl = sky%dc*(1+sky%z) ! [Mpc/h]
   sky%mstars = sam%galaxy%mstars_disk+sam%galaxy%mstars_bulge
   sky%mag = convert_absmag2appmag(convert_stellarmass2absmag(sky%mstars/para%h,1.0),dl/para%h)
   sky%vrad = sum(sam%galaxy%velocity*elos)
   mHI = (sam%galaxy%matom_disk+sam%galaxy%matom_bulge)/1.35/para%h ! [Msun] HI mass
   sky%SHI = convert_luminosity2flux(real(mHI,8)*real(L2MHI,8)*Lsun,dl/para%h)
   
end function convert_properties


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
   character(len=255)         :: filename
   type(type_sky),allocatable :: sky(:)
   integer*8                  :: n,i
   
   ! load auxilary data
   call load_parameters
   call load_box_list
   call load_snapshot_list
   
   ! load data
   filename = trim(para%path_output)//'mocksky.bin'
   open(1,file=trim(filename),action='read',form='unformatted',access='stream')
   read(1) n
   allocate(sky(n))
   read(1) sky
   close(1)
   
   ! create HDF5 file
   write(filename,'(A,A)') trim(para%path_output),'mocksky.hdf5'
   call hdf5_create(filename)
   
   ! open HDF5 file
   call hdf5_open(filename)
   
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
   call hdf5_add_group('Galaxies')
   call hdf5_write_data('Galaxies/id_galaxy_sky',sky%id_galaxy_sky,'unique galaxy ID in mock sky')
   call hdf5_write_data('Galaxies/id_galaxy_sam',sky%id_galaxy_sam,'galaxy ID in SAM')
   call hdf5_write_data('Galaxies/id_halo_sky',sky%id_galaxy_sky,'unique parent halo ID in mock sky')
   call hdf5_write_data('Galaxies/id_halo_sam',sky%id_galaxy_sam,'parent halo ID in SAM')
   call hdf5_write_data('Galaxies/snapshot',sky%snapshot,'snapshot ID')
   call hdf5_write_data('Galaxies/subsnapshot',sky%subsnapshot,'subsnapshot ID')
   call hdf5_write_data('Galaxies/box',sky%tile,'tile ID in tiling array')
   call hdf5_write_data('Galaxies/type',sky%typ,'galaxy type (0=central, 1=satellite in halo, 2=orphan)')
   call hdf5_write_data('Galaxies/z',sky%z, &
   & 'apparent redshift (Hubble flow + peculiar motion of galaxy and observer)')
   call hdf5_write_data('Galaxies/dc',sky%dc,'[Mpc/h] comoving distance')
   call hdf5_write_data('Galaxies/RA',sky%ra/degree,'[deg] right ascension')
   call hdf5_write_data('Galaxies/Dec',sky%dec/degree,'[deg] declination')
   call hdf5_write_data('Galaxies/inclination',sky%inclination/degree, &
   & '[deg] inclination = angle between line-of-sight and spin axis')
   call hdf5_write_data('Galaxies/pa',sky%pa/degree,'[deg] position angle from north to east')
   call hdf5_write_data('Galaxies/mag',sky%mag, &
   & 'apparent magnitude (generic: M/L ratio of 1, no k-correction)')
   call hdf5_write_data('Galaxies/SHI',sky%SHI,'[W/m^2] integrated HI line flux')
   call hdf5_write_data('Galaxies/vrad',sky%vrad,'[proper km/s] radial peculiar velocity')
   call hdf5_write_data('Galaxies/mstars',sky%mstars,'[Msun/h] stellar mass')
   
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
   
   ! finalize
   deallocate(sky)

end subroutine make_hdf5
   
end subroutine handle_custom_arguments

end module module_user