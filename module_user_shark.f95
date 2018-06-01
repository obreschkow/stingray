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
use module_hdf5_utilities


! ==============================================================================================================
! SET DEFAULT PARAMETER FILENAME
! ==============================================================================================================

character(len=255),parameter  :: parameter_filename_default = 'parameters_shark.txt'


! ==============================================================================================================
! VARIABLE TYPES
! ==============================================================================================================

! specify the galaxy properties output by the semi-analytic model
! these are mostly intrinsic galaxy properties
! they must be stored as fortran binary file exactly in this order

type type_galaxy_sam

   integer*8   :: id_galaxy      ! unique galaxy ID
   integer*8   :: id_halo        ! unique ID of parent halo
   integer*4   :: typ            ! galaxy type (0=central, 1=satellite in halo, 2=orphan)
   real*4      :: position(3)    ! [Mpc/h] position of galaxy centre in simulation box
   real*4      :: velocity(3)    ! [proper km/s] peculiar velocity
   real*4      :: J(3)           ! [proper Msun/h Mpc/h km/s ?] angular momentum
   real*4      :: mstars_disk    ! [Msun] stellar mass disk
   real*4      :: mstars_bulge   ! [Msun] stellar mass bulge
   real*4      :: matom_disk     ! [Msun] atomic gas mass disk
   real*4      :: matom_bulge    ! [Msun] atomic gas mass bulge
   real*4      :: rdisk_star     ! [cMpc/h] half-mass of stars in the disk
   real*4      :: rbulge_star    ! [cMpc/h] half-mass of stars in the bulge
   real*4      :: rdisk_gas      ! [cMpc/h] half-mass of gas in the disk
   real*4      :: rbulge_gas     ! [cMpc/h] half-mass of gas in the bulge

end type type_galaxy_sam

! specify the galaxy properties in the mock cone
! these are mostly apparent galaxy properties

type type_galaxy_cone

   integer*8   :: id             ! unique galaxy ID
   integer*4   :: snapshot       ! snapshot
   integer*4   :: subsnapshot    ! subsnapshot
   integer*4   :: typ            ! galaxy type (0=central, 1=satellite in halo, 2=orphan)
   real*4      :: z              ! apparent redshift
   real*4      :: dc             ! [simulation units = Mpc/h] comoving distance
   real*4      :: ra             ! [rad] right ascension
   real*4      :: dec            ! [rad] declination
   real*4      :: inclination    ! [rad] inclination
   real*4      :: pa             ! [rad] position angle from North to East
   real*4      :: mag            ! apparent magnitude (generic band)
   real*8      :: SHI            ! [W/m^2] integrated HI line flux
   real*4      :: vrad           ! [proper km/s] radial peculiar velocity

end type type_galaxy_cone

contains

! mapping from type_galaxy_sam onto three basic properties of type_galaxy_base,
! needed internally to place the galaxy in the cone

function extract_base(sam) result(base)
   
   implicit none
   type(type_galaxy_sam),intent(in) :: sam
   type(type_galaxy_base)           :: base
   
   base%groupid   = sam%id_halo     ! unique group or parent halo identifier
   base%xbox      = sam%position    ! [length_unit of simulation] position if the galaxy in the box

end function extract_base


! ==============================================================================================================
! SELECTION OF GALAXIES IN THE MOCK OBSERVING CONE FILE
! ==============================================================================================================

! selection function acting on intrinsic properties and position in the cone, applied when producing the intrinsic cone

logical function intrinsic_selection(sam)

   implicit none
   type(type_galaxy_sam),intent(in)   :: sam
   
   if (.false.) then; write(*) sam; end if ! dummy statement to avoid compiler warnings for unused arguments
   
   intrinsic_selection = (sam%mstars_disk>1e9)
   
end function intrinsic_selection

! selection function acting on apparent properties, applied when producing the apparent cone

logical function apparent_selection(cone)

   implicit none
   type(type_galaxy_cone),intent(in)  :: cone
   
   if (.false.) then; write(*) cone; end if ! dummy statement to avoid compiler warnings for unused arguments
   
   apparent_selection = .true.
   
   apparent_selection = (cone%ra>=34.1*degree)
   apparent_selection = apparent_selection .and. (cone%ra<=37.05*degree)
   apparent_selection = apparent_selection .and. (cone%dec>=-5.2*degree)
   apparent_selection = apparent_selection .and. (cone%dec<=-4.2*degree)

end function apparent_selection


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
   call hdf5_read_dataset('/Cosmology/h',para%h)
   call hdf5_read_dataset('/Cosmology/OmegaL',para%OmegaL)
   call hdf5_read_dataset('/Cosmology/OmegaM',para%OmegaM)
   call hdf5_read_dataset('/runInfo/lbox',para%L)
   para%length_unit = Mpc/para%h
   call hdf5_close()
   
end subroutine make_automatic_parameters

! write mock-cone galaxy into binary file

subroutine write_galaxy(cone)

   ! choose which variables of the structure 'cone' to save
   type(type_galaxy_cone),intent(in) :: cone
   write(1) cone
   
end subroutine write_galaxy

! load redshifts
! this routine must allocate the array snapshot and fill in its real*4-valued property 'redshift'

subroutine load_redshifts(snapshot)

   implicit none
   type(type_snapshot),intent(inout),allocatable   :: snapshot(:)
   character(len=255)                              :: filename
   integer*4                                       :: isnapshot
   real*8                                          :: z

   allocate(snapshot(para%snapshot_min:para%snapshot_max))
   
   do isnapshot = para%snapshot_min,para%snapshot_max
      write(filename,'(A,I0,A)') trim(para%path_input),isnapshot,'/0/galaxies.hdf5'
      call hdf5_open(filename) ! NB: this routine also checks if the file exists
      call hdf5_read_dataset('/runInfo/redshift',z)
      snapshot(isnapshot)%redshift = real(z,4)
      call hdf5_close()
   end do
   
end subroutine load_redshifts

! load SAM snapshot file

subroutine load_sam_snapshot(index,subindex,sam,snapshotname)

   ! variable declaration
   implicit none
   integer*4,intent(in)                            :: index             ! snapshot index
   integer*4,intent(in)                            :: subindex          ! subindex, if the snapshot is split into several files
   type(type_galaxy_sam),allocatable,intent(out)   :: sam(:)            ! derived type of all the relevant SAM properties
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
   call hdf5_read_dataset(g//'id_galaxy',sam%id_galaxy)
   call hdf5_read_dataset(g//'id_halo',sam%id_halo)
   call hdf5_read_dataset(g//'type',sam%typ)
   call hdf5_read_dataset(g//'position_x',sam%position(1))
   call hdf5_read_dataset(g//'position_y',sam%position(2))
   call hdf5_read_dataset(g//'position_z',sam%position(3))
   call hdf5_read_dataset(g//'velocity_x',sam%velocity(1))
   call hdf5_read_dataset(g//'velocity_y',sam%velocity(2))
   call hdf5_read_dataset(g//'velocity_z',sam%velocity(3))
   call hdf5_read_dataset(g//'L_x',sam%J(1))
   call hdf5_read_dataset(g//'L_y',sam%J(2))
   call hdf5_read_dataset(g//'L_z',sam%J(3))
   call hdf5_read_dataset(g//'mstars_disk',sam%mstars_disk)
   call hdf5_read_dataset(g//'mstars_bulge',sam%mstars_bulge)
   call hdf5_read_dataset(g//'matom_disk',sam%matom_disk)
   call hdf5_read_dataset(g//'matom_bulge',sam%matom_bulge)
   call hdf5_read_dataset(g//'rdisk_star',sam%rdisk_star)
   call hdf5_read_dataset(g//'rbulge_star',sam%rbulge_star)
   call hdf5_read_dataset(g//'rdisk_gas',sam%rdisk_gas)
   call hdf5_read_dataset(g//'rbulge_gas',sam%rbulge_gas)
   
   ! close file
   call hdf5_close()
   
   ! return snapshot name for screen output
   write(snapshotname,'(A,I0,A,I0,A,I0,A)') 'snapshot ',index,', subindex ',subindex,' (',n,' galaxies)'
   
end subroutine load_sam_snapshot


! ==============================================================================================================
! CONVERSION BETWEEN INTRINSIC AND APPARENT GALAXY PROPERTIES
! ==============================================================================================================

subroutine rotate_vectors(sam)

   ! This function rotates all the vector-properties of the galaxy, except for the position, which has already
   ! been processed when producing the intrinsic cone.

   implicit none
   type(type_galaxy_sam),intent(inout)    :: sam
   
   sam%velocity   = rotate(sam%velocity,.false.)
   sam%J          = rotate(sam%J,.true.)
   
end subroutine rotate_vectors

function convert_properties(base,sam) result(cone)

   ! NB: All vectors in sam have been rotated via rotate_vectors, before this function is called
   !     except for the galaxy position. The correct position (rotated+translated) is provided in the vector
   !     base%xcone in units of box side-lengths. Only use this position vector, not the position vector of the
   !     SAM to compute sky coordinates etc.

   implicit none
   
   type(type_galaxy_base),intent(in)      :: base     ! base properties of galaxy in the cone
   type(type_galaxy_sam),intent(in)       :: sam      ! intrinsic galaxy properties from SAM
   type(type_galaxy_cone)                 :: cone     ! apparent galaxy properties
   real*4                                 :: dl       ! [simulation length units] luminosity distance to observer
   real*4                                 :: elos(3)  ! unit vector pointing from the observer to the object in comoving space
   
   ! make sky coordinates
   cone%dc  = base%dc
   cone%ra  = base%ra
   cone%dec = base%dec
   
   ! make distances and redshift, provided the position x and galaxy-velocity in km/s
   call make_redshift(base%xcone,sam%velocity,z=cone%z)
   
   ! make inclination and position angle
   call make_inclination_and_pa(base%xcone,sam%J,inclination=cone%inclination,pa=cone%pa)
   
   ! copy basic constants
   cone%snapshot     = base%snapshot
   cone%subsnapshot  = base%subsnapshot
   cone%id           = sam%id_galaxy
   cone%typ          = sam%typ
   
   ! convert intrinsic to apparent properties
   dl = cone%dc*(1+cone%z)
   cone%mag    = convert_absmag2appmag(convert_stellarmass2absmag(sam%mstars_disk+sam%mstars_bulge,1.0),dl)
   cone%vrad   = sum(sam%velocity*elos)
   cone%SHI    = convert_luminosity2flux(real(sam%matom_disk+sam%matom_bulge,8)*real(LMratioHI,8)*Lsun,dl)
   
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
   case ('my.task')
      call out('Here, specify what to do as "my.task"')
      if (len(custom_option)>0) call out('Using the custom option: '//custom_option)
   case default
      success = .false.
   end select
   
end subroutine handle_custom_arguments

end module module_user