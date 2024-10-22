! **********************************************************************************************************************************
! This module defines the interface between the semi-analytic model (SAM) and stingray
!
! Each SAM must have its own module, named "module_user_routines_[sam].f03", where "sam" is the name of the semi-analytic model
! specified in the makefile. The development of stingray was primarily driven with the SAM "shark" in mind; hence the module
! module_user_routines_shark is particularly advanced.
! 
! Instructions of how to adapt this module to a particular SAM and/or add new intrinsic or apparent galaxy properties are given
! in the comments below. Throughout this module, the user has programmatic access to a list of global parameters,
! routines and operators. Important examples are:
! + physical constants and unit conversion factors, specified in the "shared_module_constants"
! + routines for reading/writing HDF5 data, specified in the "shared_module_hdf5"
! + vector handing routines and vector operators, specified in the "shared_module_vectors"
! + astrophysical conversion routines, specified in the "module_conversion"
! + various mathematical functions/operators, contained in the "shared_module_maths";
!   for instance, this module contains the operator .safedivide. to perform a division (/) but avoid infinities and assign 0/0=0
! + user parameters, set by user in the parameter file, accessible via the variable "para" of derived data type type_para
!   (see module_global for an overview of the variables in this type)
! **********************************************************************************************************************************


module module_user_routines

! **********************************************************************************************************************************
! MODULE INTERFACE (DO NOT EDIT)
! **********************************************************************************************************************************

! default modules, do not edit
use shared_module_core
use shared_module_constants
use shared_module_parameters
use shared_module_hdf5
use shared_module_maths
use shared_module_vectors
use shared_module_cosmology
use module_global
use module_parameters
use module_conversion
use module_emission_lines

! required objects
public :: type_sam
public :: type_sky_galaxy
public :: type_sky_group ! (can be empty)
public :: make_sky_object
public :: make_sky_galaxy
public :: make_sky_group
public :: make_auto_parameters ! (can be empty)
public :: get_redshift
public :: load_sam_snapshot
public :: write_hdf5

private


! **********************************************************************************************************************************
! DECLARATION OF CUSTOM TYPES
! **********************************************************************************************************************************

! Here, specify the class type_sam containing the SAM-properties to be loaded for making the mock sky. How these properties are
! read from files is specified in the subroutine load_sam_snapshot further below.
! The class type_sam requires three class functions:
! get_position: returns xyz-position of the galaxy in the simulation box
! get_groupid: returns the unique group id of the galaxy
! is_group_center: logical function specifying if the galaxy is the group center

integer*4,parameter  :: nco = 10 ! number of rotational CO lines

type type_sam

   ! header properties
   integer*4   :: id_galaxy         ! unique galaxy ID
   integer*8   :: id_halo           ! unique ID of parent halo
   integer*4   :: snapshot          ! snapshot ID
   integer*4   :: subvolume         ! subvolume index
   integer*4   :: typ               ! galaxy type (0=central, 1=satellite in halo, 2=orphan)
   
   ! properties of 1st generation halo (host halo)
   real*4      :: mvir_hosthalo     ! [Msun/h]
   real*4      :: vvir_hosthalo     ! [km/s]	virial velocity of hosthalo
   
   ! properties of subhalo associated with the particular galaxy
   real*4      :: mvir_subhalo      ! [Msun/h]
   real*4      :: cnfw_subhalo      ! [-] concentration of NFW fit to subhalo
   real*4      :: vvir_subhalo      ! [km/s]	virial velocity of subhalo
   real*4      :: vmax_subhalo      ! [km/s]	maximum circular velocity of subhalo
      
   ! dynamical properties   
   real*4      :: position(3)       ! [Mpc/h] position of galaxy centre in simulation box
   real*4      :: velocity(3)       ! [proper km/s] peculiar velocity
   real*4      :: J(3)              ! [proper Msun/h pMpc/h km/s] angular momentum
   real*4      :: jbulge            ! [cMpc/h km/s] specific angular momentum of the bulge
   real*4      :: jdisk             ! [cMpc/h km/s] specific angular momentum of the disk
   
   ! stellar properties
   real*4      :: mstars_disk       ! [Msun/h] stellar mass disk
   real*4      :: mstars_bulge      ! [Msun/h] stellar mass bulge
   real*4      :: rstar_disk        ! [cMpc/h] half-mass radius of stars in the disk
   real*4      :: rstar_bulge       ! [cMpc/h] half-mass radius of stars in the bulge
   real*4      :: sfr_disk          ! [Msun/h/Gyr] star formation rate in the disk
   real*4      :: sfr_burst         ! [Msun/h/Gyr] star formation rate in the bulge
   
   ! gas properties
   real*4      :: mgas_disk         ! [Msun/h] gas mass disk
   real*4      :: mgas_bulge        ! [Msun/h] gas mass bulge
   real*4      :: matom_disk        ! [Msun/h] atomic gas mass disk
   real*4      :: matom_bulge       ! [Msun/h] atomic gas mass bulge
   real*4      :: mmol_disk         ! [Msun/h] molecular gas mass disk
   real*4      :: mmol_bulge        ! [Msun/h] molecular gas mass bulge
   real*4      :: rgas_disk         ! [cMpc/h] half-mass radius of gas in the disk
   real*4      :: rgas_bulge        ! [cMpc/h] half-mass radius of gas in the bulge
   real*4      :: mgas_metals_disk  ! [Msun/h] mass of metals locked up in the disk
   real*4      :: mgas_metals_bulge ! [Msun/h] mass of metals locked up in the bulge
   
   ! black hole properties
   real*4      :: mbh               ! [Msun/h] black hole mass
   real*4      :: mbh_acc_hh        ! [Msun/h/Gyr] accretion rate in hot-halo mode
   real*4      :: mbh_acc_sb        ! [Msun/h/Gyr] accretion rate in starburst mode
   
   ! luminosities (loaded from additional file, if option "luminosities" selected)
   real*4      :: lco_disk(nco)     ! [Jy km/s Mpc^2] luminosities of the first nco transitions of CO (12/16) in the disk
   real*4      :: lco_bulge(nco)    ! [Jy km/s Mpc^2] luminosities of the first nco transitions of CO (12/16) in the bulge
   real*4      :: lum_agn_hx        ! [1e40 erg/s] hard x-ray luminosity of the AGN
   
contains

   procedure   :: get_position      => sam_get_position     ! required function of type real*4(dim=3)
   procedure   :: get_groupid       => sam_get_groupid      ! required function of type integer*8
   procedure   :: is_group_center   => sam_is_group_center  ! required function of type logical*4

end type type_sam

! Here, specify the class type_sky, which contains all the properties shared by galaxies *and* groups and the mock sky.
   
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
   real*4      :: vpec(3)        ! [proper km/s] 3D peculiar velocity
   real*4      :: vpecrad        ! [proper km/s] radial peculiar velocity
   integer*8   :: id_halo_sam    ! galaxy parent halo ID in the SAM
   integer*8   :: id_group_sky   ! unique group ID in the mock sky
   
end type type_sky

! Here, specify the class extension type_sky_galaxy, which contains the additional properties of galaxies in the mock sky,
! which are not shared with properties of groups.

type,extends(type_sky) :: type_sky_galaxy

   ! intrinsic properties copied directly from SAM
   integer*4   :: id_galaxy_sam              ! galaxy ID in the SAM
   integer*4   :: typ                        ! galaxy type (0=central, 1=satellite, 2=orphan)
   real*4      :: mstars_disk                ! [Msun/h] stellar mass of the disk
   real*4      :: mstars_bulge               ! [Msun/h] stellar mass of the bulge
   real*4      :: mgas_disk                  ! [Msun/h] gas mass of the disk
   real*4      :: mgas_bulge                 ! [Msun/h] gas mass of the bulge
   real*4      :: matom_disk                 ! [Msun/h] atomic gas mass of the disk
   real*4      :: matom_bulge                ! [Msun/h] atomic gas mass of the bulge
   real*4      :: mmol_disk                  ! [Msun/h] molecular gas mass of the disk
   real*4      :: mmol_bulge                 ! [Msun/h] molecular gas mass of the bulge
   real*4      :: jbulge                     ! [cMpc/h km/s] specific angular momentum of the bulge
   real*4      :: jdisk                      ! [cMpc/h km/s] specific angular momentum of the disk
   real*4      :: J(3)                       ! [Msun pMpc km/s] total angular momentum
   real*4      :: rstar_disk_intrinsic       ! [cMpc/h] intrinsic half-mass radius of stars in the disk
   real*4      :: rstar_bulge_intrinsic      ! [cMpc/h] intrinsic half-mass radius of stars in the bulge
   real*4      :: rgas_disk_intrinsic        ! [cMpc/h] intrinsic half-mass radius of gas in the disk
   real*4      :: rgas_bulge_intrinsic       ! [cMpc/h] intrinsic half-mass radius of gas in the bulge
   real*4      :: mvir_hosthalo              ! [Msun/h] mass of 1st generation halo (i.e. direct host of type 0 galaxies)
   real*4      :: mvir_subhalo               ! [Msun/h] subhalo mass
   real*4      :: vvir_hosthalo              ! [km/s]	virial velocity of hosthalo
   real*4      :: vvir_subhalo               ! [km/s]	virial velocity of subhalo
   real*4      :: vmax_subhalo               ! [km/s]	maximum circular velocity of subhalo
   real*4      :: cnfw_subhalo               ! [-] concentration of NFW fit to subhalo
   real*4      :: sfr_disk                   ! [Msun/Gyr/h] star formation rate in the disk
   real*4      :: sfr_burst                  ! [Msun/Gyr/h] star formation rate in the bulge
   real*4      :: mbh                        ! [Msun/h] black hole mass
   real*4      :: mbh_acc_hh                 ! [Msun/Gyr/h] accretion rate in hot-halo mode
   real*4      :: mbh_acc_sb                 ! [Msun/Gyr/h] accretion rate in starburst mode

   ! newly computed apparent properties
   integer*8   :: id_galaxy_sky              ! unique ID in the mock sky
   real*4      :: inclination                ! [rad] inclination
   real*4      :: pa                         ! [rad] position angle from North to East
   real*4      :: mag                        ! apparent magnitude (generic: M/L ratio of 1, no k-correction)
   real*4      :: rstar_disk_apparent        ! [arcsec] apparent semi-major axis of half-mass ellipse of stars in the disk
   real*4      :: rstar_bulge_apparent       ! [arcsec] apparent semi-major axis of half-mass ellipse of stars in the bulge
   real*4      :: rgas_disk_apparent         ! [arcsec] apparent semi-major axis of half-mass ellipse of gas in the disk
   real*4      :: rgas_bulge_apparent        ! [arcsec] apparent semi-major axis of half-mass ellipse of gas in the bulge
   real*4      :: zgas_disk                  ! metallicity of the gas in the disk
   real*4      :: zgas_bulge                 ! metallicity of the gas in the disk
   real*4      :: hiline_flux_int            ! [W/m^2] integrated HI line flux
   real*4      :: hiline_flux_int_vel        ! [Jy km/s] velocity-integrated HI line flux
   
   ! properties only computed if option "luminosities"
   real*4      :: lum_agn_hx                 ! [1e40 erg/s] hard x-ray luminosity of the AGN
   real*4      :: coline_flux_int(nco)       ! [W/m^2] integrated line fluxes of CO(1-0), C(2-1), ... transition
   real*4      :: coline_flux_int_vel(nco)   ! [Jy km/s] velocity-integrated line fluxes of CO(1-0), C(2-1), ... transition
   
   ! properties only computed if option "line_shapes" is selected
   type(type_line_shape) :: hiline_shape     ! shape-parameters of HI emission line (see module module_emission_lines)
   type(type_line_shape) :: h2line_shape     ! shape-parameters of molecular emission lines (see module module_emission_lines)
   
end type type_sky_galaxy

! Here, specify the class extension type_sky_group, which contains the additional properties of groups in the mock sky,
! which are not shared with properties of galaxies.

type,extends(type_sky) :: type_sky_group
   
   real*4      :: mvir                 ! [Msun/h] virial mass of group
   real*4      :: vvir                 ! [km/s]	virial velocity of group halo
   real*4      :: vmax                 ! [km/s]	maximum circular velocity of group halo
   real*4      :: cnfw                 ! [-] concentration of NFW fit to group halo
   integer*4   :: group_ntot           ! total number of galaxies in group
   integer*4   :: group_nsel           ! number of selected galaxies in group
   integer*4   :: group_flag           ! 0=group complete, >0 group truncated
   real*4      :: sigma_los_detected   ! [km/s] line-of-sight velocity dispersion of selected galaxies
   real*4      :: sigma_3D_all         ! [km/s] intrinsic 3D velocity dispersion of all galaxies in group (VR output)
   
end type type_sky_group

contains

! In order to place the objects in the mock sky, the class type_sam must have the following function
! enabling stingray to extract the position of each object.

pure function sam_get_position(self) result(position)
   class(type_sam),intent(in) :: self
   real*4 :: position(3)
   position = self%position ! [default length unit as specified in parameter file, here Mpc/h] 3D position if the galaxy in the box
end function sam_get_position

! In order to associate different galaxies with groups of galaxies, class type_sam must have the
! following two functions, enabling stingray to extract the group-id and central/satellite-status of each object

integer*8 function sam_get_groupid(self)
   class(type_sam),intent(in) :: self
   call nil(self) ! dummy statement to avoid compiler warnings if self is not used
   sam_get_groupid = self%id_halo ! unique group identifier
   ! NB: if groups are irrelevant and/or not available in the SAM,
   ! use "sam_get_groupid = 0" (as a dummy)
   ! and set "make_groups = 0" in the parameter file
end function sam_get_groupid

logical*4 function sam_is_group_center(self)
   class(type_sam),intent(in) :: self
   call nil(self) ! dummy statement to avoid compiler warnings if self is not used
   sam_is_group_center = self%typ==0
   ! NB: if groups are irrelevant and/or not available in the SAM,
   ! use "sam_is_group_center = .true." (as a dummy)
   ! and set "make_groups = 0" in the parameter file
end function sam_is_group_center


! **********************************************************************************************************************************
! CONVERSION BETWEEN INTRINSIC AND APPARENT GALAXY PROPERTIES
! **********************************************************************************************************************************

! The following polymorphic subroutine works both on galaxies and groups and makes their shared properties listed in the
! base type type_sky.

subroutine make_sky_object(sky_object,sam,base)

   ! required interface (do not edit)
   implicit none
   class(type_sky),intent(out)   :: sky_object  ! polymorphic argument works on type_sky_galaxy and type_sky_group
   type(type_sam),intent(in)     :: sam         ! SAM properties of the galaxy, as listed in the user defined type_sam;
                                                ! in the case of groups, this is the central galaxy
   type(type_base),intent(in)    :: base        ! basic properties about the placement in the mock sky:
                                 ! real*4       base%cartesian%x [user length unit] comoving cartesian x-coordinate
                                 ! real*4       base%cartesian%y [user length unit] comoving cartesian y-coordinate
                                 ! real*4       base%cartesian%z [user length unit] comoving cartesian z-coordinate
                                 ! real*4       base%spherical%dc [user length unit] comoving distance
                                 ! real*4       base%spherical%ra [rad] right ascension
                                 ! real*4       base%spherical%dec [rad] declination
                                 ! real*4       base%snapshot_redshift redshift of snapshot
                                 ! integer*4    base%index%tile tile index
                                 ! integer*4    base%index%shell shell index
                                 ! integer*4    base%index%snapshot snapshot index
                                 ! integer*4    base%index%subvolume subvolume index
                                 ! integer*8    base%index%group unique group id (-1 for isolated galaxies)
                                 ! integer*8    base%index%galaxy unique galaxy id (for groups this is the id of the
                                 !              central galaxy or -1, if the central galaxy is not selected)
                                 ! real*4(3)    base%transformation%translation [box sides] translation vector from SAM to sky
                                 ! real*4(3,3)  base%transformation%rotation rotation matrix from SAM to sky coordinates
                                 ! logical      base%transformation%inverted flag that is true if the three axes have been inverted
   ! end required interface
   
   ! internal variables (add here as needed)
   
   ! dummy statements 
   call nil(sky_object,sam,base) ! dummy statement to avoid compiler warnings (do not edit)
   
   ! checks
   if (base%index%snapshot/=sam%snapshot) call error('snapshot index mismatch')
   if (base%index%subvolume/=sam%subvolume) call error('subvolume index mismatch')

   ! indices
   sky_object%snapshot      = sam%snapshot
   sky_object%subvolume     = sam%subvolume
   sky_object%tile          = base%index%tile
   sky_object%id_halo_sam   = sam%id_halo
   sky_object%id_group_sky  = base%index%group
   
   ! sky coordinates
   sky_object%dc  = base%spherical%dc  ! [simulation length unit = Mpc/h]
   sky_object%ra  = base%spherical%ra  ! [rad]
   sky_object%dec = base%spherical%dec ! [rad]
   
   ! peculiar velocity
   sky_object%vpec      = rotate_vector(sam%velocity,base,ispseudovector=.false.)
   sky_object%vpecrad   = sky_object%vpec.dot.unitvector(base%cartesian)
   
   ! redshifts from the position and velocity of the object
   call make_redshift(components(base%cartesian*(para%length_unit/unit%Mpc)), &
   & sky_object%vpec,zobs=sky_object%zobs,zcmb=sky_object%zcmb,zcos=sky_object%zcos)
      
end subroutine make_sky_object

! The following subroutine makes all the galaxy properties specific to galaxies, i.e. not shared with groups.

subroutine make_sky_galaxy(sky_galaxy,sam,base)

   ! required interface (do not edit)
   implicit none
   class(type_sky_galaxy),intent(inout)   :: sky_galaxy  ! galaxy object with all shared properties (listed in the base type
                                                         ! type_sky) already computed
                                                         ! via make_sky_object
   type(type_sam),intent(in)              :: sam         ! SAM properties of the galaxy, as listed in the user defined type_sam
   type(type_base),intent(in)             :: base        ! basic properties (see details in subroutine make_sky_object)
   ! end required interface
   
   real*4,parameter  :: L2MHI = 6.27e-9      ! (LHI/Lsun)/(MHI/Msun)
   real*4,parameter  :: sigma_gas = 10.0     ! [km/s] standard velocity dispersion of cold gas
   real*4,parameter  :: fCO = 115.2712018    ! [GHz] rest-frame frequency of CO(1-0) transition
   real*4            :: dl                   ! [simulation length units] luminosity distance to observer
   real*4            :: mstars               ! [Msun] total stellar mass
   real*4            :: mHI                  ! [Msun/h] total HI mass
   real*4            :: LCO                  ! [Jy km/s Mpc^2] total velocity-integrated flux of CO lines
   real*4            :: wavelength           ! [m]
   integer*4         :: j
   
   call nil(sky_galaxy,sam,base) ! dummy statement to avoid compiler warnings (do not edit)
   
   ! INTRINSIC PROPERTIES (most of these properties are simply copied from the type_sam object to the type_sky_galaxy object)
   
   ! header properties
   sky_galaxy%id_galaxy_sam         = sam%id_galaxy
   sky_galaxy%typ                   = sam%typ
   
   ! properties of 1st generation halo (host halo)
   sky_galaxy%mvir_hosthalo         = sam%mvir_hosthalo
   sky_galaxy%vvir_hosthalo         = sam%vvir_hosthalo
   
   ! properties of subhalo associated with the particular galaxy
   sky_galaxy%mvir_subhalo          = sam%mvir_subhalo
   sky_galaxy%cnfw_subhalo          = sam%cnfw_subhalo
   sky_galaxy%vvir_subhalo          = sam%vvir_subhalo
   sky_galaxy%vmax_subhalo          = sam%vmax_subhalo
   
   ! intrinsic angular momentum
   sky_galaxy%J                     = rotate_vector(sam%J,base,ispseudovector=.true.)
   sky_galaxy%jbulge                = sam%jbulge
   sky_galaxy%jdisk                 = sam%jdisk
   sky_galaxy%sfr_disk              = sam%sfr_disk
   sky_galaxy%sfr_burst             = sam%sfr_burst
   
   ! stellar properties
   sky_galaxy%mstars_disk           = sam%mstars_disk
   sky_galaxy%mstars_bulge          = sam%mstars_bulge
   sky_galaxy%rstar_disk_intrinsic  = sam%rstar_disk
   sky_galaxy%rstar_bulge_intrinsic = sam%rstar_bulge
   
   ! gas properties
   sky_galaxy%mgas_disk             = sam%mgas_disk
   sky_galaxy%mgas_bulge            = sam%mgas_bulge
   sky_galaxy%matom_disk            = sam%matom_disk
   sky_galaxy%matom_bulge           = sam%matom_bulge
   sky_galaxy%mmol_disk             = sam%mmol_disk
   sky_galaxy%mmol_bulge            = sam%mmol_bulge
   sky_galaxy%rgas_disk_intrinsic   = sam%rgas_disk
   sky_galaxy%rgas_bulge_intrinsic  = sam%rgas_bulge
   
   ! derived gas metal fractions
   sky_galaxy%zgas_disk             = sam%mgas_metals_disk.safedivide.sam%mgas_disk
   sky_galaxy%zgas_bulge            = sam%mgas_metals_bulge.safedivide.sam%mgas_bulge
   
   ! black hole properties
   sky_galaxy%mbh                   = sam%mbh
   sky_galaxy%mbh_acc_sb            = sam%mbh_acc_sb
   sky_galaxy%mbh_acc_hh            = sam%mbh_acc_hh
   
      
   ! APPARENT PROPERTIES
   
   sky_galaxy%id_galaxy_sky   = base%index%galaxy ! unique galaxy id
      
   ! inclination and position angle
   call make_inclination_and_pa(components(base%cartesian),sky_galaxy%J,inclination=sky_galaxy%inclination,pa=sky_galaxy%pa)
   
   ! apparent radii (note division by comoving distance, because intrinsic radii given in comoving units)
   sky_galaxy%rstar_disk_apparent = sam%rstar_disk/sky_galaxy%dc/unit%arcsec ! [arcsec]
   sky_galaxy%rstar_bulge_apparent = sam%rstar_bulge/sky_galaxy%dc/unit%arcsec ! [arcsec]
   sky_galaxy%rgas_disk_apparent = sam%rgas_disk/sky_galaxy%dc/unit%arcsec ! [arcsec]
   sky_galaxy%rgas_bulge_apparent = sam%rgas_bulge/sky_galaxy%dc/unit%arcsec ! [arcsec]
   
   ! rough generic optical magnitude
   dl = sky_galaxy%dc*(1+sky_galaxy%zobs)/para%h ! [Mpc]
   mstars = (sam%mstars_disk+sam%mstars_bulge)/para%h ! [Msun]
   sky_galaxy%mag = convert_absmag2appmag(convert_stellarmass2absmag(mstars,1.0),dl)
      
   ! integrated 21cm flux
   mHI = (sam%matom_disk+sam%matom_bulge)/1.35 ! [Msun/h] HI mass (without Helium)
   sky_galaxy%hiline_flux_int = real(convert_luminosity2flux(real(mHI/para%h,8)*real(L2MHI,8)*unit%Lsun,dl),4) ! [W/m^2]
   sky_galaxy%hiline_flux_int_vel = convert_intflux2velintflux(sky_galaxy%hiline_flux_int,0.21106114,sky_galaxy%zobs)
   
   ! option "luminosities"
   if (option('luminosities')) then
   
      ! x-rays
      sky_galaxy%lum_agn_hx = sam%lum_agn_hx
   
      ! integrated fluxes of rotational CO transitions
      do j = 1,nco
         LCO = sam%lco_disk(j)+sam%lco_bulge(j) ! [Jy km/s Mpc^2] velocity-integrated luminosity of CO(j-[j-1]) transition
         wavelength = const%c/(fCO*1e9)/j ! [m] rest-frame wavelength of CO(j-[j-1]) transition
         sky_galaxy%coline_flux_int(j) = LCO/(4.0*pi*dl**2)*1e-23/wavelength ! [W/m^2] integrated flux
         sky_galaxy%coline_flux_int_vel(j) = LCO/(4.0*pi*dl**2)*(1+sky_galaxy%zobs) ! [Jy km/s] velocity-integrated flux
      end do
   end if
   
   ! option "line_shapes"
   if (option('line_shapes')) then
      ! shape parameters of atomic and molecular emission lines
      call make_line_profiles 
   end if
   
contains
   
   subroutine make_line_profiles
   
      implicit none
      type(type_line_input)   :: line_input
      type(type_line_shape)   :: line_shape(2)
      real*4                  :: a,c_halo,cMpch2pkpc,rvir
      real,parameter          :: G = 4.30241e-6     ! [(km/s)^2 kpc/Msun] gravitational constant
      
      ! scale factor for comoving-physical conversion
      a = 1.0/(1.0+base%snapshot_redshift) ! scale factor at redshift of snapshot
      cMpch2pkpc = 1e3*a/para%h
   
      ! set halo mass, radius and concentration
      if (sam%mvir_subhalo>0.0) then
         c_halo             = sam%cnfw_subhalo
         line_input%mhalo   = 0.1931472*sam%mvir_subhalo/((log(1+c_halo)-c_halo/(1+c_halo))) ! [Msun/h] Halo mass inside characteristic radius
         if (sam%vvir_subhalo>0.0) then
            rvir            = G*sam%mvir_subhalo/para%h/sam%vvir_subhalo**2 ! [pkpc]
         else
            rvir            = (G*sam%mvir_subhalo)**(1.0/3.0)/para%h
         end if
         line_input%rhalo   = rvir/c_halo ! [pkpc]
      else
         line_input%mhalo   = 0.0
         line_input%rhalo   = 1.0 ! just a dummy value that is irrelevant (but must be different from 0)
      end if
   
      ! set disk mass and radius
      line_input%mdisk      = sam%mstars_disk+sam%mgas_disk   ! [Msun/h] disk mass (all baryons)
      if (line_input%mdisk>0.0) then
         line_input%rdisk   = (sam%mstars_disk*sam%rstar_disk+sam%mgas_disk*sam%rgas_disk)/line_input%mdisk*cMpch2pkpc ! [pkpc] half-mass radius
      else
         line_input%rdisk   = 1.0 ! just a dummy value that is irrelevant (but must be different from 0)
      end if
   
      ! set bulge mass and radius
      line_input%mbulg      = sam%mstars_bulge+sam%mgas_bulge   ! [Msun/h] bulge mass
      if (line_input%mbulg>0.0) then
         line_input%rbulg   = (sam%mstars_bulge*sam%rstar_bulge+sam%mgas_bulge*sam%rgas_bulge)/line_input%mbulg*cMpch2pkpc ! [pkpc] half-mass radius
      else
         line_input%rbulg   = 1.0 ! just a dummy value that is irrelevant (but must be different from 0)
      end if
      
      ! correct halo mass for disk and bulge
      line_input%mhalo      = max(0.0,line_input%mhalo-line_input%mdisk-line_input%mbulg) 
   
      ! set cold gas properties
      line_input%mHI  = mHI ! [Msun/h] HI mass
      line_input%mH2  = ((sam%mgas_disk+sam%mgas_bulge)-(sam%matom_disk+sam%matom_bulge))/1.35 ! [Msun/h] H2 mass !
      line_input%rgas = (sam%mgas_disk*sam%rgas_disk+sam%mgas_bulge*sam%rgas_bulge)/(sam%mgas_disk+sam%mgas_bulge)*cMpch2pkpc ! [pkpc] half-mass radius of cold gas (CHECK)
      line_input%dispersionHI = sigma_gas ! [km/s] velocity dispersion of HI
      line_input%dispersionH2 = sigma_gas ! [km/s] velocity dispersion of H2
      
      ! set other properties
      line_input%z      = sky_galaxy%zobs    ! [-] redshift
      line_input%incl   = sky_galaxy%inclination    ! [rad] inclination (0 = face-on, pi/2 = edge-on)
      line_input%dc     = sky_galaxy%dc     ! [Mpc/h] comoving distance

      ! make line profile
      if ((line_input%mHI>0.0).or.(line_input%mH2>0.0)) call make_emission_line(line_input,line_shape)
   
      ! save line parameters
      if (line_input%mHI>0.0) then
         sky_galaxy%hiline_shape = line_shape(1)
      else
         sky_galaxy%hiline_shape = empty_line_shape
      end if
      if (line_input%mH2>0.0) then
         sky_galaxy%h2line_shape = line_shape(2)
      else
         sky_galaxy%h2line_shape = empty_line_shape
      end if
         
   end subroutine make_line_profiles
   
end subroutine make_sky_galaxy

! The following subroutine makes all the properties specific to groups, i.e. not shared with galaxies.
! Note that the input arrays sam(:), sky_galaxy(:) and selected(:) are ordered such that the central galaxy of the group
! is the first object, sam(1), sky_galaxy(1), selected(1).

subroutine make_sky_group(sky_group,sam,sky_galaxy,selected,base,groupflag)

   ! required interface (do not edit)
   implicit none
   class(type_sky_group),intent(inout) :: sky_group      ! group object with all shared properties (type_sky) already computed
                                                         ! via make_sky_object
   type(type_sam),intent(in)           :: sam(:)         ! array with the SAM-properties of *all* the galaxies in the group,
                                                         ! including those not selected
   type(type_sky_galaxy),intent(in)    :: sky_galaxy(:)  ! array with the properties in type_sky_galaxy (incl. type_sk) of *all* the
                                                         ! galaxies in the group, including those not selected; however only the
                                                         ! values of the selected galaxies, i.e. those where selected(i)==.true.,
                                                         ! have been computed. DO NOT use the uninitialised sky properties of
                                                         ! unselected galaxies.
   logical*4,intent(in)                :: selected(:)    ! logical array showing which galaxies in the group have passed the full
                                                         ! selection of the mock sky, as defined in the user module
                                                         ! module_user_selection_[...]
   type(type_base),intent(in)          :: base           ! basic properties (see details in subroutine make_sky_object)
   integer*4                           :: groupflag      ! indicates of the group is clipped: 0 if group unclipped, >0 if clipped by
                                                         ! + survey edge (+1)
                                                         ! + snapshot limit (+2)
                                                         ! + tile limit (+4)
                                                         ! + shell limit (+8)
   ! end required interface
   
   integer*4                           :: i,n,nsel
   real*4                              :: dv,v0(3),vr
   
   call nil(sky_group,sam(1),sky_galaxy(1),selected(1),base) ! dummy statement to avoid compiler warnings (do not edit)
      
   ! basic halo properties
   sky_group%mvir          = sam(1)%mvir_hosthalo
   sky_group%vvir          = sam(1)%vvir_hosthalo
   sky_group%vmax          = sam(1)%vmax_subhalo
   sky_group%cnfw          = sam(1)%cnfw_subhalo
   
   ! group properties
   sky_group%group_ntot    = size(selected)
   sky_group%group_flag    = groupflag
   sky_group%group_nsel    = count(selected)
   
   ! 3D velocity dispersion of all galaxies in group, whether selected or not
   n = size(selected)
   v0 = 0
   do i = 1,n
      v0 = v0+sam(i)%velocity
   end do
   v0 = v0/n
   dv = 0
   do i = 1,n
      dv = dv+sum((sam(i)%velocity-v0)**2)
   end do
   sky_group%sigma_3d_all = sqrt(dv/(n-1)) ! [km/s] velocity dispersion (std)
   
   ! 1D velocity dispersion along the line-of-sight of selected objects only
   nsel = count(selected)
   if (nsel>=2) then
      vr = 0
      do i = 1,n
         if (selected(i)) vr = vr+sky_galaxy(i)%vpecrad
      end do
      vr = vr/nsel ! mean peculiar velocity of observed galaxies
      dv = 0
      do i = 1,n
         if (selected(i)) dv = dv+(sky_galaxy(i)%vpecrad-vr)**2
      end do
      sky_group%sigma_los_detected = sqrt(dv/(nsel-1)) ! standard deviation of peculiar velocity of observed galaxies
   else
      sky_group%sigma_los_detected = 0
   end if
   
end subroutine make_sky_group


! **********************************************************************************************************************************
! I/O routines
! **********************************************************************************************************************************

! Function returns filename of the SAM galaxies, given a snapshot and subvolume index
! This function is private to this module. It is only used by other user-defined routines in this module.

function filename_sam(isnapshot,isubvolume,set) result(fn)
   
   implicit none
   integer*4,intent(in)       :: isnapshot,isubvolume
   integer*4,intent(in)       :: set
   character(:),allocatable   :: fn
   
   if (set==1) then
      fn = dir(para%path_input,isnapshot,isubvolume,'galaxies.hdf5')
   else if (set==2) then
      fn = dir(para%path_input,isnapshot,isubvolume,'CO_SLED.hdf5')
   else
      call deverror('unknown set')
   end if
   
end function filename_sam


! The following routine sets parameters automatically (e.g. from information provided with the snapshot files).
! These automatic parameters are adopted if and only if the values in the parameter files are set to 'auto'.
! Otherwise the values provided in the parameter file take priority.
! When the routine make_auto_parameters is called, only the parameter para%path_input is already available.
! Note: You cannot set parameter directly using para%x = ..., but only using the routine
! set_auto_parameter(name,value)

subroutine make_auto_parameters
   
   implicit none
   
   character(len=255)   :: filename
   integer*4            :: isnapshot,isubvolume
   integer*4            :: snapshot_min, snapshot_max
   real*4               :: x,h
   
   ! assume that first subvolume always exists
   call set_auto_parameter('subvolume_min',0)
   
   ! find minimal snapshot, assuming that the first subvolume exists
   isnapshot = -1
   do while (.not.exists(filename_sam(isnapshot,0,1)).and.(isnapshot<limit%n_snapshots_max))
      isnapshot = isnapshot+1
   end do
   if (isnapshot==limit%n_snapshots_max) call error('no snapshots found in '//trim(para%path_input))
   snapshot_min = isnapshot
   call set_auto_parameter('snapshot_min',snapshot_min)
   
   ! find maximal snapshot, assuming that the first subvolume exists
   do while (exists(filename_sam(isnapshot,0,1)))
      isnapshot = isnapshot+1
   end do
   snapshot_max = isnapshot-1
   call set_auto_parameter('snapshot_max',snapshot_max)
   
   ! find maximal subvolume, assuming that the last snapshot is representative
   isubvolume = 0
   do while (exists(filename_sam(snapshot_max,isubvolume,1)))
      isubvolume = isubvolume+1
   end do
   call set_auto_parameter('subvolume_max',isubvolume-1)
   
   ! other parameters
   filename = filename_sam(snapshot_min,0,1)
   call out('File of automatic parameters: '//trim(filename))
   call hdf5_open(filename)
   call hdf5_read_data('/cosmology/h',h); call set_auto_parameter('h',h)
   call hdf5_read_data('/cosmology/omega_l',x); call set_auto_parameter('omega_l',x)
   call hdf5_read_data('/cosmology/omega_m',x); call set_auto_parameter('omega_m',x)
   call hdf5_read_data('/cosmology/omega_b',x); call set_auto_parameter('omega_b',x)
   call hdf5_read_data('/run_info/lbox',x); call set_auto_parameter('box_side',x)
   call set_auto_parameter('length_unit',unit%Mpc/h) ! [m] length unit expressed in SI
   call hdf5_close()
   
end subroutine make_auto_parameters

! The following public function returns the redshift of a particular snapshot.

function get_redshift(isnapshot) result(z)

   ! required interface (do not edit)
   implicit none
   integer*4,intent(in) :: isnapshot   ! snapshot index
   real*4               :: z           ! redshift
   ! end required interface
   
   call hdf5_open(filename_sam(isnapshot,para%subvolume_min,1))
   call hdf5_read_data('/run_info/redshift',z,convert=.true.)
   call hdf5_close()
   
end function get_redshift

! The following public routine loads the galaxies of a specific snapshot and subvolume into the allocatable array sam(:) of
! type type_sam, which contains the user-defined SAM properties to be considered by stingray.

subroutine load_sam_snapshot(isnapshot,isubvolume,sam)

   ! required interface (do not edit)
   implicit none
   integer*4,intent(in)                            :: isnapshot         ! snapshot index
   integer*4,intent(in)                            :: isubvolume        ! subvolume index
   type(type_sam),allocatable,intent(out)          :: sam(:)            ! class of SAM properties
   ! end required interface
   
   integer*4                                       :: n                 ! number of galaxies
   character(*),parameter                          :: g = '/galaxies/'  ! group name
   integer*4                                       :: i,j,k,nopt
   type(type_sam),allocatable                      :: tmp(:)
   real*4,allocatable                              :: array4(:,:)
   
   ! read main galaxy properties
   
   ! open file
   call hdf5_open(filename_sam(isnapshot,isubvolume,1))
   
   ! determine number of galaxies in this (sub)snapshot
   call hdf5_get_dataset_size(g//'id_galaxy',n_elements=n)
   
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
   call hdf5_read_data(g//'vvir_subhalo',sam%vvir_subhalo)
   call hdf5_read_data(g//'vmax_subhalo',sam%vmax_subhalo)
   call hdf5_read_data(g//'vvir_hosthalo',sam%vvir_hosthalo)
   call hdf5_read_data(g//'m_bh',sam%mbh)
   call hdf5_read_data(g//'bh_accretion_rate_hh',sam%mbh_acc_hh)
   call hdf5_read_data(g//'bh_accretion_rate_sb',sam%mbh_acc_sb)
   call hdf5_read_data(g//'specific_angular_momentum_disk_star',sam%jdisk)
   call hdf5_read_data(g//'specific_angular_momentum_bulge_star',sam%jbulge)
   call hdf5_read_data(g//'sfr_disk',sam%sfr_disk)
   call hdf5_read_data(g//'sfr_burst',sam%sfr_burst)
   call hdf5_read_data(g//'mgas_metals_disk',sam%mgas_metals_disk)
   call hdf5_read_data(g//'mgas_metals_bulge',sam%mgas_metals_bulge)
   
   ! close file
   call hdf5_close()
   
   
   ! read optional galaxy properties
   
   ! default values
   do k = 1,nco
      sam%lco_bulge(k)  = 0
      sam%lco_disk(k)   = 0
   end do
   sam%lum_agn_hx = 0
   
   if (option('luminosities')) then
   
      ! open file
      call hdf5_open(filename_sam(isnapshot,isubvolume,2)) ! NB: this routine also checks if the file exists
   
      ! check number of galaxies
      call hdf5_get_dataset_size(g//'id_galaxy',n_elements=nopt)
      if (nopt>n) call error('more galaxies in luminosity file than in main galaxy file')
      
      ! read galaxy data
      allocate(tmp(nopt))
      allocate(array4(nco,nopt))
      call hdf5_read_data(g//'id_galaxy',tmp%id_galaxy,convert=.true.)
      call hdf5_read_data(g//'LCO_bulge',array4,convert=.true.)
      do k = 1,nco; tmp%lco_bulge(k) = array4(k,:); end do
      call hdf5_read_data(g//'LCO_disk',array4,convert=.true.)
      do k = 1,nco; tmp%lco_disk(k) = array4(k,:); end do
      call hdf5_read_data(g//'Lum_AGN_HardXray',tmp%lum_agn_hx,convert=.true.)
   
      ! close file
      call hdf5_close()

      ! match IDs and write tmp => sam
      j = 1
      do i = 1,nopt
         do while (sam(j)%id_galaxy.ne.tmp(i)%id_galaxy)
            j = j+1
            if (j>n) call error('galaxy IDs of luminosity data cannot be matched to main galaxy file ('// &
            & val2str(isnapshot)//','//val2str(isubvolume)//','//val2str(tmp(i)%id_galaxy)//')')
         end do
         do k = 1,nco
            sam(j)%lco_bulge(k)  = tmp(i)%lco_bulge(k)
            sam(j)%lco_disk(k)   = tmp(i)%lco_disk(k)
         end do
         sam(j)%lum_agn_hx = tmp(i)%lum_agn_hx
      end do
      
   end if
   
   ! assign other properties
   sam%snapshot = isnapshot
   sam%subvolume = isubvolume
   
end subroutine load_sam_snapshot

! The following public routine writes the arrays sky_galaxy(:) and sky_group(:) into a custom HDF5 file.

subroutine write_hdf5(filename_hdf5,sky_galaxy,sky_group)
   
   implicit none
   character(*),intent(in)                      :: filename_hdf5  ! output filename
   type(type_sky_galaxy),intent(in),allocatable :: sky_galaxy(:)
   type(type_sky_group),intent(in),allocatable  :: sky_group(:)
   character(:),allocatable                     :: name
   character(len=255)                           :: shark_version,shark_git_revision,shark_timestamp
   integer*4                                    :: i,j
   
   allocate(character(1)::name) ! empty allocation to avoid compiler flags
   
   ! load shark version information
   call hdf5_open(filename_sam(para%snapshot_max,0,1))
   call hdf5_read_data('/run_info/shark_version',shark_version)
   call hdf5_read_data('/run_info/shark_git_revision',shark_git_revision)
   call hdf5_read_data('/run_info/timestamp',shark_timestamp)
   call hdf5_close()
   
   ! open mock sky HDF5 file
   call hdf5_open(filename_hdf5,.true.)
   
   ! write shark version properties into existing group "run_info"
   name = 'run_info/'
   call hdf5_write_data(name//'shark_version',trim(shark_version),'Version of Shark used to produce this mock sky')
   call hdf5_write_data(name//'shark_git_revision',trim(shark_git_revision),'Git revision of shark used to produce this data')
   call hdf5_write_data(name//'shark_timestamp',trim(shark_timestamp),'Time at which this shark execution started')
   
   ! write group "galaxies"
   name = 'galaxies/'
   call hdf5_add_group(name)
   call hdf5_write_data(name//'snapshot',sky_galaxy%snapshot,'snapshot index')
   call hdf5_write_data(name//'subvolume',sky_galaxy%subvolume,'subvolume index')
   call hdf5_write_data(name//'tile',sky_galaxy%tile,'tile index in tiling array')
   call hdf5_write_data(name//'zobs',sky_galaxy%zobs,'redshift in observer-frame')
   call hdf5_write_data(name//'zcmb',sky_galaxy%zcmb,'redshift in CMB frame')
   call hdf5_write_data(name//'zcos',sky_galaxy%zcos,'cosmological redshift without peculiar motions')
   call hdf5_write_data(name//'dc',sky_galaxy%dc,'[Mpc/h] comoving distance')
   call hdf5_write_data(name//'ra',sky_galaxy%ra/unit%degree,'[deg] right ascension')
   call hdf5_write_data(name//'dec',sky_galaxy%dec/unit%degree,'[deg] declination')
   call hdf5_write_data(name//'id_galaxy_sky',(/(i,i=lbound(sky_galaxy,1),ubound(sky_galaxy,1))/), &
   & 'unique galaxy ID in mock sky, from 1 to n')
   call hdf5_write_data(name//'id_galaxy_sky_smart',sky_galaxy%id_galaxy_sky,'unique galaxy ID in mock sky, smart format')
   if (para%make_groups) call hdf5_write_data(name//'id_group_sky',sky_galaxy%id_group_sky, &
   & 'unique group ID if galaxy is in a group, -1 otherwise')
   call hdf5_write_data(name//'id_galaxy_sam',sky_galaxy%id_galaxy_sam,'galaxy ID in SAM')
   call hdf5_write_data(name//'id_halo_sam',sky_galaxy%id_halo_sam,'host halo ID in SAM')
   call hdf5_write_data(name//'type',sky_galaxy%typ,'galaxy type (0=central, 1=satellite in halo, 2=orphan)')
   call hdf5_write_data(name//'inclination',sky_galaxy%inclination/unit%degree, &
   & '[deg] inclination = angle between line-of-sight and spin axis')
   call hdf5_write_data(name//'pa',sky_galaxy%pa/unit%degree,'[deg] position angle from north to east')
   call hdf5_write_data(name//'mag',sky_galaxy%mag, &
   & 'apparent magnitude (generic: M/L ratio of 1, no k-correction)')
   call hdf5_write_data(name//'vpec_x',sky_galaxy%vpec(1),'[proper km/s] x-component of peculiar velocity')
   call hdf5_write_data(name//'vpec_y',sky_galaxy%vpec(2),'[proper km/s] y-component of peculiar velocity')
   call hdf5_write_data(name//'vpec_z',sky_galaxy%vpec(3),'[proper km/s] z-component of peculiar velocity')
   call hdf5_write_data(name//'vpec_r',sky_galaxy%vpecrad,'[proper km/s] line-of-sight peculiar velocity')
   call hdf5_write_data(name//'mstars_disk',sky_galaxy%mstars_disk,'[Msun/h] stellar mass in the disk')
   call hdf5_write_data(name//'mstars_bulge',sky_galaxy%mstars_bulge,'[Msun/h] stellar mass in the bulge')
   call hdf5_write_data(name//'mgas_disk',sky_galaxy%mgas_disk,'[Msun/h] gas mass in the disk')
   call hdf5_write_data(name//'mgas_bulge',sky_galaxy%mgas_bulge,'[Msun/h] gas mass in the bulge')
   call hdf5_write_data(name//'matom_disk',sky_galaxy%matom_disk,'[Msun/h] atomic gas mass in the disk')
   call hdf5_write_data(name//'matom_bulge',sky_galaxy%matom_bulge,'[Msun/h] atomic mass in the bulge')
   call hdf5_write_data(name//'mmol_disk',sky_galaxy%mmol_disk,'[Msun/h] molecular gas mass in the disk')
   call hdf5_write_data(name//'mmol_bulge',sky_galaxy%mmol_bulge,'[Msun/h] molecular gas mass in the bulge')
   call hdf5_write_data(name//'l_x',sky_galaxy%J(1),'[Msun pMpc km/s] x-component of total angular momentum')
   call hdf5_write_data(name//'l_y',sky_galaxy%J(2),'[Msun pMpc km/s] y-component of total angular momentum')
   call hdf5_write_data(name//'l_z',sky_galaxy%J(3),'[Msun pMpc km/s] z-component of total angular momentum')
   call hdf5_write_data(name//'jdisk',sky_galaxy%jdisk, '[cMpc/h km/s] specific angular momentum of the disk')
   call hdf5_write_data(name//'jbulge',sky_galaxy%jbulge, '[cMpc/h km/s] specific angular momentum of the bulge')
   call hdf5_write_data(name//'mvir_hosthalo',sky_galaxy%mvir_hosthalo,'[Msun/h] host halo mass')
   call hdf5_write_data(name//'mvir_subhalo',sky_galaxy%mvir_subhalo,'[Msun/h] subhalo mass')
   call hdf5_write_data(name//'zgas_disk',sky_galaxy%zgas_disk,'metallicity of the gas in the disk')
   call hdf5_write_data(name//'zgas_bulge',sky_galaxy%zgas_bulge,'metallicity of the gas in the bulge')
   call hdf5_write_data(name//'sfr_disk',sky_galaxy%sfr_disk,'[Msun/h/Gyr] star formation rate in the disk')
   call hdf5_write_data(name//'sfr_burst',sky_galaxy%sfr_burst,'[Msun/h/Gyr] star formation rate in the bulge')
   call hdf5_write_data(name//'mbh',sky_galaxy%mbh,'[Msun/h] black hole mass')
   call hdf5_write_data(name//'mbh_acc_hh',sky_galaxy%mbh_acc_hh,'[Msun/h/Gyr] black hole accretion rate in the hot halo mode')
   call hdf5_write_data(name//'mbh_acc_sb',sky_galaxy%mbh_acc_sb,'[Msun/h/Gyr] black hole accretion rate in the starburst mode')
   call hdf5_write_data(trim(name)//'/vvir_hosthalo',sky_galaxy%vvir_hosthalo,'[km/s] host halo virial velocity')
   call hdf5_write_data(trim(name)//'/vvir_subhalo',sky_galaxy%vvir_subhalo,'[km/s] subhalo virial velocity')
   call hdf5_write_data(trim(name)//'/cnfw_subhalo',sky_galaxy%cnfw_subhalo,'NFW concentration parameter of subhalo')
   call hdf5_write_data(name//'rstar_disk_apparent',sky_galaxy%rstar_disk_apparent,&
   &'[arcsec] apparent semi-major axis of half-mass ellipse of stellar disk')
   call hdf5_write_data(name//'rstar_bulge_apparent',sky_galaxy%rstar_bulge_apparent,&
   &'[arcsec] apparent semi-major axis of half-mass ellipse of stellar bulge')
   call hdf5_write_data(name//'rgas_disk_apparent',sky_galaxy%rgas_disk_apparent,&
   &'[arcsec] apparent semi-major axis of half-mass ellipse of gas disk')
   call hdf5_write_data(name//'rgas_bulge_apparent',sky_galaxy%rgas_bulge_apparent,&
   &'[arcsec] apparent semi-major axis of half-mass ellipse of gas bulge')
   call hdf5_write_data(name//'rstar_disk_intrinsic',sky_galaxy%rstar_disk_intrinsic,&
   &'[cMpc/h] intrinsic half-mass radius of stellar disk')
   call hdf5_write_data(name//'rstar_bulge_intrinsic',sky_galaxy%rstar_bulge_intrinsic,&
   &'[cMpc/h] intrinsic half-mass radius of stellar bulge')
   call hdf5_write_data(name//'rgas_disk_intrinsic',sky_galaxy%rgas_disk_intrinsic,&
   &'[cMpc/h] intrinsic half-mass radius of gas disk')
   call hdf5_write_data(name//'rgas_bulge_intrinsic',sky_galaxy%rgas_bulge_intrinsic,&
   &'[cMpc/h] intrinsic half-mass radius of gas bulge')
   call hdf5_write_data(name//'hiline_flux_int',sky_galaxy%hiline_flux_int,'[W/m^2] integrated HI line flux')
   call hdf5_write_data(name//'hiline_flux_int_vel',sky_galaxy%hiline_flux_int_vel, &
   & '[Jy km/s] velocity-integrated HI line flux')
   
   if (option('luminosities')) then
   
      call hdf5_write_data(name//'lum_agn_hx',sky_galaxy%lum_agn_hx,'[1e40 erg/s] hard x-ray luminosity of the AGN')
      do j = 1,nco ! CO transition
         call hdf5_write_data(name//'coline_flux_int_'//val2str(j),sky_galaxy%coline_flux_int(j),&
         &'[W/m^2] integrated CO('//val2str(j)//'-'//val2str(j-1)//') line flux')
         call hdf5_write_data(name//'coline_flux_int_vel_'//val2str(j),sky_galaxy%coline_flux_int_vel(j),&
         &'[Jy km/s] velocity-integrated CO('//val2str(j)//'-'//val2str(j-1)//') line flux')
      end do
      
   end if
   
   if (option('line_shapes')) then
   
      call hdf5_write_data(name//'hiline_flux_peak',sky_galaxy%hiline_shape%speak, &
      & '[s/km] normalised peak HI line flux density of inclined galaxy (multiply by hiline_flux_int_vel to get Jy values)')
      call hdf5_write_data(name//'hiline_flux_central',sky_galaxy%hiline_shape%scentral, &
      & '[s/km] normalised central HI line flux density of inclined galaxy (multiply by hiline_flux_int_vel to get Jy values)')
      call hdf5_write_data(name//'hiline_width_peak',sky_galaxy%hiline_shape%wpeak, &
      & '[km/s] HI line-width between flux peaks of inclined galaxy in rest-frame velocity units')
      call hdf5_write_data(name//'hiline_width_50',sky_galaxy%hiline_shape%w50, &
      & '[km/s] HI line-width at 50% of the peak flux of inclined galaxy in rest-frame velocity units')
      call hdf5_write_data(name//'hiline_width_20',sky_galaxy%hiline_shape%w20, &
      & '[km/s] HI line-width at 20% of the peak flux of inclined galaxy in rest-frame velocity units')
      
      call hdf5_write_data(name//'hiline_flux_peak_eo',sky_galaxy%hiline_shape%speak_eo, &
      & '[s/km] normalised peak HI line flux density of edge-on galaxy (multiply by hiline_flux_int_vel to get Jy values)')
      call hdf5_write_data(name//'hiline_flux_central_eo',sky_galaxy%hiline_shape%scentral_eo, &
      & '[s/km] normalised central HI line flux density of edge-on galaxy (multiply by hiline_flux_int_vel to get Jy values)')
      call hdf5_write_data(name//'hiline_width_peak_eo',sky_galaxy%hiline_shape%wpeak_eo, &
      & '[km/s] HI line-width between flux peaks of edge-on galaxy in rest-frame velocity units')
      call hdf5_write_data(name//'hiline_width_50_eo',sky_galaxy%hiline_shape%w50_eo, &
      & '[km/s] HI line-width at 50% of the peak flux of edge-on galaxy in rest-frame velocity units')
      call hdf5_write_data(name//'hiline_width_20_eo',sky_galaxy%hiline_shape%w20_eo, &
      & '[km/s] HI line-width at 20% of the peak flux of edge-on galaxy in rest-frame velocity units')
      
      call hdf5_write_data(name//'h2line_flux_peak',sky_galaxy%h2line_shape%speak, &
      & '[s/km] normalised peak molecular line flux density of inclined galaxy (multiply by velocity-integrated flux to get Jy)')
      call hdf5_write_data(name//'h2line_flux_central',sky_galaxy%h2line_shape%scentral, &
      & '[s/km] normalised central molecular line flux density of inclined galaxy (multiply by velocity-integrated flux to get Jy)')
      call hdf5_write_data(name//'h2line_width_peak',sky_galaxy%h2line_shape%wpeak, &
      & '[km/s] molecular line-width between flux peaks of inclined galaxy in rest-frame velocity units')
      call hdf5_write_data(name//'h2line_width_50',sky_galaxy%h2line_shape%w50, &
      & '[km/s] molecular line-width at 50% of the peak flux of inclined galaxy in rest-frame velocity units')
      call hdf5_write_data(name//'h2line_width_20',sky_galaxy%h2line_shape%w20, &
      & '[km/s] molecular line-width at 20% of the peak flux of inclined galaxy in rest-frame velocity units')
      
      call hdf5_write_data(name//'h2line_flux_peak_eo',sky_galaxy%h2line_shape%speak_eo, &
      & '[s/km] normalised peak molecular line flux density of edge-on galaxy (multiply by velocity-integrated flux to get Jy)')
      call hdf5_write_data(name//'h2line_flux_central_eo',sky_galaxy%h2line_shape%scentral_eo, &
      & '[s/km] normalised central molecular line flux density of edge-on galaxy (multiply by velocity-integrated flux to get Jy)')
      call hdf5_write_data(name//'h2line_width_peak_eo',sky_galaxy%h2line_shape%wpeak_eo, &
      & '[km/s] molecular line-width between flux peaks of edge-on galaxy in rest-frame velocity units')
      call hdf5_write_data(name//'h2line_width_50_eo',sky_galaxy%h2line_shape%w50_eo, &
      & '[km/s] molecular line-width at 50% of the peak flux of edge-on galaxy in rest-frame velocity units')
      call hdf5_write_data(name//'h2line_width_20_eo',sky_galaxy%h2line_shape%w20_eo, &
      & '[km/s] molecular line-width at 20% of the peak flux of edge-on galaxy in rest-frame velocity units')
      
   end if
   
   ! Group "groups"
   if (para%make_groups) then
   
      name = 'groups/'
      call hdf5_add_group(name)
      call hdf5_write_data(name//'snapshot',sky_group%snapshot,'snapshot index')
      call hdf5_write_data(name//'subvolume',sky_group%subvolume,'subvolume index')
      call hdf5_write_data(name//'tile',sky_group%tile,'tile index in tiling array')
      call hdf5_write_data(name//'zobs',sky_group%zobs,'redshift in observer-frame')
      call hdf5_write_data(name//'zcmb',sky_group%zcmb,'redshift in CMB frame')
      call hdf5_write_data(name//'zcos',sky_group%zcos,'cosmological redshift without peculiar motions')
      call hdf5_write_data(name//'dc',sky_group%dc,'[Mpc/h] comoving distance')
      call hdf5_write_data(name//'ra',sky_group%ra/unit%degree,'[deg] right ascension')
      call hdf5_write_data(name//'dec',sky_group%dec/unit%degree,'[deg] declination')
      call hdf5_write_data(name//'vpec_x',sky_group%vpec(1),'[proper km/s] x-component of peculiar velocity')
      call hdf5_write_data(name//'vpec_y',sky_group%vpec(2),'[proper km/s] y-component of peculiar velocity')
      call hdf5_write_data(name//'vpec_z',sky_group%vpec(3),'[proper km/s] z-component of peculiar velocity')
      call hdf5_write_data(name//'vpec_r',sky_group%vpecrad,'[proper km/s] line-of-sight peculiar velocity')
      call hdf5_write_data(name//'sigma_3d_all',sky_group%sigma_3d_all,&
      &'[proper km/s] 3D peculiar velocity dispersion of ALL group members, including non-detections')
      call hdf5_write_data(name//'sigma_los_detected',sky_group%sigma_los_detected,&
      &'[proper km/s] line-of-sight peculiar velocity dispersion of detected group members')
      call hdf5_write_data(name//'id_group_sky',sky_group%id_group_sky,'unique parent halo ID in mock sky')
      call hdf5_write_data(name//'id_halo_sam',sky_group%id_halo_sam,'parent halo ID in SAM')
      call hdf5_write_data(name//'mvir',sky_group%mvir,'[Msun/h] virial mass')
      call hdf5_write_data(name//'n_galaxies_total',sky_group%group_ntot, &
      & 'total number of galaxies that live in the same group (host halo)')
      call hdf5_write_data(name//'n_galaxies_selected',sky_group%group_nsel, &
      & 'number of galaxies that live in the same group (host halo) and are present in the mock sky')
      call hdf5_write_data(name//'flag',sky_group%group_flag,'group flag: 0 if group complete, '//&
      &'>0 if truncated by survey edge (+1), snapshot limit (+2), tile edge (+4), shell edge (+8)')
      
   end if
   
   ! close HDF5 file
   call hdf5_close()

end subroutine write_hdf5

end module module_user_routines