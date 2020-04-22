! **********************************************************************************************************************************
! This module defines the interface between the semi-analytic model (SAM) and stingray
!
! Each SAM must have its own module, named "module_user_routines_[sam].f03", where "sam" is the name of the semi-analytic model
! specified in the makefile.
! 
! Instructions of how to adapt this module to a particular SAM and/or add new intrinsic or apparent galaxy properties are given
! in the comments below.
! **********************************************************************************************************************************

module module_user_routines

! **********************************************************************************************************************************
! INTERFACE (DO NOT EDIT)
! **********************************************************************************************************************************

! default modules, do not edit
use shared_module_core
use shared_module_hdf5
use shared_module_maths
use shared_module_cosmology
use module_global
use module_conversion
use module_emission_lines

! required objects
public :: type_sam
public :: type_sky_galaxy
public :: type_sky_group ! (can be empty)
public :: make_automatic_parameters ! (can be empty)
public :: make_redshifts
public :: load_sam_snapshot
public :: make_hdf5

private


! **********************************************************************************************************************************
! TYPE DECLARATIONS
! **********************************************************************************************************************************

! Here, specify the class of SAM-properties to be loaded for making the mock sky. How these properties are
! read from files is specified in the subroutine load_sam_snapshot below
! The class requires three class functions:
! get_position: returns xyz-position of the galaxy
! get_groupid: returns the group id of the galaxy
! is_group_center: logical function specifying if the galaxy is the group_center

type type_sam

   integer*4   :: id_galaxy      ! unique galaxy ID
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
   real*4      :: jbulge         ! [km/s * cMpc/h] specific angular momentum of the bulge.
   real*4      :: jdisk          ! [km/s * cMpc/h] specific angular momentum of the disk.
   real*4      :: sfr_disk       ! [Msun/Gyr/h] star formation rate in the disk
   real*4      :: sfr_burst      ! [Msun/Gyr/h] star formation rate in the bulge
   real*4      :: mgas_metals_disk  ! [Msun/h] mass of metals locked up in the disk
   real*4      :: mgas_metals_bulge ! [Msun/h] mass of metals locked up in the bulge
   real*4      :: mvir_hosthalo  ! [Msun/h]
   real*4      :: mvir_subhalo   ! [Msun/h]
   real*4      :: cnfw_subhalo   ! [-] concentration of NFW fit to subhalo
   real*4      :: vvir_hosthalo  ! [km/s]	virial velocity of hosthalo
   real*4      :: vvir_subhalo   ! [km/s]	virial velocity of subhalo
   real*4      :: vmax_subhalo   ! [km/s]	maximum circular velocity of subhalo
   real*4      :: zgas_disk      ! metallicity of the gas in the disk
   real*4      :: zgas_bulge     ! metallicity of the gas in the bulge
   real*4      :: mbh            ! [Msun/h] black hole mass
   real*4      :: mbh_acc_hh     ! [Msun/Gyr/h] accretion rate in hot-halo mode
   real*4      :: mbh_acc_sb     ! [Msun/Gyr/h] accretion rate in starburst mode
   real*4      :: lco_disk
   real*4      :: lco_bulge
   
contains

   procedure   :: get_position      => sam_get_position     ! required function of type real*4(dim=3)
   procedure   :: get_groupid       => sam_get_groupid      ! required function of type integer*8
   procedure   :: is_group_center   => sam_is_group_center  ! required function of type logical*4

end type type_sam

! Here, specify the class or classes of mock object properties in the sky, such as apparent galaxy properties.
! See the existing routines below for clarifications. The class must contain the following class functions:
! make_from_sam
   
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

type,extends(type_sky) :: type_sky_galaxy ! must exist

   integer*8   :: id_galaxy_sky           ! unique ID in the mock sky
   integer*4   :: id_galaxy_sam           ! galaxy ID in the SAM
   integer*4   :: typ                     ! galaxy type (0=central, 1=satellite, 2=orphan)
   
   real*4      :: inclination             ! [rad] inclination
   real*4      :: pa                      ! [rad] position angle from North to East
   real*4      :: mag                     ! apparent magnitude (generic: M/L ratio of 1, no k-correction)
   
   real*4      :: mstars_disk             ! [Msun/h] stellar mass of the disk
   real*4      :: mstars_bulge            ! [Msun/h] stellar mass of the bulge
   real*4      :: mgas_disk               ! [Msun/h] gas mass of the disk
   real*4      :: mgas_bulge              ! [Msun/h] gas mass of the bulge
   real*4      :: matom_disk              ! [Msun/h] atomic gas mass of the disk
   real*4      :: matom_bulge             ! [Msun/h] atomic gas mass of the bulge
   real*4      :: mmol_disk               ! [Msun/h] molecular gas mass of the disk
   real*4      :: mmol_bulge              ! [Msun/h] molecular gas mass of the bulge
   real*4      :: jbulge                  ! [km/s * cMpc/h] specific angular momentum of the bulge.
   real*4      :: jdisk                   ! [km/s * cMpc/h] specific angular momentum of the disk. 
 
   real*4      :: J(3)                    ! [Msun pMpc km/s] total angular momentum
   
   real*4      :: mvir_hosthalo           ! [Msun/h] mass of 1st generation halo (i.e. direct host of type 0 galaxies)
   real*4      :: mvir_subhalo            ! [Msun/h] halo mass
   
   real*4      :: rstar_disk_apparent     ! [arcsec] apparent semi-major axis of half-mass ellipse of stars in the disk
   real*4      :: rstar_bulge_apparent    ! [arcsec] apparent semi-major axis of half-mass ellipse of stars in the bulge
   real*4      :: rgas_disk_apparent      ! [arcsec] apparent semi-major axis of half-mass ellipse of gas in the disk
   real*4      :: rgas_bulge_apparent     ! [arcsec] apparent semi-major axis of half-mass ellipse of gas in the bulge
   
   real*4      :: rstar_disk_intrinsic    ! [cMpc/h] intrinsic half-mass radius of stars in the disk
   real*4      :: rstar_bulge_intrinsic   ! [cMpc/h] intrinsic half-mass radius of stars in the bulge
   real*4      :: rgas_disk_intrinsic     ! [cMpc/h] intrinsic half-mass radius of gas in the disk
   real*4      :: rgas_bulge_intrinsic    ! [cMpc/h] intrinsic half-mass radius of gas in the bulge
   
   real*4      :: zgas_disk      ! metallicity of the gas in the disk
   real*4      :: zgas_bulge     ! metallicity of the gas in the disk

   real*4      :: sfr_disk       ! [Msun/Gyr/h] star formation rate in the disk
   real*4      :: sfr_burst      ! [Msun/Gyr/h] star formation rate in the bulge

   real*4      :: mbh            ! [Msun/h] black hole mass
   real*4      :: mbh_acc_hh     ! [Msun/Gyr/h] accretion rate in hot-halo mode
   real*4      :: mbh_acc_sb     ! [Msun/Gyr/h] accretion rate in starburst mode

   !group_disp_part	?	real*4	km/s	i		?	ASGR: Intrinsic 3D dispersion of the group halo (VR output)
   !group_disp_sat_los	-	real*4	km/s	a			ASGR: Measured line-of-site dispersion of the group satellites (close to what we would measure in real life)

   real*4      :: vvir_hosthalo           ! [km/s]	virial velocity of hosthalo
   real*4      :: vvir_subhalo            ! [km/s]	virial velocity of subhalo
   real*4      :: vmax_subhalo            ! [km/s]	maximum circular velocity of subhalo
   real*4      :: cnfw_subhalo            ! [-] concentration of NFW fit to subhalo
   
   ! HI line   
   real*4      :: hiline_flux_int            ! [W/m^2] integrated HI line flux
   real*4      :: hiline_flux_int_vel        ! [Jy km/s] velocity-integrated HI line flux
   type(type_line_shape) :: hiline_shape     ! shape-parameters of HI emission line (see module module_emission_lines)
      
   ! CO lines
   real*4      :: coline_flux_int(10)        ! [W/m^2] integrated line fluxes of CO(1-0), C(2-1), ..., C(10-9) transition
   real*4      :: coline_flux_int_vel(10)    ! [Jy km/s] velocity-integrated line fluxes of CO(1-0), C(2-1), ..., C(10-9) transition
   type(type_line_shape) :: coline_shape     ! shape-parameters of CO emission lines (see module module_emission_lines)

   contains
   
   procedure   :: make_from_sam  => make_sky_galaxy   ! required subroutine
   
end type type_sky_galaxy

type,extends(type_sky) :: type_sky_group ! must exist
   
   real*4      :: mvir                 ! [Msun/h] virial mass of group
   real*4      :: vvir                 ! [km/s]	virial velocity of group halo
   real*4      :: vmax                 ! [km/s]	maximum circular velocity of group halo
   real*4      :: cnfw                 ! [-] concentration of NFW fit to group halo
   integer*4   :: group_ntot           ! total number of galaxies in group
   integer*4   :: group_nsel           ! number of selected galaxies in group
   integer*4   :: group_flag           ! 0=group complete, >0 group truncated
   real*4      :: sigma_los_detected   ! [km/s] line-of-sight velocity dispersion of selected galaxies
   real*4      :: sigma_3D_all         ! [km/s] 3D velocity dispersion of all galaxies in group
   
   contains
   
   procedure   :: make_from_sam  => make_sky_group   ! required subroutine
   
end type type_sky_group

contains

! In order to place the objects in the mock sky, the class type_sam must have the following function
! enabling stingray to extract the position of each object.

function sam_get_position(self) result(position)
   class(type_sam) :: self
   real*4 :: position(3)
   position = self%position ! [length_unit of simulation] position if the galaxy in the box
end function sam_get_position

! In order to associate different galaxies with groups of galaxies, class type_sam must have the
! following two functions, enabling stingray to extract the group-id and central/satellite-status of each object

integer*8 function sam_get_groupid(self)
   class(type_sam) :: self
   call nil(self) ! dummy statement to avoid compiler warnings
   sam_get_groupid = self%id_halo ! unique group identifier
   ! NB: if groups are irrelevant and/or not available in the SAM,
   ! use "sam_get_groupid = 0" (as a dummy)
   ! and set "make_groups = 0" in the parameter file
end function sam_get_groupid

logical*4 function sam_is_group_center(self)
   class(type_sam) :: self
   call nil(self) ! dummy statement to avoid compiler warnings
   sam_is_group_center = self%typ==0
   ! NB: if groups are irrelevant and/or not available in the SAM,
   ! use "sam_is_group_center = .true." (as a dummy)
   ! and set "make_groups = 0" in the parameter file
end function sam_is_group_center


! **********************************************************************************************************************************
! CONVERSION BETWEEN INTRINSIC AND APPARENT GALAXY PROPERTIES
! **********************************************************************************************************************************

subroutine make_sky_object(sky_object,sam,base,groupid)

   implicit none

   class(type_sky),intent(out)            :: sky_object
   type(type_sam),intent(in)              :: sam
   type(type_base),intent(in)             :: base                 ! basic properties of the position of this galaxy in the sky
   integer*8,intent(in)                   :: groupid              ! unique group index in sky
   real*4                                 :: pos(3)               ! [simulation length units] position vector of galaxy
   real*4                                 :: elos(3)              ! unit vector pointing from the observer to the object in comoving space
   
   call nil(sky_object,sam,base,groupid) ! dummy statement to avoid compiler warnings

   ! indices
   sky_object%snapshot      = sam%snapshot
   sky_object%subvolume     = sam%subvolume
   sky_object%tile          = base%tile
   sky_object%id_halo_sam   = sam%id_halo
   sky_object%id_group_sky  = groupid
   
   ! sky coordinates
   sky_object%dc  = base%pos%dc*para%box_side    ! [simulation length unit = Mpc/h]
   sky_object%ra  = base%pos%ra                  ! [rad]
   sky_object%dec = base%pos%dec                 ! [rad]
   
   ! redshifts from the position and velocity of the object
   call sph2car(sky_object%dc,sky_object%ra,sky_object%dec,pos,astro=.true.)
   elos = pos/norm(pos) ! line-of-sight vector
   call make_redshift(pos*(para%length_unit/unit%Mpc),sam%velocity,zobs=sky_object%zobs,zcmb=sky_object%zcmb,zcos=sky_object%zcos)
   
   ! peculiar velocity
   sky_object%vpec      = convert_vector(sam%velocity,tileindex=base%tile,ispseudovector=.false.)
   sky_object%vpecrad   = sum(sky_object%vpec*elos)
   
end subroutine make_sky_object

subroutine make_sky_galaxy(sky_galaxy,sam,base,groupid,galaxyid)

   implicit none
   class(type_sky_galaxy),intent(out)  :: sky_galaxy
   type(type_sam),intent(in)           :: sam
   type(type_base),intent(in)          :: base                 ! basic properties of the position of this galaxy in the sky
   integer*8,intent(in)                :: groupid              ! unique group in sky index
   integer*8,intent(in)                :: galaxyid             ! (only for galaxy types) unique galaxy in sky index
   
   real*4                              :: pos(3)               ! [simulation length units] position vector of galaxy
   real*4                              :: dl                   ! [simulation length units] luminosity distance to observer
   real*4                              :: elos(3)              ! unit vector pointing from the observer to the object in comoving space
   real*4                              :: mstars,mHI,mH2,LCO
   
   real*4,parameter        :: L2MHI = 6.27e-9 ! (LHI/Lsun)/(MHI/Msun)
   real*4,parameter        :: sigma_gas = 10.0 ! [km/s] standard velocity dispersion of cold gas
   real*4,parameter        :: X_CO = 2.12 ! [1e20/(K km/s cm^2)] H2-CO(1-0) conversion factor of the MW (alpha=3.4)
   
   
   call nil(sky_galaxy,sam,base,groupid,galaxyid) ! dummy statement to avoid compiler warnings
   
   ! basics
   call make_sky_object(sky_galaxy,sam,base,groupid)
   
   
   ! INTRINSIC PROPERTIES
   
   ! make IDs
   sky_galaxy%id_galaxy_sam   = sam%id_galaxy   ! copy other properties
   sky_galaxy%id_galaxy_sky   = galaxyid        ! unique galaxy id
      
   ! basic properties
   sky_galaxy%typ             = sam%typ
   
   ! intrinsic halo properties
   sky_galaxy%mvir_hosthalo   = sam%mvir_hosthalo
   sky_galaxy%mvir_subhalo    = sam%mvir_subhalo
   sky_galaxy%vvir_hosthalo   = sam%vvir_hosthalo
   sky_galaxy%vvir_subhalo    = sam%vvir_subhalo
   sky_galaxy%vmax_subhalo    = sam%vmax_subhalo
   sky_galaxy%cnfw_subhalo    = sam%cnfw_subhalo
   
   ! intrinsic masses
   sky_galaxy%mstars_disk     = sam%mstars_disk
   sky_galaxy%mstars_bulge    = sam%mstars_bulge
   sky_galaxy%mgas_disk       = sam%mgas_disk   
   sky_galaxy%mgas_bulge      = sam%mgas_bulge  
   sky_galaxy%matom_disk      = sam%matom_disk  
   sky_galaxy%matom_bulge     = sam%matom_bulge 
   sky_galaxy%mmol_disk       = sam%mmol_disk   
   sky_galaxy%mmol_bulge      = sam%mmol_bulge  
   sky_galaxy%sfr_disk        = sam%sfr_disk
   sky_galaxy%sfr_burst       = sam%sfr_burst

   ! intrinsic angular momentum
   sky_galaxy%J      = convert_vector(sam%J,tileindex=base%tile,ispseudovector=.true.)
   sky_galaxy%jbulge = sam%jbulge
   sky_galaxy%jdisk  = sam%jdisk

   ! intrinsic radii
   sky_galaxy%rstar_disk_intrinsic = sam%rstar_disk ! [cMpc/h]
   sky_galaxy%rstar_bulge_intrinsic = sam%rstar_bulge ! [cMpc/h]
   sky_galaxy%rgas_disk_intrinsic = sam%rgas_disk ! [cMpc/h]
   sky_galaxy%rgas_bulge_intrinsic = sam%rgas_bulge ! [cMpc/h]
 
   if(sam%mgas_disk > 0) then 
      sky_galaxy%zgas_disk   = sam%mgas_metals_disk / sam%mgas_disk
   else 
      sky_galaxy%zgas_disk   = 0
   end if 
   if(sam%mgas_bulge > 0) then 
      sky_galaxy%zgas_bulge  = sam%mgas_metals_bulge/ sam%mgas_bulge
   else 
      sky_galaxy%zgas_bulge  = 0
   end if
   

   sky_galaxy%mbh = sam%mbh
   sky_galaxy%mbh_acc_sb = sam%mbh_acc_sb
   sky_galaxy%mbh_acc_hh = sam%mbh_acc_hh

   ! APPARENT PROPERTIES
   
   ! inclination and position angle
   call sph2car(sky_galaxy%dc,sky_galaxy%ra,sky_galaxy%dec,pos,astro=.true.)
   elos = pos/norm(pos) ! position vector
   call make_inclination_and_pa(pos,sky_galaxy%J,inclination=sky_galaxy%inclination,pa=sky_galaxy%pa)
   
   ! apparent radii (note division by comoving distance, because intrinsic radii given in comoving units)
   sky_galaxy%rstar_disk_apparent = sam%rstar_disk/sky_galaxy%dc/unit%arcsec ! [arcsec]
   sky_galaxy%rstar_bulge_apparent = sam%rstar_bulge/sky_galaxy%dc/unit%arcsec ! [arcsec]
   sky_galaxy%rgas_disk_apparent = sam%rgas_disk/sky_galaxy%dc/unit%arcsec ! [arcsec]
   sky_galaxy%rgas_bulge_apparent = sam%rgas_bulge/sky_galaxy%dc/unit%arcsec ! [arcsec]
   
   ! rough generic optical magnitude
   dl = sky_galaxy%dc*(1+sky_galaxy%zobs) ! [Mpc/h]
   mstars = (sam%mstars_disk+sam%mstars_bulge)/para%h ! [Msun]
   sky_galaxy%mag = convert_absmag2appmag(convert_stellarmass2absmag(mstars,1.0),dl/para%h)
      
   ! cold gas emission lines
   mHI = (sam%matom_disk+sam%matom_bulge)/1.35 ! [Msun/h] HI mass  
   mH2 = ((sam%mgas_disk+sam%mgas_bulge)-(sam%matom_disk+sam%matom_bulge))/1.35 ! [Msun/h] H2 mass
   sky_galaxy%hiline_flux_int = real(convert_luminosity2flux(real(mHI/para%h,8)*real(L2MHI,8)*unit%Lsun,dl/para%h),4) ! [W/m^2]
   sky_galaxy%hiline_flux_int_vel = convert_intflux2velintflux(sky_galaxy%hiline_flux_int,0.21106114,sky_galaxy%zobs)
   LCO = mH2/para%h/(313*X_CO) ! [Jy km/s Mpc^2] CO(1-0) luminosity (note this equation has a J^2 dependence)
   sky_galaxy%coline_flux_int = LCO/(4.0*pi*(dl/para%h)**2)/0.00260075761e23 ! [W/m^2] integrated flux
   sky_galaxy%coline_flux_int_vel = convert_intflux2velintflux(sky_galaxy%hiline_flux_int,0.00260075761,sky_galaxy%zobs)
   if (para%line_parameters==1) call make_line_profiles
   
contains
   
   subroutine make_line_profiles
   
      implicit none
      type(type_line_input)               :: line_input
      type(type_line_shape)               :: line_shape(2)
      real*4                              :: a,c_halo,cMpch2pkpc,rvir
      real,parameter                      :: G = 4.30241e-6     ! [(km/s)^2 kpc/Msun] gravitational constant
      
      ! scale factor for comoving-physical conversion
      a = 1.0/(1.0+snapshot(sam%snapshot)%redshift) ! scale factor at redshift os snapshot
      cMpch2pkpc = 1e3*a/para%h
   
      ! halo
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
   
      ! disk
      line_input%mdisk      = sam%mstars_disk+sam%mgas_disk   ! [Msun/h] disk mass (all baryons)
      if (line_input%mdisk>0.0) then
         line_input%rdisk   = (sam%mstars_disk*sam%rstar_disk+sam%mgas_disk*sam%rgas_disk)/line_input%mdisk*cMpch2pkpc ! [pkpc] half-mass radius
      else
         line_input%rdisk   = 1.0 ! just a dummy value that is irrelevant (but must be different from 0)
      end if
   
      ! bulge
      line_input%mbulg      = sam%mstars_bulge+sam%mgas_bulge   ! [Msun/h] bulge mass
      if (line_input%mbulg>0.0) then
         line_input%rbulg   = (sam%mstars_bulge*sam%rstar_bulge+sam%mgas_bulge*sam%rgas_bulge)/line_input%mbulg*cMpch2pkpc ! [pkpc] half-mass radius
      else
         line_input%rbulg   = 1.0 ! just a dummy value that is irrelevant (but must be different from 0)
      end if
      
      ! correct halo mass for disk and bulge
      line_input%mhalo      = max(0.0,line_input%mhalo-line_input%mdisk-line_input%mbulg) 
   
      ! gas
      line_input%mHI        = mHI ! [Msun/h] HI mass
      line_input%mH2        = mH2 ! [Msun/h] H2 mass
      line_input%rgas       = (sam%mgas_disk*sam%rgas_disk+sam%mgas_bulge*sam%rgas_bulge)/(sam%mgas_disk+sam%mgas_bulge)*cMpch2pkpc ! [pkpc] half-mass radius of cold gas (CHECK)
      !line_input%rgas       = sam%rgas_disk*cMpch2pkpc
      
      ! other properties
      line_input%dispersionHI = sigma_gas ! [km/s] velocity dispersion of HI
      line_input%dispersionH2 = sigma_gas ! [km/s] velocity dispersion of H2
      line_input%z          = sky_galaxy%zobs    ! [-] redshift
      line_input%incl       = sky_galaxy%inclination    ! [rad] inclination (0 = face-on, pi/2 = edge-on)
      line_input%Dc         = sky_galaxy%dc     ! [Mpc/h] comoving distance

      ! make line profile
      if ((line_input%mHI>0.0).or.(line_input%mH2>0.0)) call make_emission_line(line_input,line_shape)
   
      ! save line parameters
      if (line_input%mHI>0.0) then
         sky_galaxy%hiline_shape = line_shape(1)
      else
         sky_galaxy%hiline_shape = empty_line_shape
      end if
      if (line_input%mH2>0.0) then
         sky_galaxy%coline_shape = line_shape(2)
      else
         sky_galaxy%coline_shape = empty_line_shape
      end if
         
   end subroutine make_line_profiles
   
end subroutine make_sky_galaxy

subroutine make_sky_group(sky_group,sam,sky_galaxy,selected,base,groupid,group_nselected)

   ! the first object in sam, sky_galaxy, selected is the central galaxy of the group;
   ! sky_galaxy only exists for selected objects

   implicit none
   class(type_sky_group),intent(out)   :: sky_group
   type(type_sam),intent(in)           :: sam(:)
   type(type_sky_galaxy),intent(in)    :: sky_galaxy(:)
   logical*4,intent(in)                :: selected(:)
   type(type_base),intent(in)          :: base                 ! basic properties of the position of this galaxy in the sky
   integer*8,intent(in)                :: groupid              ! unique group in sky index
   integer*4,intent(in)                :: group_nselected   
   integer*4                           :: i,n,count
   real*4                              :: dv,v0(3),vr
   
   call nil(sky_group,sam(1),sky_galaxy(1),selected(1),base,groupid,group_nselected) ! dummy statement to avoid compiler warnings
   
   ! base properties
   call make_sky_object(sky_group,sam(1),base,groupid)
   
   ! basic halo properties
   sky_group%mvir          = sam(1)%mvir_hosthalo
   sky_group%vvir          = sam(1)%vvir_hosthalo
   sky_group%vmax          = sam(1)%vmax_subhalo
   sky_group%cnfw          = sam(1)%cnfw_subhalo
   
   ! group properties
   sky_group%group_ntot    = base%group_ntot
   sky_group%group_flag    = base%group_flag
   sky_group%group_nsel    = group_nselected
   
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
   count = 0
   if (group_nselected>=2) then
      vr = 0
      do i = 1,n
         if (selected(i)) vr = vr+sky_galaxy(i)%vpecrad
      end do
      vr = vr/group_nselected
      dv = 0
      do i = 1,n
         if (selected(i)) then
            dv = dv+(sky_galaxy(i)%vpecrad-vr)**2
            count = count+1
         end if
      end do
      if (count.ne.group_nselected) call error('count.ne.group_nselected')
      sky_group%sigma_los_detected = sqrt(dv/(group_nselected-1)) ! standard deviation
   else
      sky_group%sigma_los_detected = 0
   end if
   
end subroutine make_sky_group


! **********************************************************************************************************************************
! I/O routines
! **********************************************************************************************************************************

! Function returns filename of the SAM galaxies, given a snapshot and subvolume index

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


! Set parameters automatically (e.g. from information provided with the snapshot files).
! These automatic parameters are only adopted if the values in the parameter files are set to 'auto'.
! Otherwise the parameter file overwrites these values.

subroutine make_automatic_parameters
   
   implicit none
   
   character(len=255)   :: filename
   integer*4            :: isnapshot,isubvolume
   
   ! assume that first subvolume always exists
   para%subvolume_min = 0
   
   ! find minimal snapshot, assuming that the first subvolume exists
   isnapshot = -1
   do while (.not.exists(filename_sam(isnapshot,0,1)).and.(isnapshot<limit%n_snapshots_max))
      isnapshot = isnapshot+1
   end do
   if (isnapshot==limit%n_snapshots_max) call error('no snapshots found in '//trim(para%path_input))
   para%snapshot_min = isnapshot
   
   ! find maximal snapshot, assuming that the first subvolume exists
   do while (exists(filename_sam(isnapshot,0,1)))
      isnapshot = isnapshot+1
   end do
   para%snapshot_max = isnapshot-1
   
   ! find maximal subvolume, assuming that the last snapshot is representative
   isubvolume = 0
   do while (exists(filename_sam(para%snapshot_max,isubvolume,1)))
      isubvolume = isubvolume+1
   end do
   para%subvolume_max = isubvolume-1
   
   filename = filename_sam(para%snapshot_min,0,1)
   call out('File of automatic parameters: '//trim(filename))
   call hdf5_open(filename)
   call hdf5_read_data('/cosmology/h',para%h)
   call hdf5_read_data('/cosmology/omega_l',para%omega_l)
   call hdf5_read_data('/cosmology/omega_m',para%omega_m)
   call hdf5_read_data('/cosmology/omega_b',para%omega_b)
   call hdf5_read_data('/run_info/lbox',para%box_side)
   para%length_unit = unit%Mpc/para%h ! [m] length unit expressed in SI
   call hdf5_close()
   
end subroutine make_automatic_parameters

! load redshifts
! this routine must write the redshift of each snapshot into the real*4-array
! snapshot(isnapshot)%redshift

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

subroutine load_sam_snapshot(isnapshot,isubvolume,sam)

   ! variable declaration
   implicit none
   integer*4,intent(in)                            :: isnapshot         ! snapshot index
   integer*4,intent(in)                            :: isubvolume        ! subvolume index
   type(type_sam),allocatable,intent(out)          :: sam(:)            ! class of SAM properties
   integer*8                                       :: n                 ! number of galaxies
   character(*),parameter                          :: g = '/galaxies/'  ! group name
   !integer*8,allocatable                           :: id_check(:)
   
   ! read main galaxy properties
   
   ! open file
   call hdf5_open(filename_sam(isnapshot,isubvolume,1)) ! NB: this routine also checks if the file exists
   
   ! determine number of galaxies in this (sub)snapshot
   n = hdf5_dataset_size(g//'id_galaxy')

   ! allocate array
   if (allocated(sam)) deallocate(sam)
   allocate(sam(n))
   
   ! read file
   call hdf5_read_data(g//'id_galaxy',sam%id_galaxy,convert=.true.)
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
   call hdf5_read_data(g//'specific_angular_momentum_disk_star',sam%jdisk)
   call hdf5_read_data(g//'specific_angular_momentum_bulge_star',sam%jbulge)
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
   call hdf5_read_data(g//'sfr_disk',sam%sfr_disk)
   call hdf5_read_data(g//'sfr_burst',sam%sfr_burst)
   call hdf5_read_data(g//'m_bh',sam%mbh)
   call hdf5_read_data(g//'bh_accretion_rate_hh',sam%mbh_acc_hh)
   call hdf5_read_data(g//'bh_accretion_rate_sb',sam%mbh_acc_sb)
   call hdf5_read_data(g//'mgas_metals_disk',sam%mgas_metals_disk)
   call hdf5_read_data(g//'mgas_metals_bulge',sam%mgas_metals_bulge)
   call hdf5_read_data(g//'mvir_hosthalo',sam%mvir_hosthalo)
   call hdf5_read_data(g//'mvir_subhalo',sam%mvir_subhalo)
   call hdf5_read_data(g//'cnfw_subhalo',sam%cnfw_subhalo)
   call hdf5_read_data(g//'vvir_subhalo',sam%vvir_subhalo)
   call hdf5_read_data(g//'vmax_subhalo',sam%vmax_subhalo)
   call hdf5_read_data(g//'vvir_hosthalo',sam%vvir_hosthalo)
   
   ! close file
   call hdf5_close()
   
   !! open file
   !call hdf5_open(filename_sam(isnapshot,isubvolume,2)) ! NB: this routine also checks if the file exists
   !
   !! check number of galaxies
   !if (hdf5_dataset_size(g//'id_galaxy').ne.n) call error('wrong number of galaxies in CO file')
   !
   !! read file
   !call hdf5_read_data(g//'LCO_bulge',sam%lco_bulge)
   !call hdf5_read_data(g//'LCO_disk',sam%lco_disk)
   !
   !! close file
   !call hdf5_close()
   
   ! assign other properties
   sam%snapshot = isnapshot
   sam%subvolume = isubvolume
   
end subroutine load_sam_snapshot

! write sky_galaxy(:) and sky_group(:) into a HDF5 file

subroutine make_hdf5(filename_hdf5,sky_galaxy,sky_group,total_stats,subvolume_stats,isubvolume)
   
   implicit none
   character(*),intent(in)                      :: filename_hdf5  ! output filename
   type(type_sky_galaxy),intent(in),allocatable :: sky_galaxy(:)
   type(type_sky_group),intent(in),allocatable  :: sky_group(:)
   type(type_skystats),intent(in)               :: total_stats
   type(type_skystats),intent(in),optional      :: subvolume_stats ! only provided if the sky corresponds to a subvolume
   integer*4,intent(in),optional                :: isubvolume
   character(:),allocatable                     :: name
   character(len=255)                           :: shark_version,shark_git_revision,shark_timestamp
   integer*4                                    :: i,j
   
   allocate(character(1)::name) ! empty allocation to avoid compiler flags
   
   ! load auxilary data from shark output
   call hdf5_open(filename_sam(para%snapshot_min,0,1))
   call hdf5_read_data('/run_info/shark_version',shark_version)
   call hdf5_read_data('/run_info/shark_git_revision',shark_git_revision)
   call hdf5_read_data('/run_info/timestamp',shark_timestamp)
   call hdf5_close()
   
   ! create and open HDF5 file
   call hdf5_create(filename_hdf5)
   call hdf5_open(filename_hdf5,.true.)
   
   ! write group "parameters"
   name = 'parameters'
   call hdf5_add_group(name)
   call hdf5_write_data(name//'/survey',para%survey,'name of simulated survey')
   call hdf5_write_data(name//'/path_output',para%path_output)
   call hdf5_write_data(name//'/path_input',para%path_input)
   call hdf5_write_data(name//'/box_side',para%box_side, &
   & '[length_unit] comoving side-length of simulation box in multiples of length_unit')
   call hdf5_write_data(name//'/length_unit',para%length_unit,'[m] SI-value of comoving length unit')
   call hdf5_write_data(name//'/snapshot_min',para%snapshot_min,'index of earliest snapshot used for the mock sky')
   call hdf5_write_data(name//'/snapshot_max',para%snapshot_max,'index of latest snapshot used for the mock sky')
   call hdf5_write_data(name//'/subvolume_min',para%subvolume_min,'index of first subvolume used for the mock sky')
   call hdf5_write_data(name//'/subvolume_max',para%subvolume_max,'index of last subvolume used for the mock sky')
   call hdf5_write_data(name//'/h',para%h,'[-] Hubble parameter H0=h*100 km/s/Mpc')
   call hdf5_write_data(name//'/omega_l',para%omega_l,'energy density of dark energy relative to closure density')
   call hdf5_write_data(name//'/omega_m',para%omega_m,'energy density of all matter relative to closure density')
   call hdf5_write_data(name//'/omega_b',para%omega_b,'energy density of baryonic matter relative to closure density')
   call hdf5_write_data(name//'/dc_min',para%dc_min,'[length_unit] minimal comoving distance of mock sky')
   call hdf5_write_data(name//'/dc_max',para%dc_max,'[length_unit] maximal comoving distance of mock sky')
   call hdf5_write_data(name//'/ra_min',para%ra_min/unit%degree,'[deg] min right ascension of mock sky')
   call hdf5_write_data(name//'/ra_max',para%ra_max/unit%degree,'[deg] max right ascension of mock sky')
   call hdf5_write_data(name//'/dec_min',para%dec_min/unit%degree,'[deg] min declination of mock sky')
   call hdf5_write_data(name//'/dec_max',para%dec_max/unit%degree,'[deg] max declination of mock sky')
   call hdf5_write_data(name//'/zaxis_ra',para%zaxis_ra/unit%degree, &
   & '[deg] RA coordinate of the SAM z-axis in spherical survey coordinates')
   call hdf5_write_data(name//'/zaxis_dec',para%zaxis_dec/unit%degree, &
   & '[deg] Dec coordinate the SAM z-axis in spherical survey coordinates')
   call hdf5_write_data(name//'/xy_angle',para%xy_angle/unit%degree,'[deg] Rotation of the SAM (x,y)-plane on the sky')
   call hdf5_write_data(name//'/seed',para%seed,'seed for the random number generator of symmetry operations')
   call hdf5_write_data(name//'/translate',para%translate, &
   & 'logical flag specifying if random translations are applied (0=false, 1=true)')
   call hdf5_write_data(name//'/rotate',para%rotate, &
   & 'logical flag specifying if random rotations are applied (0=false, 1=true)')
   call hdf5_write_data(name//'/invert',para%invert, &
   & 'logical flag specifying if random inversions are applied (0=false, 1=true)')
   call hdf5_write_data(name//'/velocity_norm',para%velocity_norm, &
   & '[km/s] observer velocity relative to CMB rest-frame')
   call hdf5_write_data(name//'/velocity_ra',para%velocity_ra, &
   & '[deg] RA coordinate to which the observer is moving relative to CMB rest-frame')
   call hdf5_write_data(name//'/velocity_dec',para%velocity_dec, &
   & '[deg] Dec coordinate to which the observer is moving relative to CMB rest-frame')
   call hdf5_write_data(name//'/search_angle',para%search_angle, &
   & '[deg] typical angle in which overlaps between survey volume and tiling grid are searched')
   call hdf5_write_data(name//'/volume_search_level',para%volume_search_level, &
   & 'specifies the number of search points (2^#)^3 inside each tile')
   call hdf5_write_data(name//'/line_parameters',para%line_parameters, &
   & 'logical flag specifying if global emission line parameters are saved (0=false, 1=true)')
   call hdf5_write_data(name//'/keep_binaries',para%keep_binaries, &
   & 'logical flag specifying if binary output files are kept in additino to this HDF5 (0=false, 1=true)')
   call hdf5_write_data(name//'/keep_log',para%keep_binaries, &
   & 'logical flag specifying if the logfile is kept after successful runs (0=false, 1=true)')
   call hdf5_write_data(name//'/skyrotation',para%sky_rotation, &
   & 'Rotation matrix to map xyz-coordinates of the tiling structure onto sky coordinates.')
   
   ! write group "galaxies"
   name = 'galaxies'
   call hdf5_add_group(name)
   call hdf5_write_data(name//'/snapshot',sky_galaxy%snapshot,'snapshot index')
   call hdf5_write_data(name//'/subvolume',sky_galaxy%subvolume,'subvolume index')
   call hdf5_write_data(name//'/tile',sky_galaxy%tile,'tile index in tiling array')
   call hdf5_write_data(name//'/zobs',sky_galaxy%zobs,'redshift in observer-frame')
   call hdf5_write_data(name//'/zcmb',sky_galaxy%zcmb,'redshift in CMB frame')
   call hdf5_write_data(name//'/zcos',sky_galaxy%zcos,'cosmological redshift without peculiar motions')
   call hdf5_write_data(name//'/dc',sky_galaxy%dc,'[Mpc/h] comoving distance')
   call hdf5_write_data(name//'/ra',sky_galaxy%ra/unit%degree,'[deg] right ascension')
   call hdf5_write_data(name//'/dec',sky_galaxy%dec/unit%degree,'[deg] declination')
   call hdf5_write_data(name//'/id_galaxy_sky',(/(i,i=lbound(sky_galaxy,1),ubound(sky_galaxy,1))/), &
   & 'unique galaxy ID in mock sky, from 1 to n')
   call hdf5_write_data(name//'/id_galaxy_sky_smart',sky_galaxy%id_galaxy_sky,'unique galaxy ID in mock sky, smart format')
   if (para%make_groups==1) call hdf5_write_data(name//'/id_group_sky',sky_galaxy%id_group_sky, &
   & 'unique group ID if galaxy is in a group, -1 otherwise')
   call hdf5_write_data(name//'/id_galaxy_sam',sky_galaxy%id_galaxy_sam,'galaxy ID in SAM')
   call hdf5_write_data(name//'/id_halo_sam',sky_galaxy%id_halo_sam,'host halo ID in SAM')
   call hdf5_write_data(name//'/type',sky_galaxy%typ,'galaxy type (0=central, 1=satellite in halo, 2=orphan)')
   call hdf5_write_data(name//'/inclination',sky_galaxy%inclination/unit%degree, &
   & '[deg] inclination = angle between line-of-sight and spin axis')
   call hdf5_write_data(name//'/pa',sky_galaxy%pa/unit%degree,'[deg] position angle from north to east')
   call hdf5_write_data(name//'/mag',sky_galaxy%mag, &
   & 'apparent magnitude (generic: M/L ratio of 1, no k-correction)')
   call hdf5_write_data(trim(name)//'/vpec_x',sky_galaxy%vpec(1),'[proper km/s] x-component of peculiar velocity')
   call hdf5_write_data(trim(name)//'/vpec_y',sky_galaxy%vpec(2),'[proper km/s] y-component of peculiar velocity')
   call hdf5_write_data(trim(name)//'/vpec_z',sky_galaxy%vpec(3),'[proper km/s] z-component of peculiar velocity')
   call hdf5_write_data(trim(name)//'/vpec_r',sky_galaxy%vpecrad,'[proper km/s] line-of-sight peculiar velocity')
   call hdf5_write_data(trim(name)//'/mstars_disk',sky_galaxy%mstars_disk,'[Msun/h] stellar mass in the disk')
   call hdf5_write_data(trim(name)//'/mstars_bulge',sky_galaxy%mstars_bulge,'[Msun/h] stellar mass in the bulge')
   call hdf5_write_data(trim(name)//'/mgas_disk',sky_galaxy%mgas_disk,'[Msun/h] gas mass in the disk')
   call hdf5_write_data(trim(name)//'/mgas_bulge',sky_galaxy%mgas_bulge,'[Msun/h] gas mass in the bulge')
   call hdf5_write_data(trim(name)//'/matom_disk',sky_galaxy%matom_disk,'[Msun/h] atomic gas mass in the disk')
   call hdf5_write_data(trim(name)//'/matom_bulge',sky_galaxy%matom_bulge,'[Msun/h] atomic mass in the bulge')
   call hdf5_write_data(trim(name)//'/mmol_disk',sky_galaxy%mmol_disk,'[Msun/h] molecular gas mass in the disk')
   call hdf5_write_data(trim(name)//'/mmol_bulge',sky_galaxy%mmol_bulge,'[Msun/h] molecular gas mass in the bulge')
   call hdf5_write_data(trim(name)//'/l_x',sky_galaxy%J(1),'[Msun pMpc km/s] x-component of total angular momentum')
   call hdf5_write_data(trim(name)//'/l_y',sky_galaxy%J(2),'[Msun pMpc km/s] y-component of total angular momentum')
   call hdf5_write_data(trim(name)//'/l_z',sky_galaxy%J(3),'[Msun pMpc km/s] z-component of total angular momentum')
   call hdf5_write_data(trim(name)//'/jdisk',sky_galaxy%jdisk, '[km/s * cMpc/h] specific angular momentum of the disk')
   call hdf5_write_data(trim(name)//'/jbulge',sky_galaxy%jbulge, '[km/s * cMpc/h] specific angular momentum of the bulge')
   call hdf5_write_data(trim(name)//'/mvir_hosthalo',sky_galaxy%mvir_hosthalo,'[Msun/h] host halo mass')
   call hdf5_write_data(trim(name)//'/mvir_subhalo',sky_galaxy%mvir_subhalo,'[Msun/h] subhalo mass')
   call hdf5_write_data(trim(name)//'/zgas_disk',sky_galaxy%zgas_disk,'metallicity of the gas in the disk')
   call hdf5_write_data(trim(name)//'/zgas_bulge',sky_galaxy%zgas_bulge,'metallicity of the gas in the bulge')
   call hdf5_write_data(trim(name)//'/sfr_disk',sky_galaxy%sfr_disk,'star formation rate in the disk [Msun/Gyr/h]')
   call hdf5_write_data(trim(name)//'/sfr_burst',sky_galaxy%sfr_burst,'star formation rate in the bulge [Msun/Gyr/h]')

   call hdf5_write_data(trim(name)//'/mbh',sky_galaxy%mbh,'[Msun/h] black hole mass')
   call hdf5_write_data(trim(name)//'/mbh_acc_hh',sky_galaxy%mbh_acc_hh,'[Msun/Gyr/h] black hole accretion rate in the hot halo mode')
   call hdf5_write_data(trim(name)//'/mbh_acc_sb',sky_galaxy%mbh_acc_sb,'[Msun/Gyr/h] black hole accretion rate in the starburst mode')

   call hdf5_write_data(trim(name)//'/vvir_hosthalo',sky_galaxy%vvir_hosthalo,'[km/s] host halo virial velocity')
   call hdf5_write_data(trim(name)//'/vvir_subhalo',sky_galaxy%vvir_subhalo,'[Msun/h] subhalo virial velocity')
   call hdf5_write_data(trim(name)//'/cnfw_subhalo',sky_galaxy%cnfw_subhalo,'narro-frenk-white concentration of subhalo')

   call hdf5_write_data(trim(name)//'/rstar_disk_apparent',sky_galaxy%rstar_disk_apparent,&
   &'[arcsec] apparent semi-major axis of half-mass ellipse of stellar disk')
   call hdf5_write_data(name//'/rstar_bulge_apparent',sky_galaxy%rstar_bulge_apparent,&
   &'[arcsec] apparent semi-major axis of half-mass ellipse of stellar bulge')
   call hdf5_write_data(name//'/rgas_disk_apparent',sky_galaxy%rgas_disk_apparent,&
   &'[arcsec] apparent semi-major axis of half-mass ellipse of gas disk')
   call hdf5_write_data(name//'/rgas_bulge_apparent',sky_galaxy%rgas_bulge_apparent,&
   &'[arcsec] apparent semi-major axis of half-mass ellipse of gas bulge')
   call hdf5_write_data(name//'/rstar_disk_intrinsic',sky_galaxy%rstar_disk_intrinsic,&
   &'[cMpc/h] intrinsic half-mass radius of stellar disk')
   call hdf5_write_data(name//'/rstar_bulge_intrinsic',sky_galaxy%rstar_bulge_intrinsic,&
   &'[cMpc/h] intrinsic half-mass radius of stellar bulge')
   call hdf5_write_data(name//'/rgas_disk_intrinsic',sky_galaxy%rgas_disk_intrinsic,&
   &'[cMpc/h] intrinsic half-mass radius of gas disk')
   call hdf5_write_data(name//'/rgas_bulge_intrinsic',sky_galaxy%rgas_bulge_intrinsic,&
   &'[cMpc/h] intrinsic half-mass radius of gas bulge')
   call hdf5_write_data(name//'/hiline_flux_int',sky_galaxy%hiline_flux_int,'[W/m^2] integrated HI line flux')
   call hdf5_write_data(name//'/hiline_flux_int_vel',sky_galaxy%hiline_flux_int_vel, &
   & '[Jy km/s] velocity-integrated HI line flux')
   do j = 1,10 ! CO transition
      call hdf5_write_data(name//'/coline_flux_int_'//val2str(j),sky_galaxy%coline_flux_int(j),&
      &'[W/m^2] integrated CO('//val2str(j)//'-'//val2str(j-1)//') line flux')
      call hdf5_write_data(name//'/coline_flux_int_vel_'//val2str(j),sky_galaxy%coline_flux_int_vel(j),&
      &'[Jy km/s] velocity-integrated CO('//val2str(j)//'-'//val2str(j-1)//') line flux')
   end do
   
   if (para%line_parameters==1) then
   
      call hdf5_write_data(name//'/hiline_flux_peak',sky_galaxy%hiline_shape%speak, &
      & '[s/km] normalised peak HI line flux density of inclined galaxy (multiply by hiline_flux_int_vel to get Jy values)')
      call hdf5_write_data(name//'/hiline_flux_central',sky_galaxy%hiline_shape%scentral, &
      & '[s/km] normalised central HI line flux density of inclined galaxy (multiply by hiline_flux_int_vel to get Jy values)')
      call hdf5_write_data(name//'/hiline_width_peak',sky_galaxy%hiline_shape%wpeak, &
      & '[km/s] HI line-width between flux peaks of inclined galaxy in rest-frame velocity units')
      call hdf5_write_data(name//'/hiline_width_50',sky_galaxy%hiline_shape%w50, &
      & '[km/s] HI line-width at 50% of the peak flux of inclined galaxy in rest-frame velocity units')
      call hdf5_write_data(name//'/hiline_width_20',sky_galaxy%hiline_shape%w20, &
      & '[km/s] HI line-width at 20% of the peak flux of inclined galaxy in rest-frame velocity units')
      
      call hdf5_write_data(name//'/hiline_flux_peak_eo',sky_galaxy%hiline_shape%speak_eo, &
      & '[s/km] normalised peak HI line flux density of edge-on galaxy (multiply by hiline_flux_int_vel to get Jy values)')
      call hdf5_write_data(name//'/hiline_flux_central_eo',sky_galaxy%hiline_shape%scentral_eo, &
      & '[s/km] normalised central HI line flux density of edge-on galaxy (multiply by hiline_flux_int_vel to get Jy values)')
      call hdf5_write_data(name//'/hiline_width_peak_eo',sky_galaxy%hiline_shape%wpeak_eo, &
      & '[km/s] HI line-width between flux peaks of edge-on galaxy in rest-frame velocity units')
      call hdf5_write_data(name//'/hiline_width_50_eo',sky_galaxy%hiline_shape%w50_eo, &
      & '[km/s] HI line-width at 50% of the peak flux of edge-on galaxy in rest-frame velocity units')
      call hdf5_write_data(name//'/hiline_width_20_eo',sky_galaxy%hiline_shape%w20_eo, &
      & '[km/s] HI line-width at 20% of the peak flux of edge-on galaxy in rest-frame velocity units')
      
      call hdf5_write_data(name//'/coline_flux_peak',sky_galaxy%coline_shape%speak, &
      & '[s/km] normalised peak CO line flux density of inclined galaxy (multiply by coline_flux_int_vel to get Jy values)')
      call hdf5_write_data(name//'/coline_flux_central',sky_galaxy%coline_shape%scentral, &
      & '[s/km] normalised central CO line flux density of inclined galaxy (multiply by coline_flux_int_vel to get Jy values)')
      call hdf5_write_data(name//'/coline_width_peak',sky_galaxy%coline_shape%wpeak, &
      & '[km/s] CO line-width between flux peaks of inclined galaxy in rest-frame velocity units')
      call hdf5_write_data(name//'/coline_width_50',sky_galaxy%coline_shape%w50, &
      & '[km/s] CO line-width at 50% of the peak flux of inclined galaxy in rest-frame velocity units')
      call hdf5_write_data(name//'/coline_width_20',sky_galaxy%coline_shape%w20, &
      & '[km/s] CO line-width at 20% of the peak flux of inclined galaxy in rest-frame velocity units')
      
      call hdf5_write_data(name//'/coline_flux_peak_eo',sky_galaxy%coline_shape%speak_eo, &
      & '[s/km] normalised peak CO line flux density of edge-on galaxy (multiply by coline_flux_int_vel to get Jy values)')
      call hdf5_write_data(name//'/coline_flux_central_eo',sky_galaxy%coline_shape%scentral_eo, &
      & '[s/km] normalised central CO line flux density of edge-on galaxy (multiply by coline_flux_int_vel to get Jy values)')
      call hdf5_write_data(name//'/coline_width_peak_eo',sky_galaxy%coline_shape%wpeak_eo, &
      & '[km/s] CO line-width between flux peaks of edge-on galaxy in rest-frame velocity units')
      call hdf5_write_data(name//'/coline_width_50_eo',sky_galaxy%coline_shape%w50_eo, &
      & '[km/s] CO line-width at 50% of the peak flux of edge-on galaxy in rest-frame velocity units')
      call hdf5_write_data(name//'/coline_width_20_eo',sky_galaxy%coline_shape%w20_eo, &
      & '[km/s] CO line-width at 20% of the peak flux of edge-on galaxy in rest-frame velocity units')
      
   end if
   
   ! Group "groups"
   if (para%make_groups==1) then
   
      name = 'groups'
      call hdf5_add_group(name)
      call hdf5_write_data(name//'/snapshot',sky_group%snapshot,'snapshot index')
      call hdf5_write_data(name//'/subvolume',sky_group%subvolume,'subvolume index')
      call hdf5_write_data(name//'/tile',sky_group%tile,'tile index in tiling array')
      call hdf5_write_data(name//'/zobs',sky_group%zobs,'redshift in observer-frame')
      call hdf5_write_data(name//'/zcmb',sky_group%zcmb,'redshift in CMB frame')
      call hdf5_write_data(name//'/zcos',sky_group%zcos,'cosmological redshift without peculiar motions')
      call hdf5_write_data(name//'/dc',sky_group%dc,'[Mpc/h] comoving distance')
      call hdf5_write_data(name//'/ra',sky_group%ra/unit%degree,'[deg] right ascension')
      call hdf5_write_data(name//'/dec',sky_group%dec/unit%degree,'[deg] declination')
      call hdf5_write_data(name//'/vpec_x',sky_group%vpec(1),'[proper km/s] x-component of peculiar velocity')
      call hdf5_write_data(name//'/vpec_y',sky_group%vpec(2),'[proper km/s] y-component of peculiar velocity')
      call hdf5_write_data(name//'/vpec_z',sky_group%vpec(3),'[proper km/s] z-component of peculiar velocity')
      call hdf5_write_data(name//'/vpec_r',sky_group%vpecrad,'[proper km/s] line-of-sight peculiar velocity')
      call hdf5_write_data(name//'/sigma_3d_all',sky_group%sigma_3d_all,&
      &'[proper km/s] 3D peculiar velocity dispersion of ALL group members, including non-detections')
      call hdf5_write_data(name//'/sigma_los_detected',sky_group%sigma_los_detected,&
      &'[proper km/s] line-of-sight peculiar velocity dispersion of detected group members')
      call hdf5_write_data(name//'/id_group_sky',sky_group%id_group_sky,'unique parent halo ID in mock sky')
      call hdf5_write_data(name//'/id_halo_sam',sky_group%id_halo_sam,'parent halo ID in SAM')
      call hdf5_write_data(name//'/mvir',sky_group%mvir,'[Msun/h] virial mass')
      call hdf5_write_data(name//'/n_galaxies_total',sky_group%group_ntot, &
      & 'total number of galaxies that live in the same group (host halo)')
      call hdf5_write_data(name//'/n_galaxies_selected',sky_group%group_nsel, &
      & 'number of galaxies that live in the same group (host halo) and are present in the mock sky')
      call hdf5_write_data(name//'/flag',sky_group%group_flag, &
      & 'group flag (0 if group complete, >0 if truncated by survey edge (+1), snapshot limit (+2), tile edge (+4))')
      
   end if
   
   ! write group "tiling"
   name = 'tiling'
   call hdf5_add_group(name)
   call hdf5_write_data(name//'/tile_id',(/(i,i=1,size(tile),1)/),'unique index of cubic tile')
   call hdf5_write_data(name//'/center_x',tile%ix(1),'x-coordinate of box-centre in units of box side length')
   call hdf5_write_data(name//'/center_y',tile%ix(2),'y-coordinate of box-centre in units of box side length')
   call hdf5_write_data(name//'/center_z',tile%ix(3),'z-coordinate of box-centre in units of box side length')
   call hdf5_write_data(name//'/dc_min',tile%dmin,'minimum comoving distance in units of box side length')
   call hdf5_write_data(name//'/dc_max',tile%dmax,'maximum comoving distance in units of box side length')
   call hdf5_write_data(name//'/rotation',tile%rotation, & 
   & 'index [1,...,6] of rotation, where 1 is the identity (negative if inversion)')
   call hdf5_write_data(name//'/translation_x',tile%translation(1), & 
   & 'x-component of translation vector in units of box side length')
   call hdf5_write_data(name//'/translation_y',tile%translation(2), & 
   & 'y-component of translation vector in units of box side length')
   call hdf5_write_data(name//'/translation_z',tile%translation(3), & 
   & 'z-component of translation vector in units of box side length')
   
   ! write group "run_info"
   name = 'run_info'
   call hdf5_add_group(name)
   call hdf5_write_data(name//'/stingray_version',version,'Version of Stingray used to produce this mock sky')
   call hdf5_write_data(name//'/shark_version',trim(shark_version),'Version of Shark used to produce this mock sky')
   call hdf5_write_data(name//'/shark_git_revision',trim(shark_git_revision),'Git revision of shark used to produce this data')
   call hdf5_write_data(name//'/shark_timestamp',trim(shark_timestamp),'Time at which this shark execution started')
   
   ! write group "sky"
   name = 'sky'
   call hdf5_add_group(name)
   call hdf5_write_data(name//'/n_galaxies',total_stats%n_galaxies,'Number of galaxies')
   call hdf5_write_data(name//'/n_groups',total_stats%n_groups,'Number of groups')
   call hdf5_write_data(name//'/n_replica_mean',real(total_stats%n_distinct,4)/total_stats%n_galaxies,&
   & 'Mean number a galaxy was replicated)')
   call hdf5_write_data(name//'/n_replica_max',total_stats%n_replica_max,'Maximum number a galaxy was replicated')
   
   ! write group "subvolume"
   if (present(subvolume_stats).and.present(isubvolume)) then
      name = 'subvolume'
      call hdf5_add_group(name)
      call hdf5_write_data(name//'/subvolume_index',isubvolume,'Index of the SAM subvolume used in this file')
      call hdf5_write_data(name//'/n_galaxies',total_stats%n_galaxies,'Number of galaxies in this subvolume')
      call hdf5_write_data(name//'/n_groups',total_stats%n_groups,'Number of groups in this subvolume')
      call hdf5_write_data(name//'/n_replica_mean',real(total_stats%n_distinct,4)/total_stats%n_galaxies,&
      & 'Mean number a galaxy was replicated in this subvolume)')
      call hdf5_write_data(name//'/n_replica_max',total_stats%n_replica_max,&
      &'Maximum number a galaxy was replicated in this subvolume')   
   end if
   
   ! write group "snapshots"
   name = 'snapshots'
   call hdf5_add_group(name)
   call hdf5_write_data(name//'/id',(/(i,i=para%snapshot_min,para%snapshot_max)/), &
   & 'snapshot number')
   call hdf5_write_data(name//'/z',snapshot%redshift, &
   & 'redshift corresponding to the cosmic time of this snapshot')
   call hdf5_write_data(name//'/dc_min',snapshot%dmin*para%box_side, &
   & '[Mpc/h] minimal comoving distance at which this snapshot is used')
   call hdf5_write_data(name//'/dc_max',snapshot%dmax*para%box_side, &
   & '[Mpc/h] maximal comoving distance at which this snapshot is used')
   call hdf5_write_data(name//'/n_tiles',snapshot%n_tiles, &
   & 'Number of tiles this snapshot has been considered for, irrespective of whether a galaxy was selected')
   
   ! close HDF5 file
   call hdf5_close()

end subroutine make_hdf5

end module module_user_routines
