module module_user_routines

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
use module_emission_lines

! custom modules
use module_hdf5

private

public :: type_sam   ! SAM class, requires procedures get_position, get_groupid, is_selected
public :: type_sky_galaxy
public :: type_sky_group
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
   & '/Users/do/Dropbox/Code/Fortran/stingray/stingray/parameters.txt'


! ==============================================================================================================
! TYPE DECLARATIONS
! ==============================================================================================================

! Here, specify the class of SAM-properties to be loaded for making the mock sky. How these properties are
! read from files is specified in the subroutine load_sam_snapshot below
! The class requires four class functions:
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
   real*4      :: mvir_hosthalo  ! [Msun/h]
   real*4      :: mvir_subhalo   ! [Msun/h]
   real*4      :: cnfw_subhalo   ! [-] concentration of NFW fit to subhalo
   real*4      :: vvir_hosthalo  ! [km/s]	virial velocity of hosthalo
   real*4      :: vvir_subhalo   ! [km/s]	virial velocity of subhalo
   real*4      :: vmax_subhalo   ! [km/s]	maximum circular velocity of subhalo
   
contains

   procedure   :: get_position      => sam_get_position     ! required function
   procedure   :: get_groupid       => sam_get_groupid      ! required function               
   procedure   :: is_group_center   => sam_is_group_center  ! required function

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
   real*4      :: vpec(3)        ! [proper km/s] 3D peculiar velocity
   real*4      :: vpecrad        ! [proper km/s] radial peculiar velocity
   integer*8   :: id_halo_sam    ! galaxy parent halo ID in the SAM
   integer*8   :: id_group_sky   ! unique group ID in the mock sky
   
   contains

   procedure   :: write_to_file  => sky_write_to_file    ! required subroutine
   
end type type_sky

type,extends(type_sky) :: type_sky_galaxy ! must exist

   integer*8   :: id_galaxy_sky           ! unique ID in the mock sky
   integer*8   :: id_galaxy_sam           ! galaxy ID in the SAM
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
   
   real*4      :: J(3)                    ! [Msun pMpc km/s] total angular momentum
   
   real*4      :: mvir_hosthalo           ! [Msun/h] mass of 1st generation halo (i.e. direct host of type 0 galaxies)
   real*4      :: mvir_subhalo            ! [Msun/h] halo mass
   
   real*4      :: rstar_disk_apparent     ! [arcsec] apparent half-mass radius of stars in the disk
   real*4      :: rstar_bulge_apparent    ! [arcsec] apparent half-mass radius of stars in the bulge
   real*4      :: rgas_disk_apparent      ! [arcsec] apparent half-mass radius of gas in the disk
   real*4      :: rgas_bulge_apparent     ! [arcsec] apparent half-mass radius of gas in the bulge
   
   real*4      :: rstar_disk_intrinsic    ! [cMpc/h] intrinsic half-mass radius of stars in the disk
   real*4      :: rstar_bulge_intrinsic   ! [cMpc/h] intrinsic half-mass radius of stars in the bulge
   real*4      :: rgas_disk_intrinsic     ! [cMpc/h] intrinsic half-mass radius of gas in the disk
   real*4      :: rgas_bulge_intrinsic    ! [cMpc/h] intrinsic half-mass radius of gas in the bulge
   
   real*4      :: vvir_hosthalo           ! [km/s]	virial velocity of hosthalo
   real*4      :: vvir_subhalo            ! [km/s]	virial velocity of subhalo
   real*4      :: vmax_subhalo            ! [km/s]	maximum circular velocity of subhalo
   real*4      :: cnfw_subhalo            ! [-] concentration of NFW fit to subhalo
   
   ! HI line
   real*4      :: hiline_flux_int         ! [W/m^2] integrated HI line flux
   real*4      :: hiline_flux_int_vel     ! [Jy km/s] velocity-integrated HI line flux
   real*4      :: hiline_flux_peak        ! [s/km] normalised peak HI line flux density of inclined galaxy (multiply by hiline_flux_int_vel to get Jy)
   real*4      :: hiline_flux_central     ! [s/km] normalised central HI line flux density of inclined galaxy
   real*4      :: hiline_width_peak       ! [km/s] line width of inclined galaxy in rest-frame velocity at peak flux
   real*4      :: hiline_width_50         ! [km/s] line width of inclined galaxy in rest-frame velocity at 50% or peak flux
   real*4      :: hiline_width_20         ! [km/s] line width of inclined galaxy in rest-frame velocity at 20% of peak flux
   real*4      :: hiline_flux_peak_eo     ! [s/km] peak HI line flux density of edge-on galaxy
   real*4      :: hiline_flux_central_eo  ! [s/km] central HI line flux density of edge-on galaxy
   real*4      :: hiline_width_peak_eo    ! [km/s] line width of edge-on galaxy in rest-frame velocity at peak flux
   real*4      :: hiline_width_50_eo      ! [km/s] line width of edge-on galaxy in rest-frame velocity at 50% or peak flux
   real*4      :: hiline_width_20_eo      ! [km/s] line width of edge-on galaxy in rest-frame velocity at 20% of peak flux
   
   ! CO (1-0) line
   real*4      :: coline_flux_int         ! [W/m^2] integrated CO(1-0) line flux
   real*4      :: coline_flux_int_vel     ! [Jy km/s] velocity-integrated CO(1-0) line flux
   real*4      :: coline_flux_peak        ! [s/km] peak CO line flux density of inclined galaxy (multiply by coline_flux_int_vel to get Jy)
   real*4      :: coline_flux_central     ! [s/km] central CO line flux density of inclined galaxy
   real*4      :: coline_width_peak       ! [km/s] line width of inclined galaxy in rest-frame velocity at peak flux
   real*4      :: coline_width_50         ! [km/s] line width of inclined galaxy in rest-frame velocity at 50% or peak flux
   real*4      :: coline_width_20         ! [km/s] line width of inclined galaxy in rest-frame velocity at 20% of peak flux
   real*4      :: coline_flux_peak_eo     ! [s/km] peak CO line flux density of edge-on galaxy
   real*4      :: coline_flux_central_eo  ! [s/km] central CO line flux density of edge-on galaxy
   real*4      :: coline_width_peak_eo    ! [km/s] line width of edge-on galaxy in rest-frame velocity at peak flux
   real*4      :: coline_width_50_eo      ! [km/s] line width of edge-on galaxy in rest-frame velocity at 50% or peak flux
   real*4      :: coline_width_20_eo      ! [km/s] line width of edge-on galaxy in rest-frame velocity at 20% of peak flux
   
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

! In order to place the objects in the mock sky, the class type_sam must have the following functions
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
! CONVERSION BETWEEN INTRINSIC AND APPARENT GALAXY PROPERTIES
! ==============================================================================================================

subroutine make_sky_object(sky_object,sam,base,groupid)

   class(type_sky),intent(out)            :: sky_object
   type(type_sam),intent(in)              :: sam
   type(type_base),intent(in)             :: base                 ! basic properties of the position of this galaxy in the sky
   integer*8,intent(in)                   :: groupid              ! unique group in sky index
   real*4                                 :: vector_rotation(3,3) ! rotation matrix for vectors
   real*4                                 :: pos(3)               ! [simulation length units] position vector of galaxy
   real*4                                 :: elos(3)              ! unit vector pointing from the observer to the object in comoving space
   
   call nil(sky_object,sam,base,groupid) ! dummy statement to avoid compiler warnings

   ! sky coordinates
   sky_object%dc  = base%dc*para%box_side  ! [Mpc/h]
   sky_object%ra  = base%ra         ! [rad]
   sky_object%dec = base%dec        ! [rad]
   
   ! make redshift, provided the galaxy position [simulation units] galaxy-velocity [km/s]
   call sph2car(sky_object%dc,sky_object%ra,sky_object%dec,pos)
   elos = pos/norm(pos) ! line-of-sight vector
   call make_redshift(pos*(para%length_unit/Mpc),sam%velocity,&
   &zobs=sky_object%zobs,zcmb=sky_object%zcmb,zcos=sky_object%zcos)
   
   ! make IDs
   sky_object%tile          = base%tile
   sky_object%snapshot      = sam%snapshot
   sky_object%subvolume     = sam%subvolume
   sky_object%id_halo_sam   = sam%id_halo
   sky_object%id_group_sky  = groupid         ! unique group id
   
   ! peculiar velocity
   vector_rotation      = tile(base%tile)%Rvector
   sky_object%vpec      = rotate(vector_rotation,sam%velocity)
   sky_object%vpecrad   = sum(sky_object%vpec*elos)
   
end subroutine make_sky_object

subroutine make_sky_galaxy(sky_galaxy,sam,base,groupid,galaxyid)

   implicit none
   class(type_sky_galaxy),intent(out)  :: sky_galaxy
   type(type_sam),intent(in)           :: sam
   type(type_base),intent(in)          :: base                 ! basic properties of the position of this galaxy in the sky
   integer*8,intent(in)                :: groupid              ! unique group in sky index
   integer*8,intent(in)                :: galaxyid             ! (only for galaxy types) unique galaxy in sky index
   real*4                              :: pseudo_rotation(3,3) ! rotation matrix for pseudo-vectors (axial vectors)
   real*4                              :: pos(3)               ! [simulation length units] position vector of galaxy
   real*4                              :: dl                   ! [simulation length units] luminosity distance to observer
   real*4                              :: elos(3)              ! unit vector pointing from the observer to the object in comoving space
   real*4                              :: mstars,mHI,mH2,LCO
   
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
   
   ! intrinsic angular momentum
   pseudo_rotation   = tile(base%tile)%Rpseudo
   sky_galaxy%J      = rotate(pseudo_rotation,sam%J)
      
   ! intrinsic radii
   sky_galaxy%rstar_disk_intrinsic = sam%rstar_disk ! [cMpc/h]
   sky_galaxy%rstar_bulge_intrinsic = sam%rstar_bulge ! [cMpc/h]
   sky_galaxy%rgas_disk_intrinsic = sam%rgas_disk ! [cMpc/h]
   sky_galaxy%rgas_bulge_intrinsic = sam%rgas_bulge ! [cMpc/h]
   
      
   ! APPARENT PROPERTIES
   
   ! inclination and position angle
   call sph2car(sky_galaxy%dc,sky_galaxy%ra,sky_galaxy%dec,pos)
   elos = pos/norm(pos) ! position vector
   call make_inclination_and_pa(pos,sky_galaxy%J,inclination=sky_galaxy%inclination,pa=sky_galaxy%pa)
   
   ! apparent radii (note division by comoving distance, because intrinsic radii given in comoving units)
   sky_galaxy%rstar_disk_apparent = sam%rstar_disk/sky_galaxy%dc/degree*3600.0 ! [arcsec]
   sky_galaxy%rstar_bulge_apparent = sam%rstar_bulge/sky_galaxy%dc/degree*3600.0 ! [arcsec]
   sky_galaxy%rgas_disk_apparent = sam%rgas_disk/sky_galaxy%dc/degree*3600.0 ! [arcsec]
   sky_galaxy%rgas_bulge_apparent = sam%rgas_bulge/sky_galaxy%dc/degree*3600.0 ! [arcsec]
   
   ! rough generic optical magnitude
   dl = sky_galaxy%dc*(1+sky_galaxy%zobs) ! [Mpc/h]
   mstars = (sam%mstars_disk+sam%mstars_bulge)/para%h ! [Msun]
   sky_galaxy%mag = convert_absmag2appmag(convert_stellarmass2absmag(mstars,1.0),dl/para%h)
      
   ! cold gas emission lines
   mHI = (sam%matom_disk+sam%matom_bulge)/1.35 ! [Msun/h] HI mass  
   mH2 = ((sam%mgas_disk+sam%mgas_bulge)-(sam%matom_disk+sam%matom_bulge))/1.35 ! [Msun/h] H2 mass
   sky_galaxy%hiline_flux_int = real(convert_luminosity2flux(real(mHI/para%h,8)*real(L2MHI,8)*Lsun,dl/para%h),4) ! [W/m^2]
   sky_galaxy%hiline_flux_int_vel = convert_intflux2velintflux(sky_galaxy%hiline_flux_int,0.21106114,sky_galaxy%zobs)
   LCO = mH2/para%h/(313*X_CO) ! [Jy km/s Mpc^2] CO(1-0) luminosity (note this equation has a J^2 dependence)
   sky_galaxy%coline_flux_int = LCO/(4.0*pi*(dl/para%h)**2)/0.00260075761e23 ! [W/m^2] integrated flux
   sky_galaxy%coline_flux_int_vel = convert_intflux2velintflux(sky_galaxy%hiline_flux_int,0.00260075761,sky_galaxy%zobs)
   if (para%line_parameters==1) call make_line_profiles
   
contains
   
   subroutine make_line_profiles
   
      implicit none
      type(type_lineinput)                :: line_input
      type(type_line)                     :: line(4)
      real*4                              :: a,c_halo,cMpch2pkpc,rvir
      real,parameter                      :: G = 4.3022682e-6     ! [(km/s)^2 kpc/Msun] gravitational constant
      
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
      if ((line_input%mHI>0.0).or.(line_input%mH2>0.0)) call line_profile(line_input,para%h,line)
   
      ! save line parameters
      if (line_input%mHI>0.0) then
         sky_galaxy%hiline_flux_peak_eo = line(1)%speak
         sky_galaxy%hiline_flux_central_eo = line(1)%scentral
         sky_galaxy%hiline_width_peak_eo = line(1)%wpeak
         sky_galaxy%hiline_width_50_eo = line(1)%w50
         sky_galaxy%hiline_width_20_eo = line(1)%w20
      
         sky_galaxy%hiline_flux_peak = line(2)%speak
         sky_galaxy%hiline_flux_central = line(2)%scentral
         sky_galaxy%hiline_width_peak = line(2)%wpeak
         sky_galaxy%hiline_width_50 = line(2)%w50
         sky_galaxy%hiline_width_20 = line(2)%w20
      else
         sky_galaxy%hiline_flux_peak_eo = 0
         sky_galaxy%hiline_flux_central_eo = 0
         sky_galaxy%hiline_width_peak_eo = 0
         sky_galaxy%hiline_width_50_eo = 0
         sky_galaxy%hiline_width_20_eo = 0
      
         sky_galaxy%hiline_flux_peak = 0
         sky_galaxy%hiline_flux_central = 0
         sky_galaxy%hiline_width_peak = 0
         sky_galaxy%hiline_width_50 = 0
         sky_galaxy%hiline_width_20 = 0
      end if
   
      if (line_input%mH2>0.0) then
         sky_galaxy%coline_flux_peak_eo = line(3)%speak
         sky_galaxy%coline_flux_central_eo = line(3)%scentral
         sky_galaxy%coline_width_peak_eo = line(3)%wpeak
         sky_galaxy%coline_width_50_eo = line(3)%w50
         sky_galaxy%coline_width_20_eo = line(3)%w20
      
         sky_galaxy%coline_flux_peak = line(4)%speak
         sky_galaxy%coline_flux_central = line(4)%scentral
         sky_galaxy%coline_width_peak = line(4)%wpeak
         sky_galaxy%coline_width_50 = line(4)%w50
         sky_galaxy%coline_width_20 = line(4)%w20
      else
         sky_galaxy%coline_flux_peak_eo = 0
         sky_galaxy%coline_flux_central_eo = 0
         sky_galaxy%coline_width_peak_eo = 0
         sky_galaxy%coline_width_50_eo = 0
         sky_galaxy%coline_width_20_eo = 0
      
         sky_galaxy%coline_flux_peak = 0
         sky_galaxy%coline_flux_central = 0
         sky_galaxy%coline_width_peak = 0
         sky_galaxy%coline_width_50 = 0
         sky_galaxy%coline_width_20 = 0
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
   logical,intent(in)                  :: selected(:)
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
   call out('File of automatic parameters: '//trim(filename))
   call hdf5_open(filename)
   call hdf5_read_data('/run_info/tot_n_subvolumes',nsub)
   para%subvolume_min = 0
   para%subvolume_max = nsub-1
   call hdf5_read_data('/cosmology/h',para%h)
   call hdf5_read_data('/cosmology/omega_l',para%omega_l)
   call hdf5_read_data('/cosmology/omega_m',para%omega_m)
   call hdf5_read_data('/cosmology/omega_b',para%omega_b)
   call hdf5_read_data('/run_info/lbox',para%box_side)
   para%length_unit = Mpc/para%h
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
   call hdf5_read_data(g//'vvir_subhalo',sam%vvir_subhalo)
   call hdf5_read_data(g//'vmax_subhalo',sam%vmax_subhalo)
   call hdf5_read_data(g//'vvir_hosthalo',sam%vvir_hosthalo)
   
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
   logical,parameter                   :: show_test_sum = .true.
   character(len=255)                  :: filename_bin
   character(len=255)                  :: filename_hdf5
   type(type_sky_galaxy),allocatable   :: sky_galaxy(:)
   type(type_sky_group),allocatable    :: sky_group(:)
   integer*8                           :: n,i,n_galaxies,n_groups
   character(len=255)                  :: name,str
   character(len=255)                  :: filename
   character(len=255)                  :: shark_version,shark_git_revision,shark_timestamp
   real*8                              :: test(0:10),expected_sum
   real*4                              :: n_replica_mean
   integer*4                           :: n_replica_max
   
   call tic
   call out('CONVERT MOCK SKY FROM BINARY TO HDF5')
   
   ! load auxilary data
   call load_parameters
   call load_tile_list
   call load_snapshot_list
   
   ! load auxilary data from shark output
   write(filename,'(A,I0,A)') trim(para%path_input),para%snapshot_min,'/0/galaxies.hdf5'
   call hdf5_open(filename)
   call hdf5_read_data('/run_info/shark_version',shark_version)
   call hdf5_read_data('/run_info/shark_git_revision',shark_git_revision)
   call hdf5_read_data('/run_info/timestamp',shark_timestamp)
   call hdf5_close()
   
   ! create HDF5 file
   filename_hdf5 = trim(para%path_output)//'mocksky.hdf5'
   call hdf5_create(filename_hdf5)
   
   ! open HDF5 file
   call hdf5_open(filename_hdf5,.true.)
   
   ! Group "Parameters"
   call hdf5_add_group('parameters')
   call hdf5_write_data('parameters/survey',para%survey,'name of simulated survey')
   call hdf5_write_data('parameters/path_output',para%path_output)
   call hdf5_write_data('parameters/path_input',para%path_input)
   call hdf5_write_data('parameters/box_side',para%box_side, &
   & '[length_unit] comoving side-length of simulation box in multiples of length_unit')
   call hdf5_write_data('parameters/length_unit',para%length_unit,'[m] SI-value of comoving length unit')
   call hdf5_write_data('parameters/snapshot_min',para%snapshot_min, &
   & 'index of earliest snapshot used for the mock sky')
   call hdf5_write_data('parameters/snapshot_max',para%snapshot_max, &
   & 'index of latest snapshot used for the mock sky')
   call hdf5_write_data('parameters/subvolume_min',para%subvolume_min, &
   & 'index of first subvolume used for the mock sky')
   call hdf5_write_data('parameters/subvolume_max',para%subvolume_max, &
   & 'index of last subvolume used for the mock sky')
   call hdf5_write_data('parameters/h',para%h,'[-] Hubble parameter H0=h*100 km/s/Mpc')
   call hdf5_write_data('parameters/omega_l',para%omega_l, &
   & 'energy density of dark energy relative to closure density')
   call hdf5_write_data('parameters/omega_m',para%omega_m, &
   & 'energy density of all matter relative to closure density')
   call hdf5_write_data('parameters/omega_b',para%omega_b, &
   & 'energy density of baryonic matter relative to closure density')
   call hdf5_write_data('parameters/dc_min',para%dc_min, &
   & '[length_unit] minimal comoving distance of mock sky')
   call hdf5_write_data('parameters/dc_max',para%dc_max, &
   & '[length_unit] maximal comoving distance of mock sky')
   call hdf5_write_data('parameters/ra_min',para%ra_min/degree, &
   & '[deg] min right ascension of mock sky')
   call hdf5_write_data('parameters/ra_max',para%ra_max/degree, &
   & '[deg] max right ascension of mock sky')
   call hdf5_write_data('parameters/dec_min',para%dec_min/degree, &
   & '[deg] min declination of mock sky')
   call hdf5_write_data('parameters/dec_max',para%dec_max/degree, &
   & '[deg] max declination of mock sky')
   call hdf5_write_data('parameters/zaxis_ra',para%zaxis_ra/degree, &
   & '[deg] RA coordinate of the SAM z-axis in spherical survey coordinates')
   call hdf5_write_data('parameters/zaxis_dec',para%zaxis_dec/degree, &
   & '[deg] Dec coordinate the SAM z-axis in spherical survey coordinates')
   call hdf5_write_data('parameters/xy_angle',para%xy_angle/degree, &
   & '[deg] Rotation of the SAM (x,y)-plane on the sky')
   call hdf5_write_data('parameters/seed',para%seed, &
   & 'seed for the random number generator of symmetry operations')
   call hdf5_write_data('parameters/translate',para%translate, &
   & 'logical flag specifying if random translations are applied (0=false, 1=true)')
   call hdf5_write_data('parameters/rotate',para%rotate, &
   & 'logical flag specifying if random rotations are applied (0=false, 1=true)')
   call hdf5_write_data('parameters/invert',para%invert, &
   & 'logical flag specifying if random inversions are applied (0=false, 1=true)')
   call hdf5_write_data('parameters/velocity_norm',para%velocity_norm, &
   & '[km/s] observer velocity relative to CMB rest-frame')
   call hdf5_write_data('parameters/velocity_ra',para%velocity_ra, &
   & '[deg] RA coordinate to which the observer is moving relative to CMB rest-frame')
   call hdf5_write_data('parameters/velocity_dec',para%velocity_dec, &
   & '[deg] Dec coordinate to which the observer is moving relative to CMB rest-frame')
   call hdf5_write_data('parameters/search_angle',para%search_angle, &
   & '[deg] typical angle in which overlaps between survey volume and tiling grid are searched')
   call hdf5_write_data('parameters/volume_search_level',para%volume_search_level, &
   & 'specifies the number of search points (2^#)^3 inside each tile')
   call hdf5_write_data('parameters/line_parameters',para%line_parameters, &
   & 'logical flag specifying if global emission line parameters are saved (0=false, 1=true)')
   call hdf5_write_data('parameters/skyrotation',para%sky_rotation, &
   & 'Rotation matrix to map xyz-coordinates of the tiling structure onto sky coordinates.')
   
   ! Group "galaxies"
   name = 'galaxies'
   filename_bin = trim(para%path_output)//'mocksky_galaxies.bin'
   open(1,file=trim(filename_bin),action='read',form='unformatted',access='stream')
   read(1) n,n_replica_mean,n_replica_max
   n_galaxies = n
   allocate(sky_galaxy(n))
   read(1) sky_galaxy
   close(1)
   call hdf5_add_group(trim(name))
   call hdf5_write_data(trim(name)//'/snapshot',sky_galaxy%snapshot,'snapshot index')
   call hdf5_write_data(trim(name)//'/subvolume',sky_galaxy%subvolume,'subvolume index')
   call hdf5_write_data(trim(name)//'/tile',sky_galaxy%tile,'tile index in tiling array')
   call hdf5_write_data(trim(name)//'/zobs',sky_galaxy%zobs,'redshift in observer-frame')
   call hdf5_write_data(trim(name)//'/zcmb',sky_galaxy%zcmb,'redshift in CMB frame')
   call hdf5_write_data(trim(name)//'/zcos',sky_galaxy%zcos,'cosmological redshift without peculiar motions')
   call hdf5_write_data(trim(name)//'/dc',sky_galaxy%dc,'[Mpc/h] comoving distance')
   call hdf5_write_data(trim(name)//'/ra',sky_galaxy%ra/degree,'[deg] right ascension')
   call hdf5_write_data(trim(name)//'/dec',sky_galaxy%dec/degree,'[deg] declination')
   call hdf5_write_data(trim(name)//'/id_galaxy_sky',sky_galaxy%id_galaxy_sky,'unique galaxy ID in mock sky')
   call hdf5_write_data(trim(name)//'/id_group_sky',sky_galaxy%id_group_sky,'unique group ID if galaxy is in a group, -1 otherwise')
   call hdf5_write_data(trim(name)//'/id_galaxy_sam',sky_galaxy%id_galaxy_sam,'galaxy ID in SAM')
   call hdf5_write_data(trim(name)//'/id_halo_sam',sky_galaxy%id_halo_sam,'host halo ID in SAM')
   call hdf5_write_data(trim(name)//'/type',sky_galaxy%typ,'galaxy type (0=central, 1=satellite in halo, 2=orphan)')
   call hdf5_write_data(trim(name)//'/inclination',sky_galaxy%inclination/degree, &
   & '[deg] inclination = angle between line-of-sight and spin axis')
   call hdf5_write_data(trim(name)//'/pa',sky_galaxy%pa/degree,'[deg] position angle from north to east')
   call hdf5_write_data(trim(name)//'/mag',sky_galaxy%mag, &
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
   call hdf5_write_data(trim(name)//'/mvir_hosthalo',sky_galaxy%mvir_hosthalo,'[Msun/h] host halo mass')
   call hdf5_write_data(trim(name)//'/mvir_subhalo',sky_galaxy%mvir_subhalo,'[Msun/h] subhalo mass')
   call hdf5_write_data(trim(name)//'/rstar_disk_apparent',sky_galaxy%rstar_disk_apparent,&
   &'[arcsec] apparent half-mass radius of stellar disk')
   call hdf5_write_data(trim(name)//'/rstar_bulge_apparent',sky_galaxy%rstar_bulge_apparent,&
   &'[arcsec] apparent half-mass radius of stellar bulge')
   call hdf5_write_data(trim(name)//'/rgas_disk_apparent',sky_galaxy%rgas_disk_apparent,&
   &'[arcsec] apparent half-mass radius of gas disk')
   call hdf5_write_data(trim(name)//'/rgas_bulge_apparent',sky_galaxy%rgas_bulge_apparent,&
   &'[arcsec] apparent half-mass radius of gas bulge')
   call hdf5_write_data(trim(name)//'/rstar_disk_intrinsic',sky_galaxy%rstar_disk_intrinsic,&
   &'[cMpc/h] intrinsic half-mass radius of stellar disk')
   call hdf5_write_data(trim(name)//'/rstar_bulge_intrinsic',sky_galaxy%rstar_bulge_intrinsic,&
   &'[cMpc/h] intrinsic half-mass radius of stellar bulge')
   call hdf5_write_data(trim(name)//'/rgas_disk_intrinsic',sky_galaxy%rgas_disk_intrinsic,&
   &'[cMpc/h] intrinsic half-mass radius of gas disk')
   call hdf5_write_data(trim(name)//'/rgas_bulge_intrinsic',sky_galaxy%rgas_bulge_intrinsic,&
   &'[cMpc/h] intrinsic half-mass radius of gas bulge')
   call hdf5_write_data(trim(name)//'/hiline_flux_int',sky_galaxy%hiline_flux_int,'[W/m^2] integrated HI line flux')
   call hdf5_write_data(trim(name)//'/hiline_flux_int_vel',sky_galaxy%hiline_flux_int_vel, &
   & '[Jy km/s] velocity-integrated HI line flux')
   call hdf5_write_data(trim(name)//'/coline_flux_int',sky_galaxy%coline_flux_int,'[W/m^2] integrated CO(1-0) line flux')
   call hdf5_write_data(trim(name)//'/coline_flux_int_vel',sky_galaxy%coline_flux_int_vel, &
   & '[Jy km/s] velocity-integrated CO(1-0) line flux')
   
   if (para%line_parameters==1) then
   
      call hdf5_write_data(trim(name)//'/hiline_flux_peak',sky_galaxy%hiline_flux_peak, &
      & '[s/km] normalised peak HI line flux density of inclined galaxy (multiply by hiline_flux_int_vel to get Jy values)')
      call hdf5_write_data(trim(name)//'/hiline_flux_central',sky_galaxy%hiline_flux_central, &
      & '[s/km] normalised central HI line flux density of inclined galaxy (multiply by hiline_flux_int_vel to get Jy values)')
      call hdf5_write_data(trim(name)//'/hiline_width_peak',sky_galaxy%hiline_width_peak, &
      & '[km/s] HI line-width between flux peaks of inclined galaxy in rest-frame velocity units')
      call hdf5_write_data(trim(name)//'/hiline_width_50',sky_galaxy%hiline_width_50, &
      & '[km/s] HI line-width at 50% of the peak flux of inclined galaxy in rest-frame velocity units')
      call hdf5_write_data(trim(name)//'/hiline_width_20',sky_galaxy%hiline_width_20, &
      & '[km/s] HI line-width at 20% of the peak flux of inclined galaxy in rest-frame velocity units')
      
      call hdf5_write_data(trim(name)//'/hiline_flux_peak_eo',sky_galaxy%hiline_flux_peak_eo, &
      & '[s/km] normalised peak HI line flux density of edge-on galaxy (multiply by hiline_flux_int_vel to get Jy values)')
      call hdf5_write_data(trim(name)//'/hiline_flux_central_eo',sky_galaxy%hiline_flux_central_eo, &
      & '[s/km] normalised central HI line flux density of edge-on galaxy (multiply by hiline_flux_int_vel to get Jy values)')
      call hdf5_write_data(trim(name)//'/hiline_width_peak_eo',sky_galaxy%hiline_width_peak_eo, &
      & '[km/s] HI line-width between flux peaks of edge-on galaxy in rest-frame velocity units')
      call hdf5_write_data(trim(name)//'/hiline_width_50_eo',sky_galaxy%hiline_width_50_eo, &
      & '[km/s] HI line-width at 50% of the peak flux of edge-on galaxy in rest-frame velocity units')
      call hdf5_write_data(trim(name)//'/hiline_width_20_eo',sky_galaxy%hiline_width_20_eo, &
      & '[km/s] HI line-width at 20% of the peak flux of edge-on galaxy in rest-frame velocity units')
      
      call hdf5_write_data(trim(name)//'/coline_flux_peak',sky_galaxy%coline_flux_peak, &
      & '[s/km] normalised peak CO line flux density of inclined galaxy (multiply by coline_flux_int_vel to get Jy values)')
      call hdf5_write_data(trim(name)//'/coline_flux_central',sky_galaxy%coline_flux_central, &
      & '[s/km] normalised central CO line flux density of inclined galaxy (multiply by coline_flux_int_vel to get Jy values)')
      call hdf5_write_data(trim(name)//'/coline_width_peak',sky_galaxy%coline_width_peak, &
      & '[km/s] CO line-width between flux peaks of inclined galaxy in rest-frame velocity units')
      call hdf5_write_data(trim(name)//'/coline_width_50',sky_galaxy%coline_width_50, &
      & '[km/s] CO line-width at 50% of the peak flux of inclined galaxy in rest-frame velocity units')
      call hdf5_write_data(trim(name)//'/coline_width_20',sky_galaxy%coline_width_20, &
      & '[km/s] CO line-width at 20% of the peak flux of inclined galaxy in rest-frame velocity units')
      
      call hdf5_write_data(trim(name)//'/coline_flux_peak_eo',sky_galaxy%coline_flux_peak_eo, &
      & '[s/km] normalised peak CO line flux density of edge-on galaxy (multiply by coline_flux_int_vel to get Jy values)')
      call hdf5_write_data(trim(name)//'/coline_flux_central_eo',sky_galaxy%coline_flux_central_eo, &
      & '[s/km] normalised central CO line flux density of edge-on galaxy (multiply by coline_flux_int_vel to get Jy values)')
      call hdf5_write_data(trim(name)//'/coline_width_peak_eo',sky_galaxy%coline_width_peak_eo, &
      & '[km/s] CO line-width between flux peaks of edge-on galaxy in rest-frame velocity units')
      call hdf5_write_data(trim(name)//'/coline_width_50_eo',sky_galaxy%coline_width_50_eo, &
      & '[km/s] CO line-width at 50% of the peak flux of edge-on galaxy in rest-frame velocity units')
      call hdf5_write_data(trim(name)//'/coline_width_20_eo',sky_galaxy%coline_width_20_eo, &
      & '[km/s] CO line-width at 20% of the peak flux of edge-on galaxy in rest-frame velocity units')
      
      test(0) = sum(sky_galaxy%hiline_flux_peak*sky_galaxy%hiline_width_50_eo)+ &
              & sum(sky_galaxy%coline_flux_central_eo*sky_galaxy%coline_width_20)
   
   else
   
      test(0) = 0
      
   end if
   
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
   n_groups = n
   allocate(sky_group(n))
   read(1) sky_group
   close(1)
   call hdf5_add_group(trim(name))
   call hdf5_write_data(trim(name)//'/snapshot',sky_group%snapshot,'snapshot index')
   call hdf5_write_data(trim(name)//'/subvolume',sky_group%subvolume,'subvolume index')
   call hdf5_write_data(trim(name)//'/tile',sky_group%tile,'tile index in tiling array')
   call hdf5_write_data(trim(name)//'/zobs',sky_group%zobs,'redshift in observer-frame')
   call hdf5_write_data(trim(name)//'/zcmb',sky_group%zcmb,'redshift in CMB frame')
   call hdf5_write_data(trim(name)//'/zcos',sky_group%zcos,'cosmological redshift without peculiar motions')
   call hdf5_write_data(trim(name)//'/dc',sky_group%dc,'[Mpc/h] comoving distance')
   call hdf5_write_data(trim(name)//'/ra',sky_group%ra/degree,'[deg] right ascension')
   call hdf5_write_data(trim(name)//'/dec',sky_group%dec/degree,'[deg] declination')
   call hdf5_write_data(trim(name)//'/vpec_x',sky_group%vpec(1),'[proper km/s] x-component of peculiar velocity')
   call hdf5_write_data(trim(name)//'/vpec_y',sky_group%vpec(2),'[proper km/s] y-component of peculiar velocity')
   call hdf5_write_data(trim(name)//'/vpec_z',sky_group%vpec(3),'[proper km/s] z-component of peculiar velocity')
   call hdf5_write_data(trim(name)//'/vpec_r',sky_group%vpecrad,'[proper km/s] line-of-sight peculiar velocity')
   call hdf5_write_data(trim(name)//'/sigma_3d_all',sky_group%sigma_3d_all,&
   &'[proper km/s] 3D peculiar velocity dispersion of ALL group members, including non-detections')
   call hdf5_write_data(trim(name)//'/sigma_los_detected',sky_group%sigma_los_detected,&
   &'[proper km/s] line-of-sight peculiar velocity dispersion of detected group members')
   call hdf5_write_data(trim(name)//'/id_group_sky',sky_group%id_group_sky,'unique parent halo ID in mock sky')
   call hdf5_write_data(trim(name)//'/id_halo_sam',sky_group%id_halo_sam,'parent halo ID in SAM')
   call hdf5_write_data(trim(name)//'/mvir',sky_group%mvir,'[Msun/h] virial mass')
   call hdf5_write_data(trim(name)//'/n_galaxies_total',sky_group%group_ntot, &
   & 'total number of galaxies that live in the same group (host halo)')
   call hdf5_write_data(trim(name)//'/n_galaxies_selected',sky_group%group_nsel, &
   & 'number of galaxies that live in the same group (host halo) and are present in the mock sky')
   call hdf5_write_data(trim(name)//'/flag',sky_group%group_flag, &
   & 'group flag (0 if group complete, >0 if truncated by survey edge (+1), snapshot limit (+2), tile edge (+4))')
   test(6) = sum(sky_group%tile)
   test(7) = sum(sky_group%zcmb)
   test(2) = sum(sky_group%group_flag)+sum(sky_group%group_flag**2)
   deallocate(sky_group)
   
   ! Group "Tiling"
   call hdf5_add_group('tiling')
   call hdf5_write_data('tiling/tile_id',(/(i,i=1,size(tile),1)/),'unique index of cubic tile')
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
   call hdf5_write_data('run_info/stingray_version',version,'Version of Stingray used to produce this mock sky')
   call hdf5_write_data('run_info/shark_version',trim(shark_version),'Version of Shark used to produce this mock sky')
   call hdf5_write_data('run_info/shark_git_revision',trim(shark_git_revision),'Git revision of shark used to produce this data')
   call hdf5_write_data('run_info/shark_timestamp',trim(shark_timestamp),'Time at which this shark execution started')
   call hdf5_write_data('run_info/n_replica_mean',n_replica_mean,'Mean number a galaxy was replicated in the mock sky)')
   call hdf5_write_data('run_info/n_replica_max',n_replica_max,'Maximum number a galaxy was replicated in the mock sky')
   call hdf5_write_data('run_info/n_galaxies',n_galaxies,'Total number of galaxies in the mock sky')
   call hdf5_write_data('run_info/n_groups',n_groups,'Total number of groups in the mock sky')
   
   ! Group "Snapshots"
   call hdf5_add_group('snapshots')
   call hdf5_write_data('snapshots/id',(/(i,i=para%snapshot_min,para%snapshot_max,1)/), &
   & 'snapshot number')
   call hdf5_write_data('snapshots/z',snapshot%redshift, &
   & 'redshift corresponding to the cosmic time of this snapshot')
   call hdf5_write_data('snapshots/dc_min',snapshot%dmin*para%box_side, &
   & '[Mpc/h] minimal comoving distance at which this snapshot is used')
   call hdf5_write_data('snapshots/dc_max',snapshot%dmax*para%box_side, &
   & '[Mpc/h] maximal comoving distance at which this snapshot is used')
   call hdf5_write_data('snapshots/n_tiles',snapshot%n_tiles, &
   & 'Number of tiles this snapshot has been considered for, irrespective of whether a galaxy was selected')
   
   ! close HDF5 file
   call hdf5_close()
   
   ! check test sum
   if (show_test_sum) then
      expected_sum = 831651.793945312
      if (abs(sum(test)/expected_sum-1)>1e-5) then
         write(str,'(F20.9)') sum(test)
         call out('Test sum FAILED, found = '//trim(str))
      else
         call out('Test sum OK.')
      end if
   end if
   
   call toc

end subroutine make_hdf5
   
end subroutine custom_routines

end module module_user_routines