module module_user

use module_constants
use module_types
use module_system
use module_parameters
use module_cosmology

! ==============================================================================================================
! VARIABLE TYPES
! ==============================================================================================================

! specify the galaxy properties output by the semi-analytic model
! these are mostly intrinsic galaxy properties
! they must be stored as fortran binary file exactly in this order

type type_galaxy_sam

   integer*8   :: id             ! unique galaxy ID
   integer*8   :: haloid         ! unique ID of parent halo
   real*4      :: position(3)    ! [Mpc/h] position of galaxy centre in simulation box
   real*4      :: magR           ! absolute r-band magnitude
   real*4      :: v(3)           ! [proper km/s] peculiar velocity 
   real*4      :: j(3)           ! [proper kpc km/s] specific angular momentum
   real*4      :: MHI            ! [Msun] HI mass

end type type_galaxy_sam

! specify the galaxy properties in the mock cone
! these are mostly apparent galaxy properties

type type_galaxy_cone

   integer*8   :: id             ! unique galaxy ID
   integer*4   :: snapshot       ! snapshot
   real*4      :: z              ! apparent redshift
   real*4      :: dc             ! [Mpc] comoving distance
   real*4      :: ra             ! [rad] right ascension
   real*4      :: decl           ! [rad] declination
   real*4      :: magR           ! apparent r-band magnitude
   real*4      :: v(3)           ! [proper km/s] peculiar velocity
   real*4      :: j(3)           ! [proper kpc km/s] specific angular momentum
   real*8      :: SHI            ! [W/m^2] integrated HI line flux

end type type_galaxy_cone

! all galaxy properties ine one type; do not modify this type

type type_galaxy_all
   
   type(type_galaxy_base)  :: base
   type(type_galaxy_sam)   :: sam
   type(type_galaxy_cone)  :: cone

end type type_galaxy_all

contains

! mapping from type_galaxy_sam onto three basic properties of type_galaxy_base

function extract_base(sam) result(base)
   
   implicit none
   type(type_galaxy_sam)   :: sam
   type(type_galaxy_base)  :: base
   
   base%id        = sam%id
   base%groupid   = sam%haloid
   base%xbox      = sam%position

end function extract_base


! ==============================================================================================================
! SELECTION OF GALAXIES IN THE MOCK OBSERVING CONE FILE
! ==============================================================================================================

! selection function applied to SAM-properties, *before* producing intrinsic cone

logical function pre_selection(sam)

   implicit none
   type(type_galaxy_sam)   :: sam
   
   pre_selection = sam%id>1e5

end function pre_selection

! selection function applied to SAM-properties, when converting the intrinsic cone into the apparent cone

logical function post_selection(cone)

   implicit none
   type(type_galaxy_cone)  :: cone
   
   post_selection = cone%id>1e5

end function post_selection


! ==============================================================================================================
! IO routines
! ==============================================================================================================

! write mock-cone galaxy into binary file
! this function allows the user to select only certain properties, if needed, and/or add intrinsic properties

subroutine write_galaxy(galaxy)

   ! choose which variables of the structure 'galaxy' to save. this structure has the substructures
   ! galaxy%base  = basic cone geometry properties
   ! galaxy%sam   = SAM output, mainly intrinsic galaxy properties
   ! galaxy%cone  = apparent galaxy properties in the cone

   implicit none
   type(type_galaxy_all),intent(in) :: galaxy
   write(1) galaxy%cone
   
end subroutine write_galaxy

! load SAM snapshot file

subroutine load_sam_snapshot(index,sam)

   ! variable declaration
   implicit none
   integer*4,intent(in)                            :: index
   type(type_galaxy_sam),allocatable,intent(inout) :: sam(:)
   character(len=255)                              :: filename
   integer*8                                       :: i,n
   integer*4                                       :: bytesgal
   integer*8                                       :: bytestot
   
   ! determine number of galaxies from file size
   bytesgal = bytes_per_galaxy_sam()
   filename = trim(snapshot_filename(index))
   inquire(file=trim(filename), size=bytestot)
   if (modulo(bytestot,bytesgal).ne.0) then
      call out('ERROR: Size of snapshot file inconsistent with type_galaxy_sam.')
      stop
   end if
   n = bytestot/bytesgal
   
   ! allocate
   if (allocated(sam)) deallocate(sam)
   allocate(sam(n))
   
   ! read file
   open(1,file=trim(filename),action='read',form='unformatted',status='old',access='stream')
   do i = 1,n
      read(1) sam(i)
   end do
   close(1)
   
   contains
   
   integer*4 function bytes_per_galaxy_sam()
      implicit none
      character(len=255)      :: filename
      type(type_galaxy_sam)   :: sam
      filename = trim(para%path_output)//'.tmpsizeof'
      open(1,file=trim(filename),action='write',form="unformatted",status='replace',access='stream')
      write(1) sam
      close(1)
      inquire(file=trim(filename), size=bytes_per_galaxy_sam)
   end function bytes_per_galaxy_sam
   
end subroutine load_sam_snapshot


! ==============================================================================================================
! CONVERSION BETWEEN INTRINSIC AND APPARENT GALAXY PROPERTIES
! ==============================================================================================================

function convert_properties(base,sam) result(cone)

   implicit none
   type(type_galaxy_base),intent(in)      :: base
   type(type_galaxy_sam),intent(in)       :: sam      ! intrinsic galaxy properties from SAM
   type(type_galaxy_cone)                 :: cone     ! apparent galaxy properties
   real*4                                 :: zcos     ! cosmological redshift due to Hubble flow
   real*4                                 :: zpec     ! redshift due to the peculiar motion of object relative to Hubble flow
   real*4                                 :: zobs     ! redshift due to the peculiar motion of the observer relative to the Hubble flow
   real*4                                 :: z        ! total redshift
   real*4                                 :: dc       ! [simulation length units] comoving distance to observer
   real*4                                 :: dl       ! [simulation length units] luminosity distance to observer
   real*4                                 :: da       ! [simulation length units] angular diameter distance to observer
   real*4                                 :: norm     ! [box side-length] comoving distance to observer
   real*4                                 :: elos(3)  ! unit vector pointing from the observer to the object in comoving space
   
   ! compute basic redshift and distance
   norm = sqrt(sum(base%xcone**2))
   elos = base%xcone/norm
   dc = norm*para%L*(para%length_unit/Mpc) ! [Mpc]
   zcos = dc_to_redshift(dc)
   zpec = min(0.1,max(-0.1,sum(sam%v*base%xcone)/c*1e3))          ! limited to 0.1 to avoid relativistic regime, ok for all practical purposes
   zobs = min(0.1,max(-0.1,-sum(para%velocity*base%xcone)/c*1e3)) ! limited to 0.1 to avoid relativistic regime, ok for all practical purposes
   z = (1+zcos)*(1+zpec)*(1+zobs)-1 ! following Davis & Scrimgeour 2014
   dl = dc*(1+z)
   da = dc/(1+z)
   
   ! copy basic variables
   cone%id        = sam%id
   cone%snapshot  = base%snapshot
   
   ! make geometric properties
   cone%z      = z
   cone%dc     = norm*para%L ! [simulation units] comoving distance
   cone%ra     = atan2(base%xcone(1),base%xcone(3))
   cone%decl   = asin(base%xcone(2)/norm)
   
   ! convert intrinsic to apparent properties
   cone%magR   = convert_abs2appmag(sam%magR,dl)
   cone%v      = convert_vector(sam%v,base%rotation)
   cone%j      = convert_pseudovector(sam%j,base%rotation)
   cone%SHI    = convert_luminosity2flux(real(sam%MHI,8)*real(LMratioHI,8)*Lsun,dl)
   
   contains
   
   ! Hereafter follows the library of basic conversion functions.
   ! Make sure to check these functions before writing your own.

   function convert_luminosity2flux(L,dl) result(S)

      implicit none
      real*8,intent(in) :: L     ! [W] Luminosity
      real*4,intent(in) :: dl    ! [Mpc] luminosity distance
      real*8            :: S     ! [W/m^2] Flux
   
      S = L/real(dl,8)**2/ASphereMpc
   
   end function
   
   function convert_abs2appmag(M,dl) result(mag)

      implicit none
      real*4,intent(in) :: M     ! absolute magnitude
      real*4,intent(in) :: dl    ! [Mpc] luminosity distance
      real*4            :: mag   ! apparent magnitude
   
      mag = M+5*log10(dl)+25
   
   end function convert_abs2appmag

   function convert_vector(x,rotation) result(y)

      implicit none
      real*4,intent(in)    :: x(3)
      integer*4,intent(in) :: rotation
      real*4               :: y(3)
   
      select case (abs(rotation))
         case (1) ! identity
            y = x
         case(2) ! invert x-axis, while permuting y and z
            y = (/-x(1),x(3),x(2)/)
         case(3) ! invert y-axis, while permuting z and x
            y = (/x(3),-x(2),x(1)/)
         case(4) ! invert z-axis, while permuting x and y
            y = (/x(2),x(1),-x(3)/)
         case(5) ! permute (x,y,z) -> (y,z,x)
            y = (/x(2),x(3),x(1)/)
         case(6) ! permute (x,y,z) -> (z,x,y)
            y = (/x(3),x(1),x(2)/)
      end select
   
      if (rotation<0) y = -y ! inversion
   
   end function convert_vector

   function convert_pseudovector(x,rotation) result(y)

      implicit none
      real*4,intent(in)    :: x(3)
      integer*4,intent(in) :: rotation
      real*4               :: y(3)
   
      y = convert_vector(x,abs(rotation))
   
   end function convert_pseudovector
   
end function convert_properties

end module module_user