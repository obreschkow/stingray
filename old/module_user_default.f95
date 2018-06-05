module module_user

! ==============================================================================================================
! IMPORT MODULES
! ==============================================================================================================

! default modules, do not edit
use module_constants
use module_types
use module_system
use module_cosmology
use module_conversion

! custom modules
! ...


! ==============================================================================================================
! SET DEFAULT PARAMETER FILENAME
! ==============================================================================================================

character(len=255),parameter  :: parameter_filename_default = 'parameters_default.txt'


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
   real*4      :: mag            ! absolute magnitude (generic band)
   real*4      :: MHI            ! [Msun] HI mass
   real*4      :: v(3)           ! [proper km/s] peculiar velocity 
   real*4      :: j(3)           ! [proper kpc km/s] specific angular momentum

end type type_galaxy_sam

! specify the galaxy properties in the mock cone
! these are mostly apparent galaxy properties

type type_galaxy_cone

   integer*8   :: id             ! unique galaxy ID
   integer*4   :: snapshot       ! snapshot
   real*4      :: z              ! apparent redshift
   real*4      :: dc             ! [simulation units = Mpc/h] comoving distance
   real*4      :: ra             ! [rad] right ascension
   real*4      :: decl           ! [rad] declination
   real*4      :: mag            ! apparent magnitude (generic band)
   real*8      :: SHI            ! [W/m^2] integrated HI line flux
   real*4      :: v(3)           ! [proper km/s] peculiar velocity
   real*4      :: j(3)           ! [proper kpc km/s] specific angular momentum

end type type_galaxy_cone

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

! set parameters automatically (e.g. from information provided with the snapshot files)
! this overwrites parameters specified in the parameter textfile

subroutine make_automatic_parameters
   
   implicit none
   
   ! Here, specify how certain parameters, if any, are to be set automatically
   ! ...
   
end subroutine make_automatic_parameters

! write mock-cone galaxy into binary file
! this function allows the user to select only certain properties, if needed, and/or add intrinsic properties

subroutine write_galaxy(cone)

   ! choose which variables of the structure 'cone' to save
   implicit none
   type(type_galaxy_cone),intent(in) :: cone
   write(1) cone
   
end subroutine write_galaxy

! load redshifts

subroutine load_redshifts(snapshot)

   implicit none
   type(type_snapshot),intent(inout),allocatable   :: snapshot(:)
   integer*4                                       :: isnapshot,i
   character(len=255)                              :: filename
   real*4                                          :: redshift
   
   ! make filename
   write(filename,'(A)') trim(para%path_input)//'redshifts.txt'
   call check_exists(trim(filename))
   
   ! read redshifts into snapshot structure
   allocate(snapshot(para%snapshot_min:para%snapshot_max))
   open(1,file=trim(filename),action='read',form="formatted")
   do isnapshot = para%snapshot_min,para%snapshot_max
      read(1,*) i,redshift
      if (isnapshot.ne.i) then
         call out('ERROR: snapshot_min and snapshot_max inconsistent with snapshots in redshift file.')
         close(1)
         stop
      else
         snapshot(isnapshot)%redshift = redshift
      end if
   end do
   close(1)
    
end subroutine load_redshifts

! load SAM snapshot file

subroutine load_sam_snapshot(index,subindex,sam,snapshotname)

   ! variable declaration
   implicit none
   integer*4,intent(in)                            :: index          ! snapshot index
   integer*4,intent(in)                            :: subindex       ! subindex, if the snapshot is split into several files
   type(type_galaxy_sam),allocatable,intent(out)   :: sam(:)         ! derived type of all the relevant SAM properties
   character(len=100),intent(out)                  :: snapshotname   ! snapshot name to be returned for user display
   character(len=255)                              :: filename
   integer*8                                       :: i,n
   integer*4                                       :: bytesgal
   integer*8                                       :: bytestot
   
   if (.false.) write(*,*) subindex ! avoids compiler warning if argument unused
   
   ! make filename
   write(filename,'(A,A,I0.3)') trim(para%path_input),'snapshot_',index
   call check_exists(trim(filename))
   
   ! determine number of galaxies from file size
   bytesgal = bytes_per_galaxy_sam()
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
   
   ! return snapshot name for screen output
   write(snapshotname,'(A,I0,A,I0,A)') 'snapshot ',index,' (',n,' galaxies)'
   
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
   cone%mag    = convert_absmag2appmag(sam%mag,dl)
   cone%v      = convert_vector(sam%v,base%rotation)
   cone%j      = convert_pseudovector(sam%j,base%rotation)
   cone%SHI    = convert_luminosity2flux(real(sam%MHI,8)*real(LMratioHI,8)*Lsun,dl)
   
end function convert_properties


! ==============================================================================================================
! HANDLE CUSTOM TASKS (= optional subroutines and arguments)
! ==============================================================================================================

subroutine handle_custom_arguments(task,custom_option,success)

   implicit none
   character(*),intent(in) :: task
   character(*),intent(in) :: custom_option
   logical,intent(out)     :: success
   integer*4               :: ngalaxies
   
   if (.false.) write(*,*) custom_option ! avoids compiler warning if argument unused
   
   ! custom task handler
   success = .true.
   select case (trim(task))
   case ('my.task')
      call out('Here, specify what to do as "my.task"')
      if (len(custom_option)>0) call out('Using the custom option: '//custom_option)
   case ('make.fake.data')
      if (len(custom_option)>0) then
         read(custom_option,*) ngalaxies
      else
         ngalaxies = 100
      end if
      call make_fake_data(ngalaxies)
   case default
      success = .false.
   end select
   
end subroutine handle_custom_arguments

subroutine make_fake_data(ngalaxies)

   implicit none
   integer*4,intent(in)                :: ngalaxies
   integer*4                           :: snapshot
   integer*4                           :: i
   type(type_galaxy_sam),allocatable   :: galaxy(:)
   
   call tic
   call out('MAKE FAKE DATA')

   allocate(galaxy(ngalaxies))
   
   do snapshot = para%snapshot_min,para%snapshot_max
   
      do i = 1,ngalaxies
         call random_number(galaxy(i)%position)
         galaxy(i)%id = i+int(1e5)*snapshot
         galaxy(i)%haloid = galaxy(i)%haloid
         galaxy(i)%position = galaxy(i)%position*para%L
         galaxy(i)%mag  = get_normal_random_number(-23.0,1.0)
         galaxy(i)%MHI  = 10.0**(get_normal_random_number(9.0,1.0))
         galaxy(i)%v(1) = get_normal_random_number(0.0,1e2)
         galaxy(i)%v(2) = get_normal_random_number(0.0,1e2)
         galaxy(i)%v(3) = get_normal_random_number(0.0,1e2)
         galaxy(i)%j(1) = get_normal_random_number(0.0,1e3)
         galaxy(i)%j(2) = get_normal_random_number(0.0,1e3)
         galaxy(i)%j(3) = get_normal_random_number(0.0,1e3)
      end do
      call save_snapshot(snapshot)
      
   end do
   
   call save_redshifts
   
   call toc
   
   contains
   
   subroutine save_redshifts

      implicit none
      integer*4   :: isnapshot
   
      call out('Save redshifts')
      
      open(1,file=trim(para%path_input)//'redshifts.txt',action='write',form='formatted',status='replace')
      do isnapshot = para%snapshot_min,para%snapshot_max
         write(1,*) isnapshot,(para%snapshot_max-isnapshot)*0.1
      end do
      close(1)
    
   end subroutine save_redshifts

   subroutine save_snapshot(index)

      ! variable declaration
      implicit none
      integer*4,intent(in) :: index
      character(len=255)   :: fn,txt
      integer*8            :: i,n
   
      ! write user info
      write(fn,'(A,A,I0.3)') trim(para%path_input),'snapshot_',index
      call out('Save snapshot '//trim(fn))
   
      ! write header
      open(1,file=trim(fn),action='write',form='unformatted',status='replace',access='stream')
   
      ! write IDs and positions
      n = size(galaxy)
      do i = 1,n
         write(1) galaxy(i)
      end do
   
      ! read IDs
      close(1)
   
      ! output basic statistics
      call out('Number of galaxies:',n)
      write(txt,'(A,F0.4,A,F0.4)') 'Position range: ',minval((/minval(galaxy(:)%position(1)),minval(galaxy(:)%position(2)), &
      & minval(galaxy(:)%position(3))/)) &
      &,' to ',maxval((/maxval(galaxy(:)%position(1)),maxval(galaxy(:)%position(2)),maxval(galaxy(:)%position(3))/))
      call out(txt)
      write(txt,'(A,I0,A,I0)')     'Identifier range: ',minval(galaxy%id),' to ',maxval(galaxy%id)
      call out(txt)
   
   end subroutine save_snapshot
   
end subroutine make_fake_data

end module module_user