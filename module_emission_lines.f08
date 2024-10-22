! **********************************************************************************************************************************
! This Fortran module calculates emission line profiles of inclined, rotating galactic disks with velocity dispersion,
! following the prescription of Obreschkow et al. 2009, ApJ 698, 1467. A the moment, only HI lines and molecular are generated
! by the module, but the implementation of other lines is straight forward, given a prescription for their spatial distribution.
!
! Note that flux densities are given in normalised units (s/km), such that the total line integral along the rest-frame velocity
! axis is equal to unity.
!
! Developed by Danail Obreschkow, danail.obreschkow@icrar.org
! **********************************************************************************************************************************

module module_emission_lines

   use shared_module_core
   use shared_module_constants
   use shared_module_cosmology

   public   :: make_emission_line
   public   :: type_line_input
   public   :: type_line_shape
   public   :: type_line_profile
   public   :: empty_line_shape
   
   private
   
   type type_line_input
   
      real*4   :: mhalo          ! [Msun/h] Halo mass without disk and bulge
      real*4   :: mdisk          ! [Msun/h] Disk mass (all baryons)
      real*4   :: mbulg          ! [Msun/h] bulge mass
      real*4   :: mHI            ! [Msun/h] HI mass
      real*4   :: mH2            ! [Msun/h] H2 mass
      real*4   :: rhalo          ! [kpc, physical] characteristic NFW radius
      real*4   :: rdisk          ! [kpc, physical] half-mass radius of disk (stars+HI)
      real*4   :: rbulg          ! [kpc, physical] half-mass radius of bulge (stars+HI)
      real*4   :: rgas           ! [kpc, physical] half-mass radius of cold gas (HI+H2)
      real*4   :: dispersionHI   ! [km/s] isotropic velocity dispersion of HI
      real*4   :: dispersionH2   ! [km/s] isotropic velocity dispersion of H2
      real*4   :: z              ! [-] redshift
      real*4   :: incl           ! [rad] inclination (0 = face-on, pi/2 = edge-on)
      real*4   :: Dc             ! [Mpc/h] comoving distance
      
   end type type_line_input
   
   type type_line_shape
   
      ! shape parameters in observed galaxy-projection
      real*4   :: speak = 0      ! [s/km] normalised peak flux density (multiply by velocity-integrated flux to get Jy)
      real*4   :: scentral =0    ! [s/km] normalised central flux density (multiply by velocity-integrated flux to get Jy)
      real*4   :: w20 = 0        ! [km/s] line with at 20% of peak flux (in rest-frame, measured from the outside)
      real*4   :: w50 = 0        ! [km/s] line with at 50% of peak flux (in rest-frame, measured from the outside)
      real*4   :: wpeak = 0      ! [km/s] line with between the peaks (in rest-frame)
      
      ! shape parameters of galaxy seen edge-on
      real*4   :: speak_eo = 0   ! [s/km] normalised peak flux density (multiply by velocity-integrated flux to get Jy)
      real*4   :: scentral_eo =0 ! [s/km] normalised central flux density (multiply by velocity-integrated flux to get Jy)
      real*4   :: w20_eo = 0     ! [km/s] line with at 20% of peak flux (in rest-frame, measured from the outside)
      real*4   :: w50_eo = 0     ! [km/s] line with at 50% of peak flux (in rest-frame, measured from the outside)
      real*4   :: wpeak_eo = 0   ! [km/s] line with between the peaks (in rest-frame)
   
   end type type_line_shape
   
   type type_line_profile
   
      real*4,allocatable  :: v(:)      ! [km/s] central velocities of velocity channels
      real*4,allocatable  :: s(:)      ! [s/km] normalised flux density in velocity channels
      real*4,allocatable  :: s_eo(:)   ! [s/km] normalised flux density in velocity channels when seen edge-on
      
   end type type_line_profile
   
   type(type_line_shape),protected   :: empty_line_shape
   
   contains
  
   subroutine make_emission_line(input,line_shape,line_profile)
   
      implicit none
    
      ! arguments
      type(type_line_input),intent(in)             :: input ! Input galaxy parameters
      type(type_line_shape),intent(out)            :: line_shape(2) ! Line shape parameters of HI and H2 lines
      type(type_line_profile),intent(out),optional :: line_profile(2) ! Line profiles of HI and H2 lines
    
      ! parameters (probably no need to change these)
      logical*4,parameter  :: error_check = .true.
      logical*4,parameter  :: adaptive_dv = .false.      ! if true, the width of the channels is adapted to keep their number const.
      real*4,parameter     :: G = 64462.3408*6.67430e-11 ! [(km/s)^2 kpc/Msun] gravitational constant
      real*4,parameter     :: xmax = 10                  ! [-] maximal integration radius in units of gas exponential scale
      real*4,parameter     :: dx = 0.02                  ! [-] annulus thickness in units of gas exponential scale
      integer*4,parameter  :: nx = nint(xmax/dx)         ! [-] number of annuli for numerical integration
      integer*4,parameter  :: nv_fixed = 200             ! [-] maximal number of velocity channels
      real*4,parameter     :: dv_fixed = 2               ! [km/s] fixed width of velocity channels
            
      ! computation variables      
      real*4               :: v_halo_sqr(nx)             ! [km^2/s^2] square of circular velocity due to halo
      real*4               :: v_disk_sqr(nx)             ! [km^2/s^2] square of circular velocity due to disk
      real*4               :: v_bulg_sqr(nx)             ! [km^2/s^2] square of circular velocity due to bluge
      real*4               :: v_circ(nx)                 ! [km/s] total circular velocity
      real*4               :: Sigma(2,nx)                ! [-] normalised surface density of HI/H2
      real*4               :: x(nx)                      ! [-] central radii of annuli relative to input%rgas
      real*4,allocatable   :: s(:)                       ! [s/km] normalised flux density
      real*4               :: sini                       ! [-] sine of inclination
      real*4               :: dispHI                     ! [km/s] velocity dispersion of HI
      real*4               :: dispH2                     ! [km/s] velocity dispersion of H2
      real*4               :: vmax                       ! [km/s] velocity of highest velocity channel
      real*4               :: dv                         ! [km/s] width of velocity channels
      real*4,allocatable   :: filter(:,:)                ! [-] normalised smoothing filter to add velocity dispersion
      integer*4            :: nfilter                    ! [-] index range for filter evaluation
      integer*4            :: nv                         ! [-] number of velocity channels
      integer*4            :: i,j,k,ij                   ! [-] integration indices
      real*4               :: Rmol                       ! [-] mH2/mHI ratio
      real*4               :: Rcmol                      ! [-] central SigmaH2/SigmaHI ratio
      real*4               :: rgasexp                    ! [pkpc] exponential scale-length or cold gas (HI+H2)
      integer*4            :: id                         ! index distinguishing HI and H2
      real*4               :: h                          ! Hubble parameter
      type(type_line_profile) :: profile(2)
      
      ! auxiliary variables, mainly used for fast computation
      real*4               :: f,dy,starget,yp,ym
      real*4,allocatable   :: y(:)
      
      ! Check input arguments
      if (error_check) then
         if (input%mhalo<0) call error('subroutine make_emission_line: All masses must be positive.')
         if (input%mdisk<0) call error('subroutine make_emission_line: All masses must be positive.')
         if (input%mbulg<0) call error('subroutine make_emission_line: All masses must be positive.')
         if (input%mHI<0) call error('subroutine make_emission_line: All masses must be positive.')
         if (input%mH2<0) call error('subroutine make_emission_line: All masses must be positive.')
         if (input%rdisk<=0) call error('subroutine make_emission_line: All radii must be positive.')
         if (input%rbulg<=0) call error('subroutine make_emission_line: All radii must be positive.')
         if (input%rgas<=0) call error('subroutine make_emission_line: All radii must be positive.')
         if (input%mhalo>1e16) call error('subroutine make_emission_line: Halo mass larger than 1e16 Msun.')
      end if
      
      ! Get cosmology
      h = cosmology%h
      
      ! Make vector of annuli for numerical integration
      x = ((/(i, i = 1, nx)/)-0.5)*dx ! radius in multiples of hydrogen exponential scale radius
      
      ! Surface density profile of HI and H2 (normalised to 1)
      Rmol = max(1e-4,min(1e4,input%mH2/(input%mHI+1.0))) ! limits imposed to avoid numerical issues
      Rcmol = (17*Rmol**2.11+8.2*Rmol**1.02)*(1-0.934*exp(-(0.102*log(Rmol)-1.1)**2.0)) ! central SigmaH2/SigmaHI ratio (approximation accurate <1% for Rmol=1e-4..1e4)
      Sigma(1,:) = exp(-x)/(1+Rcmol*exp(-1.6*x))
      Sigma(1,:) = Sigma(1,:)/sum(x*Sigma(1,:))
      Sigma(2,:) = Rcmol*exp(-2.6*x)/(1+Rcmol*exp(-1.6*x))
      Sigma(2,:) = Sigma(2,:)/sum(x*Sigma(2,:))
      
      ! Circular velocity profiles
      rgasexp = 0.5958233*input%rgas
      v_halo_sqr    = G*input%mhalo/h/input%rhalo*rotation_model(1,x*rgasexp/input%rhalo)
      v_disk_sqr    = G*input%mdisk/h/input%rdisk*rotation_model(9,x*rgasexp/input%rdisk)
      v_bulg_sqr    = G*input%mbulg/h/input%rbulg*rotation_model(4,x*rgasexp/input%rbulg)
      v_circ        = sqrt(v_halo_sqr+v_disk_sqr+v_bulg_sqr)
      
      ! Make velocity channels
      vmax = min(5e3,(maxval(v_circ)+max(input%dispersionHI,input%dispersionH2)*2)*1.2) ! [km/s] maximal velocity channel
      !if (error_check) then
      !   if (vmax>1e4) call error('subroutine make_emission_line: Circular velocity > 10,000 km/s.')
      !end if
      if (adaptive_dv) then
         nv = nv_fixed
         dv = vmax/nv
      else
         dv = dv_fixed
         nv = nint(vmax/dv)
      end if
      
      ! Make smoothing filter
      dispHI = max(1e-2,input%dispersionHI) ! [km/s] velocity dispersion, min val to avoid division by 0
      dispH2 = max(1e-2,input%dispersionH2) ! [km/s] velocity dispersion, min val to avoid division by 0
      nfilter = max(1,nint(3*dispHI/dv),nint(3*dispH2/dv))
      allocate(filter(2,-nfilter:nfilter))
      do i = -nfilter,nfilter
         filter(1,i) = exp(-i**2/(dispHI/dv)**2/2)
         filter(2,i) = exp(-i**2/(dispH2/dv)**2/2)
      end do
      filter(1,:) = filter(1,:)/sum(filter(1,:))
      filter(2,:) = filter(2,:)/sum(filter(2,:))
      
      ! allocate emission line
      allocate(y(nv),s(nv))
      do id = 1,2
         allocate(profile(id)%v(nv))
         allocate(profile(id)%s(nv))
         allocate(profile(id)%s_eo(nv))
         profile(id)%v = ((/(i, i = 1, nv)/)-0.5)*dv ! [km/s] velocity channels
      end do
      
      do k = 1,4 ! 1 = HI edge-on, 2 = HI projected, 3 = H2 edge-on, 4 = H2 projected
      
         ! Set line index
         if ((k==1).or.(k==2)) then
            id = 1 ! HI-line
         else
            id = 2 ! H2-line
         end if
         
         ! Set galaxy inclination
         if ((k==1).or.(k==3)) then
            sini = 1.0 ! edge-on
         else
            sini = max(1e-4,abs(sin(input%incl))) ! avoids division by 0 for incl = 0
         end if
      
         ! Convolve v_circ(r) and Sigma(r) into a symmetric, normalised emission line
         s = 0
         do j = 1,nx
            ! The annulus j generates an emission line 1/pi/sqrt(sini^2-(v/v_circ)^2),
            ! which analytically integrates to arcsin(y)/pi, where y = v/v_circ/sini
            y = profile(id)%v/v_circ(j)/sini
            dy = dv/v_circ(j)/sini
            f = x(j)*Sigma(id,j)/pi/dv
            do i = 1,nv
               ym = y(i)-dy/2
               yp = y(i)+dy/2
               if (yp>1) then
                  yp = 1
                  if (ym>1) exit
               end if
               s(i) = s(i)+f*(asin(yp)-asin(ym))
            end do
         end do
         
         ! Smooth lines by velocity dispersion
         profile(id)%s = 0
         do i = 1,nv
            do j = -nfilter,nfilter
               ij = i+j
               if (ij<1) ij = 1-ij
               if (ij<=nv) profile(id)%s(i) = profile(id)%s(i)+s(ij)*filter(id,j)
            end do
         end do
         
         ! Check normalisation
         if (error_check) then
            if (abs(sum(profile(id)%s*dv)-0.5)>0.01) then
               write(*,*) sum(profile(id)%s*dv)
               write(*,*) k,vmax,Rmol
               call error('subroutine make_emission_line: check arguments.')
            end if
         end if
         
         ! Measure characteristic line parameters
         line_shape(id)%speak = maxval(profile(id)%s) ! [s/km]
         line_shape(id)%scentral = profile(id)%s(1) ! [s/km]
         do i = 1,nv
            if (profile(id)%s(i)==line_shape(id)%speak) then
               line_shape(id)%wpeak = profile(id)%v(i)*2
               exit
            end if
         end do
         starget = 0.5*line_shape(id)%speak
         do i = nv-1,1,-1
            if (profile(id)%s(i)>starget) then
               f = (profile(id)%s(i)-starget)/(profile(id)%s(i)-profile(id)%s(i+1))
               line_shape(id)%w50 = ((1-f)*profile(id)%v(i)+f*profile(id)%v(i+1))*2
               exit
            end if
         end do
         starget = 0.2*line_shape(id)%speak
         do i = nv-1,1,-1
            if (profile(id)%s(i)>starget) then
               f = (profile(id)%s(i)-starget)/(profile(id)%s(i)-profile(id)%s(i+1))
               line_shape(id)%w20 = ((1-f)*profile(id)%v(i)+f*profile(id)%v(i+1))*2
               exit
            end if
         end do
         
         ! copy into _eo
         if ((k==1).or.(k==3)) then
            profile(id)%s_eo = profile(id)%s
            line_shape(id)%speak_eo = line_shape(id)%speak
            line_shape(id)%scentral_eo = line_shape(id)%scentral
            line_shape(id)%w20_eo = line_shape(id)%w20
            line_shape(id)%w50_eo = line_shape(id)%w50
            line_shape(id)%wpeak_eo = line_shape(id)%wpeak
         end if
         
      end do
      
      if (present(line_profile)) line_profile = profile
      
      !call save_profiles('/Users/do/Desktop/testprofile.txt')
      !
      !contains
      !
      !subroutine save_profiles(filename)
      !
      !   ! produces inputs for R-routine "show_emission_line.R"
      !
      !   implicit none
      !   character(*),intent(in) :: filename
      !   integer*4               :: i
      !
      !   open(1,file=trim(filename),action='write',form="formatted",status='replace')
      !   write(1,'(1Es15.6)') input%incl
      !   do k = 1,2
      !      write(1,'(10Es15.6)') line_shape(id)%speak,line_shape(id)%scentral,line_shape(id)%wpeak,line_shape(id)%w50,line_shape(id)%w20,&
      !      &line_shape(id)%speak_eo,line_shape(id)%scentral_eo,line_shape(id)%wpeak_eo,line_shape(id)%w50_eo,line_shape(id)%w20_eo
      !   end do
      !   do i = 1,size(x)
      !      write(1,'(7Es15.6)') x(i),Sigma(1,i),Sigma(2,i),sqrt(v_halo_sqr(i)),sqrt(v_disk_sqr(i)),sqrt(v_bulg_sqr(i)),v_circ(i)
      !   end do
      !   close(1)
      !
      !end subroutine save_profiles
       
   end subroutine make_emission_line
  
   function rotation_model(model,x) result(f)
      
      ! characteristic halo velocity profile, such that v = sqrt(GM/R * f),
      ! where M is the total mass and R is the half-mass radius
      ! except in the case of the NFW profile, where a total mass is ill-defined and M and R a specified below
      
      implicit none
      integer,intent(in)   :: model ! index of physical model for rotation curve
      real,intent(in)      :: x(:)  ! radius normalised to the half-mass radius, except for NFW, where it is normalised to R_s
      real,allocatable     :: f(:)
      real,allocatable     :: s(:)
      real                 :: a
      
      allocate(f(size(x)))
      allocate(s(size(x)))
      
      select case (model)
      
         case (1)
            ! NFW profile, describing CDM haloes
            ! M = mass enclosed in side R
            ! R = characteristic NFW halo
            a = 5.177399
            f = real(a*(log(1.0_8+1e-10_8+x)/(x+1e-10)-1/(1.0_8+x)),4) ! real*8 avoids negatives for values of x that are too small
         
         case (2)
            ! Exponential flat disk (1%-approximation of Freeman solution, avoiding four modified Bessel functions)
            a = 1.67835 ! R/Rexp
            s = a*x+1e-10 ! radius normalised to Rexp, 1e-10 avoids floating point exception if x=0
            f = a*(1+4.8*exp(max(-80.0,-0.35*s-3.5/s)))/(s+1/s**2+2/sqrt(s)) ! max(-80,...) avoids floating point exception for small x
      
         case (3)
            ! Plummer model for spheroids (bulges)
            a = 0.766421 ! (Plummer radius)/(half-mass radius)
            f = x**2/(a**2+x**2)**1.5
            
         case (4)
            ! Hernquist model for spheroids (bulges)
            ! NB: The Hernquist model is particularly appealing as it arises from the numerical simulation of the
            !     merger of two equal mass disk galaxies, each embedded within a dark matter halo.
            a = 0.4142137 ! (Hernquist radius)/(half-mass radius)
            f = x/(a+x)**2
            
         case (5)
            ! Kepler potential of a point mass M
            ! (radius and velocity normalised consistently with the spherical potentials above)
            f = 1.0/(x+1e-10) ! 1e-10 to avoid singularity at x=0
            
         case (6)
            ! Spherical Jaffe potential
            f = 1/(1+x)
            
         case (7)
            ! Truncated singular isothermal potential (2R contains the full mass)
            f = 0.25+0.5/(x+1e-10)-abs(0.25-0.5/(x+1e-10)) ! 1e-10 to avoid singularity at x=0
            
         case (8)
            ! Uniform sphere
            f = (x**2+2/(x+1e-10)-abs(2.0**(2.0/3.0)-2/(x+1e-10))-abs(2.0**(2.0/3.0)-x**2))/4.0 ! 1e-10 to avoid singularity at x=0
            
         case (9)
            ! Miyamoto-Nagai potential for a thick disk (half-mass axis ratio of 1/10)
            a = 0.7057713 ! (characteristic radius a+b)/(half-mass radius)
            f = x**2/(a**2+x**2)**1.5
            
         case default
            call error('Rotation curve model unknown.')
            
      end select
         
   end function rotation_model

end module module_emission_lines