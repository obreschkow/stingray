! =================================================================================================================
! WHAT IS THIS CODE?
! This module calculates emission line profiles of inclined, rotating galactic disks with velocity dispersion,
! following the prescription of Obreschkow et al. 2009, ApJ 698, 1467. A the moment, only HI lines are generated
! by the module, but the implementation of other lines (e.g. CO) is trivial, given a prescription for their
! integrated flux and spatial distribution.
! Note that the line spectra are given in units of Jy on a velocity grid in rest-frame physical km/s. In this
! standard convention an extra z-dependence comes into the relationship between mass and integrated flux Sv ("v"
! stands for velocity-integrated). The correct relationship is MHI/Msun = 2.36e-5*(Dl/Mpc)^2*Sv/(1+z).
!
! COMPILER
! The code was developed in fortran 90 standard and tested on the GNU fortran compiler version 4.8.2.
! The compiler option -O3 is recommended for accelerated performance.
!
! VERSION HISTORY
! February 2016: Original version - HI line only (contact author for other lines)
! (Code written by Danail Obreschkow, danail.obreschkow@icrar.org)
! =================================================================================================================

module module_emission_lines

use module_system

   type type_lineinput
   
      real  :: mhalo          ! [Msun/h] Halo mass without disk and bulge
      real  :: mdisk          ! [Msun/h] Disk mass (all baryons)
      real  :: mbulg          ! [Msun/h] bulge mass
      real  :: mHI            ! [Msun/h] HI mass
      real  :: mH2            ! [Msun/h] H2 mass
      real  :: rhalo          ! [kpc, physical] characteristic NFW radius
      real  :: rdisk          ! [kpc, physical] half-mass radius of disk (stars+HI)
      real  :: rbulg          ! [kpc, physical] half-mass radius of bulge (stars+HI)
      real  :: rgas           ! [kpc, physical] half-mass radius of cold gas (HI+H2)
      real  :: dispersionHI   ! [km/s] isotropic velocity dispersion of HI
      real  :: dispersionH2   ! [km/s] isotropic velocity dispersion of H2
      real  :: z              ! [-] redshift
      real  :: incl           ! [rad] inclination (0 = face-on, pi/2 = edge-on)
      real  :: Dc             ! [Mpc/h] comoving distance
      
   end type type_lineinput
   
   type type_line
   
      real,allocatable  :: v(:)     ! [km/s] central velocities of velocity channels
      real,allocatable  :: s(:)     ! [Jy] flux density in velocity channels
      real              :: sintvel  ! [Jy km/s] velocity-integrated flux
      real              :: speak    ! [Jy] peak flux density
      real              :: scentral ! [Jy] central flux density
      real              :: w20      ! [km/s] line with at 20% of peak flux (in rest-frame, measured from the outside)
      real              :: w50      ! [km/s] line with at 50% of peak flux (in rest-frame, measured from the outside)
      real              :: wpeak    ! [km/s] line with between the peaks (in rest-frame)
   
   end type type_line
   
   contains
  
   subroutine line_profile(input,h,line)
   
    implicit none
    
      ! arguments
      real,intent(in)                  :: h            ! [-] Hubble parameter H0/[100 km/s/Mpc]
      type(type_lineinput),intent(in)  :: input        ! [*] Input galaxy parameters
      type(type_line),intent(out)      :: line(4)      ! [*] Emission line properties: 1/2 (3/4) = projected/edge-on HI line (H2 line)
    
      ! parameters (probably no need to change these)
      logical,parameter    :: error_check = .true.
      logical,parameter    :: adaptive_dv = .false. ! if true, the width of the channels is adapted to keep their number const.
      real,parameter       :: G = 4.3022682e-6     ! [(km/s)^2 kpc/Msun] gravitational constant
      real,parameter       :: xmax = 10            ! [-] maximal integration radius in units of gas exponential scale
      real,parameter       :: dx = 0.02            ! [-] annulus thickness in units of gas exponential scale
      integer,parameter    :: nx = nint(xmax/dx)   ! [-] number of annuli for numerical integration
      integer,parameter    :: nv_fixed = 200       ! [-] maximal number of velocity channels
      real,parameter       :: dv_fixed = 2         ! [km/s] fixed width of velocity channels
      real,parameter       :: pi = 3.14159265359

      ! computation variables
      real                 :: v_halo_sqr(nx)       ! [km^2/s^2] square of circular velocity due to halo
      real                 :: v_disk_sqr(nx)       ! [km^2/s^2] square of circular velocity due to disk
      real                 :: v_bulg_sqr(nx)       ! [km^2/s^2] square of circular velocity due to bluge
      real                 :: v_circ(nx)           ! [km/s] total circular velocity
      real                 :: Sigma(2,nx)          ! [-] normalised surface density of HI/H2
      real                 :: x(nx)                ! [-] central radii of annuli relative to input%rgas
      real,allocatable     :: s(:)                 ! [s/km] normalised flux density
      real                 :: sini                 ! [-] sine of inclination
      real                 :: dispHI               ! [km/s] velocity dispersion of HI
      real                 :: dispH2               ! [km/s] velocity dispersion of H2
      real                 :: vmax                 ! [km/s] velocity of highest velocity channel
      real                 :: dv                   ! [km/s] width of velocity channels
      real,allocatable     :: filter(:,:)          ! [-] normalised smoothing filter to add velocity dispersion
      integer              :: nfilter              ! [-] index range for filter evaluation
      integer              :: nv                   ! [-] number of velocity channels
      integer              :: i,j,k,ij             ! [-] integration indices
      real                 :: Rmol                 ! [-] mH2/mHI ratio
      real                 :: Rcmol                ! [-] central SigmaH2/SigmaHI ratio
      real                 :: rgasexp              ! [pkpc] exponential scale-length or cold gas (HI+H2)
      integer              :: id                   ! index distinguishing HI and H2
      
      ! auxiliary variables, mainly used for fast computation
      real                 :: f,dy,starget,yp,ym
      real,allocatable     :: y(:) 
      
      ! Check input arguments
      if (error_check) then
         if (input%mhalo<0) call error('Error in subroutine line_profile. All masses must be positive.')
         if (input%mdisk<0) call error('Error in subroutine line_profile. All masses must be positive.')
         if (input%mbulg<0) call error('Error in subroutine line_profile. All masses must be positive.')
         if (input%mHI<0) call error('Error in subroutine line_profile. All masses must be positive.')
         if (input%mH2<0) call error('Error in subroutine line_profile. All masses must be positive.')
         if (input%rdisk<=0) call error('Error in subroutine line_profile. All radii must be positive.')
         if (input%rbulg<=0) call error('Error in subroutine line_profile. All radii must be positive.')
         if (input%rgas<=0) call error('Error in subroutine line_profile. All radii must be positive.')
         if (input%mhalo>1e16) call error('Error in subroutine line_profile. Halo mass excessively high.')
      end if
      
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
      !   if (vmax>1e4) call error('Error in subroutine line_profile. Circular velocity > 10,000 km/s.')
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
      
      allocate(y(nv),s(nv))
      
      do k = 1,4
      
         ! allocate emission line
         allocate(line(k)%v(nv))
         allocate(line(k)%s(nv))
         line(k)%v = ((/(i, i = 1, nv)/)-0.5)*dv
         
         ! Set galaxy inclination
         if ((k==1).or.(k==3)) then
            sini = 1.0 ! edge-on
         else
            sini = max(1e-4,abs(sin(input%incl))) ! avoids division by 0 for incl = 0
         end if
         
         ! Set line index
         if ((k==1).or.(k==2)) then
            id = 1 ! HI-line
         else
            id = 2 ! H2-line
         end if
      
         ! Convolve v_circ(r) and Sigma(r) into a symmetric, normalised emission line
         s = 0
         do j = 1,nx
            ! The annulus j generates an emission line 1/pi/sqrt(sini^2-(v/v_circ)^2),
            ! which analytically integrates to arcsin(y)/pi, where y = v/v_circ/sini
            y = line(k)%v/v_circ(j)/sini
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
         line(k)%s = 0
         do i = 1,nv
            do j = -nfilter,nfilter
               ij = i+j
               if (ij<1) ij = 1-ij
               if (ij<=nv) line(k)%s(i) = line(k)%s(i)+s(ij)*filter(id,j)
            end do
         end do
         
         ! Check normalisation
         if (error_check) then
            if (abs(sum(line(k)%s*dv)-0.5)>0.01) then
               write(*,*) sum(line(k)%s*dv)
               write(*,*) k,vmax,Rmol
               call error('Error in subroutine line_profile. Check arguments.')
            end if
         end if
         
         ! Measure characteristic line parameters
         line(k)%speak = maxval(line(k)%s) ! [s/km]
         line(k)%scentral = line(k)%s(1) ! [s/km]
         do i = 1,nv
            if (line(k)%s(i)==line(k)%speak) then
               line(k)%wpeak = line(k)%v(i)*2
               exit
            end if
         end do
         starget = 0.5*line(k)%speak
         do i = nv-1,1,-1
            if (line(k)%s(i)>starget) then
               f = (line(k)%s(i)-starget)/(line(k)%s(i)-line(k)%s(i+1))
               line(k)%w50 = ((1-f)*line(k)%v(i)+f*line(k)%v(i+1))*2
               exit
            end if
         end do
         starget = 0.2*line(k)%speak
         do i = nv-1,1,-1
            if (line(k)%s(i)>starget) then
               f = (line(k)%s(i)-starget)/(line(k)%s(i)-line(k)%s(i+1))
               line(k)%w20 = ((1-f)*line(k)%v(i)+f*line(k)%v(i+1))*2
               exit
            end if
         end do
         
      end do
      
      !if (abs((input%mdisk/1.92865059e+11+input%mbulg/1.36288600e+12)-2.0)<1e-5) then
      !if (input%mdisk<1.1e10) then
      !   call save_profiles('/Users/do/Desktop/testprofile.txt')
      !   stop
      !end if
      
      contains
      
      subroutine save_profiles(filename)
      
         ! produces inputs for R-routine "show_emission_line.R"
      
         implicit none
         character(*),intent(in) :: filename
         integer*4               :: i
      
         open(1,file=trim(filename),action='write',form="formatted",status='replace')
         write(1,'(1Es15.6)') input%incl
         do k = 1,4
            write(1,'(5Es15.6)') line(k)%speak,line(k)%scentral,line(k)%wpeak,line(k)%w50,line(k)%w20
         end do
         do i = 1,size(x)
            write(1,'(7Es15.6)') x(i),Sigma(1,i),Sigma(2,i),sqrt(v_halo_sqr(i)),sqrt(v_disk_sqr(i)),sqrt(v_bulg_sqr(i)),v_circ(i)
         end do
         close(1)
      
      end subroutine save_profiles
       
   end subroutine line_profile
  
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