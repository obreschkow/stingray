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
   
      real  :: mhalo       ! [Msun/h] Halo mass without disk and bulge
      real  :: mdisk       ! [Msun/h] Disk mass (all baryons)
      real  :: mbulg       ! [Msun/h] bulge mass
      real  :: mHI         ! [Msun/h] HI mass
      real  :: rvir        ! [kpc, physical] virial radius
      real  :: rdisk       ! [kpc, physical] exponential disk scale (stars+HI)
      real  :: rHI         ! [kpc, physical] exponential scale of HI
      real  :: dispersion  ! [km/s] velocity dispersion of the line emitting component
      real  :: z           ! [-] redshift
      real  :: incl        ! [rad] inclination (0 = face-on, pi/2 = edge-on)
      real  :: Dc          ! [Mpc/h] comoving distance
      
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
  
   subroutine line_profile(input,h,line,line_edge_on)
    
    implicit none
    
      ! arguments
      real,intent(in)                  :: h            ! [-] Hubble parameter H0/[100 km/s/Mpc]
      type(type_lineinput),intent(in)  :: input        ! [*] Input galaxy parameters
      type(type_line),intent(out)      :: line         ! [*] Observed emission line properties 
      type(type_line),intent(out)      :: line_edge_on ! [*] Edge-on emission line properties (same as observed, if incl = pi/2)
    
      ! parameters (probably no need to change these)
      logical,parameter    :: error_check = .true.
      logical,parameter    :: adaptive_dv = .false. ! if true, the width of the channels is adapted to keep their number const.
      real,parameter       :: G = 4.3022682e-6     ! [(km/s)^2 kpc/Msun] gravitational constant
      real,parameter       :: xmax = 10            ! [-] maximal integration radius relative to input%rHI
      real,parameter       :: dx = 0.02            ! [-] annulus thickness relative to input%rHI
      integer,parameter    :: nx = nint(xmax/dx)   ! [-] number of annuli for numerical integration
      integer,parameter    :: nv_fixed = 200       ! [-] maximal number of velocity channels
      real,parameter       :: dv_fixed = 2         ! [km/s] fixed width of velocity channels
      real,parameter       :: pi = 3.14159265359

      ! computation variables
      real                 :: c_halo               ! [-] halo concentration
      real                 :: c_disk               ! [-] disk concentration
      real                 :: c_bulg               ! [-] bulge concentration (in 3D)
      real                 :: v_halo_sqr(nx)       ! [km^2/s^2] square of circular velocity due to halo
      real                 :: v_disk_sqr(nx)       ! [km^2/s^2] square of circular velocity due to disk
      real                 :: v_bulg_sqr(nx)       ! [km^2/s^2] square of circular velocity due to bluge
      real                 :: v_circ(nx)           ! [km/s] total circular velocity
      real                 :: Sigma(nx)            ! [-] normalised surface density
      real                 :: x(nx)                ! [-] central radii of annuli relative to input%rdisk
      real,allocatable     :: s(:)                 ! [s/km] normalised flux density
      real                 :: sini                 ! [-] sine of inclination
      real                 :: disp                 ! [km/s] velocity dispersion
      real                 :: vmax                 ! [km/s] velocity of highest velocity channel
      real                 :: dv                   ! [km/s] width of velocity channels
      real,allocatable     :: filter(:)            ! [-] normalised smoothing filter to add velocity dispersion
      integer              :: nfilter              ! [-] index range for filter evaluation
      integer              :: nv                   ! [-] number of velocity channels
      integer              :: i,j,k,ij             ! [-] integration indices
      
      ! auxiliary variables, mainly used for fast computation
      real                 :: f,dy,starget,yp,ym
      real,allocatable     :: y(:) 
      
      ! Check input arguments
      if (error_check) then
         if (input%mhalo<=0) call error('Error in subroutine line_profile. All masses must be positive.')
         if (input%mdisk<0) call error('Error in subroutine line_profile. All masses must be positive.')
         if (input%mbulg<0) call error('Error in subroutine line_profile. All masses must be positive.')
         if (input%mHI<0)   call error('Error in subroutine line_profile. All masses must be positive.')
         if (input%rvir<=0)  call error('Error in subroutine line_profile. All radii must be positive.')
         if (input%rdisk<=0) call error('Error in subroutine line_profile. All radii must be positive.')
         if (input%rHI<=0)   call error('Error in subroutine line_profile. All radii must be positive.')
         if (input%mhalo>1e16) call error('Error in subroutine line_profile. Halo mass excessively high.')
         if (input%rdisk>input%rvir) call error('Error in subroutine line_profile: rdisk>rvir seems unphysical.')
         if (input%rHI>input%rvir) call error('Error in subroutine line_profile: rHI>rvir seems unphysical.')
      end if
 
      ! Make vector of annuli for numerical integration
      x = ((/(i, i = 1, nx)/)-0.5)*dx
      
      ! Make surface density profile of HI (normalised to 1)
      Sigma         = exp(-x)
      Sigma         = Sigma/sum(x*Sigma)
 
      ! Estimate concentration parameters of the three mass components
      c_halo        = 12.3/(1+input%z)*(input%mhalo/1.3e13)**(-0.13)
      c_disk        = input%rvir/input%rdisk
      c_bulg        = 12.0*c_disk
write(*,*) 'c'
      ! Make circular velocity profile
      v_halo_sqr    = G*input%mhalo/h/input%rvir*f_halo(c_halo,x*input%rHI/input%rvir)
      write(*,*) 'c1'
      v_disk_sqr    = G*input%mdisk/h/input%rvir*f_disk(c_disk,x*input%rHI/input%rvir)
      write(*,*) 'c2'
      v_bulg_sqr    = G*input%mbulg/h/input%rvir*f_bulg(c_bulg,x*input%rHI/input%rvir)
      write(*,*) 'c3'
      write(*,*) v_halo_sqr(1),G*input%mhalo/h/input%rvir,c_halo,input%rHI/input%rvir
      write(*,*) 'xxx'
      write(*,*) x
      v_circ        = sqrt(v_halo_sqr+v_disk_sqr+v_bulg_sqr)
      
      ! Make velocity channels
      vmax = min(5e3,(maxval(v_circ)+input%dispersion*2)*1.2) ! [km/s] maximal velocity channel
      write(*,*) 'd'
      if (error_check) then
         if (vmax/1.2>5e3) call error('Error in subroutine line_profile. Circular velocity > 5000 km/s.')
      end if
      if (adaptive_dv) then
         nv = nv_fixed
         dv = vmax/nv
      else
         dv = dv_fixed
         nv = nint(vmax/dv)
      end if
      
      allocate(line%v(nv),y(nv))
      allocate(line%s(nv),s(nv))
      line%v = ((/(i, i = 1, nv)/)-0.5)*dv
      
      ! Make smoothing filter
      disp = max(1e-2,input%dispersion) ! [km/s] velocity dispersion, min val to avoid division by 0
      nfilter = max(1,nint(3*disp/dv))
      allocate(filter(-nfilter:nfilter))
      do i = -nfilter,nfilter
         filter(i) = exp(-i**2/(disp/dv)**2/2)
      end do
      filter = filter/sum(filter)
      
      do k = 1,2
       
         ! Set galaxy inclination
         if (k==1) then
            sini = 1.0 ! edge-on
         else
            sini = max(1e-4,abs(sin(input%incl))) ! avoids division by 0 for incl = 0
         end if
      
         ! Convolve v_circ(r) and Sigma(r) into a symmetric, normalised emission line
         s = 0
         do j = 1,nx
            ! The annulus j generates an emission line 1/pi/sqrt(sini^2-(v/v_circ)^2),
            ! which analytically integrates to arcsin(y)/pi, where y = v/v_circ/sini
            y = line%v/v_circ(j)/sini
            dy = dv/v_circ(j)/sini
            f = x(j)*Sigma(j)/pi/dv
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
         line%s = 0
         do i = 1,nv
            do j = -nfilter,nfilter
               ij = i+j
               if (ij<1) ij = 1-ij
               if (ij<=nv) line%s(i) = line%s(i)+s(ij)*filter(j)
            end do
         end do
         
         ! Check normalisation
         if (error_check) then
            !if (abs(sum(line%s*dv)-0.5)>0.01) then
            if (abs((input%mhalo/1.83906473e+14+input%mdisk/1.92865059e+11+input%mbulg/1.36288600e+12)-3.0)<1e-5) then
               call save_profiles('/Users/do/Desktop/right.txt')
               write(*,*) input
               write(*,*) sum(line%s*dv)
               call error('Error in subroutine line_profile. Check arguments.')
            end if
         end if
         
         ! Re-normalise emission lines to correct flux
         line%sintvel = flux(input%mHI/h,input%Dc/h,input%z)  ! [Jy km/s] velocity-integrated flux
         line%s = line%s*line%sintvel
         
         ! Measure characteristic line parameters
         line%speak = maxval(line%s)
         line%scentral = line%s(1)
         do i = 1,nv
            if (line%s(i)==line%speak) then
               line%wpeak = line%v(i)*2
               exit
            end if
         end do
         starget = 0.5*line%speak
         do i = nv-1,1,-1
            if (line%s(i)>starget) then
               f = (line%s(i)-starget)/(line%s(i)-line%s(i+1))
               line%w50 = ((1-f)*line%v(i)+f*line%v(i+1))*2
               exit
            end if
         end do
         starget = 0.2*line%speak
         do i = nv-1,1,-1
            if (line%s(i)>starget) then
               f = (line%s(i)-starget)/(line%s(i)-line%s(i+1))
               line%w20 = ((1-f)*line%v(i)+f*line%v(i+1))*2
               exit
            end if
         end do
         
         ! Make return
         if (k==1) line_edge_on = line
         
      end do
      
      contains
      
      subroutine save_profiles(filename)
      
         implicit none
         character(*),intent(in) :: filename
         integer*4               :: i
      
         open(1,file=trim(filename),action='write',form="formatted",status='replace')
         do i = 1,size(x)
            write(1,'(6Es12.3)') x(i),Sigma(i),sqrt(v_halo_sqr(i)),sqrt(v_disk_sqr(i)),sqrt(v_bulg_sqr(i)),v_circ(i)
         end do
         close(1)
      
      end subroutine save_profiles
       
   end subroutine line_profile
  
   function flux(mass,Dc,z)
   
      ! returns the flux of the emission line as a function of the mass and distance of the emitting material
      ! currently, this is coded up for the HI line, but can be changed to any other line
      
      implicit none
      real,intent(in)   :: mass           ! [Msun WITHOUT h] HI mass
      real,intent(in)   :: Dc             ! [Mpc WITHOUT h] comoving distance
      real,intent(in)   :: z              ! [-] redshift
      real              :: flux           ! [Jy km/s] velocity-integrated line flux (different from frequency-integrated units)
      real              :: L              ! [Lsun] luminosity
      real,parameter    :: fHI = 1.4204   ! [GHz] HI rest-frame frequency
      
      L     = 6.27e-9*mass                ! [Lsun = 3.839e26 W] luminosity = total power of the HI line
      flux  = L/1.04e-3/Dc**2/(1+z)/fHI   ! [Jy km/s] (see eq. A16 in Obreschkow 2009, ApJ 702, 1321)
   
   end function flux
  
   function f_halo(c,x)
      
      ! characteristic halo velocity profile, such that v = sqrt(G*Mhalo/rvir * f_halo)
      
      implicit none
      real,intent(in)   :: c     ! concentration parameter
      real,intent(in)   :: x(:)  ! radius normalised to rvir
      real,allocatable  :: f_halo(:)
      real,allocatable  :: cx(:)
      
      allocate(f_halo(size(x)),cx(size(x)))
      cx = c*x
      f_halo = (log(1+cx)-cx/(1+cx))/x/(log(1+c)-c/(1+c))
      
   end function f_halo
   
   function f_disk(c,x)
      
      ! characteristic disk velocity profile, such that v = sqrt(G*Mdisk/rvir * f_disk)
      
      implicit none
      real,intent(in)   :: c     ! concentration parameter
      real,intent(in)   :: x(:)  ! radius normalised to rvir
      real,allocatable  :: f_disk(:)
      real,allocatable  :: cx(:)
      
      allocate(f_disk(size(x)),cx(size(x)))
      cx = c*x
      f_disk = (c+4.8*c*exp(-0.35*cx-3.5/cx))/(cx+1/cx**2+2/sqrt(cx))
      
   end function f_disk
   
   function f_bulg(c,x)
      
      ! characteristic bulge velocity profile, such that v = sqrt(G*Mbulg/rvir * f_disk)
      
      implicit none
      real,intent(in)   :: c     ! concentration parameter
      real,intent(in)   :: x(:)  ! radius normalised to rvir
      real,allocatable  :: f_bulg(:)
      real,allocatable  :: cx(:)
      
      allocate(f_bulg(size(x)),cx(size(x)))
      cx = c*x
      f_bulg = cx**2*c/(1+cx**2)**1.5
      
   end function f_bulg

end module module_emission_lines