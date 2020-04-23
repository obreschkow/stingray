module module_tiling

   use shared_module_core
   use shared_module_parameters
   use shared_module_maths
   use shared_module_constants
   use module_global
   use module_user_selection
   
   public   :: make_tiling
   public   :: is_in_fov
   
   private
   integer*4,allocatable   :: intersection(:,:,:)   ! ==  0 : not checked
                                                    ! == -1 : does not intersect
                                                    ! >= +1 : id of intersecting tile
   integer*4               :: imax,ntile,counter
   
contains

subroutine make_tiling

   implicit none
   integer*4   :: starting_point(3)
   
   call tic
   call out('MAKE 3D TILING')
   
   call set_seed(para%seed)
   call get_position_range
   
   ntile = 0
   counter = 0
   imax = ceiling(para%dc_max/para%box_side)
   if (abs(imax)>limit%n_tiles_max) call error('maximum comoving distance too large for this box side length')
   allocate(intersection(-imax:imax,-imax:imax,-imax:imax))
   intersection = 0
   ntile = 0
   call make_starting_point(starting_point)
   call check_tile(starting_point,0)
   call make_tile_list
   
   call out('Number of tiles = ',size(tile))
   call out('Number of points checked = ',counter)
   
   call toc

end subroutine make_tiling

logical function is_in_fov(pos)

   implicit none
   
   type(type_pos),intent(in)  :: pos
   
   is_in_fov = (pos%dc>=para%dc_min).and.(pos%dc<=para%dc_max).and. &
             & (pos%ra>=para%ra_min).and.(pos%ra<=para%ra_max).and. &
             & (pos%dec>=para%dec_min).and.(pos%dec<=para%dec_max)
             
end function is_in_fov

subroutine get_position_range

   implicit none
   type(type_range)  :: range
   real*4,parameter  :: x = huge(0.0_4)
   
   ! set default values
   range%dc = x
   range%ra = x
   range%dec = x
   
   ! extract ranges from selection function
   call selection_function(range=range)
   
   ! check initialization
   if (any(range%dc==x)) call error('range of comoving distances (range%dc) not provided in selection function')
   if (any(range%ra==x)) call error('range of right ascension (range%ra) not provided in selection function')
   if (any(range%dec==x)) call error('range of declination (range%dec) not provided in selection function')
   
   ! cast into parameters
   para%dc_min = range%dc(1)
   para%dc_max = range%dc(2)
   para%ra_min = range%ra(1)
   para%ra_max = range%ra(2)
   para%dec_min = range%dec(1)
   para%dec_max = range%dec(2)
   
   ! check values
   if (para%dc_min<0.0) call error('dc_min must be >=0')
   if (para%dc_max<=0.0) call error('dc_max must be >0')
   if (para%dc_max<=para%dc_min) call error('dc_min must be smaller than dc_max')
   
   if (para%ra_min<0.0) call error('ra_min must be >=0')
   if (para%ra_max>360.0) call error('ra_max must be <=360')
   if (para%ra_min>=para%ra_max) call error('ra_max must be larger than ra_min')
   
   if (para%dec_min<-90.0) call error('ra_min must be >=-90')
   if (para%dec_max>+90.0) call error('ra_max must be <=90')
   if (para%dec_min>=para%dec_max) call error('dec_max must be larger than dec_min')
   
   ! convert degrees to radian
   para%ra_min = para%ra_min*unit%degree
   para%ra_max = para%ra_max*unit%degree
   para%dec_min = para%dec_min*unit%degree
   para%dec_max = para%dec_max*unit%degree

end subroutine get_position_range

subroutine make_starting_point(ix)

   implicit none
   integer*4,intent(out)   :: ix(3)
   integer*4               :: i,j,k
   
   ix = (/0,0,0/)
   if (is_tile_in_survey(ix,.false.).and.is_tile_in_survey(ix,.true.)) return
   
   do i = -imax,imax
      do j = -imax,imax
         do k = -imax,imax
            ix = (/i,j,k/)
            if (is_tile_in_survey(ix,.false.).and.is_tile_in_survey(ix,.true.)) return
         end do
      end do
   end do
      
   call error('maximum comoving distance too large for this box side length')
   
end subroutine make_starting_point

recursive subroutine check_tile(ix,nrecursion)

   implicit none
   integer,intent(in) :: ix(3),nrecursion
   
   if (nrecursion>1e4) call error('too many recursions in tiling; please restrict the ranges of dc, ra, dec in selection_function')
   
   if (maxval(abs(ix))>imax) return
   if (intersection(ix(1),ix(2),ix(3)).ne.0) return
   
   if (is_tile_in_survey(ix,.false.)) then
      
      ! tile intersects with survey volume specified by parameter file
   
      if (is_tile_in_survey(ix,.true.)) then
         ! tile also intersects with survey volume specified by user file
         ntile = ntile+1
         if (ntile>limit%n_tiles_max) call error ('number of required tiles exceeds ',limit%n_tiles_max)
         intersection(ix(1),ix(2),ix(3)) = ntile
      else
         ! tile does not intersect with survey volume specified by user file
         intersection(ix(1),ix(2),ix(3)) = -1
      end if
      
      ! check neighboring tiles
      call check_tile(ix+(/+1,0,0/),nrecursion+1)
      call check_tile(ix+(/-1,0,0/),nrecursion+1)
      call check_tile(ix+(/0,+1,0/),nrecursion+1)
      call check_tile(ix+(/0,-1,0/),nrecursion+1)
      call check_tile(ix+(/0,0,+1/),nrecursion+1)
      call check_tile(ix+(/0,0,-1/),nrecursion+1)
      
   else
   
      ! tile does not intersect with survey volume specified by parameter file
      intersection(ix(1),ix(2),ix(3)) = -1
      
   end if

end subroutine check_tile

subroutine make_tile_list

   ! NB: random_number() uses a better prng than rand(), but depends on the system and fortran version
   !     we therefore chose to use rand()

   implicit none
   real*4                     :: dmin,dmax
   integer*4                  :: i,j,k,ntile,itile
   
   ntile = count(intersection>0)
   if (ntile==0) call error('no galaxy passes the selection criterion.')
   
   allocate(tile(ntile))
   do i = -imax,imax
      do j = -imax,imax
         do k = -imax,imax
            if (intersection(i,j,k)>0) then
      
               itile = intersection(i,j,k)
         
               ! position
               tile(itile)%ix = (/i,j,k/)
      
               ! distance range of tile
               call get_distance_range_of_cube(tile(itile)%ix,dmin,dmax)
               tile(itile)%dmin = max(para%dc_min/para%box_side,dmin)
               tile(itile)%dmax = min(para%dc_max/para%box_side,dmax)
               
               ! choose random proper rotation
               if (para%rotate==1) then
                  !call random_number(rnd)
                  tile(itile)%rotation = max(1,min(6,ceiling(rand()*6.0)))
               else
                  tile(itile)%rotation = 1
               end if
               tile(itile)%Rvector = matmul(para%sky_rotation,rot(:,:,tile(itile)%rotation))

               ! choose random inversion
               if (para%invert==1) then
                  !call random_number(rnd)
                  if (rand()<0.5) tile(itile)%rotation = -tile(itile)%rotation
               end if
               tile(itile)%Rpseudo = matmul(para%sky_rotation,rot(:,:,tile(itile)%rotation))

               ! choose random translation
               if (para%translate==1) then
                  !call random_number(tile(itile)%translation)
                  tile(itile)%translation = rand()
               else
                  tile(itile)%translation = (/0,0,0/)
               end if
      
            end if
         end do
      end do
   end do

end subroutine make_tile_list

subroutine get_distance_range_of_cube(ix,dmin,dmax)

   implicit none
   integer*4,intent(in) :: ix(3)
   real*4,intent(out)   :: dmin,dmax
   real*4               :: dmin_sq(6)
   integer*4            :: i,j,k
   
   if (maxval(abs(ix))==0) then
      dmin = 0
   else
      dmin_sq(1) = get_min_distance_to_square(ix+(/-0.5,0.0,0.0/),(/0.0,1.0,0.0/),(/0.0,0.0,1.0/))
      dmin_sq(2) = get_min_distance_to_square(ix+(/+0.5,0.0,0.0/),(/0.0,1.0,0.0/),(/0.0,0.0,1.0/))
      dmin_sq(3) = get_min_distance_to_square(ix+(/0.5,-0.5,0.0/),(/1.0,0.0,0.0/),(/0.0,0.0,1.0/))
      dmin_sq(4) = get_min_distance_to_square(ix+(/0.5,+0.5,0.0/),(/1.0,0.0,0.0/),(/0.0,0.0,1.0/))
      dmin_sq(5) = get_min_distance_to_square(ix+(/0.0,0.0,-0.5/),(/1.0,0.0,0.0/),(/0.0,1.0,0.0/))
      dmin_sq(6) = get_min_distance_to_square(ix+(/0.0,0.0,+0.5/),(/1.0,0.0,0.0/),(/0.0,1.0,0.0/))
      dmin = minval(dmin_sq)
   end if
   
   dmax = 0
   do i = -1,1,2
      do j = -1,1,2
         do k = -1,1,2
            dmax = max(dmax,norm(real(ix)+(/i*0.5,j*0.5,k*0.5/)))
         end do
      end do
   end do
               
end subroutine get_distance_range_of_cube

function get_min_distance_to_square(p,e1,e2) result(dmin)
   
   implicit none
   real*4,intent(in)    :: p(3)        ! center of square
   real*4,intent(in)    :: e1(3),e2(3) ! orthonormal basis in the plane of the square, parallel to the edges
   real*4               :: dmin        ! minimal distance from the origin (0,0,0) to the square
   real*4               :: p1,p2       ! projections of p onto square
   real*4               :: s1,s2
   
   ! project p onto square
   p1 = sum(p*e1)
   p2 = sum(p*e2)
   s1 = sign(1.0,p1) ! sign of p1
   s2 = sign(1.0,p2) ! sign of p2
   
   ! handle the situation where projection lands inside the square => closest point lies inside square
   if (max(abs(p1),abs(p2))<=0.5) then
      dmin = sqrt(sum(p**2)-sum(p*e1)**2-sum(p*e2)**2) ! equation for the minimum distance to a plane
      return
   end if
   
   ! handle the situation where the closes point lies along an edge
   if ((min(abs(p1),abs(p2))<=0.5).and.(max(abs(p1),abs(p2))>0.5)) then
      if (abs(p1)>0.5) then
         dmin = get_min_distance_to_line(p-s1*e1/2,e2)
         return
      else if (abs(p2)>0.5) then
         dmin = get_min_distance_to_line(p-s2*e2/2,e1)
         return
      end if
   end if
   
   ! handle the situation where the closest point is a corner
   if (min(abs(p1),abs(p2))>0.5) then
      dmin = norm(p-s1*e1/2-s2*e2/2)
      return
   else
      call deverror('unknown erro in get_min_distance_to_square')
   end if
   
   ! code never gets here, this is just to avoid compiler warning of uninitialised variable
   dmin = 0.0 
   return
   
end function get_min_distance_to_square

function get_min_distance_to_line(p,e) result(dmin)

   implicit none
   real*4,intent(in)    :: p(3)  ! position on the line
   real*4,intent(in)    :: e(3)  ! unit vector parallel to the line
   real*4               :: dmin  ! minimal distance from the origin (0,0,0) to the line
   
   dmin = sqrt(sum(p**2)-sum(p*e)**2)
   
end function get_min_distance_to_line

logical function is_tile_in_survey(ix,user)

   implicit none
   integer*4,intent(in)    :: ix(3) ! [box side-length] tile center in tiling coordinates
   logical,intent(in)      :: user
   integer*4               :: n2D,n3D
   real*4                  :: x(3),dx(3),dy(3),dz(3),sx(3),sy(3),sz(3),fi,fj,d
   integer*4               :: i,j,k
   integer*4,allocatable   :: index(:)
   
   sx = matmul(para%sky_rotation,(/1.0,0.0,0.0/))*para%box_side
   sy = matmul(para%sky_rotation,(/0.0,1.0,0.0/))*para%box_side
   sz = matmul(para%sky_rotation,(/0.0,0.0,1.0/))*para%box_side
   x  = matmul(para%sky_rotation,real(ix,4))*para%box_side

   d = max(0.5,sqrt(real(sum(ix**2),4))) ! [box side-length] approximate distance from origin to nearest tile face
   n2D = max(10,2*nint(0.5/(d*para%search_angle)))
   allocate(index(0:n2D))
   do i = 0,n2D/2-1
      index(i*2) = i
      index(i*2+1) = n2D-i
   end do
   index(n2D) = n2D/2

   do i = 0,n2D
      fi = real(index(i),4)/n2D-0.5
      do j = 0,n2D
         fj = real(index(j),4)/n2D-0.5
         if (is_point_in_survey(x+fi*sx+fj*sy+0.5*sz,user)) then; is_tile_in_survey = .true.; return; end if
         if (is_point_in_survey(x+fi*sx+fj*sy-0.5*sz,user)) then; is_tile_in_survey = .true.; return; end if
         if (is_point_in_survey(x+fi*sx+0.5*sy+fj*sz,user)) then; is_tile_in_survey = .true.; return; end if
         if (is_point_in_survey(x+fi*sx-0.5*sy+fj*sz,user)) then; is_tile_in_survey = .true.; return; end if
         if (is_point_in_survey(x+0.5*sx+fi*sy+fj*sz,user)) then; is_tile_in_survey = .true.; return; end if
         if (is_point_in_survey(x-0.5*sx+fi*sy+fj*sz,user)) then; is_tile_in_survey = .true.; return; end if
      end do
   end do
   
   n3D = int(2**para%volume_search_level-1,4)
   do i = 0,n3D
      dx = ((real(i,4)+0.5)/(n3D+1)-0.5)*sx
      do j = 0,n3D
         dy = ((real(j,4)+0.5)/(n3D+1)-0.5)*sy
         do k = 0,n3D
            dz = ((real(k,4)+0.5)/(n3D+1)-0.5)*sz
            if (is_point_in_survey(x+dx+dy+dz,user)) then; is_tile_in_survey = .true.; return; end if
         end do
      end do
   end do
   
   is_tile_in_survey = .false. 

end function is_tile_in_survey

logical function is_point_in_survey(x,user)

   ! Check if the point x lies inside the selected survey volume. If user = true, this is the survey volume
   ! specified by the user function selection_function, otherwise it is the survey volume specified in the
   ! parameter file.

   implicit none
   real*4,intent(in)    :: x(3)     ! [simulation unit] position of point in tiling coordinates
   logical,intent(in)   :: user
   real*4               :: ra,dec   ! [rad] sky coordinates
   real*4               :: dc       ! [simulation unit] comoving distance
   type(type_pos)       :: pos
   
   counter = counter+1
   
   call car2sph(x,dc,ra,dec,astro=.true.)
   
   if (user) then
   
      pos%dc = dc
      pos%ra = ra/unit%degree
      pos%dec = dec/unit%degree
      is_point_in_survey = .true.
      call selection_function(pos=pos,selected=is_point_in_survey)

   else
   
      if (dc<=epsilon(dc)) then
         is_point_in_survey = (para%dc_min<=0)
         return
      end if
      pos%dc = dc
      pos%ra = ra
      pos%dec = dec
      is_point_in_survey = is_in_fov(pos)
      
   end if
   
end function is_point_in_survey       
             
end module module_tiling