module module_tiling

   use module_constants
   use module_system
   use module_types
   use module_io
   use module_linalg
   use module_user_selection
   use module_parameters
   
   private
   public   :: make_tiling, is_in_fov
   
   integer*4,allocatable   :: intersection(:,:,:)   ! ==  0 : not checked
                                                    ! == -1 : does not intersect
                                                    ! >= +1 : id of intersecting tile
   integer*4               :: imax,ntile,counter
   
contains

subroutine make_tiling

   implicit none
   real*4   :: starting_point(3)
   
   call tic
   call out('MAKE 3D TILING')
   
   call load_parameters
   call set_seed(para%seed)
   ntile = 0
   counter = 0
   imax = ceiling(para%dc_max/para%L)
   if (imax>100) call error('Maximum comoving distance too large for this box side length.')
   allocate(intersection(-imax:imax,-imax:imax,-imax:imax))
   intersection = 0
   ntile = 0
   call sph2car(para%dc_min/para%L,(para%ra_min+para%ra_max)/2.0,(para%dec_min+para%dec_max)/2.0,starting_point)
   starting_point = matmul(starting_point,para%sky_rotation)
   call check_tile(nint(starting_point))
   call make_tile_list
   call save_tile_list
   
   call out('Number of tiles = ',size(tile)*1_8)
   call out('Number of points checked = ',counter*1_8)
   if (allocated(tile)) deallocate(tile)
   call toc

end subroutine make_tiling

recursive subroutine check_tile(ix)

   implicit none
   integer,intent(in) :: ix(3)
   
   if (maxval(abs(ix))>imax) return
   if (intersection(ix(1),ix(2),ix(3)).ne.0) return
   
   if (is_tile_in_survey(ix,.false.)) then
      
      ! tile intersects with survey volume specified by parameter file
   
      if (is_tile_in_survey(ix,.true.)) then
         ! tile also intersects with survey volume specified by user file
         ntile = ntile+1
         intersection(ix(1),ix(2),ix(3)) = ntile
      else
         ! tile does not intersect with survey volume specified by user file
         intersection(ix(1),ix(2),ix(3)) = -1
      end if
   
      ! check neighboring tiles
      call check_tile(ix+(/+1,0,0/))
      call check_tile(ix+(/-1,0,0/))
      call check_tile(ix+(/0,+1,0/))
      call check_tile(ix+(/0,-1,0/))
      call check_tile(ix+(/0,0,+1/))
      call check_tile(ix+(/0,0,-1/))
      
   else
   
      ! tile does not intersect with survey volume specified by parameter file
      intersection(ix(1),ix(2),ix(3)) = -1
      
   end if

end subroutine check_tile

subroutine make_tile_list

   implicit none
   real*4                     :: rand,dmin,dmax
   integer*4                  :: i,j,k,ntile,itile
   
   ntile = count(intersection>0)
   if (ntile==0) call error('No galaxy passes the selection criterion.')
   
   if (allocated(tile)) deallocate(tile)
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
               tile(itile)%dmin = max(para%dc_min/para%L,dmin)
               tile(itile)%dmax = min(para%dc_max/para%L,dmax)
               
               ! choose random proper rotation
               if (para%rotate==1) then
                  call random_number(rand)
                  tile(itile)%rotation = max(1,min(6,ceiling(rand*6.0)))
               else
                  tile(itile)%rotation = 1
               end if
               tile(itile)%Rvector = matmul(para%sky_rotation,rot(:,:,tile(itile)%rotation))

               ! choose random inversion
               if (para%invert==1) then
                  call random_number(rand)
                  if (rand<0.5) tile(itile)%rotation = -tile(itile)%rotation
               end if
               tile(itile)%Rpseudo = matmul(para%sky_rotation,rot(:,:,tile(itile)%rotation))

               ! choose random translation
               if (para%translate==1) then
                  call random_number(tile(itile)%translation)
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
      else if (abs(p2)>0.5) then
         dmin = get_min_distance_to_line(p-s2*e2/2,e1)
      end if
      return
   end if
   
   ! handle the situation where the closest point is a corner
   if (min(abs(p1),abs(p2))>0.5) then
      dmin = norm(p-s1*e1/2-s2*e2/2)
      return
   else
      write(*,*) 'Error in get_min_distance_to_square'
      stop
   end if
   
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
   
   sx = matmul(para%sky_rotation,(/1.0,0.0,0.0/))*para%L
   sy = matmul(para%sky_rotation,(/0.0,1.0,0.0/))*para%L
   sz = matmul(para%sky_rotation,(/0.0,0.0,1.0/))*para%L
   x  = matmul(para%sky_rotation,real(ix,4))*para%L

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
   
   if (para%volume_search_level>0) then
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
   end if
   
   is_tile_in_survey = .false. 

end function is_tile_in_survey

logical function is_point_in_survey(x,user)

   ! Check if the point x lies inside the selected survey volume. If user = true, this is the survey volume
   ! specified by the user function pos_selection, otherwise it is the survey volume specified in the
   ! parameter file.

   real*4,intent(in)    :: x(3)     ! [simulation unit] position of point in tiling coordinates
   logical,intent(in)   :: user
   real*4               :: ra,dec   ! [rad] sky coordinates
   real*4               :: dc       ! [simulation unit] comoving distance
   
   counter = counter+1
   
   call car2sph(x,dc,ra,dec)
   
   if (user) then
   
      is_point_in_survey = pos_selection(dc,ra/degree,dec/degree)
      
   else
   
      if (dc<=epsilon(dc)) then
         is_point_in_survey = (para%dc_min<=0)
         return
      end if
      is_point_in_survey = is_in_fov(dc,ra,dec)
      
   end if
   
end function is_point_in_survey

logical function is_in_fov(dc,ra,dec)
   implicit none
   real*4,intent(in) :: dc       ! [simulation units]
   real*4,intent(in) :: ra,dec   ! [rad]
   is_in_fov = (dc>=para%dc_min).and.(dc<=para%dc_max).and. &
             & (ra>=para%ra_min).and.(ra<=para%ra_max).and. &
             & (dec>=para%dec_min).and.(dec<=para%dec_max)
end function is_in_fov           
             
end module module_tiling