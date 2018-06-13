module module_tiling

   use module_constants
   use module_system
   use module_types
   use module_linalg
   use module_user
   use module_parameters
   
   private
   public   :: make_tiling, is_in_fov
   
   integer*4,allocatable   :: intersection(:,:,:)   ! ==  0 : not checked
                                                    ! == -1 : does not intersect
                                                    ! >= +1 : box id of intersecting box
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
   call check_boxes(nint(starting_point))
   call make_boxes
   call save_box_list
   
   call out('Number of boxes = ',size(tile)*1_8)
   call out('Number of points checked = ',counter*1_8)
   if (allocated(tile)) deallocate(tile)
   call toc
   
end subroutine make_tiling

recursive subroutine check_boxes(ix)

   implicit none
   integer,intent(in) :: ix(3)
   
   if (maxval(abs(ix))>imax) return
   if (intersection(ix(1),ix(2),ix(3)).ne.0) return
   
   if (is_box_in_survey(ix,.false.)) then
   
      if (is_box_in_survey(ix,.true.)) then
         ntile = ntile+1
         intersection(ix(1),ix(2),ix(3)) = ntile
      else
         intersection(ix(1),ix(2),ix(3)) = -1 ! box is only in survey volume specified by parameter file
      end if
   
      ! check neighboring boxes
      call check_boxes(ix+(/+1,0,0/))
      call check_boxes(ix+(/-1,0,0/))
      call check_boxes(ix+(/0,+1,0/))
      call check_boxes(ix+(/0,-1,0/))
      call check_boxes(ix+(/0,0,+1/))
      call check_boxes(ix+(/0,0,-1/))
      
   end if

end subroutine check_boxes

subroutine make_boxes

   implicit none
   real*4                     :: d ! distance from observer to cube centre to in units of box side-lengths
   real*4,parameter           :: h3 = sqrt(3.0)/2.0 ! half the space diagonal of a unit cube
   real*4                     :: rand
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
      
               ! distance range of box
               d = norm(tile(itile)%ix*1.0)
               tile(itile)%dmin = max(para%dc_min/para%L,d-h3)
               tile(itile)%dmax = min(para%dc_max/para%L,d+h3)

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

end subroutine make_boxes

logical function is_box_in_survey(ix,user)

   implicit none
   integer*4,intent(in)    :: ix(3) ! [box side-length] box center in tiling coordinates
   logical,intent(in)      :: user
   integer*4               :: n2D,n3D
   real*4                  :: x(3),dx(3),dy(3),dz(3),sx(3),sy(3),sz(3),fi,fj,d
   integer*4               :: i,j,k
   integer*4,allocatable   :: index(:)
   
   sx = matmul(para%sky_rotation,(/1.0,0.0,0.0/))*para%L
   sy = matmul(para%sky_rotation,(/0.0,1.0,0.0/))*para%L
   sz = matmul(para%sky_rotation,(/0.0,0.0,1.0/))*para%L
   x  = matmul(para%sky_rotation,real(ix,4))*para%L

   d = max(0.5,sqrt(real(sum(ix**2),4))) ! [box side-length] approximate distance from origin to nearest box face
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
         if (is_point_in_survey(x+fi*sx+fj*sy+0.5*sz,user)) then; is_box_in_survey = .true.; return; end if
         if (is_point_in_survey(x+fi*sx+fj*sy-0.5*sz,user)) then; is_box_in_survey = .true.; return; end if
         if (is_point_in_survey(x+fi*sx+0.5*sy+fj*sz,user)) then; is_box_in_survey = .true.; return; end if
         if (is_point_in_survey(x+fi*sx-0.5*sy+fj*sz,user)) then; is_box_in_survey = .true.; return; end if
         if (is_point_in_survey(x+0.5*sx+fi*sy+fj*sz,user)) then; is_box_in_survey = .true.; return; end if
         if (is_point_in_survey(x-0.5*sx+fi*sy+fj*sz,user)) then; is_box_in_survey = .true.; return; end if
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
               if (is_point_in_survey(x+dx+dy+dz,user)) then; is_box_in_survey = .true.; return; end if
            end do
         end do
      end do
   end if
   
   is_box_in_survey = .false. 

end function is_box_in_survey

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