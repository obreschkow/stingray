module module_tiling

   use module_constants
   use module_types
   use module_system
   use module_linalg
   use module_parameters
   use module_user
   
   private
   public   :: make_tiling
   
   integer*4,allocatable  :: intersection(:,:,:)   ! ==  0 : not checked
                                                   ! == -1 : does not intersect
                                                   ! >= +1 : box id of intersecting box
   integer*4            :: imax,nbox,counter
   integer*4,parameter  :: maxlevel_1D = 6 ! edge will be sampled on intervals 1/2^(maxlevel_1D+1)
   integer*4,parameter  :: maxlevel_2D = 5 ! face will be sampled on intervals 1/2^(maxlevel_1D+1)
   integer*4,parameter  :: maxlevel_3D = 0 ! cube will be sampled on intervals 1/2^(maxlevel_3D+1)
   real*4               :: coneaxis(3)
   
contains

subroutine make_tiling

   implicit none
   
   call tic
   call out('MAKE 3D TILING')
   
   call load_parameters
   call set_seed(para%seed)
   nbox = 0
   counter = 0
   coneaxis = (/sin(para%ra)*cos(para%dec),sin(para%dec),cos(para%ra)*cos(para%dec)/)
   imax = ceiling(para%dc_max/para%L)
   if (imax>100) call error('Maximum comoving distance too large for this box side length.')
   allocate(intersection(-imax:imax,-imax:imax,-imax:imax))
   intersection = 0
   nbox = 0
   call check_boxes(nint(para%axis*para%dc_min/para%L))
   call make_boxes
   call save_box_list
   
   call out('Number of boxes = ',size(box)*1_8)
   call out('Number of points checked = ',counter*1_8)
   if (allocated(box)) deallocate(box)
   call toc
   
end subroutine make_tiling

recursive subroutine check_boxes(ix)

   implicit none
   integer,intent(in) :: ix(3)
   
   if (intersection(ix(1),ix(2),ix(3)).ne.0) return
   
   if (is_box_in_survey(ix,.false.)) then
   
      if (is_box_in_survey(ix,.true.)) then
         nbox = nbox+1
         intersection(ix(1),ix(2),ix(3)) = nbox
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
   integer*4                  :: i,j,k,nbox,ibox
   
   nbox = count(intersection>0)
   if (nbox==0) call error('No galaxy passes the selection criterion.')
   
   if (allocated(box)) deallocate(box)
   allocate(box(nbox))
   do i = -imax,imax
      do j = -imax,imax
         do k = -imax,imax
            if (intersection(i,j,k)>0) then
      
               ibox = intersection(i,j,k)
         
               ! position
               box(ibox)%ix = (/i,j,k/)
      
               ! distance range of box
               d = norm(box(ibox)%ix*1.0)
               box(ibox)%dmin = max(para%dc_min/para%L,d-h3)
               box(ibox)%dmax = min(para%dc_max/para%L,d+h3)

               ! choose random proper rotation
               if (para%rotate==1) then
                  call random_number(rand)
                  box(ibox)%rotation = max(1,min(6,ceiling(rand*6.0)))
               else
                  box(ibox)%rotation = 1
               end if
               box(ibox)%Rvector = matmul(para%sky_rotation,rot(:,:,box(ibox)%rotation))

               ! choose random inversion
               if (para%invert==1) then
                  call random_number(rand)
                  if (rand<0.5) box(ibox)%rotation = -box(ibox)%rotation
               end if
               box(ibox)%Rpseudo = matmul(para%sky_rotation,rot(:,:,box(ibox)%rotation))

               ! choose random translation
               if (para%translate==1) then
                  call random_number(box(ibox)%translation)
               else
                  box(ibox)%translation = (/0,0,0/)
               end if
      
            end if
         end do
      end do
   end do

end subroutine make_boxes

logical function is_box_in_survey(ix,user)

   implicit none
   integer*4,intent(in) :: ix(3) ! [box side-length] box center in tiling coordinates
   logical,intent(in)   :: user
   integer*4,parameter  :: n2D = 2**5
   integer*4,parameter  :: n3D = 2**3
   real*4               :: x(3),dx(3),dy(3),dz(3),sx(3),sy(3),sz(3),fi,fj
   integer*4            :: i,j,k
   integer*4            :: index(0:n2D)
   
   sx = matmul(para%sky_rotation,(/1.0,0.0,0.0/))*para%L
   sy = matmul(para%sky_rotation,(/0.0,1.0,0.0/))*para%L
   sz = matmul(para%sky_rotation,(/0.0,0.0,1.0/))*para%L
   x  = matmul(para%sky_rotation,real(ix,4))*para%L
   
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
   
   do i = 0,n3D
      dx = (real(i,4)/n3D-0.5)*sx
      do j = 0,n3D
         dy = (real(j,4)/n3D-0.5)*sy
         do k = 0,n3D
            dz = (real(k,4)/n3D-0.5)*sz
            if (is_point_in_survey(x+dx+dy+dz,user)) then; is_box_in_survey = .true.; return; end if
         end do
      end do
   end do
   
   is_box_in_survey = .false. 

end function is_box_in_survey

logical function is_point_in_survey(x,user)

   ! Check if the point x lies inside the selected survey volume. If user = true, this is the survey volume
   ! specified by the user function position_selection, otherwise it is the survey volume specified in the
   ! parameter file.

   real*4,intent(in)    :: x(3)     ! [box side-length] position of point in tiling coordinates
   logical,intent(in)   :: user
   real*4               :: ra,dec   ! [rad] sky coordinates
   real*4               :: dc       ! [simulation unit] comoving distance
   
   counter = counter+1
   
   dc = norm(x)
   
   if (user) then
   
      if (dc<=epsilon(dc)) then
         ra = 0.0
         dec = 0.0
      else
         ra = modulo(atan2(x(1),x(3)),2*pi)
         dec = asin(min(1.0,x(2)/dc))
      end if
      is_point_in_survey = position_selection(ra/degree,dec/degree,dc)
      
   else
   
      if ((dc<para%dc_min).or.(dc>para%dc_max)) then
         is_point_in_survey = .false.
         return
      end if
      if (dc<=epsilon(dc)) then
         is_point_in_survey = (para%dc_min<=0)
         return
      end if
      alpha = acos(min(1.0,sum(x*coneaxis)/dc))
      is_point_in_survey = (alpha<=para%angle)
      
   end if
   
end function is_point_in_survey
   
end module module_tiling