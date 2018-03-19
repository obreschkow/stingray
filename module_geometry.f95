module module_geometry

   use module_constants
   use module_system
   use module_parameters
   
   private
   public   :: make_geometry, load_geometry, make_cone_basis
   
   type(type_box),allocatable :: box(:)
   
contains

subroutine make_geometry(seed)

   implicit none
   integer*4,intent(in) :: seed
   
   call tic
   call out('MAKE BOX GEOMETRY')
   call load_parameters
   call set_seed(seed)
   call add_box(nint(para%axis*para%dc_min/para%L))
   call save_geometry
   call out('Number of boxes = ',size(box)*1_8)
   if (allocated(box)) deallocate(box)
   call toc
   
end subroutine make_geometry

subroutine save_geometry

   implicit none
   character(len=255)   :: filename
   integer*4            :: i
   
   ! write to ascii file
   filename = trim(para%path_output)//'geometry.txt'
   call check_exists(filename)
   open(1,file=trim(filename),action='write',form="formatted",status='replace')
   write(1,'(A)') 'Stingray Box-geometry'
   write(1,'(A)') '--------------------------------------------------------------------------------'
   write(1,'(A)') 'Col  1:  Box index, starting at 1'
   write(1,'(A)') 'Col  2:  x-position of box-centre in units of box side lengths (L)'
   write(1,'(A)') 'Col  3:  y-position of box-centre in units of L'
   write(1,'(A)') 'Col  4:  z-position of box-centre in units of L'
   write(1,'(A)') 'Col  5:  min comoving distance to be considered to fill this box in units of L'
   write(1,'(A)') 'Col  6:  max comoving distance to be considered to fill this box in units of L'
   write(1,'(A)') 'Col  7:  index [1,...,6] of rotation, where 1 is the identity'
   write(1,'(A)') 'Col  8:  x-component of translation vector in units of L'
   write(1,'(A)') 'Col  9:  y-component of translation vector in units of L'
   write(1,'(A)') 'Col 10:  z-component of translation vector in units of L'
   write(1,'(A)') '--------------------------------------------------------------------------------'
   do i = 1,size(box)
      write(1,'(I6,3I6,2F13.6,I3,3F9.5)') i,box(i)
   end do
   close(1)
   
   ! write to binary file
   filename = trim(para%path_output)//'geometry.bin'
   open(1,file=trim(filename),action='write',form='unformatted',status='replace')
   write(1) size(box)
   do i = 1,size(box)
      write(1) box(i)
   end do
   close(1)

end subroutine save_geometry

subroutine load_geometry(box)

   implicit none
   type(type_box),intent(out),allocatable :: box(:)
   character(len=255)                     :: filename
   integer*4                              :: i,nbox
   
   filename = trim(para%path_output)//'geometry.bin'
   call check_exists(filename)
   open(1,file=trim(filename),action='read',form='unformatted')
   read(1) nbox
   if (allocated(box)) deallocate(box)
   allocate(box(nbox))
   do i = 1,nbox
      read(1) box(i)
   end do
   close(1)

end subroutine load_geometry

recursive subroutine add_box(ix)

   implicit none
   integer,intent(in)   :: ix(3)
   type(type_box),allocatable :: tmpbox(:)
   integer*4            :: nbox
   real*4               :: d ! distance from observer to cube centre to in units of box side-lengths
   real*4,parameter     :: h3 = sqrt(3.0)/2.0 ! half the space diagonal of a unit cube
   real*4               :: rand
   
   if (box_exists(ix)) return
   if (.not.box_intersects(ix)) return
   
   ! add new empty element to array box(:)
   if (allocated(box)) then
      nbox = size(box)
      tmpbox = box
      deallocate(box)
      allocate(box(nbox+1))
      box(1:nbox) = tmpbox
      nbox = nbox+1
   else
      nbox = 1
      allocate(box(1))
   end if
   
   ! fill in box properties
   box(nbox)%ix = ix
   d = sqrt(real(sum(ix**2)))
   box(nbox)%dmin = max(para%dc_min/para%L,d-h3)
   box(nbox)%dmax = min(para%dc_max/para%L,d+h3)
   if (para%rotate==1) then
      call random_number(rand)
      box(nbox)%rotation = max(1,min(6,ceiling(rand*6.0)))
   else
      box(nbox)%rotation = 1
   end if
   if (para%invert==1) then
      call random_number(rand)
      if (rand<0.5) box(nbox)%rotation = -box(nbox)%rotation
   end if
   if (para%translate==1) then
      call random_number(box(nbox)%translation)
   else
      box(nbox)%translation = (/0,0,0/)
   end if
   
   ! check neighboring boxes
   call add_box(ix+(/+1,0,0/))
   call add_box(ix+(/-1,0,0/))
   call add_box(ix+(/0,+1,0/))
   call add_box(ix+(/0,-1,0/))
   call add_box(ix+(/0,0,+1/))
   call add_box(ix+(/0,0,-1/))

end subroutine add_box

function box_exists(ix) result(exists)

   implicit none
   integer*4,intent(in) :: ix(3) ! integer box position with ix=(0,0,0) being the box centered on the observer
   logical              :: exists
   integer*4            :: i
   
   if (allocated(box)) then
      exists = .false.
      do i = 1,size(box)
         if (box(i)%ix(1)==ix(1).and.box(i)%ix(2)==ix(2).and.box(i)%ix(3)==ix(3)) exists = .true.
      end do
   else
      exists = .false.
   end if

end function box_exists

function crossproduct(a,b) result(c)
   implicit none
   real*4,intent(in)    :: a(3),b(3)
   real*4               :: c(3)
   c(1) = a(2)*b(3)-a(3)*b(2)
   c(2) = a(3)*b(1)-a(1)*b(3)
   c(3) = a(1)*b(2)-a(2)*b(1)
end function crossproduct

function box_intersects(ix) result(output)

   implicit none
   integer,intent(in)   :: ix(3) ! integer box position with ix=(0,0,0) being the box centered on the observer
   logical              :: output
   real*4               :: dmin,dmax ! minimum and maximum distance in units of box side-lengths
   real*4               :: d ! distance from observer to cube centre to in units of box side-lengths
   real*4,parameter     :: h3 = sqrt(3.0)/2.0 ! half the space diagonal of a unit cube
   real*4               :: alpha,beta ! angles
   real*4               :: a(3),b(3),c(3)
   real*4               :: center(3),e1(3),e2(3) ! center position of and unit vectors on square
   real*4               :: v1,v2,v(3),a1(3),a2(3)
   integer*4            :: i,j1,j2
   integer*4,parameter  :: npoints = 50
   real*4               :: norm,vectnorm,alpha1,alpha2
   
   ! check distance range
   dmin = para%dc_min/para%L
   dmax = para%dc_max/para%L
   d = sqrt(real(sum(ix**2)))
   if ((d<dmin-h3).or.(d>dmax+h3)) then
      output = .false.
      return
   end if
   
   ! check if circumscribed sphere of the box overlaps with cone
   if (h3>=d) then
      output = .true.
      return
   end if
   alpha = acos(sum(para%axis*ix)/d) ! angle between cone axis and line-of-sight to the box center
   beta = asin(h3/d) ! apparent angle of circumscribed sphere of the box
   if (alpha-beta>para%angle*sqrt(2.0)) then ! factor sqrt(2) needed in case of a cone with square base
      output = .false.
      return
   end if
   
   ! check if the cone axis passes through any of the cube's faces
   do i = 1,6

      ! parametrize square face
      select case(i)
         case(1)
            center = (/real(ix(1))+0.5,real(ix(2)),real(ix(3))/)
            e1 = (/0,1,0/)
            e2 = (/0,0,1/)
         case(2)
            center = (/real(ix(1))-0.5,real(ix(2)),real(ix(3))/)
            e1 = (/0,1,0/)
            e2 = (/0,0,1/)
         case(3)
            center = (/real(ix(1)),real(ix(2))+0.5,real(ix(3))/)
            e1 = (/1,0,0/)
            e2 = (/0,0,1/)
         case(4)
            center = (/real(ix(1)),real(ix(2))-0.5,real(ix(3))/)
            e1 = (/1,0,0/)
            e2 = (/0,0,1/)
         case(5)
            center = (/real(ix(1)),real(ix(2)),real(ix(3))+0.5/)
            e1 = (/1,0,0/)
            e2 = (/0,1,0/)
         case(6)
            center = (/real(ix(1)),real(ix(2)),real(ix(3))-0.5/)
            e1 = (/1,0,0/)
            e2 = (/0,1,0/)
      end select
   
      ! check if cone axis passes through the particular square i 
      v = crossproduct(e1,e2)
      if (abs(sum(para%axis*v))>1e-10) then ! to avoid a,e1,e2 being coplanar
         a = crossproduct(para%axis,center)
         b = crossproduct(para%axis,e1)
         c = crossproduct(para%axis,e2)
         v1 = (c(1)*a(2)-c(2)*a(1))/(b(1)*c(2)-b(2)*c(1))
         v2 = (a(1)*b(2)-a(2)*b(1))/(b(1)*c(2)-b(2)*c(1))
         if ((abs(v1)<=0.50001).and.(abs(v2)<=0.50001)) then
            a = center+e1*v1+e2*v2
            vectnorm = sqrt(sum(a**2))
            if ((vectnorm>=dmin).and.(vectnorm<=dmax)) then
               output = .true.
               return
            end if
         end if
      end if
   
   end do
   
   ! make orthogonal unit vectors, perpendicular to the cone axis
   if (para%square_base==1) call make_cone_basis(a1,a2)
   
   ! since the cone axis did not pass through any cube face, check if any point on the face lines inside the cone
   
   do i = 1,6
   
      ! parametrize square face
      select case(i)
         case(1)
            center = (/real(ix(1))+0.5,real(ix(2)),real(ix(3))/)
            e1 = (/0,1,0/)
            e2 = (/0,0,1/)
         case(2)
            center = (/real(ix(1))-0.5,real(ix(2)),real(ix(3))/)
            e1 = (/0,1,0/)
            e2 = (/0,0,1/)
         case(3)
            center = (/real(ix(1)),real(ix(2))+0.5,real(ix(3))/)
            e1 = (/1,0,0/)
            e2 = (/0,0,1/)
         case(4)
            center = (/real(ix(1)),real(ix(2))-0.5,real(ix(3))/)
            e1 = (/1,0,0/)
            e2 = (/0,0,1/)
         case(5)
            center = (/real(ix(1)),real(ix(2)),real(ix(3))+0.5/)
            e1 = (/1,0,0/)
            e2 = (/0,1,0/)
         case(6)
            center = (/real(ix(1)),real(ix(2)),real(ix(3))-0.5/)
            e1 = (/1,0,0/)
            e2 = (/0,1,0/)
      end select
      
      ! check if any point on face i lies inside the cone
      do j1 = -npoints,npoints
         do j2 = -npoints,npoints
            a = center+0.5*real(j1)/real(npoints)*e1+0.5*real(j2)/real(npoints)*e2
            norm = sqrt(sum(a**2))
            if (para%square_base==1) then
               alpha1 = abs(atan2(sum(a*a1),sum(a*para%axis)))
               alpha2 = abs(atan2(sum(a*a2),sum(a*para%axis)))
               alpha = max(alpha1,alpha2)
            else
               alpha = acos(min(1.0,sum(a*para%axis)/norm))
            end if
            if (alpha<=para%angle) then
               if ((norm>=dmin).and.(norm<=dmax)) then
                  output = .true.
                  return
               end if
            end if
         end do
      end do
   
   end do
    
   output = .false.
   return
   
end function box_intersects

subroutine make_cone_basis(a1,a2)
   implicit none
   real*4,intent(out)   :: a1(3),a2(3) ! two orthonormal basis vectors, orthogonal to the cone axis
   if (para%axis(3)>0.5) then
      a1 = crossproduct(para%axis,(/0.0,1.0,0.0/))
      a1 = a1/sqrt(sum(a1**2))
      a2 = crossproduct(para%axis,a1)
      a2 = a2/sqrt(sum(a2**2))
   else
      a1 = crossproduct(para%axis,(/0.0,0.0,1.0/))
      a1 = a1/sqrt(sum(a1**2))
      a2 = crossproduct(para%axis,a1)
      a2 = a2/sqrt(sum(a2**2))
   end if
end subroutine make_cone_basis
   
end module module_geometry