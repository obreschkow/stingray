   ! output basic statistics
   write(txt,'(A,F0.4,A,F0.4)') 'Position range: ',minval((/minval(base(:)%xbox(1)),&
   & minval(base(:)%xbox(2)),minval(base(:)%xbox(3))/)),&
   & ' to ',maxval((/maxval(base(:)%xbox(1)),maxval(base(:)%xbox(2)),maxval(base(:)%xbox(3))/))
   call out(txt)
   write(txt,'(A,I0,A,I0)')     'Identifier range: ',minval(base%id),' to ',maxval(base%id)
   call out(txt)

function box_intersects_simple(ix) result(output)

   implicit none
   integer,intent(in)   :: ix(3) ! integer box position with ix=(0,0,0) being the box centered on the observer
   logical              :: output
   real*4               :: dmin,dmax ! minimum and maximum distance in units of box side-lengths
   real*4               :: d ! distance from observer to cube centre to in units of box side-lengths
   real*4,parameter     :: h3 = sqrt(3.0)/2.0 ! half the space diagonal of a unit cube
   real*4               :: alpha,beta ! angles
   
   ! check if this is the central cube
   if ((ix(1)==0).and.(ix(2)==0).and.(ix(3)==0)) then
      output = .true.
      return
   end if
   
   ! check distance range
   dmin = para%dcmin/para%L
   dmax = para%dcmax/para%L
   d = sqrt(real(sum(ix**2)))
   if ((d<dmin-h3).or.(d>dmax+h3)) then
      output = .false.
      return
   end if
   
   ! check if circumscribed sphere of the box overlaps with cone
   alpha = acos(sum(para%axis*ix)/d) ! angle between cone axis and line-of-sight to the box center
   beta = asin(h3/d) ! apparent angle of circumscribed sphere of the box
   output = alpha-beta<para%angle
   return
   
end function box_intersects_simple