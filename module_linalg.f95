module module_linalg

   public

contains

subroutine sph2car(dc,ra,dec,x)
   real*4,intent(in)    :: dc,ra,dec
   real*4,intent(out)   :: x(3)
   x = (/sin(ra)*cos(dec),sin(dec),cos(ra)*cos(dec)/)*dc
end subroutine sph2car

subroutine car2sph(x,dc,ra,dec)
   real*4,intent(in)    :: x(3)
   real*4,intent(out)   :: dc,ra,dec
   dc = norm(x)
   if (dc<epsilon(dc)) then
      ra = 0
      dec = 0
   else
      ra = modulo(atan2(x(1),x(3)),6.2831853071795862319959)
      dec = asin(min(1.0,x(2)/dc)) ! "min" to avoid numerical issues if dc very close to x(2)
   end if
end subroutine car2sph

recursive function determinant(mat,n) result(accum)
    integer :: n
    real    :: mat(n, n)
    real    :: submat(n-1, n-1), accum
    integer :: i, sgn
    if ( n == 1 ) then
        accum = mat(1,1)
    else
        accum = 0.0
        sgn = 1
        do i = 1, n
            submat( 1:n-1, 1:i-1 ) = mat( 2:n, 1:i-1 )
            submat( 1:n-1, i:n-1 ) = mat( 2:n, i+1:n )
            accum = accum + sgn * mat(1, i) * determinant( submat, n-1 )
            sgn = - sgn
        enddo
    endif
end function determinant

function cross_product(a,b) result (c)
   
   implicit none
   real*4,intent(in) :: a(3),b(3)
   real*4            :: c(3)
   c(1) = a(2)*b(3)-a(3)*b(2)
   c(2) = a(3)*b(1)-a(1)*b(3)
   c(3) = a(1)*b(2)-a(2)*b(1)
   
end function

function norm(x) result (n)

   implicit none
   real*4,intent(in) :: x(:)
   real*4            :: n
   n = sqrt(sum(x**2))
   
end function norm

function rotation_matrix(axis,angle) result(R)

   implicit none
   real*4,intent(in) :: axis(3)
   real*4,intent(in) :: angle    ! [rad]
   real*4            :: R(3,3)
   real*4            :: c,s
   
   c = cos(angle)
   s = sin(angle)

   R(1,1) = c+axis(1)**2*(1-c)
   R(1,2) = axis(1)*axis(2)*(1-c)-axis(3)*s
   R(1,3) = axis(1)*axis(3)*(1-c)+axis(2)*s
   R(2,1) = axis(2)*axis(1)*(1-c)+axis(3)*s
   R(2,2) = c+axis(2)**2*(1-c)
   R(2,3) = axis(2)*axis(3)*(1-c)-axis(1)*s
   R(3,1) = axis(3)*axis(1)*(1-c)-axis(2)*s
   R(3,2) = axis(3)*axis(2)*(1-c)+axis(1)*s
   R(3,3) = c+axis(3)**2*(1-c)
   
end function rotation_matrix

end module module_linalg