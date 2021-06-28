! **********************************************************************************************************************************
! Shared Fortran module facilitating the handling of vectors
! Developed by Danail Obreschkow
! See subroutine example for an illustration of key functionalities
! **********************************************************************************************************************************

module shared_module_vectors

   private
     
   public :: vector4             ! 3-vector of real*4 with components x,y,z
   public :: vector8             ! 3-vector of real*8 with components x,y,z
   public :: assignment(=)       ! assign component vector to vector class and vice versa
   public :: operator(==)        ! compare all components of vectors of identical type
   public :: operator(.ne.)      ! compare all components of vectors of identical type (opposite answer of ==)
   public :: operator(+)         ! vector addition
   public :: operator(-)         ! vector subtraction
   public :: operator(*)         ! scaling of vector elements, element-by-element multiplication
   public :: operator(/)         ! scaling of vector elements, element-by-element division
   public :: operator(.dot.)     ! scalar product of two vectors
   public :: operator(.cross.)   ! cross product of two 3-vectors
   public :: components          ! function converting a vector into its components
   public :: component           ! function extracting a single component from a vector
   public :: norm                ! function returning the norm of a vector
   public :: unitvector          ! function scaling a vector to a unit vector with the same direction
   public :: scalar_product      ! scalar_product of two vectors (can also use "u.dot.v")
   public :: cross_product       ! cross_product of two 3-vectors (can also use "u.cross.v")
   public :: matmul              ! multiply vector with a square matrix (extends on intrinsic subroutine)
   public :: example_vectors     ! subroutine showing some example vector calculations
     
   type vector4
      real*4 :: x = 0.0
      real*4 :: y = 0.0
      real*4 :: z = 0.0
   contains
      procedure :: norm => norm_vector4
      procedure :: unit => unit_vector4
   end type vector4
   
   type vector8
      real*8 :: x = 0.0_8
      real*8 :: y = 0.0_8
      real*8 :: z = 0.0_8
   contains
      procedure :: norm => norm_vector8
      procedure :: unit => unit_vector8
   end type vector8
   
   interface assignment(=)
      procedure getcomponents_vector4
      procedure getcomponents_vector8
      procedure setcomponents_vector4_from_real4
      procedure setcomponents_vector4_from_real8
      procedure setcomponents_vector4_from_int4
      procedure setcomponents_vector4_from_int8
      procedure setcomponents_vector8_from_real4
      procedure setcomponents_vector8_from_real8
      procedure setcomponents_vector8_from_int4
      procedure setcomponents_vector8_from_int8
   end interface assignment(=)
   
   interface operator(==)
      procedure compare_vectors4
      procedure compare_vectors8
   end interface operator(==)
     
   interface operator(.ne.)
      procedure ncompare_vectors4
      procedure ncompare_vectors8
   end interface operator(.ne.)
   
   interface operator(+)
      procedure add_vectors4
      procedure add_vectors8
   end interface operator(+)
     
   interface operator(-)
      procedure subtract_two_vectors4
      procedure subtract_two_vectors8
   end interface operator(-)
     
   interface operator(*)
      procedure multiply_vector4_by_vector4
      procedure multiply_vector8_by_vector8
      procedure multiply_real4_with_vector4
      procedure multiply_real8_with_vector4
      procedure multiply_int4_with_vector4
      procedure multiply_int8_with_vector4
      procedure multiply_vector4_with_real4
      procedure multiply_vector4_with_real8
      procedure multiply_vector4_with_int4
      procedure multiply_vector4_with_int8
      procedure multiply_real4_with_vector8
      procedure multiply_real8_with_vector8
      procedure multiply_int4_with_vector8
      procedure multiply_int8_with_vector8
      procedure multiply_vector8_with_real4
      procedure multiply_vector8_with_real8
      procedure multiply_vector8_with_int4
      procedure multiply_vector8_with_int8
   end interface operator(*)
     
   interface operator(/)
      procedure divide_vector4_by_vector4
      procedure divide_vector8_by_vector8
      procedure divide_vector4_by_real4
      procedure divide_vector4_by_real8
      procedure divide_vector4_by_int4
      procedure divide_vector4_by_int8
      procedure divide_vector8_by_real4
      procedure divide_vector8_by_real8
      procedure divide_vector8_by_int4
      procedure divide_vector8_by_int8
   end interface operator(/)
     
   interface operator(.dot.)
      procedure scalar_product_vector4
      procedure scalar_product_vector8
      procedure scalar_product_array4
      procedure scalar_product_array8
      procedure scalar_product_a4v4
      procedure scalar_product_v4a4
      procedure scalar_product_a8v8
      procedure scalar_product_v8a8
   end interface operator(.dot.)
     
   interface operator(.cross.)
      procedure cross_product_vector4
      procedure cross_product_vector8
      procedure cross_product_array4
      procedure cross_product_array8
      procedure cross_product_a4v4
      procedure cross_product_v4a4
      procedure cross_product_a8v8
      procedure cross_product_v8a8
   end interface operator(.cross.)
   
   interface norm
      procedure norm_vector4
      procedure norm_vector8
      procedure norm_array4
      procedure norm_array8
   end interface norm   
   
   interface unitvector
      procedure unit_vector4
      procedure unit_vector8
      procedure unit_array4
      procedure unit_array8
   end interface unitvector
   
   interface scalar_product
      procedure scalar_product_vector4
      procedure scalar_product_vector8
      procedure scalar_product_array4
      procedure scalar_product_array8
      procedure scalar_product_a4v4
      procedure scalar_product_v4a4
      procedure scalar_product_a8v8
      procedure scalar_product_v8a8
   end interface scalar_product
   
   interface cross_product
      procedure cross_product_vector4
      procedure cross_product_vector8
      procedure cross_product_array4
      procedure cross_product_array8
      procedure cross_product_a4v4
      procedure cross_product_v4a4
      procedure cross_product_a8v8
      procedure cross_product_v8a8
   end interface cross_product
     
   interface matmul
      procedure matmul_matrix4_times_vector4
      procedure matmul_matrix4_times_vector8
      procedure matmul_matrix8_times_vector4
      procedure matmul_matrix8_times_vector8
      procedure matmul_vector4_times_matrix4
      procedure matmul_vector8_times_matrix4
      procedure matmul_vector4_times_matrix8
      procedure matmul_vector8_times_matrix8
   end interface matmul

   interface components
      procedure components_vector4
      procedure components_vector8
   end interface components
   
   interface component
      procedure component_vector4
      procedure component_vector8
   end interface component
     
contains

   subroutine example_vectors
   
      implicit none
      type(vector4) :: x,y ! real*4 3-vectors
      real*4 :: m(3,3),q(3),r(3)
      
      write(*,*) 'Assign vectors and print components'
      q = (/-0.5,6.1,3.2/)    ! representation of vector as standard 3-element array
      x = q                   ! make a class vector4 object from components
      write(*,*) q
      write(*,*) x
      
      write(*,*) 'Multiply vector by a constant'
      write(*,*) 2*q
      write(*,*) 2*x
      
      write(*,*) 'Divide vector by a constant'
      write(*,*) q/2
      write(*,*) x/2
      
      write(*,*) 'Vector norm'
      write(*,*) norm(q)
      write(*,*) norm(x)
      write(*,*) x%norm()
      
      write(*,*) 'Unit vectors'
      write(*,*) unitvector(q)
      write(*,*) unitvector(x)
      write(*,*) x%unit()
      
      write(*,*) 'Component-wise multiplication'
      r = (/4.3,4.1,0.5/)
      y = r
      write(*,*) q*r
      write(*,*) x*y
      
      write(*,*) 'Component-wise division'
      write(*,*) q/r
      write(*,*) x/y
      
      write(*,*) 'Scalar product'
      write(*,*) q.dot.r
      write(*,*) x.dot.y
      
      write(*,*) 'Vector product'
      write(*,*) q.cross.r
      write(*,*) x.cross.y
      
      write(*,*) 'Multiply matrix with vector'
      m = reshape((/1.3,4.2,-2.0,5.1,-2.2,1.0,-0.6,0.8,4.0/),(/3,3/))
      write(*,*) matmul(m,q)
      write(*,*) matmul(m,x)
      
      write(*,*) 'Multiply vector with matrix'
      write(*,*) matmul(q,m)
      write(*,*) matmul(x,m)
      
      write(*,*) 'Compare vectors'
      write(*,*) all(q==r),any(q.ne.r)
      write(*,*) x==y,x.ne.y
   
   end subroutine example_vectors

   function components_vector4(v) result(s)
      type(vector4),intent(in) :: v
      real*4 :: s(3)
      s = (/v%x,v%y,v%z/)
   end function components_vector4
   
   function components_vector8(v) result(s)
      type(vector8),intent(in) :: v
      real*8 :: s(3)
      s = (/v%x,v%y,v%z/)
   end function components_vector8
   
   function component_vector4(v,index) result(s)
      type(vector4),intent(in) :: v
      integer*4,intent(in):: index
      real*4 :: s
      select case(index)
      case(1)
         s = v%x
      case(2)
         s = v%y
      case(3)
         s = v%z
      case default
         s = 0
      end select
   end function component_vector4
   
   function component_vector8(v,index) result(s)
      type(vector8),intent(in) :: v
      integer*4,intent(in):: index
      real*8 :: s
      select case(index)
      case(1)
         s = v%x
      case(2)
         s = v%y
      case(3)
         s = v%z
      case default
         s = 0
      end select
   end function component_vector8
   
  subroutine getcomponents_vector4(s,v)
      class(vector4),intent(in) :: v
      real*4,intent(out) :: s(3)
      s = (/v%x,v%y,v%z/)
   end subroutine getcomponents_vector4
   
   subroutine getcomponents_vector8(s,v)
      class(vector8),intent(in) :: v
      real*8,intent(out) :: s(3)
      s = (/v%x,v%y,v%z/)
   end subroutine getcomponents_vector8
   
   subroutine setcomponents_vector4_from_real4(v,s)
      real*4,intent(in) :: s(3)
      type(vector4),intent(out) :: v
      v = vector4(s(1),s(2),s(3))
   end subroutine setcomponents_vector4_from_real4
   
   subroutine setcomponents_vector4_from_real8(v,s)
      real*8,intent(in) :: s(3)
      type(vector4),intent(out) :: v
      v = vector4(real(s(1),4),real(s(2),4),real(s(3),4))
   end subroutine setcomponents_vector4_from_real8
   
   subroutine setcomponents_vector4_from_int4(v,s)
      integer*4,intent(in) :: s(3)
      type(vector4),intent(out) :: v
      v = vector4(real(s(1),4),real(s(2),4),real(s(3),4))
   end subroutine setcomponents_vector4_from_int4
   
   subroutine setcomponents_vector4_from_int8(v,s)
      integer*8,intent(in) :: s(3)
      type(vector4),intent(out) :: v
      v = vector4(real(s(1),4),real(s(2),4),real(s(3),4))
   end subroutine setcomponents_vector4_from_int8
   
   subroutine setcomponents_vector8_from_real4(v,s)
      real*4,intent(in) :: s(3)
      type(vector8),intent(out) :: v
      v = vector8(real(s(1),8),real(s(2),8),real(s(3),8))
   end subroutine setcomponents_vector8_from_real4
   
   subroutine setcomponents_vector8_from_real8(v,s)
      real*8,intent(in) :: s(3)
      type(vector8),intent(out) :: v
      v = vector8(s(1),s(2),s(3))
   end subroutine setcomponents_vector8_from_real8
   
   subroutine setcomponents_vector8_from_int4(v,s)
      integer*4,intent(in) :: s(3)
      type(vector8),intent(out) :: v
      v = vector8(real(s(1),8),real(s(2),8),real(s(3),8))
   end subroutine setcomponents_vector8_from_int4
   
   subroutine setcomponents_vector8_from_int8(v,s)
      integer*8,intent(in) :: s(3)
      type(vector8),intent(out) :: v
      v = vector8(real(s(1),8),real(s(2),8),real(s(3),8))
   end subroutine setcomponents_vector8_from_int8
   
   pure function compare_vectors4(u,v) result(s)
      class(vector4),intent(in) :: u,v
      logical :: s
      s = u%x==v%x .and. u%y==v%y .and. u%z==v%z
   end function compare_vectors4
   
   pure function compare_vectors8(u,v) result(s)
      class(vector8),intent(in) :: u,v
      logical :: s
      s = u%x==v%x .and. u%y==v%y .and. u%z==v%z
   end function compare_vectors8
   
   pure function ncompare_vectors4(u,v) result(s)
      class(vector4),intent(in) :: u,v
      logical :: s
      s = .not.compare_vectors4(u,v)
   end function ncompare_vectors4
   
   pure function ncompare_vectors8(u,v) result(s)
      class(vector8),intent(in) :: u,v
      logical :: s
      s = .not.compare_vectors8(u,v)
   end function ncompare_vectors8
   
   pure function norm_vector4(this) result(n)
      class(vector4),intent(in) :: this
      real*4 :: n
      n = sqrt(this%x**2+this%y**2+this%z**2)
   end function norm_vector4
   
   pure function norm_vector8(this) result(n)
      class(vector8),intent(in) :: this
      real*8 :: n
      n = sqrt(this%x**2+this%y**2+this%z**2)
   end function norm_vector8
   
   pure function norm_array4(this) result(n)
      real*4,intent(in) :: this(:)
      real*4 :: n
      n = sqrt(sum(this**2))
   end function norm_array4
   
   pure function norm_array8(this) result(n)
      real*8,intent(in) :: this(:)
      real*8 :: n
      n = sqrt(sum(this**2))
   end function norm_array8
   
   pure function unit_vector4(this) result(w)
      class(vector4),intent(in) :: this
      type(vector4) :: w
      real*4 :: n
      n = this%norm()
      if (n<=epsilon(1.0)) then
         w = this*0
      else
         w = this/n
      end if
   end function unit_vector4
   
   pure function unit_vector8(this) result(w)
      class(vector8),intent(in) :: this
      type(vector8) :: w
      real*8 :: n
      n = this%norm()
      if (n<=epsilon(1.0_8)) then
         w = this*0
      else
         w = this/n
      end if
   end function unit_vector8
   
   pure function unit_array4(v) result(w)
      real*4,intent(in) :: v(:)
      real*4,allocatable :: w(:)
      real*4 :: n
      allocate(w(size(v)))
      n = norm(v)
      if (n<=epsilon(1.0)) then
         w = v*0
      else
         w = v/n
      end if
   end function unit_array4
   
   pure function unit_array8(v) result(w)
      real*8,intent(in) :: v(:)
      real*8,allocatable :: w(:)
      real*8 :: n
      allocate(w(size(v)))
      n = norm(v)
      if (n<=epsilon(1.0_8)) then
         w = v*0
      else
         w = v/n
      end if
   end function unit_array8
   
   pure function add_vectors4(u,v) result(w)
      class(vector4),intent(in) :: u,v
      type(vector4) :: w
      w = vector4(u%x+v%x,u%y+v%y,u%z+v%z)
   end function add_vectors4
   
   pure function add_vectors8(u,v) result(w)
      class(vector8),intent(in) :: u,v
      type(vector8) :: w
      w = vector8(u%x+v%x,u%y+v%y,u%z+v%z)
   end function add_vectors8
   
   pure function subtract_two_vectors4(u,v) result(w)
      class(vector4),intent(in) :: u,v
      type(vector4) :: w
      w = vector4(u%x-v%x,u%y-v%y,u%z-v%z)
   end function subtract_two_vectors4
   
   pure function subtract_two_vectors8(u,v) result(w)
      class(vector8),intent(in) :: u,v
      type(vector8) :: w
      w = vector8(u%x-v%x,u%y-v%y,u%z-v%z)
   end function subtract_two_vectors8
   
   pure function multiply_vector4_by_vector4(u,v) result(w)
      class(vector4),intent(in) :: u,v
      type(vector4) :: w
      w = vector4(u%x*v%x,u%y*v%y,u%z*v%z)
   end function multiply_vector4_by_vector4
   
   pure function multiply_vector8_by_vector8(u,v) result(w)
      class(vector8),intent(in) :: u,v
      type(vector8) :: w
      w = vector8(u%x*v%x,u%y*v%y,u%z*v%z)
   end function multiply_vector8_by_vector8
   
   pure function multiply_real4_with_vector4(k,v) result(w)
      real*4,intent(in) :: k
      class(vector4),intent(in) :: v
      type(vector4) :: w
      w = vector4(k*v%x,k*v%y,k*v%z)
   end function multiply_real4_with_vector4
   
   pure function multiply_real8_with_vector4(k,v) result(w)
      real*8,intent(in) :: k
      class(vector4),intent(in) :: v
      type(vector4) :: w
      w = vector4(real(k,4)*v%x,real(k,4)*v%y,real(k,4)*v%z)
   end function multiply_real8_with_vector4
   
   pure function multiply_int4_with_vector4(k,v) result(w)
      integer*4,intent(in) :: k
      class(vector4),intent(in) :: v
      type(vector4) :: w
      w = vector4(real(k,4)*v%x,real(k,4)*v%y,real(k,4)*v%z)
   end function multiply_int4_with_vector4
   
   pure function multiply_int8_with_vector4(k,v) result(w)
      integer*8,intent(in) :: k
      class(vector4),intent(in) :: v
      type(vector4) :: w
      w = vector4(real(k,4)*v%x,real(k,4)*v%y,real(k,4)*v%z)
   end function multiply_int8_with_vector4
   
   pure function multiply_vector4_with_real4(v,k) result(w)
      real*4,intent(in) :: k
      class(vector4),intent(in) :: v
      type(vector4) :: w
      w = vector4(k*v%x,k*v%y,k*v%z)
   end function multiply_vector4_with_real4
   
   pure function multiply_vector4_with_real8(v,k) result(w)
      real*8,intent(in) :: k
      class(vector4),intent(in) :: v
      type(vector4) :: w
      w = vector4(real(k,4)*v%x,real(k,4)*v%y,real(k,4)*v%z)
   end function multiply_vector4_with_real8
   
   pure function multiply_vector4_with_int4(v,k) result(w)
      integer*4,intent(in) :: k
      class(vector4),intent(in) :: v
      type(vector4) :: w
      w = vector4(real(k,4)*v%x,real(k,4)*v%y,real(k,4)*v%z)
   end function multiply_vector4_with_int4
   
   pure function multiply_vector4_with_int8(v,k) result(w)
      integer*8,intent(in) :: k
      class(vector4),intent(in) :: v
      type(vector4) :: w
      w = vector4(real(k,4)*v%x,real(k,4)*v%y,real(k,4)*v%z)
   end function multiply_vector4_with_int8
   
   pure function multiply_real4_with_vector8(k,v) result(w)
      real*4,intent(in) :: k
      class(vector8),intent(in) :: v
      type(vector8) :: w
      w = vector8(k*v%x,k*v%y,k*v%z)
   end function multiply_real4_with_vector8
   
   pure function multiply_real8_with_vector8(k,v) result(w)
      real*8,intent(in) :: k
      class(vector8),intent(in) :: v
      type(vector8) :: w
      w = vector8(real(k,4)*v%x,real(k,4)*v%y,real(k,4)*v%z)
   end function multiply_real8_with_vector8
   
   pure function multiply_int4_with_vector8(k,v) result(w)
      integer*4,intent(in) :: k
      class(vector8),intent(in) :: v
      type(vector8) :: w
      w = vector8(real(k,4)*v%x,real(k,4)*v%y,real(k,4)*v%z)
   end function multiply_int4_with_vector8
   
   pure function multiply_int8_with_vector8(k,v) result(w)
      integer*8,intent(in) :: k
      class(vector8),intent(in) :: v
      type(vector8) :: w
      w = vector8(real(k,4)*v%x,real(k,4)*v%y,real(k,4)*v%z)
   end function multiply_int8_with_vector8
   
   pure function multiply_vector8_with_real4(v,k) result(w)
      real*4,intent(in) :: k
      class(vector8),intent(in) :: v
      type(vector8) :: w
      w = vector8(k*v%x,k*v%y,k*v%z)
   end function multiply_vector8_with_real4
   
   pure function multiply_vector8_with_real8(v,k) result(w)
      real*8,intent(in) :: k
      class(vector8),intent(in) :: v
      type(vector8) :: w
      w = vector8(real(k,4)*v%x,real(k,4)*v%y,real(k,4)*v%z)
   end function multiply_vector8_with_real8
   
   pure function multiply_vector8_with_int4(v,k) result(w)
      integer*4,intent(in) :: k
      class(vector8),intent(in) :: v
      type(vector8) :: w
      w = vector8(real(k,4)*v%x,real(k,4)*v%y,real(k,4)*v%z)
   end function multiply_vector8_with_int4
   
   pure function multiply_vector8_with_int8(v,k) result(w)
      integer*8,intent(in) :: k
      class(vector8),intent(in) :: v
      type(vector8) :: w
      w = vector8(real(k,4)*v%x,real(k,4)*v%y,real(k,4)*v%z)
   end function multiply_vector8_with_int8
   
   pure function divide_vector4_by_vector4(u,v) result(w)
      class(vector4),intent(in) :: u,v
      type(vector4) :: w
      w = vector4(u%x/v%x,u%y/v%y,u%z/v%z)
   end function divide_vector4_by_vector4
   
   pure function divide_vector8_by_vector8(u,v) result(w)
      class(vector8),intent(in) :: u,v
      type(vector8) :: w
      w = vector8(u%x/v%x,u%y/v%y,u%z/v%z)
   end function divide_vector8_by_vector8
   
   pure function divide_vector4_by_real4(v,k) result(w)
      class(vector4),intent(in) :: v
      real*4,intent(in) :: k
      type(vector4) :: w
      w = vector4(v%x/k,v%y/k,v%z/k)
   end function divide_vector4_by_real4
   
   pure function divide_vector8_by_real4(v,k) result(w)
      class(vector8),intent(in) :: v
      real*4,intent(in) :: k
      type(vector8) :: w
      w = vector8(v%x/k,v%y/k,v%z/k)
   end function divide_vector8_by_real4
   
   pure function divide_vector4_by_real8(v,k) result(w)
      class(vector4),intent(in) :: v
      real*8,intent(in) :: k
      type(vector4) :: w
      w = vector4(v%x/real(k,4),v%y/real(k,4),v%z/real(k,4))
   end function divide_vector4_by_real8
   
   pure function divide_vector8_by_real8(v,k) result(w)
      class(vector8),intent(in) :: v
      real*8,intent(in) :: k
      type(vector8) :: w
      w = vector8(v%x/k,v%y/k,v%z/k)
   end function divide_vector8_by_real8
   
   pure function divide_vector4_by_int4(v,k) result(w)
      class(vector4),intent(in) :: v
      integer*4,intent(in) :: k
      type(vector4) :: w
      w = vector4(v%x/real(k,4),v%y/real(k,4),v%z/real(k,4))
   end function divide_vector4_by_int4
   
   pure function divide_vector8_by_int4(v,k) result(w)
      class(vector8),intent(in) :: v
      integer*4,intent(in) :: k
      type(vector8) :: w
      w = vector8(v%x/real(k,8),v%y/real(k,8),v%z/real(k,8))
   end function divide_vector8_by_int4
   
   pure function divide_vector4_by_int8(v,k) result(w)
      class(vector4),intent(in) :: v
      integer*8,intent(in) :: k
      type(vector4) :: w
      w = vector4(v%x/real(k,4),v%y/real(k,4),v%z/real(k,4))
   end function divide_vector4_by_int8
   
   pure function divide_vector8_by_int8(v,k) result(w)
      class(vector8),intent(in) :: v
      integer*8,intent(in) :: k
      type(vector8) :: w
      w = vector8(v%x/real(k,8),v%y/real(k,8),v%z/real(k,8))
   end function divide_vector8_by_int8
   
   pure function scalar_product_vector4(u,v) result(w)
      class(vector4),intent(in) :: u,v
      real*4 :: w
      w = u%x*v%x+u%y*v%y+u%z*v%z
   end function scalar_product_vector4
   
   pure function scalar_product_vector8(u,v) result(w)
      class(vector8),intent(in) :: u,v
      real*8 :: w
      w = u%x*v%x+u%y*v%y+u%z*v%z
   end function scalar_product_vector8
   
   pure function scalar_product_array4(u,v) result(w)
      real*4,intent(in) :: u(:),v(:)
      real*4 :: w
      w = sum(u*v)
   end function scalar_product_array4
   
   pure function scalar_product_array8(u,v) result(w)
      real*8,intent(in) :: u(:),v(:)
      real*8 :: w
      w = sum(u*v)
   end function scalar_product_array8
   
   function scalar_product_a4v4(u,v) result(w)
      real*4,intent(in)          :: u(3)
      class(vector4),intent(in)  :: v
      real*4 :: w
      w = u.dot.components(v)
   end function scalar_product_a4v4
   
   function scalar_product_v4a4(u,v) result(w)
      class(vector4),intent(in)  :: u
      real*4,intent(in)          :: v(3)
      real*4 :: w
      w = components(u).dot.v
   end function scalar_product_v4a4
   
   function scalar_product_a8v8(u,v) result(w)
      real*8,intent(in)          :: u(3)
      class(vector8),intent(in)  :: v
      real*8 :: w
      w = u.dot.components(v)
   end function scalar_product_a8v8
   
   function scalar_product_v8a8(u,v) result(w)
      class(vector8),intent(in)  :: u
      real*8,intent(in)          :: v(3)
      real*8 :: w
      w = components(u).dot.v
   end function scalar_product_v8a8
   
   function cross_product_a4v4(u,v) result(w)
      real*4,intent(in)          :: u(3)
      class(vector4),intent(in)  :: v
      real*4 :: w(3)
      w = u.cross.components(v)
   end function cross_product_a4v4
   
   function cross_product_v4a4(u,v) result(w)
      class(vector4),intent(in)  :: u
      real*4,intent(in)          :: v(3)
      real*4 :: w(3)
      w = components(u).cross.v
   end function cross_product_v4a4
   
   function cross_product_a8v8(u,v) result(w)
      real*8,intent(in)          :: u(3)
      class(vector8),intent(in)  :: v
      real*8 :: w(3)
      w = u.cross.components(v)
   end function cross_product_a8v8
   
   function cross_product_v8a8(u,v) result(w)
      class(vector8),intent(in)  :: u
      real*8,intent(in)          :: v(3)
      real*8 :: w(3)
      w = components(u).cross.v
   end function cross_product_v8a8
   
   pure function cross_product_vector4(u,v) result(w)
      class(vector4),intent(in) :: u,v
      type(vector4) :: w
      w = vector4(u%y*v%z-u%z*v%y,u%z*v%x-u%x*v%z,u%x*v%y-u%y*v%x)
   end function cross_product_vector4
   
   pure function cross_product_vector8(u,v) result(w)
      class(vector8),intent(in) :: u,v
      type(vector8) :: w
      w = vector8(u%y*v%z-u%z*v%y,u%z*v%x-u%x*v%z,u%x*v%y-u%y*v%x)
   end function cross_product_vector8
   
   pure function cross_product_array4(u,v) result(w)
      real*4,intent(in) :: u(3),v(3)
      real*4 :: w(3)
      w(1) = u(2)*v(3)-u(3)*v(2)
      w(2) = u(3)*v(1)-u(1)*v(3)
      w(3) = u(1)*v(2)-u(2)*v(1)
   end function cross_product_array4
   
   pure function cross_product_array8(u,v) result(w)
      real*8,intent(in) :: u(3),v(3)
      real*8 :: w(3)
      w(1) = u(2)*v(3)-u(3)*v(2)
      w(2) = u(3)*v(1)-u(1)*v(3)
      w(3) = u(1)*v(2)-u(2)*v(1)
   end function cross_product_array8
   
   pure function matmul_matrix4_times_vector4(m,v) result(w)
      real*4,intent(in) :: m(3,3)
      class(vector4),intent(in) :: v
      type(vector4) :: w
      real*4 :: q(3)
      q = matmul(m,(/v%x,v%y,v%z/))
      w = vector4(q(1),q(2),q(3))
   end function matmul_matrix4_times_vector4
   
   pure function matmul_matrix4_times_vector8(m,v) result(w)
      real*4,intent(in) :: m(3,3)
      class(vector8),intent(in) :: v
      type(vector8) :: w
      real*8 :: q(3)
      q = matmul(real(m,8),(/v%x,v%y,v%z/))
      w = vector8(q(1),q(2),q(3))
   end function matmul_matrix4_times_vector8
   
   pure function matmul_matrix8_times_vector4(m,v) result(w)
      real*8,intent(in) :: m(3,3)
      class(vector4),intent(in) :: v
      type(vector4) :: w
      real*4 :: q(3)
      q = matmul(real(m,4),(/v%x,v%y,v%z/))
      w = vector4(q(1),q(2),q(3))
   end function matmul_matrix8_times_vector4
   
   pure function matmul_matrix8_times_vector8(m,v) result(w)
      real*8,intent(in) :: m(3,3)
      class(vector8),intent(in) :: v
      type(vector8) :: w
      real*8 :: q(3)
      q = matmul(m,(/v%x,v%y,v%z/))
      w = vector8(q(1),q(2),q(3))
   end function matmul_matrix8_times_vector8
   
   pure function matmul_vector4_times_matrix4(u,m) result(w)
      real*4,intent(in) :: m(3,3)
      class(vector4),intent(in) :: u
      type(vector4) :: w
      real*4 :: q(3)
      q = matmul((/u%x,u%y,u%z/),m)
      w = vector4(q(1),q(2),q(3))
   end function matmul_vector4_times_matrix4
   
   pure function matmul_vector8_times_matrix4(u,m) result(w)
      real*4,intent(in) :: m(3,3)
      class(vector8),intent(in) :: u
      type(vector8) :: w
      real*8 :: q(3)
      q = matmul((/u%x,u%y,u%z/),real(m,8))
      w = vector8(q(1),q(2),q(3))
   end function matmul_vector8_times_matrix4
   
   pure function matmul_vector4_times_matrix8(u,m) result(w)
      real*8,intent(in) :: m(3,3)
      class(vector4),intent(in) :: u
      type(vector4) :: w
      real*4 :: q(3)
      q = matmul((/u%x,u%y,u%z/),real(m,4))
      w = vector4(q(1),q(2),q(3))
   end function matmul_vector4_times_matrix8
   
   pure function matmul_vector8_times_matrix8(u,m) result(w)
      real*8,intent(in) :: m(3,3)
      class(vector8),intent(in) :: u
      type(vector8) :: w
      real*8 :: q(3)
      q = matmul((/u%x,u%y,u%z/),m)
      w = vector8(q(1),q(2),q(3))
   end function matmul_vector8_times_matrix8
      
end module shared_module_vectors