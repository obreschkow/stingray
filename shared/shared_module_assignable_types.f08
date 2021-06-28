module shared_module_assignable_types

   private
   
   public   :: myreal4
   public   :: myinteger4
   public   :: mylogical
   public   :: mystring
   public   :: operator(+)
   public   :: operator(-)
   public   :: operator(*)
   public   :: operator(/)
   public   :: operator(==)
   public   :: operator(/=)
   public   :: operator(>)
   public   :: operator(<)
   public   :: operator(>=)
   public   :: operator(<=)
   public   :: assignment(=)
   
   type myreal4
      real*4 :: value = 0.0
      logical :: assigned = .false.
   end type myreal4
   
   type myinteger4
      integer*4 :: value = 0
      logical :: assigned = .false.
   end type myinteger4
   
   type mylogical
      logical :: value = .false.
      logical :: assigned = .false.
   end type mylogical
   
   type mystring
      character(len=255) :: value
      logical :: assigned = .false.
   end type mystring
   
   interface assignment(=)
      procedure myreal4_from_real4
      procedure myreal4_from_integer4
      procedure myreal4_from_integer8
      procedure myinteger4_from_integer4
      procedure mylogical_from_logical
      procedure mystring_from_character
      procedure real4_from_myreal4
      procedure integer4_from_myinteger4
      procedure integer8_from_myinteger4
      procedure real4_from_myinteger4
      procedure real8_from_myinteger4
      procedure logical_from_mylogical
      procedure character_from_mystring
   end interface assignment(=)
   
   interface operator(+)
      procedure myreal4_plus_myreal4
      procedure myreal4_plus_real4
      procedure myreal4_plus_real8
      procedure myreal4_plus_integer4
      procedure myreal4_plus_integer8
      procedure real4_plus_myreal4
      procedure real8_plus_myreal4
      procedure integer4_plus_myreal4
      procedure integer8_plus_myreal4
      procedure myinteger4_plus_myinteger4
      procedure myinteger4_plus_integer4
      procedure integer4_plus_myinteger4
   end interface operator(+)
   
   interface operator(-)
      procedure myreal4_minus_myreal4
      procedure myreal4_minus_real4
      procedure myreal4_minus_real8
      procedure myreal4_minus_integer4
      procedure myreal4_minus_integer8
      procedure real4_minus_myreal4
      procedure real8_minus_myreal4
      procedure integer4_minus_myreal4
      procedure integer8_minus_myreal4
      procedure myinteger4_minus_myinteger4
      procedure myinteger4_minus_integer4
      procedure integer4_minus_myinteger4
   end interface operator(-)
   
   interface operator(*)
      procedure myreal4_times_myreal4
      procedure myreal4_times_real4
      procedure myreal4_times_real8
      procedure myreal4_times_integer4
      procedure myreal4_times_integer8
      procedure real4_times_myreal4
      procedure real8_times_myreal4
      procedure integer4_times_myreal4
      procedure integer8_times_myreal4
      procedure myinteger4_times_myinteger4
      procedure myinteger4_times_integer4
      procedure integer4_times_myinteger4
   end interface operator(*)
   
   interface operator(/)
      procedure myreal4_div_myreal4
      procedure myreal4_div_real4
      procedure myreal4_div_real8
      procedure myreal4_div_integer4
      procedure myreal4_div_integer8
      procedure real4_div_myreal4
      procedure real8_div_myreal4
      procedure integer4_div_myreal4
      procedure integer8_div_myreal4
      procedure myinteger4_div_myinteger4
      procedure myinteger4_div_integer4
      procedure integer4_div_myinteger4
   end interface operator(/)
   
   interface operator(==)
      procedure myreal4_equals_myreal4
      procedure myreal4_equals_real4
      procedure real4_equals_myreal4
      procedure myinteger4_equals_myinteger4
      procedure myinteger4_equals_integer4
      procedure integer4_equals_myinteger4
   end interface operator(==)
   
   interface operator(/=)
      procedure myreal4_notequal_myreal4
      procedure myreal4_notequal_real4
      procedure real4_notequal_myreal4
      procedure myinteger4_notequal_myinteger4
      procedure myinteger4_notequal_integer4
      procedure integer4_notequal_myinteger4
   end interface operator(/=)
   
   interface operator(<)
      procedure myreal4_smallerthan_myreal4
      procedure myreal4_smallerthan_real4
      procedure real4_smallerthan_myreal4
      procedure myinteger4_smallerthan_myinteger4
      procedure myinteger4_smallerthan_integer4
      procedure integer4_smallerthan_myinteger4
   end interface operator(<)
   
   interface operator(<=)
      procedure myreal4_smallerequalthan_myreal4
      procedure myreal4_smallerequalthan_real4
      procedure real4_smallerequalthan_myreal4
      procedure myinteger4_smallerequalthan_myinteger4
      procedure myinteger4_smallerequalthan_integer4
      procedure integer4_smallerequalthan_myinteger4
   end interface operator(<=)
   
   interface operator(>)
      procedure myreal4_greaterthan_myreal4
      procedure myreal4_greaterthan_real4
      procedure real4_greaterthan_myreal4
      procedure myinteger4_greaterthan_myinteger4
      procedure myinteger4_greaterthan_integer4
      procedure integer4_greaterthan_myinteger4
   end interface operator(>)
   
   interface operator(>=)
      procedure myreal4_greaterequalthan_myreal4
      procedure myreal4_greaterequalthan_real4
      procedure real4_greaterequalthan_myreal4
      procedure myinteger4_greaterequalthan_myinteger4
      procedure myinteger4_greaterequalthan_integer4
      procedure integer4_greaterequalthan_myinteger4
   end interface operator(>=)
   
contains

elemental function myreal4_equals_myreal4(x,y) result(s)
   type(myreal4),intent(in) :: x
   type(myreal4),intent(in) :: y
   logical :: s
   s = x%value==y%value
end function myreal4_equals_myreal4

elemental function myreal4_equals_real4(x,y) result(s)
   type(myreal4),intent(in) :: x
   real*4,intent(in) :: y
   logical :: s
   s = x%value==y
end function myreal4_equals_real4

elemental function real4_equals_myreal4(x,y) result(s)
   real*4,intent(in) :: x
   type(myreal4),intent(in) :: y
   logical :: s
   s = x==y%value
end function real4_equals_myreal4

elemental function myinteger4_equals_myinteger4(x,y) result(s)
   type(myinteger4),intent(in) :: x
   type(myinteger4),intent(in) :: y
   logical :: s
   s = x%value==y%value
end function myinteger4_equals_myinteger4

elemental function myinteger4_equals_integer4(x,y) result(s)
   type(myinteger4),intent(in) :: x
   integer*4,intent(in) :: y
   logical :: s
   s = x%value==y
end function myinteger4_equals_integer4

elemental function integer4_equals_myinteger4(x,y) result(s)
   integer*4,intent(in) :: x
   type(myinteger4),intent(in) :: y
   logical :: s
   s = x==y%value
end function integer4_equals_myinteger4

elemental function myreal4_notequal_myreal4(x,y) result(s)
   type(myreal4),intent(in) :: x
   type(myreal4),intent(in) :: y
   logical :: s
   s = x%value/=y%value
end function myreal4_notequal_myreal4

elemental function myreal4_notequal_real4(x,y) result(s)
   type(myreal4),intent(in) :: x
   real*4,intent(in) :: y
   logical :: s
   s = x%value/=y
end function myreal4_notequal_real4

elemental function real4_notequal_myreal4(x,y) result(s)
   real*4,intent(in) :: x
   type(myreal4),intent(in) :: y
   logical :: s
   s = x/=y%value
end function real4_notequal_myreal4

elemental function myinteger4_notequal_myinteger4(x,y) result(s)
   type(myinteger4),intent(in) :: x
   type(myinteger4),intent(in) :: y
   logical :: s
   s = x%value/=y%value
end function myinteger4_notequal_myinteger4

elemental function myinteger4_notequal_integer4(x,y) result(s)
   type(myinteger4),intent(in) :: x
   integer*4,intent(in) :: y
   logical :: s
   s = x%value/=y
end function myinteger4_notequal_integer4

elemental function integer4_notequal_myinteger4(x,y) result(s)
   integer*4,intent(in) :: x
   type(myinteger4),intent(in) :: y
   logical :: s
   s = x/=y%value
end function integer4_notequal_myinteger4

elemental function myreal4_smallerthan_myreal4(x,y) result(s)
   type(myreal4),intent(in) :: x
   type(myreal4),intent(in) :: y
   logical :: s
   s = x%value<y%value
end function myreal4_smallerthan_myreal4

elemental function myreal4_smallerthan_real4(x,y) result(s)
   type(myreal4),intent(in) :: x
   real*4,intent(in) :: y
   logical :: s
   s = x%value<y
end function myreal4_smallerthan_real4

elemental function real4_smallerthan_myreal4(x,y) result(s)
   real*4,intent(in) :: x
   type(myreal4),intent(in) :: y
   logical :: s
   s = x<y%value
end function real4_smallerthan_myreal4

elemental function myinteger4_smallerthan_myinteger4(x,y) result(s)
   type(myinteger4),intent(in) :: x
   type(myinteger4),intent(in) :: y
   logical :: s
   s = x%value<y%value
end function myinteger4_smallerthan_myinteger4

elemental function myinteger4_smallerthan_integer4(x,y) result(s)
   type(myinteger4),intent(in) :: x
   integer*4,intent(in) :: y
   logical :: s
   s = x%value<y
end function myinteger4_smallerthan_integer4

elemental function integer4_smallerthan_myinteger4(x,y) result(s)
   integer*4,intent(in) :: x
   type(myinteger4),intent(in) :: y
   logical :: s
   s = x<y%value
end function integer4_smallerthan_myinteger4

elemental function myreal4_smallerequalthan_myreal4(x,y) result(s)
   type(myreal4),intent(in) :: x
   type(myreal4),intent(in) :: y
   logical :: s
   s = x%value<=y%value
end function myreal4_smallerequalthan_myreal4

elemental function myreal4_smallerequalthan_real4(x,y) result(s)
   type(myreal4),intent(in) :: x
   real*4,intent(in) :: y
   logical :: s
   s = x%value<=y
end function myreal4_smallerequalthan_real4

elemental function real4_smallerequalthan_myreal4(x,y) result(s)
   real*4,intent(in) :: x
   type(myreal4),intent(in) :: y
   logical :: s
   s = x<=y%value
end function real4_smallerequalthan_myreal4

elemental function myinteger4_smallerequalthan_myinteger4(x,y) result(s)
   type(myinteger4),intent(in) :: x
   type(myinteger4),intent(in) :: y
   logical :: s
   s = x%value<=y%value
end function myinteger4_smallerequalthan_myinteger4

elemental function myinteger4_smallerequalthan_integer4(x,y) result(s)
   type(myinteger4),intent(in) :: x
   integer*4,intent(in) :: y
   logical :: s
   s = x%value<=y
end function myinteger4_smallerequalthan_integer4

elemental function integer4_smallerequalthan_myinteger4(x,y) result(s)
   integer*4,intent(in) :: x
   type(myinteger4),intent(in) :: y
   logical :: s
   s = x<=y%value
end function integer4_smallerequalthan_myinteger4

elemental function myreal4_greaterthan_myreal4(x,y) result(s)
   type(myreal4),intent(in) :: x
   type(myreal4),intent(in) :: y
   logical :: s
   s = x%value>y%value
end function myreal4_greaterthan_myreal4

elemental function myreal4_greaterthan_real4(x,y) result(s)
   type(myreal4),intent(in) :: x
   real*4,intent(in) :: y
   logical :: s
   s = x%value>y
end function myreal4_greaterthan_real4

elemental function real4_greaterthan_myreal4(x,y) result(s)
   real*4,intent(in) :: x
   type(myreal4),intent(in) :: y
   logical :: s
   s = x>y%value
end function real4_greaterthan_myreal4

elemental function myinteger4_greaterthan_myinteger4(x,y) result(s)
   type(myinteger4),intent(in) :: x
   type(myinteger4),intent(in) :: y
   logical :: s
   s = x%value>y%value
end function myinteger4_greaterthan_myinteger4

elemental function myinteger4_greaterthan_integer4(x,y) result(s)
   type(myinteger4),intent(in) :: x
   integer*4,intent(in) :: y
   logical :: s
   s = x%value>y
end function myinteger4_greaterthan_integer4

elemental function integer4_greaterthan_myinteger4(x,y) result(s)
   integer*4,intent(in) :: x
   type(myinteger4),intent(in) :: y
   logical :: s
   s = x>y%value
end function integer4_greaterthan_myinteger4

elemental function myreal4_greaterequalthan_myreal4(x,y) result(s)
   type(myreal4),intent(in) :: x
   type(myreal4),intent(in) :: y
   logical :: s
   s = x%value>=y%value
end function myreal4_greaterequalthan_myreal4

elemental function myreal4_greaterequalthan_real4(x,y) result(s)
   type(myreal4),intent(in) :: x
   real*4,intent(in) :: y
   logical :: s
   s = x%value>=y
end function myreal4_greaterequalthan_real4

elemental function real4_greaterequalthan_myreal4(x,y) result(s)
   real*4,intent(in) :: x
   type(myreal4),intent(in) :: y
   logical :: s
   s = x>=y%value
end function real4_greaterequalthan_myreal4

elemental function myinteger4_greaterequalthan_myinteger4(x,y) result(s)
   type(myinteger4),intent(in) :: x
   type(myinteger4),intent(in) :: y
   logical :: s
   s = x%value>=y%value
end function myinteger4_greaterequalthan_myinteger4

elemental function myinteger4_greaterequalthan_integer4(x,y) result(s)
   type(myinteger4),intent(in) :: x
   integer*4,intent(in) :: y
   logical :: s
   s = x%value>=y
end function myinteger4_greaterequalthan_integer4

elemental function integer4_greaterequalthan_myinteger4(x,y) result(s)
   integer*4,intent(in) :: x
   type(myinteger4),intent(in) :: y
   logical :: s
   s = x>=y%value
end function integer4_greaterequalthan_myinteger4

elemental function myreal4_plus_myreal4(x,y) result(z)
   type(myreal4),intent(in) :: x,y
   real*4 :: z
   z = x%value+y%value
end function myreal4_plus_myreal4

elemental function myreal4_minus_myreal4(x,y) result(z)
   type(myreal4),intent(in) :: x,y
   real*4 :: z
   z = x%value-y%value
end function myreal4_minus_myreal4

elemental function myreal4_times_myreal4(x,y) result(z)
   type(myreal4),intent(in) :: x,y
   real*4 :: z
   z = x%value*y%value
end function myreal4_times_myreal4

elemental function myreal4_div_myreal4(x,y) result(z)
   type(myreal4),intent(in) :: x,y
   real*4 :: z
   z = x%value/y%value
end function myreal4_div_myreal4

elemental function myreal4_plus_real4(x,y) result(z)
   type(myreal4),intent(in) :: x
   real*4,intent(in) :: y
   real*4 :: z
   z = x%value+y
end function myreal4_plus_real4

elemental function myreal4_plus_real8(x,y) result(z)
   type(myreal4),intent(in) :: x
   real*8,intent(in) :: y
   real*8 :: z
   z = x%value+y
end function myreal4_plus_real8

elemental function myreal4_plus_integer4(x,y) result(z)
   type(myreal4),intent(in) :: x
   integer*4,intent(in) :: y
   real*4 :: z
   z = x%value+real(y,4)
end function myreal4_plus_integer4

elemental function myreal4_plus_integer8(x,y) result(z)
   type(myreal4),intent(in) :: x
   integer*8,intent(in) :: y
   real*4 :: z
   z = x%value+real(y,4)
end function myreal4_plus_integer8

elemental function real4_plus_myreal4(x,y) result(z)
   type(myreal4),intent(in) :: y
   real*4,intent(in) :: x
   real*4 :: z
   z = x+y%value
end function real4_plus_myreal4

elemental function real8_plus_myreal4(x,y) result(z)
   type(myreal4),intent(in) :: y
   real*8,intent(in) :: x
   real*8 :: z
   z = x+y%value
end function real8_plus_myreal4

elemental function integer4_plus_myreal4(x,y) result(z)
   type(myreal4),intent(in) :: y
   integer*4,intent(in) :: x
   real*4 :: z
   z = real(x,4)+y%value
end function integer4_plus_myreal4

elemental function integer8_plus_myreal4(x,y) result(z)
   type(myreal4),intent(in) :: y
   integer*8,intent(in) :: x
   real*4 :: z
   z = real(x,4)+y%value
end function integer8_plus_myreal4

elemental function myreal4_minus_real4(x,y) result(z)
   type(myreal4),intent(in) :: x
   real*4,intent(in) :: y
   real*4 :: z
   z = x%value-y
end function myreal4_minus_real4

elemental function myreal4_minus_real8(x,y) result(z)
   type(myreal4),intent(in) :: x
   real*8,intent(in) :: y
   real*8 :: z
   z = x%value-y
end function myreal4_minus_real8

elemental function myreal4_minus_integer4(x,y) result(z)
   type(myreal4),intent(in) :: x
   integer*4,intent(in) :: y
   real*4 :: z
   z = x%value-real(y,4)
end function myreal4_minus_integer4

elemental function myreal4_minus_integer8(x,y) result(z)
   type(myreal4),intent(in) :: x
   integer*8,intent(in) :: y
   real*4 :: z
   z = x%value-real(y,4)
end function myreal4_minus_integer8

elemental function real4_minus_myreal4(x,y) result(z)
   type(myreal4),intent(in) :: y
   real*4,intent(in) :: x
   real*4 :: z
   z = x-y%value
end function real4_minus_myreal4

elemental function real8_minus_myreal4(x,y) result(z)
   type(myreal4),intent(in) :: y
   real*8,intent(in) :: x
   real*8 :: z
   z = x-y%value
end function real8_minus_myreal4

elemental function integer4_minus_myreal4(x,y) result(z)
   type(myreal4),intent(in) :: y
   integer*4,intent(in) :: x
   real*4 :: z
   z = real(x,4)-y%value
end function integer4_minus_myreal4

elemental function integer8_minus_myreal4(x,y) result(z)
   type(myreal4),intent(in) :: y
   integer*8,intent(in) :: x
   real*4 :: z
   z = real(x,4)-y%value
end function integer8_minus_myreal4

elemental function myreal4_times_real4(x,y) result(z)
   type(myreal4),intent(in) :: x
   real*4,intent(in) :: y
   real*4 :: z
   z = x%value*y
end function myreal4_times_real4

elemental function myreal4_times_real8(x,y) result(z)
   type(myreal4),intent(in) :: x
   real*8,intent(in) :: y
   real*8 :: z
   z = x%value*y
end function myreal4_times_real8

elemental function myreal4_times_integer4(x,y) result(z)
   type(myreal4),intent(in) :: x
   integer*4,intent(in) :: y
   real*4 :: z
   z = x%value*real(y,4)
end function myreal4_times_integer4

elemental function myreal4_times_integer8(x,y) result(z)
   type(myreal4),intent(in) :: x
   integer*8,intent(in) :: y
   real*4 :: z
   z = x%value*real(y,4)
end function myreal4_times_integer8

elemental function real4_times_myreal4(x,y) result(z)
   type(myreal4),intent(in) :: y
   real*4,intent(in) :: x
   real*4 :: z
   z = x*y%value
end function real4_times_myreal4

elemental function real8_times_myreal4(x,y) result(z)
   type(myreal4),intent(in) :: y
   real*8,intent(in) :: x
   real*8 :: z
   z = x*y%value
end function real8_times_myreal4

elemental function integer4_times_myreal4(x,y) result(z)
   type(myreal4),intent(in) :: y
   integer*4,intent(in) :: x
   real*4 :: z
   z = real(x,4)*y%value
end function integer4_times_myreal4

elemental function integer8_times_myreal4(x,y) result(z)
   type(myreal4),intent(in) :: y
   integer*8,intent(in) :: x
   real*4 :: z
   z = real(x,4)*y%value
end function integer8_times_myreal4

elemental function myreal4_div_real4(x,y) result(z)
   type(myreal4),intent(in) :: x
   real*4,intent(in) :: y
   real*4 :: z
   z = x%value/y
end function myreal4_div_real4

elemental function myreal4_div_real8(x,y) result(z)
   type(myreal4),intent(in) :: x
   real*8,intent(in) :: y
   real*8 :: z
   z = x%value/y
end function myreal4_div_real8

elemental function myreal4_div_integer4(x,y) result(z)
   type(myreal4),intent(in) :: x
   integer*4,intent(in) :: y
   real*4 :: z
   z = x%value/real(y,4)
end function myreal4_div_integer4

elemental function myreal4_div_integer8(x,y) result(z)
   type(myreal4),intent(in) :: x
   integer*8,intent(in) :: y
   real*4 :: z
   z = x%value/real(y,4)
end function myreal4_div_integer8

elemental function real4_div_myreal4(x,y) result(z)
   type(myreal4),intent(in) :: y
   real*4,intent(in) :: x
   real*4 :: z
   z = x/y%value
end function real4_div_myreal4

elemental function real8_div_myreal4(x,y) result(z)
   type(myreal4),intent(in) :: y
   real*8,intent(in) :: x
   real*8 :: z
   z = x/y%value
end function real8_div_myreal4

elemental function integer4_div_myreal4(x,y) result(z)
   type(myreal4),intent(in) :: y
   integer*4,intent(in) :: x
   real*4 :: z
   z = real(x,4)/y%value
end function integer4_div_myreal4

elemental function integer8_div_myreal4(x,y) result(z)
   type(myreal4),intent(in) :: y
   integer*8,intent(in) :: x
   real*4 :: z
   z = real(x,4)/y%value
end function integer8_div_myreal4

elemental function myinteger4_plus_myinteger4(x,y) result(z)
   type(myinteger4),intent(in) :: x,y
   integer*4 :: z
   z = x%value+y%value
end function myinteger4_plus_myinteger4

elemental function myinteger4_plus_integer4(x,y) result(z)
   type(myinteger4),intent(in) :: x
   integer*4,intent(in) :: y
   integer*4 :: z
   z = x%value+y
end function myinteger4_plus_integer4

elemental function integer4_plus_myinteger4(x,y) result(z)
   integer*4,intent(in) :: x
   type(myinteger4),intent(in) :: y
   integer*4 :: z
   z = x+y%value
end function integer4_plus_myinteger4

elemental function myinteger4_minus_myinteger4(x,y) result(z)
   type(myinteger4),intent(in) :: x,y
   integer*4 :: z
   z = x%value-y%value
end function myinteger4_minus_myinteger4

elemental function myinteger4_minus_integer4(x,y) result(z)
   type(myinteger4),intent(in) :: x
   integer*4,intent(in) :: y
   integer*4 :: z
   z = x%value+y
end function myinteger4_minus_integer4

elemental function integer4_minus_myinteger4(x,y) result(z)
   integer*4,intent(in) :: x
   type(myinteger4),intent(in) :: y
   integer*4 :: z
   z = x+y%value
end function integer4_minus_myinteger4

elemental function myinteger4_times_myinteger4(x,y) result(z)
   type(myinteger4),intent(in) :: x,y
   integer*4 :: z
   z = x%value*y%value
end function myinteger4_times_myinteger4

elemental function myinteger4_times_integer4(x,y) result(z)
   type(myinteger4),intent(in) :: x
   integer*4,intent(in) :: y
   integer*4 :: z
   z = x%value+y
end function myinteger4_times_integer4

elemental function integer4_times_myinteger4(x,y) result(z)
   integer*4,intent(in) :: x
   type(myinteger4),intent(in) :: y
   integer*4 :: z
   z = x+y%value
end function integer4_times_myinteger4

elemental function myinteger4_div_myinteger4(x,y) result(z)
   type(myinteger4),intent(in) :: x,y
   integer*4 :: z
   z = x%value/y%value
end function myinteger4_div_myinteger4

elemental function myinteger4_div_integer4(x,y) result(z)
   type(myinteger4),intent(in) :: x
   integer*4,intent(in) :: y
   integer*4 :: z
   z = x%value+y
end function myinteger4_div_integer4

elemental function integer4_div_myinteger4(x,y) result(z)
   integer*4,intent(in) :: x
   type(myinteger4),intent(in) :: y
   integer*4 :: z
   z = x+y%value
end function integer4_div_myinteger4

subroutine myreal4_from_real4(x,y)
   type(myreal4),intent(out) :: x
   real*4,intent(in) :: y
   x%value = y
   x%assigned = .true.
end subroutine myreal4_from_real4

subroutine myreal4_from_integer4(x,y)
   type(myreal4),intent(out) :: x
   integer*4,intent(in) :: y
   x%value = real(y,4)
   x%assigned = .true.
end subroutine myreal4_from_integer4

subroutine myreal4_from_integer8(x,y)
   type(myreal4),intent(out) :: x
   integer*8,intent(in) :: y
   x%value = real(y,4)
   x%assigned = .true.
end subroutine myreal4_from_integer8

subroutine myinteger4_from_integer4(x,y)
   type(myinteger4),intent(out) :: x
   integer*4,intent(in) :: y
   x%value = y
   x%assigned = .true.
end subroutine myinteger4_from_integer4

subroutine mylogical_from_logical(x,y)
   type(mylogical),intent(out) :: x
   logical,intent(in) :: y
   x%value = y
   x%assigned = .true.
end subroutine mylogical_from_logical

subroutine mystring_from_character(x,y)
   type(mystring(*)),intent(out) :: x
   character(*),intent(in) :: y
   x%value = y
   x%assigned = .true.
end subroutine mystring_from_character

subroutine real4_from_myreal4(y,x)
   real*4,intent(out) :: y
   type(myreal4),intent(in) :: x
   y = x%value
end subroutine real4_from_myreal4

subroutine integer4_from_myinteger4(y,x)
   integer*4,intent(out) :: y
   type(myinteger4),intent(in) :: x
   y = x%value
end subroutine integer4_from_myinteger4

subroutine integer8_from_myinteger4(y,x)
   integer*8,intent(out) :: y
   type(myinteger4),intent(in) :: x
   y = x%value
end subroutine integer8_from_myinteger4

subroutine real4_from_myinteger4(y,x)
   real*4,intent(out) :: y
   type(myinteger4),intent(in) :: x
   y = x%value
end subroutine real4_from_myinteger4

subroutine real8_from_myinteger4(y,x)
   real*8,intent(out) :: y
   type(myinteger4),intent(in) :: x
   y = x%value
end subroutine real8_from_myinteger4

subroutine logical_from_mylogical(y,x)
   logical,intent(out) :: y
   type(mylogical),intent(in) :: x
   y = x%value
end subroutine logical_from_mylogical

subroutine character_from_mystring(y,x)
   character(*),intent(out) :: y
   type(mystring(*)),intent(in) :: x
   y = x%value
end subroutine character_from_mystring

end module shared_module_assignable_types