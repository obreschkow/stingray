! **********************************************************************************************************************************
! Shared Fortran module to sort a vector in increasing order
! Developed by Danail Obreschkow
! **********************************************************************************************************************************

module shared_module_sort

   public   :: sort

   private
   
   interface sort
      module procedure merge_sort_list_int4
      module procedure merge_sort_list_int8
      module procedure merge_sort_list_real4
      module procedure merge_sort_list_real8
   end interface sort

contains
   
   recursive subroutine merge_sort_list_int4(vector,index,recursive_call)
      
      implicit none
      integer*4,intent(inout)             :: vector(:)
      integer*4,intent(inout)             :: index(:)
      logical*4,intent(in),optional       :: recursive_call ! only used internally
      integer*4,allocatable               :: vector1(:),vector2(:)
      integer*4,allocatable               :: index1(:),index2(:)
      integer*4                           :: n0,n1,n2,i0,i1,i2
      
      n0 = size(vector)
      
      ! make index vector
      if (.not.present(recursive_call)) then
         do i0 = 1,n0
            index(i0) = i0
         end do
      end if
      
      ! split list in 2 sublists
      n1 = n0/2
      n2 = n0-n1
      allocate(vector1(n1),vector2(n2),index1(n1),index2(n2))
      vector1(1:n1) = vector(1:n1)
      vector2(1:n2) = vector(n0/2+1:n0)
      index1(1:n1) = index(1:n1)
      index2(1:n2) = index(n0/2+1:n0)
      
      if (n1>1) call merge_sort_list_int4(vector1,index1,recursive_call=.true.)
      if (n2>1) call merge_sort_list_int4(vector2,index2,recursive_call=.true.)
      
      ! merge sorted sublists
      i0 = 1
      i1 = 1
      i2 = 1
      do while ((i1<=n1) .and. (i2<=n2))
         if (vector1(i1)<vector2(i2)) then
            vector(i0) = vector1(i1)
            index(i0) = index1(i1)
            i1 = i1+1
         else
            vector(i0) = vector2(i2)
            index(i0) = index2(i2)
            i2 = i2+1
         end if
         i0 = i0+1
      end do
      if (i1<=n1) then
         vector(i0:i0+n1-i1) = vector1(i1:n1)
         index(i0:i0+n1-i1) = index1(i1:n1)
      else if (i2<=n2) then
         vector(i0:i0+n2-i2) = vector2(i2:n2)
         index(i0:i0+n2-i2) = index2(i2:n2)
      end if
      
   end subroutine merge_sort_list_int4

   recursive subroutine merge_sort_list_int8(vector,index,recursive_call)
      
      implicit none
      integer*8,intent(inout)             :: vector(:)
      integer*4,intent(inout)             :: index(:)
      logical*4,intent(in),optional       :: recursive_call ! only used internally
      integer*8,allocatable               :: vector1(:),vector2(:)
      integer*4,allocatable               :: index1(:),index2(:)
      integer*4                           :: n0,n1,n2,i0,i1,i2
      
      n0 = size(vector)
      
      ! make index vector
      if (.not.present(recursive_call)) then
         do i0 = 1,n0
            index(i0) = i0
         end do
      end if
      
      ! split list in 2 sublists
      n1 = n0/2
      n2 = n0-n1
      allocate(vector1(n1),vector2(n2),index1(n1),index2(n2))
      vector1(1:n1) = vector(1:n1)
      vector2(1:n2) = vector(n0/2+1:n0)
      index1(1:n1) = index(1:n1)
      index2(1:n2) = index(n0/2+1:n0)
      
      if (n1>1) call merge_sort_list_int8(vector1,index1,recursive_call=.true.)
      if (n2>1) call merge_sort_list_int8(vector2,index2,recursive_call=.true.)
      
      ! merge sorted sublists
      i0 = 1
      i1 = 1
      i2 = 1
      do while ((i1<=n1) .and. (i2<=n2))
         if (vector1(i1)<vector2(i2)) then
            vector(i0) = vector1(i1)
            index(i0) = index1(i1)
            i1 = i1+1
         else
            vector(i0) = vector2(i2)
            index(i0) = index2(i2)
            i2 = i2+1
         end if
         i0 = i0+1
      end do
      if (i1<=n1) then
         vector(i0:i0+n1-i1) = vector1(i1:n1)
         index(i0:i0+n1-i1) = index1(i1:n1)
      else if (i2<=n2) then
         vector(i0:i0+n2-i2) = vector2(i2:n2)
         index(i0:i0+n2-i2) = index2(i2:n2)
      end if
      
   end subroutine merge_sort_list_int8
   
   recursive subroutine merge_sort_list_real4(vector,index,recursive_call)
      
      implicit none
      real*4,intent(inout)                :: vector(:)
      integer*4,intent(inout)             :: index(:)
      logical*4,intent(in),optional       :: recursive_call ! only used internally
      real*4,allocatable                  :: vector1(:),vector2(:)
      integer*4,allocatable               :: index1(:),index2(:)
      integer*4                           :: n0,n1,n2,i0,i1,i2
      
      n0 = size(vector)
      
      ! make index vector
      if (.not.present(recursive_call)) then
         do i0 = 1,n0
            index(i0) = i0
         end do
      end if
      
      ! split list in 2 sublists
      n1 = n0/2
      n2 = n0-n1
      allocate(vector1(n1),vector2(n2),index1(n1),index2(n2))
      vector1(1:n1) = vector(1:n1)
      vector2(1:n2) = vector(n0/2+1:n0)
      index1(1:n1) = index(1:n1)
      index2(1:n2) = index(n0/2+1:n0)
      
      if (n1>1) call merge_sort_list_real4(vector1,index1,recursive_call=.true.)
      if (n2>1) call merge_sort_list_real4(vector2,index2,recursive_call=.true.)
      
      ! merge sorted sublists
      i0 = 1
      i1 = 1
      i2 = 1
      do while ((i1<=n1) .and. (i2<=n2))
         if (vector1(i1)<vector2(i2)) then
            vector(i0) = vector1(i1)
            index(i0) = index1(i1)
            i1 = i1+1
         else
            vector(i0) = vector2(i2)
            index(i0) = index2(i2)
            i2 = i2+1
         end if
         i0 = i0+1
      end do
      if (i1<=n1) then
         vector(i0:i0+n1-i1) = vector1(i1:n1)
         index(i0:i0+n1-i1) = index1(i1:n1)
      else if (i2<=n2) then
         vector(i0:i0+n2-i2) = vector2(i2:n2)
         index(i0:i0+n2-i2) = index2(i2:n2)
      end if
      
   end subroutine merge_sort_list_real4
   
   recursive subroutine merge_sort_list_real8(vector,index,recursive_call)
      
      implicit none
      real*8,intent(inout)                :: vector(:)
      integer*4,intent(inout)             :: index(:)
      logical*4,intent(in),optional       :: recursive_call ! only used internally
      real*8,allocatable                  :: vector1(:),vector2(:)
      integer*4,allocatable               :: index1(:),index2(:)
      integer*4                           :: n0,n1,n2,i0,i1,i2
      
      n0 = size(vector)
      
      ! make index vector
      if (.not.present(recursive_call)) then
         do i0 = 1,n0
            index(i0) = i0
         end do
      end if
      
      ! split list in 2 sublists
      n1 = n0/2
      n2 = n0-n1
      allocate(vector1(n1),vector2(n2),index1(n1),index2(n2))
      vector1(1:n1) = vector(1:n1)
      vector2(1:n2) = vector(n0/2+1:n0)
      index1(1:n1) = index(1:n1)
      index2(1:n2) = index(n0/2+1:n0)
      
      if (n1>1) call merge_sort_list_real8(vector1,index1,recursive_call=.true.)
      if (n2>1) call merge_sort_list_real8(vector2,index2,recursive_call=.true.)
      
      ! merge sorted sublists
      i0 = 1
      i1 = 1
      i2 = 1
      do while ((i1<=n1) .and. (i2<=n2))
         if (vector1(i1)<vector2(i2)) then
            vector(i0) = vector1(i1)
            index(i0) = index1(i1)
            i1 = i1+1
         else
            vector(i0) = vector2(i2)
            index(i0) = index2(i2)
            i2 = i2+1
         end if
         i0 = i0+1
      end do
      if (i1<=n1) then
         vector(i0:i0+n1-i1) = vector1(i1:n1)
         index(i0:i0+n1-i1) = index1(i1:n1)
      else if (i2<=n2) then
         vector(i0:i0+n2-i2) = vector2(i2:n2)
         index(i0:i0+n2-i2) = index2(i2:n2)
      end if
      
   end subroutine merge_sort_list_real8

end module shared_module_sort