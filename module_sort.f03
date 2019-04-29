module module_sort

   public

contains

   recursive subroutine merge_sort_list(list,level)
   
      ! sorts n rows of the n-by-2 matrix list(1:n,1:2) in increasing order of the first column
      ! for nested parallelization, the variable omp_set_nested must be set to true
      
      implicit none
      integer*8,intent(inout)       :: list(:,:)
      integer,intent(in),optional   :: level
      integer*8,allocatable         :: list1(:,:),list2(:,:)
      integer*4                     :: n0,n1,n2,i0,i1,i2
      integer*4                     :: next_level
      
      if (present(level)) then
         next_level = level+1
      else
         next_level = 1
         !$ call omp_set_nested(.true.)
      end if
      
      n0 = size(list(:,1))
      
      ! split list in 2 sublists
      n1 = n0/2
      n2 = n0-n1
      allocate(list1(n1,2),list2(n2,2))
      list1(1:n1,:) = list(1:n1,:)
      list2(1:n2,:) = list(n0/2+1:n0,:)
      
      if (n1>1) call merge_sort_list(list1,next_level)
      if (n2>1) call merge_sort_list(list2,next_level)
      
      ! merge sorted sublists
      i0 = 1
      i1 = 1
      i2 = 1
      do while ((i1<=n1) .and. (i2<=n2))
         if (list1(i1,1)<list2(i2,1)) then
            list(i0,:) = list1(i1,:)
            i1 = i1+1
         else
            list(i0,:) = list2(i2,:)
            i2 = i2+1
         end if
         i0 = i0+1
      end do
      if (i1<=n1) then
         list(i0:i0+n1-i1,:) = list1(i1:n1,:)
      else if (i2<=n2) then
         list(i0:i0+n2-i2,:) = list2(i2:n2,:)
      end if
      
   end subroutine merge_sort_list

end module module_sort