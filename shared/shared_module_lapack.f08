! **********************************************************************************************************************************
! Shared Fortran module to facilitate access to LAPACK math routines, requires LAPACK library to be included at compilation
! Developed by Danail Obreschkow
! **********************************************************************************************************************************

module shared_module_lapack

   private
   
   public   :: eigen ! returns eigen values and (optionally) eigenvectors of a symmetric matrix

contains

subroutine eigen(A,values,vectors)

   ! returns the eigenvalues (from smallest to largest)
   ! and optionally the eigenvectors vectors(:,1), vectors(:,2), ...
   ! requires LAPACK
   
   implicit none
   real*4,intent(in)             :: A(:,:)   ! symmetric NxN matrix
   real*4,intent(out)            :: values(:)
   real*4,intent(out),optional   :: vectors(:,:)
   real*4,allocatable            :: W(:)     ! eigenvalues
   real*4,allocatable            :: S(:,:)
   integer*4                     :: N        ! order of the matrix
   integer*4                     :: lwork,info
   integer*4,parameter           :: lwmax = 1000
   real*4                        :: work(lwmax)
   
   allocate(S(size(A,1),size(A,2)))
   S = A
   N = size(S,1)
   allocate(W(N))
   
   ! Query the optimal workspace.
   lwork = -1
   call sSYEV( 'Vectors', 'Upper', N, S, N, W, work, lwork, info)
   lwork = min(lwmax,int(work(1),4))
   
   ! Solve eigenproblem.
   call sSYEV( 'Vectors', 'Upper', N, S, N, W, work, lwork, info)
   
   ! Return results
   values = W
   if (present(vectors)) vectors = S
   
end subroutine eigen

end module shared_module_lapack