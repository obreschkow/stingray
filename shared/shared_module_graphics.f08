! **********************************************************************************************************************************
! Shared Fortran module with graphics routines
! Developed by Danail Obreschkow
! **********************************************************************************************************************************

module shared_module_graphics

   private
   
   public   :: raster2bitmap  ! converts and n-m-3 rgb array into a bitmap file with basic header information

contains

subroutine raster2bitmap(rgb,filename,unit)

   ! converts rgb(1:height,1:width,1:3) [0...1] to a 24-bit bitmap file

   implicit none
   real*4,intent(in)             :: rgb(:,:,:)
   character(*),intent(in)       :: filename
   integer*4,intent(in),optional :: unit ! I/O unit
   integer*4                     :: width,height
   integer*4                     :: extrabytes
   integer*4                     :: paddedsize
   integer*4                     :: i,j,k
   integer*4                     :: fileunit
   
   height = size(rgb,1)
   width = size(rgb,2)
   extrabytes = 4-modulo(width*3,4) ! number of bytes of padding to add to each horizontal line
   if (extrabytes == 4) extrabytes = 0
   paddedsize = ((width*3)+extrabytes)*height

   if (present(unit)) then
      fileunit = unit
   else
      fileunit = 1
   end if

   open(fileunit,file=trim(filename),action='write',form='unformatted',status='replace',access='stream')
   
   ! write header
   write(fileunit) 'BM',paddedsize+54,0,54 ! BMP header
   write(fileunit) 40,width,height,achar(1),achar(0),achar(24),achar(0),0,paddedsize,2835,2835,0,0 ! DIB header
   
   ! write array
   do j = 1,width
      do i = 1,height
         do k = 3,1,-1
            write(fileunit) achar(nint(rgb(i,j,k)*255))
         end do
      end do
      do i = 1,extrabytes
         write(fileunit) achar(0)
      end do
   end do
   
   close(fileunit)

end subroutine raster2bitmap

end module shared_module_graphics