module module_hdf5_utilities

   use hdf5
   use module_system
   
   private 
   public   :: hdf5_open, hdf5_close, hdf5_read_dataset, hdf5_dataset_size

   interface hdf5_read_dataset
      module procedure read_dataset_0d_int4
      module procedure read_dataset_0d_int8
      module procedure read_dataset_0d_real4
      module procedure read_dataset_0d_real8
      module procedure read_dataset_1d_int4
      module procedure read_dataset_1d_int8
      module procedure read_dataset_1d_real4
      module procedure read_dataset_1d_real8
   end interface hdf5_read_dataset
   
   integer(hid_t) :: file_id
   
contains

   subroutine hdf5_open(filename)
   
      implicit none
      character(*),intent(in)    :: filename
      integer*4                  :: status

      if (exists(filename)) then
         call h5open_f(status) ! Initialize the Fortran interface
         call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,status) ! Open an hdf5 file
      end if
      
   end subroutine hdf5_open
   
   subroutine hdf5_close()
   
      implicit none
      integer*4                  :: status
      integer*4                  :: error
      call h5fclose_f(file_id,error) ! Terminate access to the file
      call h5close_f(status) ! Close the Fortran interface
      
   end subroutine hdf5_close

   integer*8 function hdf5_dataset_size(dataset)
   
      implicit none
      character(*),intent(in)    :: dataset
      integer*4                  :: error
      integer(hid_t)             :: dataset_id
      integer(hid_t)             :: type_id
      integer(hsize_t)           :: size,bytespervar
   
      ! determine number of galaxies
      call h5dopen_f(file_id, dataset, dataset_id, error) ! Open dataset
      call h5dget_storage_size_f(dataset_id, size, error)
      call h5dget_type_f(dataset_id, type_id, error)
      call h5tget_size_f(type_id, bytespervar, error)
      call h5dclose_f(dataset_id, error) ! Terminate access to the dataset
      hdf5_dataset_size = size/bytespervar
      
   end function hdf5_dataset_size
   
      subroutine read_dataset_0d_int4(dataset, dat)
   
      implicit none
      character(*),intent(in)    :: dataset
      integer*4,intent(inout)    :: dat
      integer*4                  :: error
      integer(hsize_t)           :: size(2)
      integer(hid_t)             :: dataset_id
      
      call h5dopen_f(file_id, dataset, dataset_id, error)
      call h5dread_f(dataset_id, get_mem_type_id(dataset_id,'int*4'), dat, size, error)
      call h5dclose_f(dataset_id, error)
   
   end subroutine read_dataset_0d_int4
   
   subroutine read_dataset_0d_int8(dataset, dat)
   
      implicit none
      character(*),intent(in)    :: dataset
      integer*8,intent(inout)    :: dat
      integer*4                  :: error
      integer(hsize_t)           :: size(2)
      integer(hid_t)             :: dataset_id
      
      call h5dopen_f(file_id, dataset, dataset_id, error)
      call h5dread_f(dataset_id, get_mem_type_id(dataset_id,'int*8'), dat, size, error)
      call h5dclose_f(dataset_id, error)
   
   end subroutine read_dataset_0d_int8
   
   subroutine read_dataset_0d_real4(dataset, dat)
   
      implicit none
      character(*),intent(in)    :: dataset
      real*4,intent(inout)       :: dat
      integer*4                  :: error
      integer(hsize_t)           :: size(2)
      integer(hid_t)             :: dataset_id
      
      call h5dopen_f(file_id, dataset, dataset_id, error)
      call h5dread_f(dataset_id, get_mem_type_id(dataset_id,'real*4'), dat, size, error)
      call h5dclose_f(dataset_id, error)
   
   end subroutine read_dataset_0d_real4
   
   subroutine read_dataset_0d_real8(dataset, dat)
   
      implicit none
      character(*),intent(in)    :: dataset
      real*8,intent(inout)       :: dat
      integer*4                  :: error
      integer(hsize_t)           :: size(2)
      integer(hid_t)             :: dataset_id
      
      call h5dopen_f(file_id, dataset, dataset_id, error)
      call h5dread_f(dataset_id, get_mem_type_id(dataset_id,'real*8'), dat, size, error)
      call h5dclose_f(dataset_id, error)
   
   end subroutine read_dataset_0d_real8
   
   subroutine read_dataset_1d_int4(dataset, dat)
   
      implicit none
      character(*),intent(in)    :: dataset
      integer*4,intent(inout)    :: dat(:)
      integer*4                  :: error
      integer(hsize_t)           :: size(2)
      integer(hid_t)             :: dataset_id
      
      call h5dopen_f(file_id, dataset, dataset_id, error)
      call h5dread_f(dataset_id, get_mem_type_id(dataset_id,'int*4'), dat, size, error)
      call h5dclose_f(dataset_id, error)
   
   end subroutine read_dataset_1d_int4
   
   subroutine read_dataset_1d_int8(dataset, dat)
   
      implicit none
      character(*),intent(in)    :: dataset
      integer*8,intent(inout)    :: dat(:)
      integer*4                  :: error
      integer(hsize_t)           :: size(2)
      integer(hid_t)             :: dataset_id
      
      call h5dopen_f(file_id, dataset, dataset_id, error)
      call h5dread_f(dataset_id, get_mem_type_id(dataset_id,'int*8'), dat, size, error)
      call h5dclose_f(dataset_id, error)
   
   end subroutine read_dataset_1d_int8
   
   subroutine read_dataset_1d_real4(dataset, dat)
   
      implicit none
      character(*),intent(in)    :: dataset
      real*4,intent(inout)       :: dat(:)
      integer*4                  :: error
      integer(hsize_t)           :: size(2)
      integer(hid_t)             :: dataset_id
      
      call h5dopen_f(file_id, dataset, dataset_id, error)
      call h5dread_f(dataset_id, get_mem_type_id(dataset_id,'real*4'), dat, size, error)
      call h5dclose_f(dataset_id, error)
   
   end subroutine read_dataset_1d_real4
   
   subroutine read_dataset_1d_real8(dataset, dat)
   
      implicit none
      character(*),intent(in)    :: dataset
      real*8,intent(inout)       :: dat(:)
      integer*4                  :: error
      integer(hsize_t)           :: size(2)
      integer(hid_t)             :: dataset_id
      
      call h5dopen_f(file_id, dataset, dataset_id, error)
      call h5dread_f(dataset_id, get_mem_type_id(dataset_id,'real*8'), dat, size, error)
      call h5dclose_f(dataset_id, error)
   
   end subroutine read_dataset_1d_real8
   
   integer(hid_t) function get_mem_type_id(dataset_id, expected_fortran_type)

      implicit none
   
      integer(hid_t),intent(in)        :: dataset_id
      character(*),intent(in),optional :: expected_fortran_type ! Must be 'int*4','int*8','real*4','real*8'
      character(6)                     :: fortran_type
      integer*4                        :: error, class, sign
      integer(size_t)                  :: size
      integer(hid_t)                   :: datatype_id

      call h5dget_type_f(dataset_id, datatype_id, error)
      call h5tget_class_f(datatype_id, class, error)
      call h5tget_size_f(datatype_id, size, error) 
      call h5tget_sign_f(datatype_id, sign, error)
      
      if (class == 0) then ! integer
         if (size == 4) then
            if (sign == 0) then
               get_mem_type_id = H5T_STD_U32LE
            else if (sign == 1) then
               get_mem_type_id = H5T_STD_I32LE
            else
               stop('ERROR: Unknown type.')
            end if
         else if (size == 8) then
            if (sign == 0) then
               get_mem_type_id = H5T_STD_U64LE
            else if (sign == 1) then
               get_mem_type_id = H5T_STD_I64LE
            else
               stop('ERROR: Unknown type.')
            end if
         else
            stop('ERROR: Unknown type.')
         end if
      else if (class == 1) then ! floating point
         if (size == 4) then
            get_mem_type_id = H5T_IEEE_F32LE
         else if (size == 8) then
            get_mem_type_id = H5T_IEEE_F64LE
         else
            stop('ERROR: Unknown type.')
         end if
      else
         stop('ERROR: Unknown type.')
      end if
      
      if (present(expected_fortran_type)) then
         if (class == 0) then
            write(fortran_type,'(A,I0)') 'int*',size
         else
            write(fortran_type,'(A,I0)') 'real*',size
         end if
         if (trim(fortran_type).ne.expected_fortran_type) stop('HDF5: Data type does not match expectation.')
      end if

   end function get_mem_type_id
   
end module module_hdf5_utilities