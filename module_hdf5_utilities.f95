module module_hdf5_utilities

   use hdf5
   use module_system
   
   private 
   public   :: hdf5_create, hdf5_open, hdf5_close
   public   :: hdf5_dataset_size
   public   :: hdf5_add_group
   public   :: hdf5_read_data, hdf5_write_data

   interface hdf5_read_data
      module procedure read_dataset_0d_int4
      module procedure read_dataset_0d_int8
      module procedure read_dataset_0d_real4
      module procedure read_dataset_0d_real8
      module procedure read_dataset_1d_int4
      module procedure read_dataset_1d_int8
      module procedure read_dataset_1d_real4
      module procedure read_dataset_1d_real8
   end interface hdf5_read_data
   
   interface hdf5_write_data
      module procedure write_dataset_0d_int4
      module procedure write_dataset_0d_int8
      module procedure write_dataset_0d_real4
      module procedure write_dataset_0d_real8
      module procedure write_dataset_1d_int4
      module procedure write_dataset_1d_int8
      module procedure write_dataset_1d_real4
      module procedure write_dataset_1d_real8
   end interface hdf5_write_data
   
   integer(hid_t) :: file_id
   
contains

   subroutine hdf5_open(filename) ! opens an existing HDF5 file for read and write
   
      implicit none
      character(*),intent(in)    :: filename
      integer*4                  :: status

      if (exists(filename)) then
         call h5open_f(status) ! Initialize the Fortran interface
         call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,status) ! Open an hdf5 file
      end if
      
   end subroutine hdf5_open
   
   subroutine hdf5_create(filename) ! creates a new HDF5 file and closes it
   
      implicit none
      character(*),intent(in)    :: filename
      integer*4                  :: status

      call h5open_f(status) ! Initialize the Fortran interface
      call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,status) ! Create an hdf5 file
      call hdf5_close()
      
   end subroutine hdf5_create
   
   subroutine hdf5_add_group(groupname) ! adds a new group to the open file (=> hdf5_open must be called first)
   
      implicit none
      character(*),intent(in)    :: groupname
      integer(hid_t)             :: group_id
      integer*4                  :: error
   
      call h5gcreate_f(file_id, groupname, group_id, error)
      call h5gclose_f(group_id, error)
      
   end subroutine hdf5_add_group
   
   subroutine hdf5_close()
   
      implicit none
      integer*4                  :: status
      integer*4                  :: err
      call h5fclose_f(file_id,err) ! Terminate access to the file
      call h5close_f(status) ! Close the Fortran interface
      
   end subroutine hdf5_close

   integer*8 function hdf5_dataset_size(dataset)
   
      implicit none
      character(*),intent(in)    :: dataset
      integer*4                  :: err
      integer(hid_t)             :: dataset_id
      integer(hid_t)             :: type_id
      integer(hsize_t)           :: size,bytespervar
   
      ! determine number of galaxies
      call h5dopen_f(file_id, dataset, dataset_id, err) ! Open dataset
      call h5dget_storage_size_f(dataset_id, size, err)
      call h5dget_type_f(dataset_id, type_id, err)
      call h5tget_size_f(type_id, bytespervar, err)
      call h5dclose_f(dataset_id, err) ! Terminate access to the dataset
      hdf5_dataset_size = size/bytespervar
      
   end function hdf5_dataset_size
   
   ! ===========================================================================================================
   ! subroutines used by interface hdf5_write_dataset
   ! ===========================================================================================================
   
   subroutine write_dataset_0d_int4(datasetname,dat,attribute)
   
      implicit none
      character(*),intent(in)          :: datasetname
      integer*4,intent(in)             :: dat
      character(*),intent(in),optional :: attribute ! optional explanation
      
      call write_dataset_1d_int4(datasetname,(/dat/),attribute)
      
   end subroutine write_dataset_0d_int4
   
   subroutine write_dataset_0d_int8(datasetname,dat,attribute)
   
      implicit none
      character(*),intent(in)          :: datasetname
      integer*8,intent(in)             :: dat
      character(*),intent(in),optional :: attribute ! optional explanation
      
      call write_dataset_1d_int8(datasetname,(/dat/),attribute)
      
   end subroutine write_dataset_0d_int8
   
   subroutine write_dataset_0d_real4(datasetname,dat,attribute)
   
      implicit none
      character(*),intent(in)          :: datasetname
      real*4,intent(in)                :: dat
      character(*),intent(in),optional :: attribute ! optional explanation
      
      call write_dataset_1d_real4(datasetname,(/dat/),attribute)
      
   end subroutine write_dataset_0d_real4
   
   subroutine write_dataset_0d_real8(datasetname,dat,attribute)
   
      implicit none
      character(*),intent(in)          :: datasetname
      real*8,intent(in)                :: dat
      character(*),intent(in),optional :: attribute ! optional explanation
      
      call write_dataset_1d_real8(datasetname,(/dat/),attribute)
      
   end subroutine write_dataset_0d_real8
   
   subroutine write_dataset_1d_int4(datasetname,dat,attribute)
   
      implicit none
      character(*),intent(in)          :: datasetname
      integer*4,intent(in)             :: dat(:)
      character(*),intent(in),optional :: attribute ! optional explanation
      integer(hid_t)                   :: space_id
      integer(hid_t)                   :: dataset_id
      integer*4                        :: error
      integer*4,parameter              :: rank = 1
      integer*8                        :: n(1)
      
      n(1) = size(dat,1)
      call h5screate_simple_f(rank,n,space_id,error)  ! Create the space for the dataset
      call h5dcreate_f(file_id,datasetname,H5T_STD_I32LE,space_id,dataset_id,error) ! create dataset
      call h5dwrite_f(dataset_id,H5T_STD_I32LE,dat,n,error) ! write data set
      call h5sclose_f(space_id,error) ! close data space
      if (present(attribute)) call write_attribute(attribute,dataset_id)
      call h5dclose_f(dataset_id,error) ! close data set
      
   end subroutine write_dataset_1d_int4
   
   subroutine write_dataset_1d_int8(datasetname,dat,attribute)
   
      implicit none
      character(*),intent(in)          :: datasetname
      integer*8,intent(in)             :: dat(:)
      character(*),intent(in),optional :: attribute ! optional explanation
      integer(hid_t)                   :: space_id
      integer(hid_t)                   :: dataset_id
      integer*4                        :: error
      integer*4,parameter              :: rank = 1
      integer*8                        :: n(1)
      
      n(1) = size(dat,1)
      call h5screate_simple_f(rank,n,space_id,error)  ! Create the space for the dataset
      call h5dcreate_f(file_id,datasetname,H5T_STD_I64LE,space_id,dataset_id,error) ! create dataset
      call h5dwrite_f(dataset_id,H5T_STD_I64LE,dat,n,error) ! write data set
      call h5sclose_f(space_id,error) ! close data space
      if (present(attribute)) call write_attribute(attribute,dataset_id)
      call h5dclose_f(dataset_id,error) ! close data set
      
   end subroutine write_dataset_1d_int8
           
   subroutine write_dataset_1d_real4(datasetname,dat,attribute)
   
      implicit none
      character(*),intent(in)          :: datasetname
      real*4,intent(in)                :: dat(:)
      character(*),intent(in),optional :: attribute ! optional explanation
      integer(hid_t)                   :: space_id
      integer(hid_t)                   :: dataset_id
      integer*4                        :: error
      integer*4,parameter              :: rank = 1
      integer*8                        :: n(1)
      
      n(1) = size(dat,1)
      call h5screate_simple_f(rank,n,space_id,error)  ! Create the space for the dataset
      call h5dcreate_f(file_id,datasetname,H5T_IEEE_F32LE,space_id,dataset_id,error) ! create dataset
      call h5dwrite_f(dataset_id,H5T_IEEE_F32LE,dat,n,error) ! write data set
      call h5sclose_f(space_id,error) ! close data space
      if (present(attribute)) call write_attribute(attribute,dataset_id)
      call h5dclose_f(dataset_id,error) ! close data set
      
   end subroutine write_dataset_1d_real4
   
   subroutine write_dataset_1d_real8(datasetname,dat,attribute)
   
      implicit none
      character(*),intent(in)          :: datasetname
      real*8,intent(in)                :: dat(:)
      character(*),intent(in),optional :: attribute ! optional explanation
      integer(hid_t)                   :: space_id
      integer(hid_t)                   :: dataset_id
      integer*4                        :: error
      integer*4,parameter              :: rank = 1
      integer*8                        :: n(1)
      
      n(1) = size(dat,1)
      call h5screate_simple_f(rank,n,space_id,error)  ! Create the space for the dataset
      call h5dcreate_f(file_id,datasetname,H5T_IEEE_F64LE,space_id,dataset_id,error) ! create dataset
      call h5dwrite_f(dataset_id,H5T_IEEE_F64LE,dat,n,error) ! write data set
      call h5sclose_f(space_id,error) ! close data space
      if (present(attribute)) call write_attribute(attribute,dataset_id)
      call h5dclose_f(dataset_id,error) ! close data set
      
   end subroutine write_dataset_1d_real8
   
   subroutine write_attribute(attribute,dataset_id)
   
      implicit none
      character(*),intent(in)       :: attribute
      integer(hid_t),intent(in)     :: dataset_id
      integer(hid_t)                :: attr_id ! Attribute identifier
      integer(hid_t)                :: aspace_id ! Attribute Dataspace identifier
      integer(hid_t)                :: atype_id ! Attribute Dataspace identifier
      integer(hid_t),dimension(1)   :: adims = (/1/) ! Attribute dimension
      integer*4,parameter           :: arank = 1 ! Attribure rank
      integer(size_t)               :: attrlen ! Length of the attribute string
      character(*),parameter        :: aname = "Comment" ! Attribute name
      integer(hsize_t),dimension(1) :: data_dims = (/1/)
      integer*4                     :: error
      
      attrlen = len(attribute)
      call h5screate_simple_f(arank,adims,aspace_id,error) ! Create scalar data space for the attribute.
      call h5tcopy_f(H5T_NATIVE_CHARACTER,atype_id,error) ! Create datatype for the attribute.
      call h5tset_size_f(atype_id,attrlen,error)
      call h5acreate_f(dataset_id,aname,atype_id,aspace_id,attr_id,error)
      call h5awrite_f(attr_id,atype_id,attribute,data_dims,error) ! Write the attribute data.
      call h5aclose_f(attr_id,error) ! Close the attribute.
      
   end subroutine write_attribute
   
   
   ! ===========================================================================================================
   ! subroutines used by interface hdf5_read_dataset
   ! ===========================================================================================================
   
   subroutine read_dataset_0d_int4(dataset, dat)
   
      implicit none
      character(*),intent(in)    :: dataset
      integer*4,intent(inout)    :: dat
      integer*4                  :: err
      integer(hsize_t)           :: size(2)
      integer(hid_t)             :: dataset_id
      
      call h5dopen_f(file_id, dataset, dataset_id, err)
      call h5dread_f(dataset_id, get_mem_type_id(dataset_id,'int*4'), dat, size, err)
      call h5dclose_f(dataset_id, err)
   
   end subroutine read_dataset_0d_int4
   
   subroutine read_dataset_0d_int8(dataset, dat)
   
      implicit none
      character(*),intent(in)    :: dataset
      integer*8,intent(inout)    :: dat
      integer*4                  :: err
      integer(hsize_t)           :: size(2)
      integer(hid_t)             :: dataset_id
      
      call h5dopen_f(file_id, dataset, dataset_id, err)
      call h5dread_f(dataset_id, get_mem_type_id(dataset_id,'int*8'), dat, size, err)
      call h5dclose_f(dataset_id, err)
   
   end subroutine read_dataset_0d_int8
   
   subroutine read_dataset_0d_real4(dataset, dat)
   
      implicit none
      character(*),intent(in)    :: dataset
      real*4,intent(inout)       :: dat
      integer*4                  :: err
      integer(hsize_t)           :: size(2)
      integer(hid_t)             :: dataset_id
      
      call h5dopen_f(file_id, dataset, dataset_id, err)
      call h5dread_f(dataset_id, get_mem_type_id(dataset_id,'real*4'), dat, size, err)
      call h5dclose_f(dataset_id, err)
   
   end subroutine read_dataset_0d_real4
   
   subroutine read_dataset_0d_real8(dataset, dat)
   
      implicit none
      character(*),intent(in)    :: dataset
      real*8,intent(inout)       :: dat
      integer*4                  :: err
      integer(hsize_t)           :: size(2)
      integer(hid_t)             :: dataset_id
      
      call h5dopen_f(file_id, dataset, dataset_id, err)
      call h5dread_f(dataset_id, get_mem_type_id(dataset_id,'real*8'), dat, size, err)
      call h5dclose_f(dataset_id, err)
   
   end subroutine read_dataset_0d_real8
   
   subroutine read_dataset_1d_int4(dataset, dat)
   
      implicit none
      character(*),intent(in)    :: dataset
      integer*4,intent(inout)    :: dat(:)
      integer*4                  :: err
      integer(hsize_t)           :: size(2)
      integer(hid_t)             :: dataset_id
      
      call h5dopen_f(file_id, dataset, dataset_id, err)
      call h5dread_f(dataset_id, get_mem_type_id(dataset_id,'int*4'), dat, size, err)
      call h5dclose_f(dataset_id, err)
   
   end subroutine read_dataset_1d_int4
   
   subroutine read_dataset_1d_int8(dataset, dat)
   
      implicit none
      character(*),intent(in)    :: dataset
      integer*8,intent(inout)    :: dat(:)
      integer*4                  :: err
      integer(hsize_t)           :: size(2)
      integer(hid_t)             :: dataset_id
      
      call h5dopen_f(file_id, dataset, dataset_id, err)
      call h5dread_f(dataset_id, get_mem_type_id(dataset_id,'int*8'), dat, size, err)
      call h5dclose_f(dataset_id, err)
   
   end subroutine read_dataset_1d_int8
   
   subroutine read_dataset_1d_real4(dataset, dat)
   
      implicit none
      character(*),intent(in)    :: dataset
      real*4,intent(inout)       :: dat(:)
      integer*4                  :: err
      integer(hsize_t)           :: size(2)
      integer(hid_t)             :: dataset_id
      
      call h5dopen_f(file_id, dataset, dataset_id, err)
      call h5dread_f(dataset_id, get_mem_type_id(dataset_id,'real*4'), dat, size, err)
      call h5dclose_f(dataset_id, err)
   
   end subroutine read_dataset_1d_real4
   
   subroutine read_dataset_1d_real8(dataset, dat)
   
      implicit none
      character(*),intent(in)    :: dataset
      real*8,intent(inout)       :: dat(:)
      integer*4                  :: err
      integer(hsize_t)           :: size(2)
      integer(hid_t)             :: dataset_id
      
      call h5dopen_f(file_id, dataset, dataset_id, err)
      call h5dread_f(dataset_id, get_mem_type_id(dataset_id,'real*8'), dat, size, err)
      call h5dclose_f(dataset_id, err)
   
   end subroutine read_dataset_1d_real8
   
   integer(hid_t) function get_mem_type_id(dataset_id, expected_fortran_type)

      implicit none
   
      integer(hid_t),intent(in)        :: dataset_id
      character(*),intent(in),optional :: expected_fortran_type ! Must be 'int*4','int*8','real*4','real*8'
      character(6)                     :: fortran_type
      integer*4                        :: err, class, sign
      integer(size_t)                  :: size
      integer(hid_t)                   :: datatype_id

      call h5dget_type_f(dataset_id, datatype_id, err)
      call h5tget_class_f(datatype_id, class, err)
      call h5tget_size_f(datatype_id, size, err) 
      call h5tget_sign_f(datatype_id, sign, err)
      
      if (class == 0) then ! integer
         if (size == 4) then
            if (sign == 0) then
               get_mem_type_id = H5T_STD_U32LE
            else if (sign == 1) then
               get_mem_type_id = H5T_STD_I32LE
            else
               call error('Unknown type')
            end if
         else if (size == 8) then
            if (sign == 0) then
               get_mem_type_id = H5T_STD_U64LE
            else if (sign == 1) then
               get_mem_type_id = H5T_STD_I64LE
            else
               call error('Unknown type')
            end if
         else
            call error('Unknown type')
         end if
      else if (class == 1) then ! floating point
         if (size == 4) then
            get_mem_type_id = H5T_IEEE_F32LE
         else if (size == 8) then
            get_mem_type_id = H5T_IEEE_F64LE
         else
            call error('Unknown type')
         end if
      else
         call error('Unknown type')
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