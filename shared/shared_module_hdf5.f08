! **********************************************************************************************************************************
! Shared Fortran module to read and write HDF5 data, using the hdf5 library
! Developed by Danail Obreschkow
!
! Additional useful resources:
! https://support.hdfgroup.org/HDF5/examples/api-fortran.html
! https://github.com/tiasus/HDF5_utils
! **********************************************************************************************************************************
   
module shared_module_hdf5

   use hdf5
   use shared_module_core
   
   private
   
   public   :: hdf5_create             ! routine to create a new HDF5 file
   public   :: hdf5_open               ! routine to open an existing HDF5 file for read or write access
   public   :: hdf5_close              ! routine to close the currently open HDF5 file
   public   :: hdf5_get_dataset_size   ! routine to determine rank and dimensions of a dataset
   public   :: hdf5_add_group          ! routine to add a group to to the currently open HDF5 file
   public   :: hdf5_object_exists      ! logical function to check if a group or dataset exits in the currently open HDF5 file
   public   :: hdf5_read_data          ! routine to read dataset of strings or numerical scalars/vectors/matrices
   public   :: hdf5_write_data         ! routine to write dataset of string or numerical scalar/vector/matrix, with opt attribute 

   interface hdf5_read_data
      module procedure read_dataset_0d_string
      module procedure read_dataset_0d_int4
      module procedure read_dataset_0d_int8
      module procedure read_dataset_0d_real4
      module procedure read_dataset_0d_real8
      module procedure read_dataset_1d_int4
      module procedure read_dataset_1d_int8
      module procedure read_dataset_1d_real4
      module procedure read_dataset_1d_real8
      module procedure read_dataset_2d_int4
      module procedure read_dataset_2d_int8
      module procedure read_dataset_2d_real4
      module procedure read_dataset_2d_real8
   end interface hdf5_read_data
   
   interface hdf5_write_data
      module procedure write_dataset_0d_string
      module procedure write_dataset_0d_int4
      module procedure write_dataset_0d_int8
      module procedure write_dataset_0d_real4
      module procedure write_dataset_0d_real8
      module procedure write_dataset_0d_log
      module procedure write_dataset_1d_int4
      module procedure write_dataset_1d_int8
      module procedure write_dataset_1d_real4
      module procedure write_dataset_1d_real8
      module procedure write_dataset_1d_log
      module procedure write_dataset_2d_int4
      module procedure write_dataset_2d_int8
      module procedure write_dataset_2d_real4
      module procedure write_dataset_2d_real8
      module procedure write_dataset_2d_log
   end interface hdf5_write_data
   
   integer(hid_t) :: file_id = 0
   
contains

subroutine hdf5_open(filename,write_access) ! opens an existing HDF5 file for read and write

   implicit none
   character(*),intent(in)       :: filename
   logical*4,optional,intent(in) :: write_access  ! if set to true, the file is opened for writing, otherwise for reading
   integer*4                     :: status
   logical*4                     :: wa
   
   wa = .false.
   if (present(write_access)) wa = write_access
   
   call h5open_f(status) ! Initialize the Fortran interface
   if (wa) then
      call check_file(filename,'rw')
      call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,status)
   else
      call check_file(filename,'r')
      call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,status)
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
   integer*4                  :: err
   integer*4                  :: i
   character(255)             :: parentgroup

   if (file_id==0) call error('no HDF5 file open to add group '//trim(groupname))
   
   do i = len(trim(groupname))-1,2,-1
      if (groupname(i:i)=='/') then
         parentgroup = groupname(1:i)
         if (.not.hdf5_object_exists(parentgroup)) then
            call error('could not find HDF5 parent group "'//trim(parentgroup)//'" to add group "'//trim(groupname)//'"')
         end if
         exit
      end if
   end do
   
   call h5gcreate_f(file_id, groupname, group_id, err)
   call h5gclose_f(group_id, err)
   
end subroutine hdf5_add_group

subroutine hdf5_close()

   implicit none
   integer*4                  :: status
   integer*4                  :: err
   
   if (file_id==0) call error('no HDF5 file open to close')
   call h5fclose_f(file_id,err) ! Terminate access to the file
   call h5close_f(status) ! Close the Fortran interface
   file_id = 0
   
end subroutine hdf5_close

subroutine hdf5_get_dataset_size(dataset, bytes_per_element, n_elements, rank, dimensions)

   implicit none
   character(*),intent(in)                      :: dataset
   integer*4,intent(out),optional               :: bytes_per_element
   integer*4,intent(out),optional               :: n_elements
   integer*4,intent(out),optional               :: rank
   integer*4,intent(out),allocatable,optional   :: dimensions(:)
   integer*4                                    :: my_rank
   integer*4,allocatable                        :: my_dimensions(:)
   integer(hid_t)                               :: dataset_id, dataspace_id, datatype_id
   integer(hsize_t)                             :: dset_dims(6), max_dims(6), size
   integer*4                                    :: hdferror
   
   ! checks
   if (file_id==0) call error('no HDF5 file open to read dataset '//trim(dataset))
   if (.not.hdf5_object_exists(dataset)) call error('could not find HDF5 dataset '//trim(dataset))

   ! open dataset
   call h5dopen_f(file_id, dataset, dataset_id, hdferror)

   ! get bytes per element
   if (present(bytes_per_element)) then
      call h5dget_type_f(dataset_id, datatype_id, hdferror)
      call h5tget_size_f(datatype_id, size, hdferror)
      bytes_per_element = int(size,4)
   end if

   ! get rank (ndims)
   call h5dget_space_f(dataset_id, dataspace_id, hdferror)
   call h5sget_simple_extent_ndims_f(dataspace_id, my_rank, hdferror)
   if (present(rank)) rank = my_rank

   ! get dims
   if (allocated(my_dimensions)) deallocate(my_dimensions)
   allocate(my_dimensions(my_rank))
   if (my_rank>0) then
      call h5sget_simple_extent_dims_f(dataspace_id, dset_dims(1:my_rank), max_dims(1:my_rank), hdferror)
      my_dimensions = int(dset_dims(1:my_rank),4)
   end if
   if (present(dimensions)) dimensions = my_dimensions
   
   ! get number of elements
   if (present(n_elements)) then
      if (my_rank>0) then
         n_elements = product(my_dimensions)
      else
         n_elements = 1
      end if
   end if

   ! close ids
   call h5sclose_f(dataspace_id, hdferror)  
   call h5dclose_f(dataset_id, hdferror)

end subroutine hdf5_get_dataset_size

function hdf5_object_exists(object) result(exists)

   implicit none
   character(*),intent(in) :: object ! name of group or dataset
   logical                 :: exists
   integer*4               :: hdferror
   integer*4               :: i
   character(255)          :: group
   
   do i = 1,len(trim(object))
      if (object(i:i)=='/') then
         group = object(1:i)
         call h5lexists_f(file_id, trim(group), exists, hdferror)
         if (.not.exists) return
      end if
   end do
   
   call h5lexists_f(file_id, object, exists, hdferror)

end function hdf5_object_exists

! ===========================================================================================================
! subroutines used by interface hdf5_write_dataset
! ===========================================================================================================

subroutine write_dataset_0d_string(dataset,dat,attribute)

   implicit none
   character(*),intent(in)          :: dataset
   character(*),intent(in)          :: dat
   character(*),intent(in),optional :: attribute ! optional explanation
   integer(hid_t)                   :: space_id
   integer(hid_t)                   :: dataset_id,filetype
   integer*4                        :: status,hdferr
   integer*4,parameter              :: rank = 1
   integer*8                        :: n(1)
   
   call check_file_for_writing(dataset)
   
   n(1) = max(1,len(dat))
   call h5screate_simple_f(rank,(/1_8/),space_id,status)  ! Create the space for the dataset
   call H5Tcopy_f(H5T_FORTRAN_S1, filetype, hdferr)
   call H5Tset_size_f(filetype, n(1), hdferr)
   call h5dcreate_f(file_id,dataset,filetype,space_id,dataset_id,status) ! create dataset
   call h5dwrite_f(dataset_id,filetype,dat,n,status) ! write data set
   call h5sclose_f(space_id,status) ! close data space
   if (present(attribute)) call write_attribute(attribute,dataset_id)
   call h5dclose_f(dataset_id,status) ! close data set
   
end subroutine write_dataset_0d_string

subroutine write_dataset_0d_int4(dataset,dat,attribute)

   implicit none
   character(*),intent(in)          :: dataset
   integer*4,intent(in)             :: dat
   character(*),intent(in),optional :: attribute ! optional explanation
   
   call write_dataset_1d_int4(dataset,(/dat/),attribute)
   
end subroutine write_dataset_0d_int4

subroutine write_dataset_0d_log(dataset,dat,attribute)

   implicit none
   character(*),intent(in)          :: dataset
   logical,intent(in)               :: dat
   character(*),intent(in),optional :: attribute ! optional explanation
   
   call write_dataset_1d_int4(dataset,(/log2int(dat)/),attribute)
   
end subroutine write_dataset_0d_log

subroutine write_dataset_0d_int8(dataset,dat,attribute)

   implicit none
   character(*),intent(in)          :: dataset
   integer*8,intent(in)             :: dat
   character(*),intent(in),optional :: attribute ! optional explanation
   
   call write_dataset_1d_int8(dataset,(/dat/),attribute)
   
end subroutine write_dataset_0d_int8

subroutine write_dataset_0d_real4(dataset,dat,attribute)

   implicit none
   character(*),intent(in)          :: dataset
   real*4,intent(in)                :: dat
   character(*),intent(in),optional :: attribute ! optional explanation
   
   call write_dataset_1d_real4(dataset,(/dat/),attribute)
   
end subroutine write_dataset_0d_real4

subroutine write_dataset_0d_real8(dataset,dat,attribute)

   implicit none
   character(*),intent(in)          :: dataset
   real*8,intent(in)                :: dat
   character(*),intent(in),optional :: attribute ! optional explanation
   
   call write_dataset_1d_real8(dataset,(/dat/),attribute)
   
end subroutine write_dataset_0d_real8

subroutine write_dataset_1d_int4(dataset,dat,attribute)

   implicit none
   character(*),intent(in)          :: dataset
   integer*4,intent(in)             :: dat(:)
   character(*),intent(in),optional :: attribute ! optional explanation
   integer(hid_t)                   :: space_id
   integer(hid_t)                   :: dataset_id
   integer*4                        :: status
   integer*4,parameter              :: rank = 1
   integer*8                        :: n(1)
   
   call check_file_for_writing(dataset)
   
   n(1) = size(dat,1)
   call h5screate_simple_f(rank,n,space_id,status)  ! Create the space for the dataset
   call h5dcreate_f(file_id,dataset,H5T_STD_I32LE,space_id,dataset_id,status) ! create dataset
   call h5dwrite_f(dataset_id,H5T_STD_I32LE,dat,n,status) ! write data set
   call h5sclose_f(space_id,status) ! close data space
   if (present(attribute)) call write_attribute(attribute,dataset_id)
   call h5dclose_f(dataset_id,status) ! close data set
   
end subroutine write_dataset_1d_int4

subroutine write_dataset_1d_log(dataset,dat,attribute)

   implicit none
   character(*),intent(in)          :: dataset
   logical,intent(in)               :: dat(:)
   character(*),intent(in),optional :: attribute ! optional explanation
   
   call write_dataset_1d_int4(dataset,log2int(dat),attribute)
   
end subroutine write_dataset_1d_log

subroutine write_dataset_1d_int8(dataset,dat,attribute)

   implicit none
   character(*),intent(in)          :: dataset
   integer*8,intent(in)             :: dat(:)
   character(*),intent(in),optional :: attribute ! optional explanation
   integer(hid_t)                   :: space_id
   integer(hid_t)                   :: dataset_id
   integer*4                        :: status
   integer*4,parameter              :: rank = 1
   integer*8                        :: n(1)
   
   call check_file_for_writing(dataset)
   
   n(1) = size(dat,1)
   call h5screate_simple_f(rank,n,space_id,status)  ! Create the space for the dataset
   call h5dcreate_f(file_id,dataset,H5T_STD_I64LE,space_id,dataset_id,status) ! create dataset
   call h5dwrite_f(dataset_id,H5T_STD_I64LE,dat,n,status) ! write data set
   call h5sclose_f(space_id,status) ! close data space
   if (present(attribute)) call write_attribute(attribute,dataset_id)
   call h5dclose_f(dataset_id,status) ! close data set
   
end subroutine write_dataset_1d_int8
        
subroutine write_dataset_1d_real4(dataset,dat,attribute)

   implicit none
   character(*),intent(in)          :: dataset
   real*4,intent(in)                :: dat(:)
   character(*),intent(in),optional :: attribute ! optional explanation
   integer(hid_t)                   :: space_id
   integer(hid_t)                   :: dataset_id
   integer*4                        :: status
   integer*4,parameter              :: rank = 1
   integer*8                        :: n(1)
   
   call check_file_for_writing(dataset)
   
   n(1) = size(dat,1)
   call h5screate_simple_f(rank,n,space_id,status)  ! Create the space for the dataset
   call h5dcreate_f(file_id,dataset,H5T_IEEE_F32LE,space_id,dataset_id,status) ! create dataset
   call h5dwrite_f(dataset_id,H5T_IEEE_F32LE,dat,n,status) ! write data set
   call h5sclose_f(space_id,status) ! close data space
   if (present(attribute)) call write_attribute(attribute,dataset_id)
   call h5dclose_f(dataset_id,status) ! close data set
   
end subroutine write_dataset_1d_real4

subroutine write_dataset_1d_real8(dataset,dat,attribute)

   implicit none
   character(*),intent(in)          :: dataset
   real*8,intent(in)                :: dat(:)
   character(*),intent(in),optional :: attribute ! optional explanation
   integer(hid_t)                   :: space_id
   integer(hid_t)                   :: dataset_id
   integer*4                        :: status
   integer*4,parameter              :: rank = 1
   integer*8                        :: n(1)
   
   call check_file_for_writing(dataset)
   
   n(1) = size(dat,1)
   call h5screate_simple_f(rank,n,space_id,status)  ! Create the space for the dataset
   call h5dcreate_f(file_id,dataset,H5T_IEEE_F64LE,space_id,dataset_id,status) ! create dataset
   call h5dwrite_f(dataset_id,H5T_IEEE_F64LE,dat,n,status) ! write data set
   call h5sclose_f(space_id,status) ! close data space
   if (present(attribute)) call write_attribute(attribute,dataset_id)
   call h5dclose_f(dataset_id,status) ! close data set
   
end subroutine write_dataset_1d_real8

subroutine write_dataset_2d_int4(dataset,dat,attribute)

   implicit none
   character(*),intent(in)          :: dataset
   integer*4,intent(in)             :: dat(:,:)
   character(*),intent(in),optional :: attribute ! optional explanation
   integer(hid_t)                   :: space_id
   integer(hid_t)                   :: dataset_id
   integer*4                        :: status
   integer*4,parameter              :: rank = 2
   integer*8                        :: n(2)
   
   call check_file_for_writing(dataset)
   
   n(1) = size(dat,1)
   n(2) = size(dat,2)
   call h5screate_simple_f(rank,n,space_id,status)  ! Create the space for the dataset
   call h5dcreate_f(file_id,dataset,H5T_STD_I32LE,space_id,dataset_id,status) ! create dataset
   call h5dwrite_f(dataset_id,H5T_STD_I32LE,dat,n,status) ! write data set
   call h5sclose_f(space_id,status) ! close data space
   if (present(attribute)) call write_attribute(attribute,dataset_id)
   call h5dclose_f(dataset_id,status) ! close data set
   
end subroutine write_dataset_2d_int4

subroutine write_dataset_2d_log(dataset,dat,attribute)

   implicit none
   character(*),intent(in)          :: dataset
   logical,intent(in)               :: dat(:,:)
   character(*),intent(in),optional :: attribute ! optional explanation
   
   call write_dataset_2d_int4(dataset,log2int(dat),attribute)
   
end subroutine write_dataset_2d_log


subroutine write_dataset_2d_int8(dataset,dat,attribute)

   implicit none
   character(*),intent(in)          :: dataset
   integer*8,intent(in)             :: dat(:,:)
   character(*),intent(in),optional :: attribute ! optional explanation
   integer(hid_t)                   :: space_id
   integer(hid_t)                   :: dataset_id
   integer*4                        :: status
   integer*4,parameter              :: rank = 2
   integer*8                        :: n(2)
   
   call check_file_for_writing(dataset)
   
   n(1) = size(dat,1)
   n(2) = size(dat,2)
   call h5screate_simple_f(rank,n,space_id,status)  ! Create the space for the dataset
   call h5dcreate_f(file_id,dataset,H5T_STD_I64LE,space_id,dataset_id,status) ! create dataset
   call h5dwrite_f(dataset_id,H5T_STD_I64LE,dat,n,status) ! write data set
   call h5sclose_f(space_id,status) ! close data space
   if (present(attribute)) call write_attribute(attribute,dataset_id)
   call h5dclose_f(dataset_id,status) ! close data set
   
end subroutine write_dataset_2d_int8

subroutine write_dataset_2d_real4(dataset,dat,attribute)

   implicit none
   character(*),intent(in)          :: dataset
   real*4,intent(in)                :: dat(:,:)
   character(*),intent(in),optional :: attribute ! optional explanation
   integer(hid_t)                   :: space_id
   integer(hid_t)                   :: dataset_id
   integer*4                        :: status
   integer*4,parameter              :: rank = 2
   integer*8                        :: n(2)
   
   call check_file_for_writing(dataset)
   
   n(1) = size(dat,1)
   n(2) = size(dat,2)
   call h5screate_simple_f(rank,n,space_id,status)  ! Create the space for the dataset
   call h5dcreate_f(file_id,dataset,H5T_IEEE_F32LE,space_id,dataset_id,status) ! create dataset
   call h5dwrite_f(dataset_id,H5T_IEEE_F32LE,dat,n,status) ! write data set
   call h5sclose_f(space_id,status) ! close data space
   if (present(attribute)) call write_attribute(attribute,dataset_id)
   call h5dclose_f(dataset_id,status) ! close data set
   
end subroutine write_dataset_2d_real4

subroutine write_dataset_2d_real8(dataset,dat,attribute)

   implicit none
   character(*),intent(in)          :: dataset
   real*8,intent(in)                :: dat(:,:)
   character(*),intent(in),optional :: attribute ! optional explanation
   integer(hid_t)                   :: space_id
   integer(hid_t)                   :: dataset_id
   integer*4                        :: status
   integer*4,parameter              :: rank = 2
   integer*8                        :: n(2)
   
   call check_file_for_writing(dataset)
   
   n(1) = size(dat,1)
   n(2) = size(dat,2)
   call h5screate_simple_f(rank,n,space_id,status)  ! Create the space for the dataset
   call h5dcreate_f(file_id,dataset,H5T_IEEE_F64LE,space_id,dataset_id,status) ! create dataset
   call h5dwrite_f(dataset_id,H5T_IEEE_F64LE,dat,n,status) ! write data set
   call h5sclose_f(space_id,status) ! close data space
   if (present(attribute)) call write_attribute(attribute,dataset_id)
   call h5dclose_f(dataset_id,status) ! close data set
   
end subroutine write_dataset_2d_real8

subroutine write_attribute(attribute,dataset_id)

   implicit none
   character(*),intent(in)       :: attribute
   integer(hid_t),intent(in)     :: dataset_id
   integer(hid_t)                :: attr_id ! Attribute identifier
   integer(hid_t)                :: aspace_id ! Attribute Dataspace identifier
   integer(hid_t)                :: atype_id ! Attribute Dataspace identifier
   integer*8,dimension(1)        :: adims = (/1/) ! Attribute dimension
   integer*4,parameter           :: arank = 1 ! Attribure rank
   integer(size_t)               :: attrlen ! Length of the attribute string
   character(*),parameter        :: aname = "Comment" ! Attribute name
   integer(hsize_t),dimension(1) :: data_dims = (/1/)
   integer*4                     :: status
   
   attrlen = len(attribute)
   call h5screate_simple_f(arank,adims,aspace_id,status) ! Create scalar data space for the attribute.
   call h5tcopy_f(H5T_NATIVE_CHARACTER,atype_id,status) ! Create datatype for the attribute.
   call h5tset_size_f(atype_id,attrlen,status)
   call h5acreate_f(dataset_id,aname,atype_id,aspace_id,attr_id,status)
   call h5awrite_f(attr_id,atype_id,attribute,data_dims,status) ! Write the attribute data.
   call h5aclose_f(attr_id,status) ! Close the attribute.
   
end subroutine write_attribute

subroutine check_file_for_writing(dataset)

   character(*),intent(in) :: dataset
   integer*4               :: i
   character(255)          :: group
   
   if (file_id==0) call error('no HDF5 file open to read dataset '//trim(dataset))
   
   do i = len(trim(dataset)),2,-1
      if (dataset(i:i)=='/') then
         group = dataset(1:i)
         if (.not.hdf5_object_exists(group)) then
            call error('could not find HDF5 group "'//trim(group)//'" to add dataset "'//trim(dataset)//'"')
         end if
         exit
      end if
   end do
   
   if (hdf5_object_exists(dataset)) call error('dataset "'//trim(dataset)//'" already exists')
   
end subroutine check_file_for_writing


! ===========================================================================================================
! subroutines used by interface hdf5_read_dataset
! ===========================================================================================================

subroutine read_dataset_0d_string(dataset, dat, convert)

   implicit none
   character(*),intent(in)       :: dataset
   character(255),intent(inout)  :: dat
   logical*4,intent(in),optional :: convert
   integer*4                     :: err,hdferr
   integer(hsize_t)              :: size(2)
   integer(hid_t)                :: dataset_id,filetype
   integer*8                     :: n(1)
   
   if (present(convert)) call error('HDF5: String variables cannot be converted in HDF5.')
   
   n(1) = 255
   
   call hdf5_check_dataset(dataset,0)
   call h5dopen_f(file_id, dataset, dataset_id, err) ! Open dataset
   call H5Tcopy_f(H5T_FORTRAN_S1, filetype, hdferr)
   call H5Tset_size_f(filetype, n(1), hdferr)
   call h5dread_f(dataset_id, filetype, dat, size, err)
   call h5dclose_f(dataset_id, err)

end subroutine read_dataset_0d_string

subroutine read_dataset_0d_int4(dataset, dat, convert)

   implicit none
   character(*),intent(in)          :: dataset
   integer*4,intent(inout)          :: dat
   logical*4,intent(in),optional    :: convert
   integer*4                        :: dat1(1)

   call hdf5_check_dataset(dataset,0)
   call read_numeric(dataset, i4=dat1, convert=convert)
   dat = dat1(1)
   
end subroutine read_dataset_0d_int4

subroutine read_dataset_0d_int8(dataset, dat, convert)

   implicit none
   character(*),intent(in)          :: dataset
   integer*8,intent(inout)          :: dat
   logical*4,intent(in),optional    :: convert
   integer*8                        :: dat1(1)

   call hdf5_check_dataset(dataset,0)
   call read_numeric(dataset, i8=dat1, convert=convert)
   dat = dat1(1)
   
end subroutine read_dataset_0d_int8

subroutine read_dataset_0d_real4(dataset, dat, convert)

   implicit none
   character(*),intent(in)          :: dataset
   real*4,intent(inout)             :: dat
   logical*4,intent(in),optional    :: convert
   real*4                           :: dat1(1)

   call hdf5_check_dataset(dataset,0)
   call read_numeric(dataset, r4=dat1, convert=convert)
   dat = dat1(1)
   
end subroutine read_dataset_0d_real4

subroutine read_dataset_0d_real8(dataset, dat, convert)

   implicit none
   character(*),intent(in)          :: dataset
   real*8,intent(inout)             :: dat
   logical*4,intent(in),optional    :: convert
   real*8                           :: dat1(1)

   call hdf5_check_dataset(dataset,0)
   call read_numeric(dataset, r8=dat1, convert=convert)
   dat = dat1(1)
   
end subroutine read_dataset_0d_real8

subroutine read_dataset_1d_int4(dataset, dat, convert)

   implicit none
   character(*),intent(in)          :: dataset
   integer*4,intent(inout)          :: dat(:)
   logical*4,intent(in),optional    :: convert
   integer*4                        :: dims(1)

   dims(1) = size(dat)
   call hdf5_check_dataset(dataset,1,dims)
   call read_numeric(dataset, i4=dat, convert=convert)
   
end subroutine read_dataset_1d_int4

subroutine read_dataset_1d_int8(dataset, dat, convert)

   implicit none
   character(*),intent(in)          :: dataset
   integer*8,intent(inout)          :: dat(:)
   logical*4,intent(in),optional    :: convert
   integer*4                        :: dims(1)

   dims(1) = size(dat)
   call hdf5_check_dataset(dataset,1,dims)
   call read_numeric(dataset, i8=dat, convert=convert)
   
end subroutine read_dataset_1d_int8

subroutine read_dataset_1d_real4(dataset, dat, convert)

   implicit none
   character(*),intent(in)          :: dataset
   real*4,intent(inout)             :: dat(:)
   logical*4,intent(in),optional    :: convert
   integer*4                        :: dims(1)
   
   dims(1) = size(dat)
   call hdf5_check_dataset(dataset,1,dims)
   call read_numeric(dataset, r4=dat, convert=convert)
   
end subroutine read_dataset_1d_real4

subroutine read_dataset_1d_real8(dataset, dat, convert)

   implicit none
   character(*),intent(in)          :: dataset
   real*8,intent(inout)             :: dat(:)
   logical*4,intent(in),optional    :: convert
   integer*4                        :: dims(1)
   
   dims(1) = size(dat)
   call hdf5_check_dataset(dataset,1,dims)
   call read_numeric(dataset, r8=dat, convert=convert)
   
end subroutine read_dataset_1d_real8

subroutine read_dataset_2d_int4(dataset, dat, convert)

   implicit none
   character(*),intent(in)          :: dataset
   integer*4,intent(inout)          :: dat(:,:)
   logical*4,intent(in),optional    :: convert
   integer*4,allocatable            :: dat1(:)
   integer*4                        :: dims(2)
   
   dims = (/size(dat,1),size(dat,2)/)
   call hdf5_check_dataset(dataset,2,dims)
   allocate(dat1(product(dims)))
   call read_numeric(dataset, i4=dat1, convert=convert)
   dat = reshape(dat1,dims)
   
end subroutine read_dataset_2d_int4

subroutine read_dataset_2d_int8(dataset, dat, convert)

   implicit none
   character(*),intent(in)          :: dataset
   integer*8,intent(inout)          :: dat(:,:)
   logical*4,intent(in),optional    :: convert
   integer*8,allocatable            :: dat1(:)
   integer*4                        :: dims(2)
   
   dims = (/size(dat,1),size(dat,2)/)
   call hdf5_check_dataset(dataset,2,dims)
   allocate(dat1(product(dims)))
   call read_numeric(dataset, i8=dat1, convert=convert)
   dat = reshape(dat1,dims)
   
end subroutine read_dataset_2d_int8

subroutine read_dataset_2d_real4(dataset, dat, convert)

   implicit none
   character(*),intent(in)          :: dataset
   real*4,intent(inout)             :: dat(:,:)
   logical*4,intent(in),optional    :: convert
   real*4,allocatable               :: dat1(:)
   integer*4                        :: dims(2)
   
   dims = (/size(dat,1),size(dat,2)/)
   call hdf5_check_dataset(dataset,2,dims)
   allocate(dat1(product(dims)))
   call read_numeric(dataset, r4=dat1, convert=convert)
   dat = reshape(dat1,dims)
   
end subroutine read_dataset_2d_real4

subroutine read_dataset_2d_real8(dataset, dat, convert)

   implicit none
   character(*),intent(in)          :: dataset
   real*8,intent(inout)             :: dat(:,:)
   logical*4,intent(in),optional    :: convert
   real*8,allocatable               :: dat1(:)
   integer*4                        :: dims(2)
   
   dims = (/size(dat,1),size(dat,2)/)
   call hdf5_check_dataset(dataset,2,dims)
   allocate(dat1(product(dims)))
   call read_numeric(dataset, r8=dat1, convert=convert)
   dat = reshape(dat1,dims)
   
end subroutine read_dataset_2d_real8

subroutine read_numeric(dataset, i4, i8, r4, r8, convert)

   implicit none
   character(*),intent(in)          :: dataset
   integer*4,intent(out),optional   :: i4(:)
   integer*8,intent(out),optional   :: i8(:)
   real*4,intent(out),optional      :: r4(:)
   real*8,intent(out),optional      :: r8(:)
   logical*4,intent(in),optional    :: convert
   integer*8                        :: n
   character(2)                     :: expected_type
   character(2)                     :: detected_type
   integer*4,allocatable            :: dati4(:)
   integer*8,allocatable            :: dati8(:)
   real*4,allocatable               :: datr4(:)
   real*8,allocatable               :: datr8(:)
   integer*4                        :: k
   integer*4                        :: err
   integer(hsize_t)                 :: hsize(2)
   integer(hid_t)                   :: dataset_id
   integer(hid_t)                   :: hdf_type
   character(255)                   :: msg
   
   ! determine expected data type and size
   k = 0
   if (present(i4)) then
      expected_type = 'i4'
      n = size(i4,1,4)
      k = k+1
   end if
   if (present(i8)) then
      expected_type = 'i8'
      n = size(i8,1,4)
      k = k+1
   end if
   if (present(r4)) then
      expected_type = 'r4'
      n = size(r4,1,4)
      k = k+1
   end if
   if (present(r8)) then
      expected_type = 'r8'
      n = size(r8,1,4)
      k = k+1
   end if
   if (k.ne.1) call deverror('HDF5: read_numeric must have exactly one expected data type')
   
   ! open dataset with checks
   call h5dopen_f(file_id, dataset, dataset_id, err) ! Open dataset
   
   ! determine detected data type
   call get_mem_type_id(dataset_id, hdf_type, detected_type)
   
   ! read data
   select case(detected_type)
   case('i4')
      allocate(dati4(n))
      call h5dread_f(dataset_id, hdf_type, dati4, hsize, err)
   case('i8')
      allocate(dati8(n))
      call h5dread_f(dataset_id, hdf_type, dati8, hsize, err) 
   case('r4')
      allocate(datr4(n))
      call h5dread_f(dataset_id, hdf_type, datr4, hsize, err)
   case('r8')
      allocate(datr8(n))
      call h5dread_f(dataset_id, hdf_type, datr8, hsize, err)
   case default
      call deverror('HDF5: unknown detected type: ',expected_type)
   end select
   if (err.ne.0) call error('HDF5: cannot read dataset '//trim(dataset))
   
   if (detected_type==expected_type) then
   
      select case(detected_type)
      case('i4')
         i4 = dati4
      case('i8')
         i8 = dati8
      case('r4')
         r4 = datr4
      case('r8')
         r8 = datr8
      end select
   
   else
   
      ! handle differences between expected and detected type
      write(msg,'(6A)') 'Expected type ',expected_type,' but detected type ', &
      & detected_type,' in dataset ',dataset
      if (present(convert)) then
         if (.not.convert) call error(trim(msg))
      else 
         call error(trim(msg))
      end if
      
      select case(detected_type)
      case('i4')
         if (expected_type=='i8') i8 = int(dati4,8)
         if (expected_type=='r4') r4 = real(dati4,4)
         if (expected_type=='r8') r8 = real(dati4,8)
      case('i8')
         if (expected_type=='i4') i4 = int(dati8,4)
         if (expected_type=='r4') r4 = real(dati8,4)
         if (expected_type=='r8') r8 = real(dati8,8)
      case('r4')
         if (expected_type=='i4') i4 = int(datr4,4)
         if (expected_type=='i8') i8 = int(datr4,8)
         if (expected_type=='r8') r8 = real(datr4,8)
      case('r8')
         if (expected_type=='i4') i4 = int(datr8,4)
         if (expected_type=='i8') i8 = int(datr8,8)
         if (expected_type=='r4') r4 = real(datr8,4)
      end select
      
   end if
   
   call h5dclose_f(dataset_id, err)
   
end subroutine read_numeric

subroutine get_mem_type_id(dataset_id, hdf_type, fortran_type)

   implicit none

   integer(hid_t),intent(in)        :: dataset_id
   integer(hid_t),intent(out)       :: hdf_type
   character(2),intent(out)         :: fortran_type
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
            hdf_type = H5T_STD_U32LE
         else if (sign == 1) then
            hdf_type = H5T_STD_I32LE
         else
            call error('HDF5: unknown data type')
         end if
      else if (size == 8) then
         if (sign == 0) then
            hdf_type = H5T_STD_U64LE
         else if (sign == 1) then
            hdf_type = H5T_STD_I64LE
         else
            call error('HDF5: unknown data type')
         end if
      else
         call error('HDF5: unknown data type')
      end if
   else if (class == 1) then ! floating point
      if (size == 4) then
         hdf_type = H5T_IEEE_F32LE
      else if (size == 8) then
         hdf_type = H5T_IEEE_F64LE
      else
         call error('HDF5: unknown data type')
      end if
   else
      call error('HDF5: unknown data type')
   end if
   
   if (class == 0) then
      write(fortran_type,'(A,I0)') 'i',size
   else
      write(fortran_type,'(A,I0)') 'r',size
   end if

end subroutine get_mem_type_id

subroutine hdf5_check_dataset(dataset,expected_rank,expected_dimensions)

   implicit none
   character(*),intent(in)          :: dataset
   integer*4,intent(in)             :: expected_rank
   integer*4,intent(in),optional    :: expected_dimensions(:)
   integer*4                        :: detected_rank
   integer*4,allocatable            :: detected_dimensions(:)
   integer*4                        :: i
   
   call hdf5_get_dataset_size(dataset,rank=detected_rank,dimensions=detected_dimensions)
   
   if (detected_rank.ne.expected_rank) then
      if (.not.((expected_rank==0).and.(detected_dimensions(1)==1))) then
         call hdf5_close()
         call error('Expected rank '//val2str(expected_rank)//' but detected rank '//val2str(detected_rank)//&
         &' in dataset '//trim(dataset))
      end if
   end if
   
   if ((expected_rank>0).and.(present(expected_dimensions))) then
      if (size(expected_dimensions).ne.expected_rank) call deverror('HDF5: size(expected_dimensions) differs from expected_rank')
      do i = 1,expected_rank
         if (expected_dimensions(i).ne.detected_dimensions(i)) then
            call hdf5_close()
            call error('HDF5: Expected dimensions differ from detected dimensions in dataset '//trim(dataset))
         end if
      end do
   end if
   
end subroutine hdf5_check_dataset
   
end module shared_module_hdf5