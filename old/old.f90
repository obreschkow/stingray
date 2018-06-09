! write info (ascii)
   inquire(file=filename, size=bytes)
   filename = trim(para%path_output)//'mocksurvey_info.txt'
   open(1,file=trim(filename),action='write',form="formatted",status='replace')
   write(1,'(A,I10)') 'Number.of.galaxies.in.apparent.sky        ',m
   write(1,'(A,I10)') 'Number.of.bytes.per.galaxy.in.apparent.sky',bytes/m
   close(1)
   
   ! write info (binary)
   inquire(file=filename, size=bytes)
   filename = trim(para%path_output)//'mocksurvey_info.bin'
   open(1,file=trim(filename),action='write',form="unformatted",status='replace')
   write(1) m
   write(1) bytes/m
   close(1)
   
   
   ! write info
   inquire(file=filename_sky_intrinsic, size=bytes)
   open(1,file=trim(para%path_output)//'mocksurvey_intrinsic_info.txt',action='write',form="formatted",status='replace')
   write(1,'(A,I10)') 'Number.of.galaxies.in.intrinsic.sky        ',nmockgalaxies
   write(1,'(A,I10)') 'Number.of.bytes.per.galaxy.in.intrinsic.sky',bytes/nmockgalaxies
   close(1)
   
   function extract_base(sam) result(base)
   
   implicit none
   type(type_sam),intent(in) :: sam
   type(type_base)           :: base
   
   base%groupid   = sam%galaxy%id_halo     
   base%xsam      = sam%galaxy%position    

end function extract_base

subroutine make_sky_coordinates(x,dc,ra,dec)

   implicit none
   real*4,intent(in)                :: x(3)  ! [box side-length] position vector
   real*4,intent(out),optional      :: dc    ! [length units of simulation] comoving distance
   real*4,intent(out),optional      :: ra    ! [rad] right ascension
   real*4,intent(out),optional      :: dec   ! [rad] declination
   real*4                           :: normx
   
   normx = norm(x)
   if (normx<=epsilon(normx)) call error('make_sky_coordinates: norm of x is zero.')
   
   if (present(dc))  dc  = normx*para%L
   if (present(ra))  ra  = modulo(atan2(x(1),x(3)),2*pi)
   if (present(dec)) dec = asin(min(1.0,x(2)/normx))
   
end subroutine make_sky_coordinates

logical function is_box_in_survey(ix,user)

   implicit none
   integer*4,intent(in) :: ix(3) ! [box side-length] box center in tiling coordinates
   logical,intent(in)   :: user
   
   ! check 8 corners
   if (is_intersecting_0D(real(ix,4),user)) then
      is_box_in_survey = .true.
      return
   end if
   
   ! check 12 edges
   if (is_intersecting_1D(real(ix,4),(/0.0,0.5,0.5/),0,user)) then
      is_box_in_survey = .true.
      return
   end if
   
   ! check 6 faces
   if (is_intersecting_2D(real(ix,4),(/0.0,0.0,0.5/),0,user)) then
      is_box_in_survey = .true.
      return
   end if
   
   ! check cube
   if (user) then
      if (is_intersecting_3D(real(ix,4),(/0.0,0.0,0.0/),0,user)) then
         is_box_in_survey = .true.
         return
      end if
   end if
   
   is_box_in_survey = .false.

end function is_box_in_survey

function is_intersecting_0D(x,user) result(selected)

   real*4,intent(in)       :: x(3)        ! [box side-length] central position of box
   logical,intent(in)      :: user
   logical                 :: selected
   
   if (is_point_in_survey(x+(/+0.5,+0.5,+0.5/),user)) then; selected = .true.; return; end if
   if (is_point_in_survey(x+(/+0.5,+0.5,-0.5/),user)) then; selected = .true.; return; end if
   if (is_point_in_survey(x+(/+0.5,-0.5,+0.5/),user)) then; selected = .true.; return; end if
   if (is_point_in_survey(x+(/+0.5,-0.5,-0.5/),user)) then; selected = .true.; return; end if
   if (is_point_in_survey(x-(/+0.5,+0.5,+0.5/),user)) then; selected = .true.; return; end if
   if (is_point_in_survey(x-(/+0.5,+0.5,-0.5/),user)) then; selected = .true.; return; end if
   if (is_point_in_survey(x-(/+0.5,-0.5,+0.5/),user)) then; selected = .true.; return; end if
   if (is_point_in_survey(x-(/+0.5,-0.5,-0.5/),user)) then; selected = .true.; return; end if
   
   selected = .false.

end function is_intersecting_0D

recursive function is_intersecting_1D(x,offset,level,user) result(selected)

   real*4,intent(in)       :: x(3)        ! [box side-length] central position of box
   real*4,intent(in)       :: offset(3)   ! [box side-length] position
   integer*4,intent(in)    :: level       ! level of recursion, starting from 0
   logical,intent(in)      :: user
   logical                 :: selected
   real*4                  :: dx          ! [box side-length] spacing probed at this recursion level
   integer*4               :: i
   
   if (is_point_in_survey(x+(/+offset(1),+offset(2),+offset(3)/),user)) then; selected = .true.; return; end if
   if (is_point_in_survey(x+(/+offset(1),-offset(2),+offset(3)/),user)) then; selected = .true.; return; end if
   if (is_point_in_survey(x+(/-offset(1),+offset(2),+offset(3)/),user)) then; selected = .true.; return; end if
   if (is_point_in_survey(x+(/-offset(1),-offset(2),+offset(3)/),user)) then; selected = .true.; return; end if
   if (is_point_in_survey(x+(/+offset(3),+offset(1),+offset(2)/),user)) then; selected = .true.; return; end if
   if (is_point_in_survey(x+(/+offset(3),+offset(1),-offset(2)/),user)) then; selected = .true.; return; end if
   if (is_point_in_survey(x+(/+offset(3),-offset(1),+offset(2)/),user)) then; selected = .true.; return; end if
   if (is_point_in_survey(x+(/+offset(3),-offset(1),-offset(2)/),user)) then; selected = .true.; return; end if
   if (is_point_in_survey(x+(/+offset(2),+offset(3),+offset(1)/),user)) then; selected = .true.; return; end if
   if (is_point_in_survey(x+(/-offset(2),+offset(3),+offset(1)/),user)) then; selected = .true.; return; end if
   if (is_point_in_survey(x+(/+offset(2),+offset(3),-offset(1)/),user)) then; selected = .true.; return; end if
   if (is_point_in_survey(x+(/-offset(2),+offset(3),-offset(1)/),user)) then; selected = .true.; return; end if
      
   if (level<maxlevel_1D) then
      dx = 0.25/real(2**level,4)
      do i = -1,1,2
         if (is_intersecting_1D(x,offset+dx*(/i,0,0/),level+1,user)) then
            selected = .true. 
            return
         end if
      end do
   end if
   
   selected = .false.

end function is_intersecting_1D

recursive function is_intersecting_2D(x,offset,level,user) result(selected)

   real*4,intent(in)       :: x(3)        ! [box side-length] central position of box
   real*4,intent(in)       :: offset(3)   ! [box side-length] position
   integer*4,intent(in)    :: level       ! level of recursion, starting from 0
   logical,intent(in)      :: user
   logical                 :: selected
   real*4                  :: dx          ! [box side-length] spacing probed at this recursion level
   integer*4               :: i,j
   
   if (is_point_in_survey(x+(/offset(1),offset(2),offset(3)/),user)) then; selected = .true.; return; end if
   if (is_point_in_survey(x+(/offset(2),offset(3),offset(1)/),user)) then; selected = .true.; return; end if
   if (is_point_in_survey(x+(/offset(3),offset(1),offset(2)/),user)) then; selected = .true.; return; end if
   if (is_point_in_survey(x-(/offset(1),offset(2),offset(3)/),user)) then; selected = .true.; return; end if
   if (is_point_in_survey(x-(/offset(2),offset(3),offset(1)/),user)) then; selected = .true.; return; end if
   if (is_point_in_survey(x-(/offset(3),offset(1),offset(2)/),user)) then; selected = .true.; return; end if
      
   if (level<maxlevel_2D) then
      dx = 0.25/real(2**level,4)
      do i = -1,1,2
         do j = -1,1,2
            if (is_intersecting_2D(x,offset+dx*(/i,j,0/),level+1,user)) then
               selected = .true. 
               return
            end if
         end do
      end do
   end if
   
   selected = .false.

end function is_intersecting_2D

recursive function is_intersecting_3D(x,offset,level,user) result(selected)

   real*4,intent(in)       :: x(3)        ! [box side-length] central position of box
   real*4,intent(in)       :: offset(3)   ! [box side-length] position
   integer*4,intent(in)    :: level       ! level of recursion, starting from 0
   logical,intent(in)      :: user
   logical                 :: selected
   real*4                  :: dx          ! [box side-length] spacing probed at this recursion level
   integer*4               :: i,j,k
   
   if (is_point_in_survey(x+offset,user)) then
      selected = .true.
      return
   end if
      
   if (level<maxlevel_3D) then
      dx = 0.25/real(2**level,4)
      do i = -1,1,2
         do j = -1,1,2
            do k = -1,1,2
               if (is_intersecting_3D(x,offset+dx*(/i,j,k/),level+1,user)) then
                  selected = .true. 
                  return
               end if
            end do
         end do
      end do
   end if
   
   selected = .false.

end function is_intersecting_3D

logical function is_point_in_survey(x,user)

   ! Check if the point x lies inside the selected survey volume. If user = true, this is the survey volume
   ! specified by the user function position_selection, otherwise it is the survey volume specified in the
   ! parameter file.

   real*4,intent(in)    :: x(3)     ! [box side-length] position of point in tiling coordinates
   logical,intent(in)   :: user
   real*4               :: normx    ! [box side-length]
   real*4               :: xrot(3)  ! [box side-length] position of point in rotated coordinates
   real*4               :: ra,dec   ! [rad] sky coordinates
   real*4               :: dc       ! [simulation unit] comoving distance
   
   counter = counter+1
   
   normx = norm(x)
   if (normx<=epsilon(normx)) then
      ra = 0.0
      dec = 0.0
      dc = 0.0
   else
      xrot = matmul(para%sky_rotation,x)
      dc = normx*para%L
      ra = modulo(atan2(xrot(1),xrot(3)),2*pi)
      dec = asin(min(1.0,xrot(2)/normx))
   end if
   
   if (user) then
   
      is_point_in_survey = position_selection(ra/degree,dec/degree,dc)
      return
      
   else
   
      if ((dc<para%dc_min).or.(dc>para%dc_max)) then
         is_point_in_survey = .false.
         return
      end if
      
      if (normx<=epsilon(normx)) then
         is_point_in_survey = (para%dc_min<=0)
         return
      end if
   
      alpha = acos(min(1.0,sum(x*para%axis)/normx))
      is_point_in_survey = (alpha<=para%angle)
      return
      
   end if
   
end function is_point_in_survey

subroutine make_cone_basis(a1,a2)
   implicit none
   real*4,intent(out)   :: a1(3),a2(3) ! two orthonormal basis vectors, orthogonal to the cone axis
   if (para%axis(3)>0.5) then
      a1 = cross_product(para%axis,(/0.0,1.0,0.0/))
      a1 = a1/sqrt(sum(a1**2))
      a2 = cross_product(para%axis,a1)
      a2 = a2/sqrt(sum(a2**2))
   else
      a1 = cross_product(para%axis,(/0.0,0.0,1.0/))
      a1 = a1/sqrt(sum(a1**2))
      a2 = cross_product(para%axis,a1)
      a2 = a2/sqrt(sum(a2**2))
   end if
end subroutine make_cone_basis

   ! Names (file and HDF5 objects)
      CHARACTER(LEN=15), PARAMETER :: filename = "/Users/do/Desktop/mytestfile.hdf5" ! File name
      CHARACTER(LEN=9), PARAMETER :: dsetname1 = "mydataset" ! Dataset name
     CHARACTER(LEN=8), PARAMETER :: groupname = "subgroup" ! Sub-Group 1 name
     CHARACTER(LEN=9), PARAMETER :: groupname3 = "subgroup2" ! Sub-Group 3 name
     ! Dataset 2 name
     CHARACTER(LEN=24), PARAMETER :: dsetname2 = "subgroup/another_dataset"
     ! Dataset 3 name
     CHARACTER(LEN=23), PARAMETER :: dsetname3 = "subgroup2/dataset_three"
     
     ! Identifiers
     INTEGER(HID_T) :: file_id       ! File identifier
     INTEGER(HID_T) :: group_id      ! Group identifier
     INTEGER(HID_T) :: group3_id     ! Group 3 identifier
     INTEGER(HID_T) :: dset1_id      ! Dataset 1 identifier
     INTEGER(HID_T) :: dset2_id      ! Dataset 2 identifier
     INTEGER(HID_T) :: dset3_id      ! Dataset 3 identifier
     INTEGER(HID_T) :: dspace1_id    ! Dataspace 1 identifier
     INTEGER(HID_T) :: dspace2_id    ! Dataspace 2 identifier
     INTEGER(HID_T) :: dspace3_id    ! Dataspace 3 identifier
   
     ! Integer array
     INTEGER :: rank                 ! Dataset rank
     INTEGER(HSIZE_T), DIMENSION(1) :: dims1 = (/100/) ! Dataset dimensions
     INTEGER(HSIZE_T), DIMENSION(1) :: data_dims1
     INTEGER, DIMENSION(100) :: dset_data1   ! Data buffers

     ! FP array
     INTEGER(HSIZE_T), DIMENSION(1) :: dims2 = (/50/)
     INTEGER(HSIZE_T), DIMENSION(1) :: data_dims2
     REAL, DIMENSION(50) :: dset_data2
   
     ! Integer array for dataset_three
     INTEGER(HSIZE_T), DIMENSION(1) :: dims3 = (/10/) ! Dataset dimensions
     INTEGER(HSIZE_T), DIMENSION(1) :: data_dims3     ! Dataset rank
     INTEGER, DIMENSION(10) :: dset_data3
   
     ! Misc variables (e.g. loop counters)
     INTEGER :: error ! Error flag
     INTEGER :: i,j
 ! =====================================================================

     ! Initialize the dset_data array 
     data_dims1(1) = 100
     rank = 1
     DO i = 1, 100
         dset_data1(i) = i
     END DO
   
     ! Initialize Fortran interface
     CALL h5open_f(error)   
     ! Create a new file
     CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

     ! Create dataspace 1 (the dataset is next) "dspace_id" is returned
     CALL h5screate_simple_f(rank, dims1, dspace1_id, error)
     ! Create dataset 1 with default properties "dset_id" is returned
     CALL h5dcreate_f(file_id, dsetname1, H5T_NATIVE_INTEGER, dspace1_id, &
                      dset1_id, error)
     ! Write dataset 1
     CALL h5dwrite_f(dset1_id, H5T_NATIVE_INTEGER, dset_data1, data_dims1, &
                     error)
     ! Close access to dataset 1
     CALL h5dclose_f(dset1_id, error)
     ! Close access to data space 1
     CALL h5sclose_f(dspace1_id, error)
   
     ! Create a group in the HDF5 file
     CALL h5gcreate_f(file_id, groupname, group_id, error)
  function pos_exists(ix) result(exists)

   implicit none
   integer*4,intent(in) :: ix(3) ! integer box position with ix=(0,0,0) being the box centered on the observer
   logical              :: exists
   integer*4            :: i
   
   if (allocated(pos)) then
      exists = .false.
      do i = 1,size(pos)
         if (pos(i)%ix(1)==ix(1).and.pos(i)%ix(2)==ix(2).and.pos(i)%ix(3)==ix(3)) then
            exists = .true.
            return
         end if
      end do
   else
      exists = .false.
   end if

end function pos_exists

function box_intersects(ix) result(output)

   ! check if box intersects with the global survey volume specified in the parameter file, regardless of the
   ! additional selection stated in the user function position_selection

   implicit none
   integer,intent(in)   :: ix(3) ! integer box position with ix=(0,0,0) being the box centered on the observer
   logical              :: output
   real*4               :: dmin,dmax ! minimum and maximum distance in units of box side-lengths
   real*4               :: d ! distance from observer to cube centre to in units of box side-lengths
   real*4,parameter     :: h3 = sqrt(3.0)/2.0 ! half the space diagonal of a unit cube
   real*4               :: alpha,beta ! angles
   real*4               :: a(3)
   real*4               :: center(3),e1(3),e2(3) ! center position of and unit vectors on square
   integer*4            :: i,j1,j2
   integer*4,parameter  :: npoints = 50
   real*4               :: norm,ra,dec
   
   ! check distance range
   dmin = para%dc_min/para%L
   dmax = para%dc_max/para%L
   d = sqrt(real(sum(ix**2)))
   if ((d<dmin-h3).or.(d>dmax+h3)) then
      output = .false.
      return
   end if
   
   ! check if circumscribed sphere of the box overlaps with cone
   if (h3>=d) then
      output = .true.
      return
   end if
   alpha = acos(sum(para%axis*ix)/d) ! angle between cone axis and line-of-sight to the box center
   beta = asin(h3/d) ! apparent angle of circumscribed sphere of the box
   if (alpha-beta>para%angle*sqrt(2.0)) then ! factor sqrt(2) needed in case of a cone with square base
      output = .false.
      return
   end if
   
   ! check if any point on the cube's face lies inside the global survey volume
   do i = 1,6
   
      ! parametrize square face
      select case(i)
         case(1)
            center = (/real(ix(1))+0.5,real(ix(2)),real(ix(3))/)
            e1 = (/0,1,0/)
            e2 = (/0,0,1/)
         case(2)
            center = (/real(ix(1))-0.5,real(ix(2)),real(ix(3))/)
            e1 = (/0,1,0/)
            e2 = (/0,0,1/)
         case(3)
            center = (/real(ix(1)),real(ix(2))+0.5,real(ix(3))/)
            e1 = (/1,0,0/)
            e2 = (/0,0,1/)
         case(4)
            center = (/real(ix(1)),real(ix(2))-0.5,real(ix(3))/)
            e1 = (/1,0,0/)
            e2 = (/0,0,1/)
         case(5)
            center = (/real(ix(1)),real(ix(2)),real(ix(3))+0.5/)
            e1 = (/1,0,0/)
            e2 = (/0,1,0/)
         case(6)
            center = (/real(ix(1)),real(ix(2)),real(ix(3))-0.5/)
            e1 = (/1,0,0/)
            e2 = (/0,1,0/)
      end select
      
      if (para%dec_min>pi) then
      
         ! check if any point on face i lies inside the circular cone
         do j1 = -npoints,npoints
            do j2 = -npoints,npoints
               a = center+0.5*real(j1)/real(npoints)*e1+0.5*real(j2)/real(npoints)*e2
               norm = sqrt(sum(a**2))
               alpha = acos(min(1.0,sum(a*para%axis)/norm))
               if (alpha<=para%angle) then
                  if ((norm>=dmin).and.(norm<=dmax)) then
                     output = .true.
                     return
                  end if
               end if
            end do
         end do
         
      else
      
         ! check if any point on the face i lies inside rectangular cone
         do j1 = -npoints,npoints
            do j2 = -npoints,npoints
               a = center+0.5*real(j1)/real(npoints)*e1+0.5*real(j2)/real(npoints)*e2
               a = matmul(para%sky_rotation,a)
               norm = sqrt(sum(a**2))
               ra = modulo(atan2(a(1),a(3)),2*pi)
               if ((ra>=para%ra_min).and.(ra<=para%ra_max)) then
                  dec = asin(min(1.0,a(2)/norm))
                  if ((dec>=para%dec_min).and.(dec<=para%dec_max)) then
                     if ((norm>=dmin).and.(norm<=dmax)) then
                        output = .true.
                        return
                     end if
                  end if
               end if
            end do
         end do
      
      end if
   
   end do
    
   output = .false.
   return
   
end function box_intersects
   
   ! Group "Global"
   call hdf5_add_group('Global')
   call hdf5_write_data('Global/mstars',sum(cone%mstars),'Total stellar mass')
   call hdf5_write_data('Global/SolidAngle',solidangle,'[deg^2] solid angle')
   call hdf5_write_data('Global/RAmin',ra_min,'[deg] minimum right ascension')
   call hdf5_write_data('Global/RAmax',ra_max,'[deg] maximum right ascension')
   call hdf5_write_data('Global/DECmin',dec_min,'[deg] minimum declination')
   call hdf5_write_data('Global/DECmax',dec_max,'[deg] maximum declination')
  
     ! Close the group
     CALL h5gclose_f(group_id, error)

     ! Create dataspace 2 (the dataset is next)
     data_dims2(1) = 50
     DO i = 1, 50
         dset_data2(i) = 1.0
     END DO
     ! Create dataspace 2
     CALL h5screate_simple_f(rank, dims2, dspace2_id, error)
     ! Create dataset 2 with default properties
     CALL h5dcreate_f(file_id, dsetname2, H5T_NATIVE_REAL, dspace2_id, &
                      dset2_id, error)
     ! Write dataset 2
     CALL h5dwrite_f(dset2_id, H5T_NATIVE_REAL, dset_data2, data_dims2, &
                     error)
     ! Close access to dataset 2
     CALL h5dclose_f(dset2_id, error)
     ! Close access to data space 2
     CALL h5sclose_f(dspace2_id, error)
   
     ! Create a group in the HDF5 file
     CALL h5gcreate_f(file_id, groupname3, group3_id, error)
     ! Close the group
    CALL h5gclose_f(group3_id, error)
   
    ! Create dataspace 3
    data_dims3(1) = 10
    DO i = 1, 10
        dset_data3(i) = i + 3
    END DO
    ! Create dataspace 3
    CALL h5screate_simple_f(rank, dims3, dspace3_id, error)
    ! Create dataset 3 with default properties
    CALL h5dcreate_f(file_id, dsetname3, H5T_NATIVE_INTEGER,  &
                     dspace3_id, dset3_id, error)
    ! Write dataset 3
    CALL h5dwrite_f(dset3_id, H5T_NATIVE_INTEGER, dset_data3, data_dims3, &
                    error)
    ! Close access to dataset 3
    CALL h5dclose_f(dset3_id, error)
    ! Close access to data space 3
    CALL h5sclose_f(dspace3_id, error)
   
    ! Close the file
    CALL h5fclose_f(file_id, error)
    ! Close FORTRAN interface
    CALL h5close_f(error)


   ! output basic statistics
   write(txt,'(A,F0.4,A,F0.4)') 'Position range: ',minval((/minval(base(:)%xbox(1)),&
   & minval(base(:)%xbox(2)),minval(base(:)%xbox(3))/)),&
   & ' to ',maxval((/maxval(base(:)%xbox(1)),maxval(base(:)%xbox(2)),maxval(base(:)%xbox(3))/))
   call out(txt)
   write(txt,'(A,I0,A,I0)')     'Identifier range: ',minval(base%id),' to ',maxval(base%id)
   call out(txt)

function box_intersects_simple(ix) result(output)

   implicit none
   integer,intent(in)   :: ix(3) ! integer box position with ix=(0,0,0) being the box centered on the observer
   logical              :: output
   real*4               :: dmin,dmax ! minimum and maximum distance in units of box side-lengths
   real*4               :: d ! distance from observer to cube centre to in units of box side-lengths
   real*4,parameter     :: h3 = sqrt(3.0)/2.0 ! half the space diagonal of a unit cube
   real*4               :: alpha,beta ! angles
   
   ! check if this is the central cube
   if ((ix(1)==0).and.(ix(2)==0).and.(ix(3)==0)) then
      output = .true.
      return
   end if
   
   ! check distance range
   dmin = para%dcmin/para%L
   dmax = para%dcmax/para%L
   d = sqrt(real(sum(ix**2)))
   if ((d<dmin-h3).or.(d>dmax+h3)) then
      output = .false.
      return
   end if
   
   ! check if circumscribed sphere of the box overlaps with cone
   alpha = acos(sum(para%axis*ix)/d) ! angle between cone axis and line-of-sight to the box center
   beta = asin(h3/d) ! apparent angle of circumscribed sphere of the box
   output = alpha-beta<para%angle
   return
   
end function box_intersects_simple