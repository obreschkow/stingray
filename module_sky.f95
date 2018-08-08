! Module converts intrinsic into apparent properties

module module_sky

   use module_constants
   use module_system
   use module_types
   use module_io
   use module_cosmology
   use module_conversion
   use module_user
   use module_parameters
   use module_tiling
   
   private
   public   :: make_sky
   
contains

subroutine make_sky

   implicit none
   character(len=255)      :: filename
   integer*8               :: n,i,n100
   integer*8               :: n_galaxies
   integer*8               :: n_groups
   type(type_base)         :: base
   type(type_sam)          :: sam   ! intrinsic galaxy properties from SAM
   character(len=3)        :: str
   type(type_sky_galaxy)   :: sky_galaxy
   type(type_sky_group)    :: sky_group
   integer*8               :: groupid,groupid_old
   type(type_sam)          :: sam_group_center
   type(type_base)         :: base_group_center
   integer*4               :: group_nsel
   real*4                  :: empty_real4
   integer*4               :: empty_int4
   
   ! write user info
   call tic
   call out('CONVERT INTRINSIC TO APPARENT PROPERTIES IN MOCK SKY')
   
   ! load previous steps
   call load_parameters
   call load_box_list
   
   ! set random seed (in case the assignment of sky-properties involves randomness)
   call set_seed(para%seed)
   
   ! open intrinsic sky
   filename = trim(para%path_output)//'mocksky_intrinsic.bin'
   open(1,file=trim(filename),action='read',form='unformatted',status='old',access='stream')
   read(1) n,empty_real4,empty_int4 ! number of galaxies
   
   ! initialize apparent sky (galaxies)
   write(filename,'(A,A,A,A)') trim(para%path_output),'mocksky_galaxies.bin'
   open(2,file=trim(filename),action='write',form="unformatted",status='replace',access='stream')
   write(2) 0_8 ! place holder for number of objects in mock sky
   
   ! initialize apparent sky (groups)
   write(filename,'(A,A,A,A)') trim(para%path_output),'mocksky_groups.bin'
   open(3,file=trim(filename),action='write',form="unformatted",status='replace',access='stream')
   write(3) 0_8 ! place holder for number of objects in mock sky
   
   ! convert galaxy properties and write master sky
   n100 = n/100
   groupid_old = -1
   group_nsel = 0
   n_galaxies = 0
   n_groups = 0
   do i = 1,n
   
      if (modulo(i-1,n100)==0) then
         write(str,'(I3)') nint(real(i-1,8)/n*100)
         call out('Progress: '//str//'%')
      end if
      
      read(1) base,sam
      
      Rvector = tile(base%tile)%Rvector
      Rpseudo = tile(base%tile)%Rpseudo
      call sam%rotate_vectors
      
      groupid = sam%get_groupid()
      if (groupid.ne.groupid_old) then
         
         ! make last group
         if ((groupid_old>=0).and.(group_nsel>0).and.(base_group_center%group_ntot>1)) then
            n_groups = n_groups+1
            call sky_group%make_from_sam(sam_group_center,base_group_center,groupid=n_groups,group_nselected=group_nsel)
            call sky_group%write_to_file(3)
         end if
         
         ! define new group
         sam_group_center = sam
         base_group_center = base
         group_nsel = 0
         groupid_old = groupid
         
      end if
      
      if (base%sam_selected) then
         if (base%group_ntot==1) then
            call sky_galaxy%make_from_sam(sam,base,galaxyid=n_galaxies,groupid=-1_8)
         else
            call sky_galaxy%make_from_sam(sam,base,galaxyid=n_galaxies,groupid=n_groups+1)
         end if
         if (sky_galaxy%is_selected(sam)) then
            n_galaxies = n_galaxies+1
            group_nsel = group_nsel+1
            call sky_galaxy%write_to_file(2)
         end if
      end if
      
   end do
   call out('Progress: 100%')
   
   ! add number of objects to beginning of file & close files
   write(2,pos=1) n_galaxies
   write(3,pos=1) n_groups
   close(1)
   close(2)
   close(3)
   
   ! check number of objects
   if (n_galaxies==0) call error('No galaxies in the apparent sky. Consider changing selection function.')
   !if (n_groups==0) call error('No groups in the apparent sky. Consider changing selection function.')
   
   ! user output
   call out('Number of galaxies in intrinsic sky:',n)
   call out('Number of galaxies in apparent sky:',n_galaxies)
   call out('Number of groups in apparent sky:',n_groups)
   call toc
   
end subroutine make_sky

end module module_sky