module module_sky_apparent

   use module_constants
   use module_types
   use module_system
   use module_cosmology
   use module_user
   use module_parameters
   use module_tiling
   
   private
   public   :: make_sky_apparent
   
contains

subroutine make_sky_apparent

   implicit none
   character(len=255)   :: filename
   integer*8            :: n,i,m
   integer*4            :: bytespergalaxy
   integer*8            :: bytes
   type(type_base)      :: base
   type(type_sam)       :: sam   ! intrinsic galaxy properties from SAM
   type(type_sky)       :: sky  ! apparent galaxy properties
   
   ! write user info
   call tic
   call out('CONVERT INTRINSIC SKY TO APPARENT sky')
   
   ! load previous steps
   call load_parameters
   call load_box_list
   
   ! determine number of bytes per galaxy in intrinsic sky
   filename = trim(para%path_output)//'.tmpsizeof'
   open(1,file=trim(filename),action='write',form="unformatted",status='replace',access='stream')
   write(1) base,sam
   close(1)
   inquire(file=trim(filename), size=bytespergalaxy)
   
   ! determine number of galaxies in intrinsic sky
   filename = trim(para%path_output)//'mocksurvey_intrinsic.bin'
   inquire(file=filename, size=bytes)
   if (modulo(bytes,bytespergalaxy).ne.0) then
      call error('Size of intrinsic sky file inconsistent with type_base and/or type_sam.')
   end if
   n = bytes/bytespergalaxy
   
   ! open intrinsic sky
   open(2,file=trim(filename),action='read',form='unformatted',status='old',access='stream')
   
   ! initialize master one
   filename = trim(para%path_output)//'mocksurvey.bin'
   open(1,file=trim(filename),action='write',form="unformatted",status='replace',access='stream')
   
   ! write user info
   call out('Number of galaxies in intrinsic sky:',n)
   
   ! convert galaxy properties and write master sky
   m = 0
   do i = 1,n
      read(2) base,sam
      Rvector = tile(base%tile)%Rvector
      Rpseudo = tile(base%tile)%Rpseudo
      call rotate_vectors(sam)
      sky = convert_properties(base,sam,m+1)
      if (apparent_selection(sky)) then
         m = m+1
         call write_galaxy(sky)
      end if
   end do
   
   ! close files
   close(1)
   close(2)
   
   if (m==0) call error('No galaxies in the apparent sky. Consider changing selection function.')
     
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
   
   ! user output
   call out('Number of galaxies in apparent sky:',m)
   call out('Number of bytes per galaxy in apparent sky:',bytes/m)
   call toc
   
end subroutine make_sky_apparent

end module module_sky_apparent