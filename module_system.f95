module module_system

   use module_types
   use module_constants

   public
   
   ! screen output
   logical                 :: opt_logfile
   logical                 :: opt_logscreen

   ! computation time evaluation
   integer*8               :: tstart_total
   integer*8               :: tstart
   integer*8               :: trate

contains

subroutine set_seed(seed)
   integer,intent(in)      :: seed
   integer                 :: rnd_size
   integer*4,allocatable   :: seed_array(:)
   call random_seed(size=rnd_size)
   allocate(seed_array(rnd_size))
   seed_array = seed
   call random_seed(put=seed_array)
end subroutine set_seed

function get_normal_random_number(mu,sigma)
   implicit none
   real,intent(in) :: mu,sigma
   real,parameter  :: pi = 3.14159265359
   real :: get_normal_random_number
   real :: z,u1,u2
   call random_number(u1)
   call random_number(u2)
   z = sqrt(-2.0*log(1.0-u1))*cos(2*pi*u2)
   get_normal_random_number = mu+z*sigma
end function get_normal_random_number

subroutine out_open(version)
   implicit none
   character(*),intent(in) :: version
   if (opt_logfile) then
      open(999,file=trim(para%path_output)//fn_log,action='write',status='replace',form='formatted')
      close(999)
   end if
   call hline
   call out('RUNNING stingray '//version)
   call tic_total
end subroutine out_open

subroutine out_close
   implicit none
   call toc_total
end subroutine out_close

subroutine out(txt,i)
   implicit none
   character(*),intent(in)       :: txt
   integer*8,intent(in),optional :: i
   if (opt_logfile) then
      open(999,file=trim(para%path_output)//fn_log,action='write',status='old',position='append',form='formatted')
      if (present(i)) then
         write(999,'(A,I0)') trim(txt)//' ',i
      else
         write(999,'(A)') trim(txt)
      end if
      close(999)
   end if
   if (opt_logscreen) then
      if (present(i)) then
         write(*,'(A,I0)') trim(txt)//' ',i
      else
         write(*,'(A)') trim(txt)
      end if
   end if
end subroutine out

subroutine error(txt)
   implicit none
   character(*),intent(in) :: txt
   call out('ERROR: '//txt)
   stop
end subroutine error

subroutine warning(txt)
   implicit none
   character(*),intent(in) :: txt
   call out('WARNING: '//txt)
end subroutine warning

subroutine hline
   implicit none
   call out('----------------------------------------------------------------------------------')
end subroutine hline

subroutine tic
   implicit none
   call system_clock(tstart,trate)
end subroutine tic  

subroutine toc
   implicit none
   integer*8 :: tstop
   call system_clock(tstop)
   call out('Time taken: '//trim(time_string(real(tstop-tstart,8)/trate)))
   call hline
end subroutine toc

subroutine tic_total
   implicit none
   call system_clock(tstart_total,trate)
   call hline
end subroutine tic_total

subroutine toc_total
   implicit none
   integer*8 :: tstop
   call system_clock(tstop)
   call out('TOTAL TIME TAKEN: '//trim(time_string(real(tstop-tstart_total,8)/trate)))
   call hline
end subroutine toc_total

function time_string(secs) result(str)
   implicit none
   real*8,intent(in)          :: secs
   real*4                     :: seconds
   integer*4                  :: minutes,hours,days
   character(100)             :: str
   character(1)               :: zero
   minutes = int(secs/60,4)
   hours = int(secs/3600,4)
   days = int(secs/86400,4)
   seconds = real(secs-minutes*60,4)
   minutes = minutes-hours*60
   hours = hours-days*24
   if (seconds<1) then
      zero = '0'
   else
      zero = ''
   end if
   if (days>0) then
      write(str,'(I0,A,I0,A,I0,A,F0.2,A)') days,'d ',hours,'h ', minutes,'m '//trim(zero),seconds,'s'
   else if (hours>0) then
      write(str,'(I0,A,I0,A,F0.2,A)') hours,'h ',minutes,'m '//trim(zero),seconds,'s'
   else if (minutes>0) then
      write(str,'(I0,A,F0.2,A)') minutes,'m '//trim(zero),seconds,'s'
   else
      write(str,'(A,F0.2,A)') trim(zero),seconds,'s'
   end if
end function time_string

function nodot(strin) result(strout)
   ! removed dot from on the RHS of string
   implicit none
   character(*),intent(in) :: strin
   character(len=255)      :: strout
   integer*4               :: l
   l = len(trim(strin))
   if (strin(l:l)=='.') then
      strout = strin(1:l-1)
   else
      strout = strin
   end if
end function nodot

function noslash(strin) result(strout)
   ! removed slash from on the RHS of string
   implicit none
   character(*),intent(in) :: strin
   character(len=255)      :: strout
   integer*4               :: l
   l = len(trim(strin))
   if (strin(l:l)=='/') then
      strout = strin(1:l-1)
   else
      strout = strin
   end if
end function noslash

function exists(filename,do_not_stop) result(res)
   implicit none
   character(len=*),intent(in)   :: filename
   logical,intent(in),optional   :: do_not_stop
   logical                       :: res
   inquire(file=trim(filename), exist=res)
   if ((.not.res).and.(.not.present(do_not_stop))) then
      call out('ERROR: File does not exist: '//trim(filename))
      stop
   end if
end function exists

subroutine check_exists(filename)
   implicit none
   character(*),intent(in) :: filename
   if (.not.exists(filename)) stop
end subroutine check_exists

subroutine nil(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20)
   class(*),optional :: x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20
   if (.false.) then
      select type(x1);  end select;    select type(x2);  end select
      select type(x3);  end select;    select type(x4);  end select
      select type(x5);  end select;    select type(x6);  end select
      select type(x7);  end select;    select type(x8);  end select
      select type(x9);  end select;    select type(x10); end select
      select type(x11); end select;    select type(x12); end select
      select type(x13); end select;    select type(x14); end select
      select type(x15); end select;    select type(x16); end select
      select type(x17); end select;    select type(x18); end select
      select type(x19); end select;    select type(x20); end select
   end if
end subroutine nil

end module module_system