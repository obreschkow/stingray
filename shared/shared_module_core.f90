! **********************************************************************************************************************************
! Shared Fortran module with essential routines also used by other shared modules
! Developed by Danail Obreschkow
!
! This this is the only shared module, also used by other shared modules. It guarantees consistency between the other shared
! modules, e.g. in the way screen/logfile output is handled and variables are converted.
! **********************************************************************************************************************************

module shared_module_core

   private

   ! normally accessed routines/variables ******************************************************************************************
   
   ! initialise screen/logfile output
   public   :: set_logfile_name ! routine to set logfile_name programmatically
   public   :: set_verbose ! routine to set verbose programmatically
   public   :: set_version ! routine to set version programmatically
   public   :: set_help ! routine to set help programmatically
   public   :: set_copyright ! routine to set copyright programmatically
   
   ! write to screen/logfile
   public   :: start_output ! starts screen/log output, writes title
   public   :: stop_output ! stops screen/log output, writes total wall time
   public   :: out ! write genering text with optional numeric value
   public   :: hline ! write horizontal line
   public   :: error ! write error message and stop
   public   :: deverror ! write developer error message and stop
   public   :: warning ! write warning message, but continue program
   public   :: progress ! write progress in percent
   public   :: tic ! draw horizontal line and start timer (normally used at the beginning of a program section)
   public   :: toc ! write wall time since last call of tic and draw horizontal line (used at the end of a program section)

   ! string handling
   public   :: isempty ! logical function checking if a string is empty
   public   :: replace_text ! replace all occurrences of a pattern in a string
   public   :: lowercase ! turn a string into a lower case string
   public   :: remove_tabs ! removes tabs from a string
   public   :: tabs2spaces ! replaces all tabs by spaces in a character string
   public   :: timestamp ! character function returning the current date+time as a character string
   public   :: last_character ! function returning the last non-empty character of a string
   public   :: execname ! character function returning the name of the executable
   
   ! file handling
   public   :: exists ! logical function checks if file or path exists
   public   :: check_file ! subroutine checks if file or path exists, if not produce error (optionally checks permissions)
   public   :: make_path ! makes new directory and produces error if not allowed
   public   :: remove_path ! deletes directory (inclusive subdirectories) and produces error if not allowed
   public   :: dir ! make file path with separators; do not call by multiple threads simultaneously
   public   :: delete_file ! deletes a file, if the user has the rights to do so; other wise produce an error
   
   ! conversion functions
   public   :: val2str ! converts any numeric type into character string; do not call by multiple threads simultaneously
   public   :: log2int ! converts logical type into int*4 of value 0 or 1
   public   :: int2log ! converts an integer (0/1) into a logical*4, produces error message if not 0/1
   public   :: sec2time ! converts seconds into human-readable time string
   
   ! miscellaneous
   public   :: version ! read-only character string containing the version
   public   :: help ! read-only character string containing help text
   public   :: copyright ! read-only character string containing the copyright text
   public   :: nil ! subroutine with no effect, but 20 optional arguments; used to suppress compiler warnings for unused arguments
   public   :: separator ! single character '/' or '\' used to separate (sub)directories and files
   
   ! rarely accessed routines/variables ********************************************************************************************
   public   :: set_logfile_unit ! routine to set the logfile_unit programmatically (default, see below)
   public   :: start_logfile ! starts a logfile; only called once per program (e.g. via start_output)
   public   :: verbose ! read-only logical flag specifying if an output is displayed on screen; normally specified using -verbose
   public   :: logfile_name ! read-only string of the logfile name, set using set_logfile_name
   public   :: logfile_open ! read-only logical flag specifying if the logfile has been initialised
   public   :: logfile_unit ! read-only integer giving the I/O-unit used for the logfile, set using set_logfile_unit
   public   :: time_start ! wall time when start_output was called
   
   ! handling screen/logfile outputs
   logical*4,parameter           :: stop_at_warnings = .false.
   integer*4,protected           :: logfile_unit = 999 ! file unit for log-file
   logical*4,protected           :: verbose = .true. ! logical flag to produce screen output; set programmatically via set_verbose
   character(len=255),protected  :: logfile_name = '' ! logfile name; can be set programmatically using set_logfile
   logical*4,protected           :: logfile_open = .false.
   integer*8,protected           :: time_start = -1
   integer*8,protected           :: time_tic = -1
   character(*),parameter        :: hline_str = repeat('-',60)
   logical*4,protected           :: just_made_hline = .false.
   logical*4,protected           :: once_made_hline = .false.
   character(len=255),protected  :: version = '0.0'
   character(len=255),protected  :: help = 'Consult the README file for additional information.'
   character(len=255),protected  :: copyright = 'Developed by Danail Obreschkow (danail.obreschkow@icrar.org).'
   logical,protected             :: has_stopped = .false.
   
   ! constants
   character(1),parameter        :: separator = '/'
   
contains

! screen/logfile output ************************************************************************************************************

subroutine start_output(title)

   implicit none
   character(*),intent(in),optional :: title
   
   if (time_start.ne.-1) call deverror('start_output cannot be called more than once')
   
   if (.not.isempty(logfile_name)) call start_logfile
   
   call hline
   if (present(title)) then
      if (.not.isempty(title)) call out(title)
   else
      call out('Running '//execname()//' version '//trim(version)//'.')
   end if
   call system_clock(time_start)

end subroutine start_output

subroutine stop_output(delete_logfile)

   ! close output on screen/logfile and write total time taken

   implicit none
   logical*4,intent(in),optional :: delete_logfile ! if true and if a logfile exists, it will be deleted
   integer*8                     :: time_close
   integer*8                     :: time_rate
   
   if (time_start==-1) call deverror('stop_output cannot be called without first calling start_output')
   
   call system_clock(time_close,time_rate)
   
   call out('TOTAL WALL TIME: '//sec2time(real(time_close-time_start,8)/time_rate))
   call hline
   
   if (present(delete_logfile)) then
      if (delete_logfile) then
         if (logfile_open) then
            call delete_file(logfile_name)
         end if
      end if
   end if
   
   logfile_open = .false.

end subroutine stop_output

subroutine set_verbose(value)

   implicit none
   logical*4,intent(in) :: value
   verbose = value

end subroutine set_verbose

subroutine set_logfile_name(txt)

   implicit none
   character(*),intent(in) :: txt
   logfile_name = txt

end subroutine set_logfile_name

subroutine set_logfile_unit(unit)

   implicit none
   integer*4,intent(in) :: unit
   logfile_unit = unit

end subroutine set_logfile_unit

subroutine set_version(txt)

   implicit none
   character(*),intent(in) :: txt
   version = trim(txt)
   
end subroutine set_version

subroutine set_help(txt)

   implicit none
   character(*),intent(in) :: txt
   help = trim(txt)
   
end subroutine set_help

subroutine set_copyright(txt)

   implicit none
   character(*),intent(in) :: txt
   copyright = trim(txt)
   
end subroutine set_copyright

subroutine start_logfile

   implicit none
   integer*4   :: i
   
   ! basic checks
   if (logfile_open) call deverror('do not call start_logfile more than once')
   if (isempty(logfile_name)) call deverror('do not call start_logfile before defining a non-empty logfile_name')
   
   ! check if path exists and has read+write permission
   do i = len(trim(logfile_name)),2,-1
      if (logfile_name(i:i)==separator) then
         if (exists(logfile_name(1:i))) then
            call check_file(logfile_name(1:i),'rw')
         else
            call make_path(logfile_name(1:i))
            call check_file(logfile_name(1:i),'rw')
         end if
         exit
      end if
   end do
   
   ! open log file
   !$OMP CRITICAL (subroutine_start_logfile)
   open(logfile_unit,file=trim(logfile_name),action='write',status='replace',form='formatted')
   close(logfile_unit)
   logfile_open = .true.
   !$OMP END CRITICAL (subroutine_start_logfile)

end subroutine start_logfile

subroutine out(txt,value)

   ! write text to screen and/or logfile

   implicit none
   character(*),intent(in)       :: txt
   class(*),intent(in),optional  :: value    ! optional numeric value to be appended to txt
   character(len=255)            :: string
   
   !$OMP CRITICAL (subroutine_out)
   if (.not.has_stopped) then ! to avoid multiple error messages in multi-threading
   
      ! convert input arguments into a single string
      if (present(value)) then
         string = txt//val2str(value)
      else
         string = txt
      end if
      
      ! output on screen
      if (verbose) write(*,'(A)') trim(string)

      ! output to logfile      
      if (logfile_open) then
         open(logfile_unit,file=trim(logfile_name),action='write',status='old',position='append',form='formatted')
         write(logfile_unit,'(A)') trim(string)
         close(logfile_unit)
      end if
   
      just_made_hline = .false.
      
   end if
   !$OMP END CRITICAL (subroutine_out)
         
end subroutine out

subroutine hline

   ! outputs a horizontal line, unless such a line has just been written

   implicit none
   if (.not.just_made_hline) then
      call out(hline_str)
      just_made_hline = .true.
      once_made_hline = .true.
   end if
   
end subroutine hline

subroutine error(txt,value,showhelp)

   ! write error message and stop code

   implicit none
   character(*),intent(in)       :: txt
   class(*),intent(in),optional  :: value
   logical*4,intent(in),optional   :: showhelp
   
   call out('ERROR: '//txt,value)
   if (present(showhelp)) then
      if ((showhelp).and.(.not.isempty(help))) call out(help)
   end if
   if (once_made_hline) call hline
   has_stopped = .true.
   stop
   
end subroutine error

subroutine deverror(txt,value)

   ! write error message and stop code

   implicit none
   character(*),intent(in)       :: txt
   class(*),intent(in),optional  :: value
   
   call out('DEVELOPEPR ERROR: '//txt,value)
   if (once_made_hline) call hline
   has_stopped = .true.
   stop
   
end subroutine deverror

subroutine warning(txt,value)

   ! write warning; the code is stopped if stop_at_warnings==.true.

   implicit none
   character(*),intent(in)       :: txt
   class(*),intent(in),optional  :: value
   
   if (.not.once_made_hline) call hline
   call out('WARNING: '//trim(txt),value)
   if (stop_at_warnings) then
      call hline
      has_stopped = .true.
      stop
   end if
   
end subroutine warning

subroutine progress(fraction,fmt)

   ! write progress in percent

   implicit none
   real*4,intent(in)                :: fraction
   character(*),intent(in),optional :: fmt ! optional number format
   character(len=20)                :: fmt_
   character(len=30)                :: str
   
   fmt_ = 'F6.2'
   if (present(fmt)) fmt_=fmt
   write(str,'(A,'//fmt_//',A)') 'Progress: ',fraction*100,'%'
   call out(trim(str))

end subroutine progress


subroutine tic(title,value)

   ! start timer for program part, denoted with a horizontal line in the output
   
   implicit none
   character(*),optional,intent(in) :: title
   class(*),intent(in),optional     :: value
   
   call hline
   if (present(title)) call out(trim(title),value)
   call system_clock(time_tic)
   
end subroutine tic  

subroutine toc

   ! stop timer for program part, denoted with a horizontal line in the output

   implicit none
   integer*8 :: toc_time
   integer*8 :: time_rate
   
   if (time_tic==-1) call deverror('every call of toc must be preceded by a call of tic.')
   
   call system_clock(toc_time,time_rate)
   call out('Wall time: '//sec2time(real(toc_time-time_tic,8)/time_rate))
   call hline
   time_tic = -1
   
end subroutine toc


! **********************************************************************************************************************************
! conversion functions *************************************************************************************************************
! **********************************************************************************************************************************

function val2str(value) result(str)

   ! IMPORTANT: must be called without trim or adjust, e.g. use val2str(char), not val2str(trim(char))

   implicit none
   class(*),intent(in)        :: value ! value in any intrinsic numeric type, passed as unlimited polymorphic variable
   character(len=255)         :: txt
   character(:),allocatable   :: str

   select type (value)
   type is (integer(kind=4));    write(txt,'(i0)') value
   type is (integer(kind=8));    write(txt,'(i0)') value
   type is (integer(kind=16));   write(txt,'(i0)') value
   type is (real(kind=4));       write(txt,'(1pg0)') value
   type is (real(kind=8));       write(txt,'(1pg0)') value
   !type is (real(kind=16));      write(txt,'(1pg0)') value ! not supported by all compilers
   type is (logical);            write(txt,'(1l)') value
   type is (character(*));       write(txt,'(a)') value
   type is (complex(kind=4));    write(txt,'(1pg0,sp,1pg0,"i")') value
   type is (complex(kind=8));    write(txt,'(1pg0,sp,1pg0,"i")') value
   !type is (complex(kind=16));   write(txt,'(1pg0,sp,1pg0,"i")') value ! not supported by all compilers
   class default
      call deverror('unknown variable type')
   end select

   str = trim(txt)
   
end function val2str

pure elemental integer*4 function log2int(a)

   logical*4,intent(in) :: a
   
   if (a) then
      log2int = 1
   else
      log2int = 0
   end if
   
end function log2int

pure elemental logical*4 function int2log(a)

   integer*4,intent(in)  :: a
   if (a==0) then
      int2log = .false.
   else
      int2log = .true.
   end if
   
end function int2log

function sec2time(secs) result(strout)

   ! convertes seconds into human-readable text format

   implicit none
   real*8,intent(in)          :: secs
   real*4                     :: seconds
   integer*4                  :: minutes,hours,days
   character(100)             :: str
   character(1)               :: zero
   character(:),allocatable   :: strout
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
   
   strout = trim(str)
   
end function sec2time


! **********************************************************************************************************************************
! file handling ********************************************************************************************************************
! **********************************************************************************************************************************

logical*4 function exists(filename)
   
   implicit none
   character(len=*),intent(in)      :: filename ! path of file name
   
   if (isempty(filename)) call error('attempting to check existence of empty filename')
   inquire(file=trim(filename), exist=exists)
   
end function exists

subroutine check_file(filename,permission)

   implicit none
   character(*),intent(in)          :: filename ! path or file name
   character(*),intent(in),optional :: permission ! requested permission, e.g, 'r', 'w', 'rw', 'rwx'
   integer*4                        :: status
   character(4)                     :: kind
   
   ! decide if filename is a file or a path
   if (last_character(filename)==separator) then
      kind = 'path'
   else
      kind = 'file'
   end if
   
   ! check existence
   if (.not.exists(filename)) call error(kind//' does not exist: '//trim(filename))
   
   ! optional check of permissions
   if (present(permission)) then
      status = access(trim(filename),trim(permission))
      if (status.ne.0) call error('you do not have '//trim(permission)//' rights for '//kind//': '//trim(filename))
   end if
   
end subroutine check_file

subroutine make_path(path)

   implicit none
   character(*),intent(in) :: path
   integer*4               :: status
   
   status = system('mkdir -p '//trim(path))
   if (status.ne.0) call error('you do not have permission to create the directory: '//trim(path))
   
end subroutine make_path

subroutine remove_path(path)

   implicit none
   character(*),intent(in) :: path
   integer*4               :: status
   
   status = system('rm -rf '//trim(path))
   if (status.ne.0) call error('you do not have permission to remove the directory: '//trim(path))
   
end subroutine remove_path

function dir(s1,s2,s3,s4,s5,s6,s7,s8,s9,ispath) result(out)

   ! IMPORTANT: must be called without trim or adjust, e.g. use dir(str1,str2), not dir(trim(str1),str2)

   implicit none
   class(*),intent(in)           :: s1
   class(*),intent(in),optional  :: s2,s3,s4,s5,s6,s7,s8,s9
   logical,intent(in),optional   :: ispath
   character(255)                :: str
   character(:),allocatable      :: out
   
   ! concatenate pieces
   str = val2str(s1)
   if (present(s2)) str = trim(str)//separator//trim(val2str(s2))
   if (present(s3)) str = trim(str)//separator//trim(val2str(s3))
   if (present(s4)) str = trim(str)//separator//trim(val2str(s4))
   if (present(s5)) str = trim(str)//separator//trim(val2str(s5))
   if (present(s6)) str = trim(str)//separator//trim(val2str(s6))
   if (present(s7)) str = trim(str)//separator//trim(val2str(s7))
   if (present(s8)) str = trim(str)//separator//trim(val2str(s8))
   if (present(s9)) str = trim(str)//separator//trim(val2str(s9))
   
   ! terminate paths on a single [separator]
   if (present(ispath)) then
      if (ispath) str = trim(str)//separator
   end if
   
   ! avoid double separator
   str = replace_text(str,separator//separator,separator)
   
   ! produce output of correct size
   out = trim(str)
      
end function dir

subroutine delete_file(filename)

   implicit none
   character(*),intent(in) :: filename
   integer*4               :: status
   
   call check_file(filename,'rw')
   status = system('rm -r '//trim(filename))
   if (status.ne.0) call error('could not delete file ',trim(filename))
   
end subroutine delete_file


! **********************************************************************************************************************************
! character string handling ********************************************************************************************************
! **********************************************************************************************************************************

logical*4 function isempty(str)

   implicit none
   character(*),intent(in) :: str
   isempty = len(trim(str))==0
   
end function isempty

function replace_text (string,pattern,replacement)  result(out)
   
   implicit none
   character(*),intent(in)    :: string      ! input string
   character(*),intent(in)    :: pattern     ! pattern to be replaced
   character(*),intent(in)    :: replacement ! replacement text
   character(:),allocatable   :: out
   integer*4                  :: i,np,nr
   
   out = string
   np = len(pattern)
   nr = len(replacement)
   
   do
      i = index(out,pattern)
      if (i==0) exit
      out = out(:i-1)//replacement//out(i+np:)
   end do
   
end function replace_text

function lowercase(str_in) result(str_out)

   ! changes a strong to lower case

   implicit none
   character(*), intent(in) :: str_in
   character(len(str_in))      :: str_out
   integer :: ic, i

   character(26), parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
   character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

   str_out = str_in
   do i = 1, len_trim(str_in)
      ic = index(cap, str_in(i:i))
      if (ic > 0) str_out(i:i) = low(ic:ic)
   end do

end function lowercase

function remove_tabs(str) result(out)

   implicit none
   character(*),intent(in)    :: str      ! input string
   character(:),allocatable   :: out
   
   out = replace_text(str,achar(9),'')

end function remove_tabs

function tabs2spaces(str) result(out)

   implicit none
   character(*),intent(in)    :: str      ! input string
   character(:),allocatable   :: out
   
   out = replace_text(str,achar(9),' ')

end function tabs2spaces

character(1) function last_character(str)

   ! returns the right-most non empty character of a string
   implicit none
   character(*),intent(in)    :: str
   integer*4                  :: l
   l = len(trim(str))
   last_character = str(l:l)

end function

function timestamp() result(str)

   ! human readable time stamp
   implicit none
   integer*4         :: t(8)
   character(len=23) :: str
   call date_and_time(VALUES=t)
   write(str,'(I0,A,I0.2,A,I0.2,A,I0.2,A,I0.2,A,I0.2,A,I0.3)') t(1),'/',t(2),'/',t(3),'-',t(5),':',t(6),':',t(7),'.',t(8)
   
end function timestamp

function execname() result(str)

   implicit none
   character(:),allocatable   :: str
   character(255)             :: name
   integer*4                  :: i,j
   
   call get_command_argument(0,name)
   j = 0
   do i = len(trim(name)),1,-1
      if (name(i:i)==separator) then
         j = i+1
         exit
      end if
   end do
   
   if (j==0) then
      str = trim(name)
   else
      str = name(j:len(trim(name)))
   end if

end function execname


! **********************************************************************************************************************************
! miscellaneous ********************************************************************************************************************
! **********************************************************************************************************************************

subroutine nil(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16,x17,x18,x19,x20)
   ! call this routine with unused function arguments to avoid compiler warnings
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

end module shared_module_core