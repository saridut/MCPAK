program main

use m_precision
use m_strings
use m_logger
use m_globals, only: job_tag
use m_control_io, only: read_control, write_control
use m_setup, only: setup, finish
use m_mc_solver, only: mcs_run

implicit none

!*******************************************************************************

character(len=64) :: cla
character(len=:), allocatable :: key
character(len=:), allocatable :: val
character(len=:), allocatable :: fn_control
character(len=512) :: msg
real(rp):: t_start, t_end
integer :: ierr
integer :: icla
integer :: ncla ! number of command line arguments, without the command name

!Two command line arguments may be provided -- (i) fn_control=val and 
!(ii) job_tag=val

fn_control = 'control.txt'

ncla = command_argument_count()
do icla = 1, ncla
    call get_command_argument(icla, value=cla, status=ierr)
    if (ierr > 0) then
        write(*, '(a,1x,i0)') "read failure for command argument", icla
        stop
    else if (ierr == -1) then
        write(*, '(a,1x,i0)') "command argument truncated", icla
        stop
    end if

    call str_get_keyval(cla, key, val, '=')
    if (key == 'fn_control') then
        fn_control = val
    else if (key == 'job_tag') then
        job_tag = '.'//val
    end if
end do

call logger%init('mcpak.log'//job_tag, .true.)

call cpu_time(t_start)

call read_control(fn_control)

call setup()

call mcs_run(ierr)

call cpu_time(t_end)

write(msg,'(f0.6)') (t_end-t_start)

call logger%log_msg('Total execution time (s): ' // trim(adjustl(msg)))

call finish()

!*******************************************************************************

end program
