program traj2txt

use m_precision
use m_constants_math
use m_strings
use m_trajectory
use m_globals
use m_config_io

implicit none

character(len=256) :: cla_buf
character(len=:), allocatable :: fn_out
character(len=3) :: ft_out
integer :: traj_ver
integer :: ifrm_beg, ifrm_end, ifrm_stp
integer :: iframe, oframe, nofrms, iframe_step
integer :: fu_traj
integer :: i, ierr, istage

traj_ver = 2

call get_command_argument(1, cla_buf)
fn_cfg = trim(cla_buf)

call get_command_argument(2, cla_buf)
fn_traj = trim(cla_buf)

call get_command_argument(3, cla_buf)
ifrm_beg = str_to_i(cla_buf)

call get_command_argument(4, cla_buf)
ifrm_end = str_to_i(cla_buf)

call get_command_argument(5, cla_buf)
ifrm_stp = str_to_i(cla_buf)

call get_command_argument(6, cla_buf)
ft_out = trim(cla_buf)

!Read cfg file
call read_config(fn_cfg)

!Open file for reading trajectory
call traj%open(fn_traj, 'r')
write(*,*) 'numframes ', traj%num_frames

if (ifrm_beg == -1) ifrm_beg = traj%num_frames
if (ifrm_end == -1) ifrm_end = traj%num_frames

!Number of output frames
nofrms = (ifrm_end - (ifrm_beg-1))/ifrm_stp
if (nofrms > 1000) then
    write(*,*) 'More than 1000 frames will be written'
    write(*,'(a)', advance='no') 'Press "y" if ok: '
    read(*,*) cla_buf
    !Will fall through if not 'y' 
    if ( trim(cla_buf) /= 'y' ) ifrm_end = ifrm_beg - 1
end if

do iframe = ifrm_beg, ifrm_end, ifrm_stp
    oframe = (iframe - ifrm_beg + 1)/ifrm_stp
    write(*,*) oframe, iframe
    call traj%read(traj_ver, iframe, istage, nblks, nts, coordinates, ierr)

    if (ft_out == 'ldf') then
        fn_out = 'frame.txt.'//str_from_num(oframe)
        call write_ldf(fn_out, 'Frame-'//str_from_num(iframe))
    else if (ft_out == 'xyz') then
        fn_out = 'frame.xyz.'//str_from_num(oframe)
        call write_xyz(fn_out, 'Frame-'//str_from_num(iframe))
    else if (ft_out == 'cfg') then
        fn_out = 'frame.cfg.'//str_from_num(oframe)
        call write_config(fn_out, 'frame-'//str_from_num(iframe))
    else
        write(*,'(a,a)') 'unknown file type ', ft_out
        exit 
    end if

end do

!Close trajectory
call traj%close()

!*******************************************************************************

end program
