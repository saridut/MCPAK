program calc_rh

use m_precision
use m_constants_math
use m_strings
use m_trajectory
use m_vector
use m_globals
use m_config_io

implicit none

character(len=256) :: cla_buf
character(len=:), allocatable :: traj_dir
character(len=:), allocatable :: fn
integer :: num_traj, nts_beg
integer :: traj_ver, istage
integer :: i, j, itraj, iframe, ofrm, ierr, fu
type(ivector_t) :: counter
type(dvector_t) :: rh
real(rp) :: dkirk, tr, rnb, rh_ifrm
real(rp) :: rim2, rjm2
real(rp) :: rijm, irijm, irijm2
real(rp) :: ri_dot_rj
real(rp) :: C1, C2, consij
real(rp), dimension(3) :: ri, rj, rij

call get_command_argument(1, cla_buf)
fn_cfg = trim(cla_buf)

call get_command_argument(2, cla_buf)
fn_traj = trim(cla_buf)

call get_command_argument(3, cla_buf)
traj_dir = trim(cla_buf)

call get_command_argument(4, cla_buf)
num_traj = str_to_i(cla_buf)

call get_command_argument(5, cla_buf)
nts_beg = str_to_i(cla_buf)

traj_ver = 2 !Trajectory version

!Read cfg file
call read_config(fn_cfg)

call ivector_init(counter)
call dvector_init(rh)

rnb = 1.0_rp/num_atoms

do itraj = 1, num_traj
    !Open file for reading trajectory
    fn = traj_dir // '/' // fn_traj //'.' // str_from_num(itraj)
    write(*,*) fn
    call traj%open(fn, 'r')
    write(*,*) '  numframes ', traj%num_frames
    
    do iframe = 1, traj%num_frames
        call traj%read(traj_ver, iframe, istage, nblks, nts, coordinates, ierr)
        if (nblks < nts_beg) cycle
        ofrm = iframe - nts_beg + 1

        do j = 2, num_atoms
            rj = coordinates(:,j)
            rjm2 = dot_product(rj,rj)
            do i = 1, (j-1)
                ri = coordinates(:,i)
                rim2 = dot_product(ri,ri)
                ri_dot_rj = dot_product(ri,rj)
                rij = rj - ri
                rijm = norm2(rij)
                irijm = 1.0_rp/rijm
                irijm2 = irijm*irijm

                if (rijm >= 2.0_rp) then
                  C1 =  1.0_rp + (2.0_rp/3.0_rp)*irijm2
                  C2 =  1.0_rp - 2.0_rp*irijm2
                  consij = 0.75_rp*irijm
                else
                  C1 = 1.0_rp - 9.0_rp*rijm/(32.0_rp)
                  C2 = 3.0_rp*rijm/(32.0_rp)
                  consij = 1.0_rp
                end if

                tr = consij*( C1 + C2*rij(1)*rij(1)*irijm2 ) &
                   + consij*( C1 + C2*rij(2)*rij(2)*irijm2 ) &
                   + consij*( C1 + C2*rij(3)*rij(3)*irijm2 )

                dkirk = dkirk + 2*tr !Twice since symmetric
            end do
        end do

        dkirk = dkirk + 3*num_atoms !Diagonal components
        dkirk = dkirk*rnb*rnb/3
        rh_ifrm = 1.0_rp/dkirk

        if (ofrm > rh%len) then
            call rh%append( rh_ifrm )
            call counter%append( 1 )
        else
            call rh%set_val( ofrm, rh%get_val(ofrm)+rh_ifrm )
            call counter%set_val( ofrm, counter%get_val(ofrm)+1 )
        end if
    end do

    !Close trajectory
    call traj%close()
end do

open(newunit=fu, file='tmp.txt', action='write', status='replace')
write(fu,*) 'rh'
do i = 1, rh%len
    if (counter%get_val(i) > 0) write(fu, *) rh%get_val(i)/counter%get_val(i)
end do
close(fu)

call rh%delete()
call counter%delete()

!*******************************************************************************

end program
