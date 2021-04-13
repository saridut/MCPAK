module m_control_io

use m_precision
use m_strings
use m_globals

implicit none

contains

!******************************************************************************

subroutine read_control(fn)
    !! Reads simulation parameters from file

    character(len=*), intent(in) :: fn !! Name of parameters file.
    character(len=:), allocatable :: key
    character(len=:), allocatable :: val
    character(len=mxrdln) :: line
    character(len=1) :: cstr = '#' !Comment string
    integer :: fu
    integer :: ios

    open (newunit=fu, file = fn, action = 'read', status = 'old')
    
    do 
        call readline(fu, line, cstr, ios)
        if (ios /= 0) return

        call str_get_keyval(line, key, val)

        if (key=='use_mc_pvt') read(val, *) use_mc_pvt
        if (key=='use_mc_res_pvt') read(val, *) use_mc_res_pvt
        if (key=='use_mc_dspl') read(val, *) use_mc_dspl
        if (key=='use_mc_crnk') read(val, *) use_mc_crnk
        if (key=='use_mc_sc_pvt') read(val, *) use_mc_sc_pvt
        if (key=='use_verlet_tab') read(val, *) use_verlet_tab
        if (key=='mcm_mxdspl') mcm_mxdspl = str_to_d(val)
        if (key=='rcutoff') rcutoff = str_to_d(val)
        if (key=='tskin') tskin = str_to_d(val)
        if (key=='excluded_atoms') excluded_atoms = str_to_i(val)
        if (key=='lvdw') read(val,*) lvdw
        if (key=='lelectrostatics') read(val,*) lelectrostatics

        if (key=='nts_log')    nts_log    = int(str_to_d(val))
        if (key=='nts_dump')   nts_dump   = int(str_to_d(val))
        if (key=='nts_eql')    nts_eql    = int(str_to_d(val))
        if (key=='nts_eql_samp') nts_eql_samp = int(str_to_d(val))
        if (key=='block_size') block_size = int(str_to_d(val))
        if (key=='nblks_sim')  nblks_sim  = int(str_to_d(val))

        if (key == 'fn_cfg')    fn_cfg    = val
        if (key == 'fn_revive') fn_revive = val
        if (key == 'fn_stats')  fn_stats  = val
        if (key == 'fn_traj')   fn_traj   = val

        if (key == 'lrevive') read(val,*) lrevive
        if (key == 'read_seed') read(val,*) read_seed
        if (key == 'write_seed') read(val,*) write_seed
        if (key == 'write_eql_stats') read(val,*) write_eql_stats
        if (key == 'write_traj') read(val,*) write_traj
    end do

    close (fu)

    end subroutine

!******************************************************************************

subroutine write_control(fn)
    !! Write simulation parameters to file

    character(len=*), intent(in) :: fn
        !! File name
    integer :: fu

    open(newunit=fu, file=fn, action='write', status='unknown')

    write(fu, '(a,t20,l1)'  ) 'use_mc_pvt', use_mc_pvt
    write(fu, '(a,t20,l1)'  ) 'use_mc_res_pvt', use_mc_res_pvt
    write(fu, '(a,t20,l1)'  ) 'use_mc_dspl', use_mc_dspl
    write(fu, '(a,t20,l1)'  ) 'use_mc_crnk', use_mc_crnk
    write(fu, '(a,t20,l1)'  ) 'use_mc_sc_pvt', use_mc_sc_pvt
    write(fu, '(a,t20,l1)'  ) 'use_verlet_tab', use_verlet_tab
    write(fu, '(a,t20,g0.6)') 'mcm_mxdspl', mcm_mxdspl
    write(fu, '(a,t20,g0.6)') 'rcutoff', rcutoff
    write(fu, '(a,t20,g0.6)') 'tskin', tskin
    write(fu, '(a,t20,i0)'  ) 'excluded_atoms', excluded_atoms
    write(fu, '(a,t20,l1)'  ) 'lvdw', lvdw
    write(fu, '(a,t20,l1)'  ) 'lelectrostatics', lelectrostatics

    write(fu, *)
    write(fu, '(a,t20,i0)')   'nts_log',    nts_log
    write(fu, '(a,t20,i0)')   'nts_dump',   nts_dump
    write(fu, '(a,t20,i0)')   'nts_eql',    nts_eql
    write(fu, '(a,t20,i0)')   'nts_eql_samp', nts_eql_samp
    write(fu, '(a,t20,i0)')   'block_size', block_size
    write(fu, '(a,t20,i0)')   'nblks_sim',  nblks_sim

    write(fu, *)
    write(fu, '(a,t20,a)') 'fn_cfg', fn_cfg
    write(fu, '(a,t20,a)') 'fn_revive', fn_revive
    write(fu, '(a,t20,a)') 'fn_stats', fn_stats
    write(fu, '(a,t20,a)') 'fn_traj', fn_traj

    write(fu, *)
    write(fu, '(a,t20,l1)') 'lrevive', lrevive
    write(fu, '(a,t20,l1)') 'read_seed', read_seed
    write(fu, '(a,t20,l1)') 'write_seed', write_seed
    write(fu, '(a,t20,l1)') 'write_eql_stats', write_eql_stats
    write(fu, '(a,t20,l1)') 'write_traj', write_traj

    close(fu)
    
    end subroutine

!******************************************************************************

end module m_control_io
