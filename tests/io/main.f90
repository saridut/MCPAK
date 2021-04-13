program main

use m_precision
use m_strings
use m_globals
use m_config_io
use m_setup, only: config_clear

implicit none

!*******************************************************************************

character(len=:), allocatable :: fn_cfg_in
character(len=:), allocatable :: fn_cfg_out
character(len=:), allocatable :: fn_ldf
character(len=:), allocatable :: fn_xyz

fn_cfg_in = 'tcfg.cfg'
fn_cfg_out = 'tcfg-out.cfg'
fn_ldf  = 'ldf.txt'
fn_xyz  = 'test.xyz'

call read_config(fn_cfg_in)
call write_config(fn_cfg_out, 'test config')
call write_ldf(fn_ldf, 'test data')
call write_xyz(fn_xyz, 'test data')

call config_clear()

!*******************************************************************************

end program
