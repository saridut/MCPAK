module m_ia_types
!! **Type definitions for various interactions.**

use m_precision

implicit none

!******************************************************************************

type atm_specs_t
    character(len=8) :: name
    real(rp) :: mass
    integer  :: style !Style 0: Point particles, Style 1: Finite size particles.
end type atm_specs_t

!******************************************************************************

type ia_specs_t
    character(len=8) :: style
    real(rp), dimension(:), allocatable :: params
    real(rp), dimension(:), allocatable :: cache
    integer  :: tab_size = 0
    !real(rp) :: tmin = 0.0_rp
    !real(rp) :: tmax = 0.0_rp
    !real(rp) :: vtltmin = 0.0_rp
    !real(rp) :: vtgtmax = 0.0_rp
    !integer  :: intrp_meth = 0
    !    !! Interpolation method. 0: Linear, 1: Cubic spline
    !integer  :: flg_bc = 0
    !    !! Boundary condition for cubic splines. 0: Natural (second derivative
    !    !! equals zero at the two ends), 1: Spline is quadratic over the first
    !    !! and last interval, 2: Not-a-knot condition
    real(rp), dimension(:), allocatable :: tab_t
        !! (*tab_size*,) array. Independent variable.
    real(rp), dimension(:), allocatable :: tab_v
        !! (*tab_size*,) array. Tabulated values of the potential.
    real(rp), dimension(:), allocatable :: tab_vpp
        !! (*tab_size*,) array. Used only for cubic spline interpolation.
end type ia_specs_t

!******************************************************************************

type tether_t
    character(len=8) :: style
    real(rp), dimension(:), allocatable :: params
    real(rp), dimension(:), allocatable :: cache
    real(rp), dimension(3) :: point
    integer :: atm
end type tether_t

!******************************************************************************

type extrn_field_t
    character(len=8) :: style
    real(rp), dimension(:), allocatable :: params
    real(rp), dimension(:), allocatable :: cache
end type extrn_field_t


!******************************************************************************

end module m_ia_types
