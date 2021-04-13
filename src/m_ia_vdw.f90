module m_ia_vdw
!! **Vdw potentials**
!!
!! *Style _none_: None
!! *Style _lj_: Lennard-Jones (with cut-off)
!! *Style _hs_: Hard sphere
!! *Style _sqw_: Square-well potential
!! *Style _gauss_: Gaussian potential
!! *Style _tab_: Tabulated potential

use m_precision
use m_constants_math
use m_ia_types, only: ia_specs_t
use m_spline

implicit none

contains

!******************************************************************************

subroutine ia_vdw_setup(specs)
    !! Sets up parameters for vdw potential.

    type(ia_specs_t), intent(in out) :: specs

    select case(specs%style)
    case('none')
        continue
    case('lj')
        call vdw_lj_set(specs)
    case('hs')
        continue
    case('sqw')
        continue
    case('gauss')
        call vdw_gaussian_set(specs)
    case('tab')
        call vdw_tab_set(specs)
    case default
        continue
    end select

    end subroutine

!******************************************************************************

subroutine ia_get_vdw_energy(ri, rj, specs, enrg, ierr)
    !! Calculates the energy due to a single interacting pair of atoms.

    real(rp), dimension(3), intent(in) :: ri
    !!Coordinate of atom *i*
    real(rp), dimension(3), intent(in) :: rj
    !!Coordinate of atom *j*
    type(ia_specs_t), intent(in) :: specs
    real(rp), intent(out) :: enrg
    integer, intent(out) :: ierr
    real(rp), dimension(3) :: rij
    real(rp) :: rij_mag
    integer :: styl

    ierr = 0; enrg = 0.0_rp
    rij = rj - ri
    rij_mag = norm2(rij)

    select case(specs%style)
    case('none')
        enrg = 0.0_rp
    case('lj')
        call vdw_lj(rij_mag, specs, enrg, ierr)
    case('hs')
        call vdw_hs(rij_mag, specs, enrg, ierr)
    case('sqw')
        call vdw_sqwell(rij_mag, specs, enrg, ierr)
    case('gauss')
        call vdw_gaussian(rij_mag, specs, enrg, ierr)
    case('tab')
        call vdw_tab(rij_mag, specs, enrg, ierr)
    end select

    end subroutine

!******************************************************************************

subroutine vdw_lj_set(specs)
    !! Setter for LJ interaction.
    !! params(1) = eps, params(2) = sigma, params(3) = rcut
    !! cache(1) = sigma**2, cache(2) = U(rcut)

    type(ia_specs_t), intent(in out) :: specs
    real(rp) :: eps
    real(rp) :: sigma
    real(rp) :: rcut
    real(rp) :: sigmasq
    real(rp) :: pot_rcut
    real(rp) :: sir2, sir6, sir12

    eps = specs%params(1); sigma = specs%params(2); rcut = specs%params(3)
    sigmasq = sigma*sigma

    sir2 = sigmasq/(rcut*rcut); sir6 = sir2*sir2*sir2; sir12 = sir6*sir6
    pot_rcut = 4*eps*(sir12 - sir6)

    allocate( specs%cache(2) )

    specs%cache(1) = sigmasq
    specs%cache(2) = pot_rcut

    end subroutine

!******************************************************************************

pure subroutine vdw_lj(r, specs, enrg, ierr)
    !! Calculates energy for LJ interaction.
    !! params(1) = eps, params(2) = sigma, params(3) = rcut
    !! cache(1) = sigma**2, cache(2) = U(rcut)

    real(rp), intent(in) :: r
    type(ia_specs_t), intent(in) :: specs
    real(rp), intent(out) :: enrg
    integer, intent(out) :: ierr
    real(rp) :: eps, sigma
    real(rp) :: rcut, sigmasq
    real(rp) :: pot_rcut
    real(rp) :: sir2, sir6, sir12

    ierr = 0; enrg = 0.0_rp

    eps = specs%params(1); sigma = specs%params(2); rcut = specs%params(3)
    sigmasq = specs%cache(1)
    pot_rcut = specs%cache(2)

    if ( r >= rcut ) then
        enrg = 0.0_rp
    else
        sir2 = sigmasq/(r*r); sir6 = sir2*sir2*sir2; sir12 = sir6*sir6
        enrg = 4*eps*(sir12 - sir6) - pot_rcut
    end if

    end subroutine

!******************************************************************************

pure subroutine vdw_hs(r, specs, enrg, ierr)
    !! Calculates energy for hard sphere interaction. Overlap will be indicated
    !! with `ierr = 1`.
    !! params(1) = r1, params(2) = r2, params(3) = tol

    real(rp), intent(in) :: r
    type(ia_specs_t), intent(in) :: specs
    real(rp), intent(out) :: enrg
    integer, intent(out) :: ierr

    ierr = 0; enrg = 0.0_rp
    if ( r < (specs%params(1)+specs%params(2)-specs%params(3)) ) ierr = 1

    end subroutine

!******************************************************************************

pure subroutine vdw_sqwell(r, specs, enrg, ierr)
    !! Calculates energy for square well interaction. Overlap will be indicated
    !! with `ierr = 1`.
    !! params(1) = eps, params(2) = sigma, params(3) = width

    real(rp), intent(in) :: r
    type(ia_specs_t), intent(in) :: specs
    real(rp), intent(out) :: enrg
    integer, intent(out) :: ierr

    if ( r <= specs%params(2) ) then
        ierr = 1; enrg = 0.0_rp
    else if ( r >= (specs%params(2)+specs%params(3)) ) then
        ierr = 0; enrg = 0.0_rp
    else
        ierr = 0; enrg = -specs%params(1)
    end if

    end subroutine

!******************************************************************************

subroutine vdw_gaussian_set(specs)
    !! Setter for Gaussian interaction.
    !! params(1) = A, params(2) = B, params(3) = rcut
    !! cache(1) = V(rcut)

    type(ia_specs_t), intent(in out) :: specs
    real(rp) :: A, B, rcut
    real(rp) :: pot_rcut

    A = specs%params(1); B = specs%params(2); rcut = specs%params(3)

    pot_rcut = A*exp(-B*rcut**2)

    allocate( specs%cache(1) )

    specs%cache(1) = pot_rcut

    end subroutine

!******************************************************************************

pure subroutine vdw_gaussian(r, specs, enrg, ierr)
    !! Calculates energy for Gaussian interaction.
    !! params(1) = A, params(2) = B, params(3) = rcut
    !! cache(1) = V(rcut)

    real(rp), intent(in) :: r
    type(ia_specs_t), intent(in) :: specs
    real(rp), intent(out) :: enrg
    integer, intent(out) :: ierr
    real(rp) :: A, B, rcut
    real(rp) :: pot_rcut
    real(rp) :: exrs

    ierr = 0; enrg = 0.0_rp

    A = specs%params(1); B = specs%params(2); rcut = specs%params(3)
    pot_rcut = specs%cache(1)

    if ( r < rcut ) then
        exrs = exp(-B*r*r)
        enrg = A*exrs - pot_rcut
    end if

    end subroutine

!******************************************************************************

subroutine vdw_tab_set(specs)
    !! Setter for tabulated interaction.
    !! params(1) = vtltmin, params(2) = vtgtmax, params(3) = intrp_meth,
    !! params(4) = flag_bc
    !! cache(1) = tmin, cache(2) = tmax

    type(ia_specs_t), intent(in out) :: specs
    integer :: intrp_meth, flag_bc
    integer :: n

    intrp_meth = int(specs%params(3))
    n = specs%tab_size

    allocate( specs%cache(2) ) 
    specs%cache(1) = minval(specs%tab_t)
    specs%cache(2) = maxval(specs%tab_t)

    if (intrp_meth == 1) then
        allocate( specs%tab_vpp(n) )
        flag_bc = int(specs%params(4))
        if (flag_bc == 0) then
            !Natural BC: second derivative equals zero at the two ends
            call spline_cubic_set(n, specs%tab_t, specs%tab_v, &
                2, 0.0_rp, 2, 0.0_rp, specs%tab_vpp)
        else if (flag_bc == 1) then
            !Spline is quadratic over the first and last interval
            call spline_cubic_set(n, specs%tab_t, specs%tab_v, &
                0, 0.0_rp, 0, 0.0_rp, specs%tab_vpp)
        else if (flag_bc == 2) then
            !Not-a-knot condition
            call spline_cubic_set(n, specs%tab_t, specs%tab_v, &
                3, 0.0_rp, 3, 0.0_rp, specs%tab_vpp)
        end if
    end if

    end subroutine

!******************************************************************************

subroutine vdw_tab(r, specs, enrg, ierr)
    !! Calculates energy for tabulated interaction.
     
    !  params(1) = vtltmin, params(2) = vtgtmax, params(3) = intrp_meth,
    !  params(4) = flag_bc
    !  cache(1) = tmin, cache(2) = tmax

    real(rp), intent(in) :: r
    type(ia_specs_t), intent(in) :: specs
    real(rp), intent(out) :: enrg
    integer, intent(out) :: ierr
    real(rp) :: tmin, tmax, vtltmin, vtgtmax
    real(rp) :: vp, vpp
    integer :: intrp_meth, n

    ierr = 0; enrg = 0.0_rp

    n = specs%tab_size
    vtltmin = specs%params(1); vtgtmax = specs%params(2)
    tmin = specs%cache(1); tmax = specs%cache(2)
    intrp_meth = int(specs%params(3))

    if ( r < tmin ) then
        enrg = vtltmin
    else if ( r > tmax ) then
        enrg = vtgtmax
    else
        if (intrp_meth == 0) then
            call spline_linear_val(n, specs%tab_t, specs%tab_v, r, enrg, vp)
        else if (intrp_meth == 1) then
            call spline_cubic_val(n, specs%tab_t, specs%tab_v, specs%tab_vpp, &
                r, enrg, vp, vpp)
        end if
    end if

    end subroutine

!******************************************************************************

end module m_ia_vdw
