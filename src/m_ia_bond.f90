module m_ia_bond
!! **Bond potentials**
!!
!! *Style _none_: None
!! *Style _harm_: Harmonic
!! *Style _fene_: Shifted-FENE
!! *Style _kg_ : Kremer-Grest
!! *Style _tab_ : Tabulated
!! 

use m_precision
use m_constants_math
use m_ia_types, only: ia_specs_t
use m_spline

implicit none

contains

!******************************************************************************

subroutine ia_bond_setup(specs)
    !! Sets up parameters for bond potentials 

    type(ia_specs_t), intent(in out) :: specs

    select case(specs%style)
    case('none')
        continue
    case('harm')
        continue
    case('fene')
        continue
    case('kg')
        call bond_kg_set(specs)
    case('tab')
        call bond_tab_set(specs)
    case default
        continue
    end select

    end subroutine

!******************************************************************************

subroutine ia_get_bond_energy(ri, rj, specs, enrg, ierr)
    !! Calculates the energy due to a bond.

    real(rp), dimension(3), intent(in) :: ri
    !!Coordinate of atom *i*
    real(rp), dimension(3), intent(in) :: rj
    !!Coordinate of atom *j*
    type(ia_specs_t), intent(in) :: specs
    real(rp), intent(out) :: enrg
    integer, intent(out)  :: ierr
    real(rp), dimension(3) :: rij
    real(rp) :: rij_mag

    ierr = 0; enrg = 0.0_rp
    rij = rj - ri
    rij_mag = norm2(rij)

    select case (specs%style)
    case ('none')
        enrg = 0.0_rp
    case ('harm')
        call bond_harm(rij_mag, specs, enrg) 
    case ('fene')
        call bond_fene(rij_mag, specs, enrg, ierr) 
    case ('kg')
        call bond_kg(rij_mag, specs, enrg, ierr) 
    case ('tab')
        call bond_tab(rij_mag, specs, enrg, ierr) 
    end select

    end subroutine

!********************************************************************************

subroutine bond_harm(r, specs, enrg)
    !! Calculates energy for harmonic bond interaction.
    !params(1) = kspring, params(2) = requil (equilibrium distance)

    real(rp), intent(in) :: r
    type(ia_specs_t), intent(in) :: specs
    real(rp), intent(out) :: enrg

    enrg = 0.5_rp*specs%params(1)*(r-specs%params(2))*(r-specs%params(2))

    end subroutine

!********************************************************************************

pure subroutine bond_fene(r, specs, enrg, ierr)
    !! Calculates energy for FENE bond interaction. If bond length exceeds
    !! maximum extensible spring length, an error will be reported as `ierr = 1`.
    !! params(1) = kspring, params(2) = r0, params(3) = delta

    real(rp), intent(in) :: r
    type(ia_specs_t), intent(in) :: specs
    real(rp), intent(out) :: enrg
    integer, intent(out)  :: ierr
    real(rp) :: kspring
    real(rp) :: r0
    real(rp) :: delta
    real(rp) :: extn
    real(rp) :: extnsq
    real(rp) :: r0sq

    ierr = 0; enrg = 0.0_rp

    kspring = specs%params(1); r0 = specs%params(2); delta = specs%params(3)
    r0sq = r0*r0
    extn = r - delta
    extnsq = extn*extn

    if ( r >= (r0+delta) ) then
        ierr = 1; return
    else
        enrg = -0.5_rp*kspring*r0sq*log(1.0_rp-extnsq/r0sq)
    end if

    end subroutine

!******************************************************************************

subroutine bond_kg_set(specs)
    !! Setter for FENE bond interaction.
    !! params(1) = kspring, params(2) = r0, params(3) = eps, params(4) = sigma
    !! cache(1) = r0^2, cache(2) = sigma^2, cache(3) = 2^(1/6)*sigma

    type(ia_specs_t), intent(in out) :: specs
    real(rp) :: kspring, r0, eps, sigma

    kspring = specs%params(1); r0 = specs%params(2)
    eps = specs%params(3); sigma = specs%params(4)

    allocate(specs%cache(3))
    specs%cache(1) = r0*r0
    specs%cache(2) = sigma*sigma
    specs%cache(3) = math_sxrt2*sigma

    end subroutine

!********************************************************************************

pure subroutine bond_kg(r, specs, enrg, ierr)
    !! Calculates energy for Kremer-Grest bond interaction. If bond length exceeds
    !! maximum extensible spring length, an error will be reported as `ierr = 1`.

    real(rp), intent(in) :: r
    type(ia_specs_t), intent(in) :: specs
    real(rp), intent(out) :: enrg
    integer, intent(out)  :: ierr
    real(rp) :: kspring, r0, eps, sigma
    real(rp) :: r0sq, sigmasq, rcut
    real(rp) :: rsq, sir2, sir12, sir6

    ierr = 0; enrg = 0.0_rp

    kspring = specs%params(1); r0 = specs%params(2)
    eps = specs%params(3); sigma = specs%params(4)
    r0sq = specs%cache(1); sigmasq = specs%cache(2); rcut = specs%cache(3)
    rsq = r*r

    if ( r >= r0 ) then
        ierr = 1; return
    else if ( r <= 0.25_rp ) then
        !If atoms are less than quarter-radius away, energy is too high.
        ierr = 1; return
    else if ( (r >= rcut) .and. (r < r0) ) then
        enrg = -0.5_rp*kspring*r0sq*log(1.0_rp - rsq/r0sq)
    else
        sir2 = sigmasq/rsq; sir6 = sir2*sir2*sir2; sir12 = sir6*sir6
        enrg = -0.5_rp*kspring*r0sq*log(1.0_rp - rsq/r0sq) &
                + 4*eps*(sir12 - sir6) + eps
    end if

    end subroutine

!******************************************************************************

subroutine bond_tab_set(specs)
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

subroutine bond_tab(r, specs, enrg, ierr)
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

end module m_ia_bond
