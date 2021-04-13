module m_ia_angle
!! **Angle potentials**
!!
!! *Style _none_: None
!! *Style _cos_: Cosine. \[U(\theta) = k(1 - \cos \theta) \], \[\theta\] is
!!      the complementary bond angle.
!! *Style _harm_: Harmonic. \[U(\theta) = k \theta^2 \], \[\theta\] is
!!      the complementary bond angle.
!! *Style _tab_ : Tabulated potential

use m_precision
use m_constants_math
use m_ia_types, only: ia_specs_t
use m_spline

implicit none

contains

!******************************************************************************

subroutine ia_angle_setup(specs)
    !! Sets up parameters for angle potentials 

    type(ia_specs_t), intent(in out) :: specs

    select case(specs%style)
    case('none')
        continue
    case('cos')
        continue
    case('harm')
        continue
    case('tab')
        call ang_tab_set(specs)
    case default
        continue
    end select

    end subroutine

!******************************************************************************

subroutine ia_get_angle_energy(rim1, ri, rip1, specs, enrg)
    !! Calculates the energy due to an angle potential.

    real(rp), dimension(3), intent(in) :: rim1
    real(rp), dimension(3), intent(in) :: ri
    real(rp), dimension(3), intent(in) :: rip1
    type(ia_specs_t), intent(in) :: specs
    real(rp), intent(out) :: enrg
    real(rp), dimension(3) :: q1
    real(rp), dimension(3) :: q2
    real(rp), dimension(3) :: q1hat
    real(rp), dimension(3) :: q2hat
    real(rp) :: q1mag
    real(rp) :: q2mag
    real(rp) :: ctheta
    integer  :: ierr

    q1 = ri - rim1; q2 = rip1 - ri
    q1mag = norm2(q1); q2mag = norm2(q2)
    q1hat = q1/q1mag; q2hat = q2/q2mag
    ctheta = dot_product(q1hat, q2hat)
    !Floating point correction
    if (ctheta > 1.0_rp) ctheta = 1.0_rp 
    if (ctheta < -1.0_rp) ctheta = -1.0_rp

    select case(specs%style)
    case('none')
        enrg = 0.0_rp
    case('cos')
        call ang_cos(ctheta, specs, enrg)
    case('harm')
        call ang_harm(ctheta, specs, enrg)
    case('tab')
        call ang_tab(ctheta, specs, enrg, ierr)
    end select

    end subroutine

!******************************************************************************

subroutine ang_cos(ctheta, specs, enrg)
    !! Calculates energy for angular cosine interaction.
    !! params(1) = k

    real(rp), intent(in) :: ctheta
    type(ia_specs_t), intent(in) :: specs
    real(rp), intent(out) :: enrg
    real(rp) :: kang

    kang = specs%params(1)
    enrg = kang*(1.0_rp - ctheta)

    end subroutine

!******************************************************************************

subroutine ang_harm(ctheta, specs, enrg)
    !! Calculates energy for harmonic angular interaction.
    !! params(1) = k

    real(rp), intent(in) :: ctheta
    type(ia_specs_t), intent(in) :: specs
    real(rp), intent(out) :: enrg
    real(rp) :: kang
    real(rp) :: theta

    kang = specs%params(1)
    theta = acos(ctheta)
    enrg = kang*theta*theta

    end subroutine

!******************************************************************************

subroutine ang_tab_set(specs)
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

subroutine ang_tab(r, specs, enrg, ierr)
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

end module m_ia_angle
