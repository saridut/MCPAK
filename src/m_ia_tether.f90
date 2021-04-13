module m_ia_tether
!! **Tether potentials**
!!
!! *Style _none_: None
!! *Style _harm_: Harmonic tether
!! *Style _fene_: Fene spring tether

use m_precision
use m_constants_math
use m_ia_types, only: tether_t

implicit none

contains

!******************************************************************************

subroutine ia_tether_setup(teth)
    !! Sets up parameters for tether potentials 

    type(tether_t), intent(in out) :: teth

    select case(teth%style)
    case('none')
        continue
    case('harm')
        continue
    case('fene')
        continue
    case default
        continue
    end select

    end subroutine

!******************************************************************************

subroutine ia_get_tether_energy(ri, teth, enrg, ierr)
    !! Calculates the energy due to a tether.

    real(rp), dimension(3), intent(in) :: ri
        !!Coordinate of atom *i*
    type(tether_t), intent(in) :: teth
    real(rp), intent(out) :: enrg
    integer, intent(out)  :: ierr
    real(rp), dimension(3):: rij
    real(rp) :: rij_mag

    ierr = 0; enrg = 0.0_rp

    rij = ri - teth%point
    rij_mag = norm2(rij)

    select case (teth%style)
    case ('none')
        enrg = 0.0_rp
    case ('harm')
        call teth_harm(rij_mag, teth, enrg) 
    case ('fene')
        call teth_fene(rij_mag, teth, enrg, ierr) 
    end select

    end subroutine

!********************************************************************************

subroutine teth_harm(r, teth, enrg)
    !! Calculates energy for harmonic tether interaction.

    !params(1) = kspring, params(2) = requil (equilibrium distance)

    real(rp), intent(in) :: r
    type(tether_t), intent(in) :: teth
    real(rp), intent(out) :: enrg

    enrg = 0.5_rp*teth%params(1)*(r-teth%params(2))*(r-teth%params(2))

    end subroutine

!********************************************************************************

subroutine teth_fene(r, teth, enrg, ierr)
    !! Calculates energy for fene tether interaction.

    !params(1) = kspring, params(2) = r0 (maximum distance)

    real(rp), intent(in) :: r
    type(tether_t), intent(in) :: teth
    real(rp), intent(out) :: enrg
    integer, intent(out) :: ierr
    real(rp) :: kspring, r0, r0sq, extnsq

    ierr = 0; enrg = 0.0_rp

    kspring = teth%params(1); r0 = teth%params(2)
    r0sq = r0*r0; extnsq = r*r

    if ( r >= r0 ) then
        ierr = 1; return
    else
        enrg = -0.5_rp*kspring*r0sq*log(1.0_rp-extnsq/r0sq)
    end if

    end subroutine

!******************************************************************************

end module m_ia_tether
