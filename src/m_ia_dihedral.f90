module m_ia_dihedral
!! **Dihedral potentials**
!!
!! *Style 0: None
!!

use m_precision
use m_constants_math
use m_ia_types, only: ia_specs_t

implicit none

contains

!******************************************************************************

subroutine ia_dihedral_setup(specs)
    !! Sets up parameters for dihedral potentials 

    type(ia_specs_t), intent(in out) :: specs

    select case(specs%style)
    case('none')
        continue
    case default
        continue
    end select

    end subroutine

!******************************************************************************

subroutine ia_get_dihedral_energy(ri, rj, rk, rl, specs, enrg)
    !! Calculates the energy due to a dihedral.

    real(rp), dimension(3), intent(in) :: ri
    real(rp), dimension(3), intent(in) :: rj
    real(rp), dimension(3), intent(in) :: rk
    real(rp), dimension(3), intent(in) :: rl
    type(ia_specs_t), intent(in) :: specs
    real(rp), intent(out) :: enrg

    enrg = 0.0_rp

    select case(specs%style)
    case('none')
        enrg = 0.0_rp
    end select

    end subroutine

!******************************************************************************

end module m_ia_dihedral
