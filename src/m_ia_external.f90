module m_ia_external
!! **External potentials**
!!
!! *Style _none_: None
!! *Style _xfrc_: Pulling force along x-axis
!! *Style _hwall_: Hard wall

use m_precision
use m_constants_math
use m_ia_types, only: extrn_field_t

implicit none

private

public :: ia_external_setup, ia_get_external_energy

contains

!******************************************************************************

subroutine ia_external_setup(ef)
    !! Sets up parameters for external potentials. Placeholder.

    type(extrn_field_t), intent(in out) :: ef

    select case(ef%style)
    case ('none')
        continue
    case default
        continue
    end select

    end subroutine

!******************************************************************************

subroutine ia_get_external_energy(ef, coordinates, enrg, ierr)
    !! Calculates the energy due to an external field.

    type(extrn_field_t), intent(in) :: ef
    real(rp), dimension(:,:), intent(in) :: coordinates
    real(rp), intent(out)  :: enrg
    integer, intent(out)   :: ierr
    integer :: iatm, jatm
    integer :: m
    real(rp) :: frcx, v, sn

    ierr = 0

    select case(ef%style)
    case('none')
        continue
    case('xfrc')
        iatm = int(ef%params(1))
        jatm = int(ef%params(2))
        frcx = ef%params(3)
        enrg = -frcx*(coordinates(1,jatm) - coordinates(1,iatm))
    case('hwall')
        enrg = 0.0_rp
        m = int(ef%params(1)) 
        v = ef%params(2)
        sn = ef%params(3)
        if (sn > 0.0_rp) then
            if ( any(coordinates(m,:) < v) ) ierr = 1
        else
            if ( any(coordinates(m,:) > v) ) ierr = 1
        end if
    case default
        enrg = 0.0_rp
    end select

    end subroutine

!******************************************************************************

end module m_ia_external
