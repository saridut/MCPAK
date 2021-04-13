module m_stats_io

use m_precision
use m_constants_math
use m_strings
use m_globals

implicit none

private

public :: stats_init, stats_finish, stats_collect, stats_write, stats_write_eql

!Overall properties
real(rp), dimension(3) :: gt_ev_avg = 0.0_rp
    !! Gyration tensor eigen values (in descending order)
real(rp) :: rgsq_avg = 0.0_rp
    !! Gyration radius squared
real(rp) :: reedsq_avg = 0.0_rp
    !! End-to-end distance squared
real(rp) :: reedx_avg = 0.0_rp
    !! Extension along x-axis
real(rp) :: asph_avg = 0.0_rp
    !! Asphericity
real(rp) :: prol_avg = 0.0_rp
    !! Prolateness
real(rp) :: dkirk_avg = 0.0_rp
    !! Kirkwood diffusivity
real(rp) :: rh_avg = 0.0_rp
    !! Hydrodynamic radius
real(rp) :: iv_avg = 0.0_rp
    !! Intrinsic viscosity

!Backbone properties
real(rp) :: rgsq_bbone_avg = 0.0_rp
    !! Backbone gyration radius squared

!Side chain properties
real(rp), dimension(3) :: gt_ev_sc_avg = 0.0_rp
    !! Side chain gyration tensor eigen values (in descending order)
real(rp) :: rgsq_sc_avg = 0.0_rp
    !! Side chain gyration radius squared
real(rp) :: reedsq_sc_avg = 0.0_rp
    !! Side chain end-to-end distance squared
real(rp) :: asph_sc_avg = 0.0_rp
    !! Side chain asphericity
real(rp) :: prol_sc_avg = 0.0_rp
    !! Side chain prolateness

!Miscellaneous stats
integer  :: avg_cnt = 0
real(rp) :: bndlen = 0.0_rp
    !! Average bond length at a single instant.
real(rp) :: bndlen_min = huge(0.0_rp)
    !! Maximum bond length
real(rp) :: bndlen_max = 0.0_rp
    !! Minimum bond length
integer :: fu_stats
    !! Unit number of fn_stats file

contains

!******************************************************************************

subroutine stats_write_eql(lhdr)
    !! Writes equilibration statistics

    logical, intent(in) :: lhdr
    real(rp) :: energy_tot
    real(rp) :: rgsq, reedsq
    real(rp), dimension(3) :: com, reev
    integer  :: i, na_bbone

    if (lhdr) then
        write(fu_stats,'(a12,2x,*(a16,2x))') 'nts',                     &
            'rgsq', 'reedsq',                                           &
            'energy_bond', 'energy_angle', 'energy_dihedral', 'energy_vdw',&
            'energy_tether', 'energy_external', 'energy_tot'
        return
    end if

    energy_tot = energy_bond + energy_angle + energy_dihedral + energy_vdw &
        + energy_tether + energy_external

    if (num_branches == 0) then
        na_bbone = num_atoms
    else
        na_bbone = branches(2,1)
    end if

    com = sum(coordinates,2)/num_atoms
    rgsq = 0.0_rp
    do i = 1, num_atoms
        rgsq = rgsq + sum((coordinates(:,i)-com)**2)
    end do
    rgsq = rgsq/num_atoms

    reev = coordinates(:,na_bbone) - coordinates(:,1)
    reedsq = sum(reev**2)

    write(fu_stats,'(i12,2x,*(es16.7,2x))') nts,         &
        rgsq, reedsq,                                    &
        energy_bond, energy_angle, energy_dihedral, energy_vdw, &
        energy_tether, energy_external, energy_tot

    end subroutine

!******************************************************************************

subroutine stats_write(lhdr)
    !! Writes production statistics and zeros out accumulators

    logical, intent(in) :: lhdr
    real(rp), dimension(mx_mcmtyp) :: ar

    if (lhdr) then
        write(fu_stats,'(a12,2x,*(a16,2x))') 'nblks',                   &
            'rgsq', 'reedsq', 'asph', 'prol',                           &
            'rgsq_bbone', 'rgsq_sc', 'reedsq_sc', 'asph_sc', 'prol_sc', &
            'gt_ev1', 'gt_ev2', 'gt_ev3',                               &
            'gt_ev_sc1', 'gt_ev_sc2', 'gt_ev_sc3',                      &
            'mc_bb_pvt_ar', 'mc_crnk_ar', 'mc_sc_pvt_ar', 'mc_res_pvt_ar', &
            'mc_dspl_ar'
        return
    end if

    ar = 0.0_rp
    where(mc_rec(1,:) > 0) ar = real(mc_rec(2,:),rp)/real(mc_rec(1,:),rp)

    !Note: gt_ev(2) and gt_ev(3) are flipped so that the three eigen values are
    !in descending order
    write(fu_stats,'(i12,2x,*(es16.7,2x))') nblks,                            &
        rgsq_avg, reedsq_avg, asph_avg, prol_avg,                             &
        rgsq_bbone_avg, rgsq_sc_avg, reedsq_sc_avg, asph_sc_avg, prol_sc_avg, &
        gt_ev_avg(1), gt_ev_avg(3), gt_ev_avg(2),                             &
        gt_ev_sc_avg(1), gt_ev_sc_avg(3), gt_ev_sc_avg(2),                    &
        ar 

    rgsq_avg=0.0_rp; reedsq_avg=0.0_rp
    asph_avg=0.0_rp; prol_avg=0.0_rp
    rgsq_bbone_avg=0.0_rp; rgsq_sc_avg=0.0_rp; reedsq_sc_avg=0.0_rp
    asph_sc_avg=0.0_rp; prol_sc_avg=0.0_rp
    gt_ev_avg = 0.0_rp
    gt_ev_sc_avg = 0.0_rp

    avg_cnt = 0
    mc_rec = 0

    end subroutine

!******************************************************************************

subroutine stats_init()
    !! Set up for stats collection

    logical :: lexists

    if (lrevive) then
        !Restarting simulation
        if (leql) then
            !Equilibration
            if (write_eql_stats) then
                !Check if file exists
                inquire(file=fn_stats//'.eq'//trim(adjustl(job_tag)), exist=lexists)
                if (lexists) then
                    !Open existing file for appending equilibration statistics
                    open(newunit=fu_stats, file=fn_stats//'.eq'//trim(adjustl(job_tag)), &
                        action='write', position='append', status='old')
                else
                    !Open new file for writing equilibration statistics
                    open(newunit=fu_stats, file=fn_stats//'.eq'//trim(adjustl(job_tag)), &
                        action='write', status='new')
                    !Write header
                    call stats_write_eql(.true.)
                end if
            end if
        else
            !Production
            !Check if file exists
            inquire(file=fn_stats//trim(adjustl(job_tag)), exist=lexists)
            if (lexists) then
                !Open existing file for appending statistics
                open(newunit=fu_stats, file=fn_stats//trim(adjustl(job_tag)), &
                    action='write', position='append', status='old')
            else
                !Open new file for writing statistics
                open(newunit=fu_stats, file=fn_stats//trim(adjustl(job_tag)), &
                    action='write', status='new')
                !Write header
                call stats_write(.true.)
            end if
        end if
    else
        !New simulation
        if (leql) then
            !Open file for writing equilibration statistics
            if (write_eql_stats) then
                open(newunit=fu_stats, file=fn_stats//'.eq'//trim(adjustl(job_tag)), &
                    action='write', status='replace')
                call stats_write_eql(.true.)
            end if
        else
            !Open file for writing production statistics
            open(newunit=fu_stats, file=fn_stats//trim(adjustl(job_tag)), &
                action='write', status='replace')
            call stats_write(.true.)
        end if
    end if

    end subroutine

!******************************************************************************

subroutine stats_finish()

    integer :: fu
    logical :: lopened

    if (leql) then
        inquire(file=fn_stats//'.eq'//trim(adjustl(job_tag)), number=fu, &
            opened=lopened)
        if (lopened) close(fu)
    end if

    if (.not. leql) then
        inquire(file=fn_stats//trim(adjustl(job_tag)), number=fu, &
            opened=lopened)
        if (lopened) close(fu)
    end if

    end subroutine

!******************************************************************************

subroutine stats_collect()
    !! Calculates and accumulates production statistics

    real(rp), dimension(3) :: gt_ev, gt_ev_sc, gt_ev_isc
    real(rp), dimension(3,3) :: S
    real(rp), dimension(3) :: ri, rj, rij
    real(rp), dimension(3) :: com, com_bbone, com_sc
    real(rp), dimension(3) :: reev
    real(rp), dimension(3) :: r_tp
    real(rp) :: rgsq, reedsq, reedx, asph, prol
    real(rp) :: rgsq_bbone
    real(rp) :: rgsq_sc, reedsq_sc, asph_sc, prol_sc
    real(rp) :: rgsq_isc, reedsq_isc, asph_isc, prol_isc
    real(rp) :: dkirk, rh, iv
    real(rp) :: rnb
    real(rp) :: radsq_tot
    real(rp) :: tr
    real(rp) :: rim2, rjm2
    real(rp) :: rijm, irijm, irijm2
    real(rp) :: ri_dot_rj
    real(rp) :: C1, C2, consij
    integer :: n_sc
    integer :: na_bbone 
    integer :: na_sc
    integer :: ibr, ibeg, iend, ia_beg, ia_end
    integer :: i, j, k

    rnb = 1.0_rp/num_atoms

    !Number of backbone atoms
    if (num_branches == 0) then
        na_bbone = num_atoms
        n_sc = 0 !Number of side chains
    else
        na_bbone = branches(2,1)
        n_sc = num_branches - 1  !Number of side chains
    end if

    rgsq = 0.0_rp; reedsq = 0.0_rp; reedx = 0.0_rp
    asph = 0.0_rp; prol = 0.0_rp
    rgsq_bbone = 0.0_rp 
    rgsq_sc = 0.0_rp; reedsq_sc = 0.0_rp; asph_sc = 0.0_rp; prol_sc = 0.0_rp
    gt_ev = 0.0_rp;
    gt_ev_sc = 0.0_rp;

    dkirk = 0.0_rp; rh = 0.0_rp; iv = 0.0_rp

    S = 0.0_rp; reev = 0.0_rp; 

    com = rnb*sum(coordinates,2)
    do i = 1, num_atoms
        ri = coordinates(:,i) - com
        S(1,1) = S(1,1) + ri(1)*ri(1)
        S(1,2) = S(1,2) + ri(1)*ri(2)
        S(1,3) = S(1,3) + ri(1)*ri(3)
        S(2,2) = S(2,2) + ri(2)*ri(2)
        S(2,3) = S(2,3) + ri(2)*ri(3)
        S(3,3) = S(3,3) + ri(3)*ri(3)
    end do
    S = S*rnb
    rgsq = S(1,1) + S(2,2) + S(3,3)
    
    call dsyevc3(S, gt_ev)
    call calc_shape(gt_ev(1), gt_ev(3), gt_ev(2), asph, prol) 

    reev = coordinates(:,na_bbone) - coordinates(:,1)
    reedsq = reev(1)*reev(1) + reev(2)*reev(2) + reev(3)*reev(3)
    reedx = reev(1)

    com_bbone = sum(coordinates(:,1:na_bbone),2)/na_bbone
    do i = 1, na_bbone
        rgsq_bbone = rgsq_bbone + sum((coordinates(:,i)-com_bbone)**2)
    end do
    rgsq_bbone = rgsq_bbone/na_bbone

    do ibr = 2, num_branches
        ia_beg = branches(3,ibr)
        ia_end = ia_beg + branches(2,ibr) - 1
        r_tp = coordinates(:,branches(1,ibr))

        !Add end-to-end distance squared of this side chain
        reedsq_sc = reedsq_sc + sum((coordinates(:,ia_end) - r_tp)**2)

        !Number of side chain atoms. Adding 1 for the tether point
        na_sc = branches(2,ibr) + 1

        !Get side chain c.o.m.
        com_sc = r_tp + sum(coordinates(:,ia_beg:ia_end),2)
        com_sc = com_sc/na_sc

        !Get radius of gyration squared of this side chain
        S = 0.0_rp
        ri = r_tp - com_sc
        S(1,1) = S(1,1) + ri(1)*ri(1)
        S(1,2) = S(1,2) + ri(1)*ri(2)
        S(1,3) = S(1,3) + ri(1)*ri(3)
        S(2,2) = S(2,2) + ri(2)*ri(2)
        S(2,3) = S(2,3) + ri(2)*ri(3)
        S(3,3) = S(3,3) + ri(3)*ri(3)

        do i = ia_beg, ia_end
            ri = coordinates(:,i) - com_sc
            S(1,1) = S(1,1) + ri(1)*ri(1)
            S(1,2) = S(1,2) + ri(1)*ri(2)
            S(1,3) = S(1,3) + ri(1)*ri(3)
            S(2,2) = S(2,2) + ri(2)*ri(2)
            S(2,3) = S(2,3) + ri(2)*ri(3)
            S(3,3) = S(3,3) + ri(3)*ri(3)
        end do
        S = S/na_sc
        rgsq_isc = S(1,1) + S(2,2) + S(3,3)

        call dsyevc3(S, gt_ev_isc)
        call calc_shape(gt_ev_isc(1), gt_ev_isc(3), gt_ev_isc(2), &
            asph_isc, prol_isc) 
        !Add radius of gyration, etc. of this side chain
        rgsq_sc = rgsq_sc + rgsq_isc
        gt_ev_sc = gt_ev_sc + gt_ev_isc
        asph_sc = asph_sc + asph_isc
        prol_sc = prol_sc + prol_isc
    end do

    if (n_sc > 0) then
        rgsq_sc = rgsq_sc/n_sc
        reedsq_sc = reedsq_sc/n_sc
        gt_ev_sc = gt_ev_sc/n_sc
        asph_sc = asph_sc/n_sc
        prol_sc = prol_sc/n_sc
    end if

    !Update averages
    avg_cnt = avg_cnt + 1
    rgsq_avg = rgsq_avg + (rgsq - rgsq_avg)/avg_cnt
    reedsq_avg = reedsq_avg + (reedsq - reedsq_avg)/avg_cnt
    reedx_avg = reedx_avg + (reedx - reedx_avg)/avg_cnt
    asph_avg = asph_avg + (asph - asph_avg)/avg_cnt
    prol_avg = prol_avg + (prol - prol_avg)/avg_cnt
    rgsq_bbone_avg = rgsq_bbone_avg + (rgsq_bbone - rgsq_bbone_avg)/avg_cnt
    rgsq_sc_avg = rgsq_sc_avg + (rgsq_sc - rgsq_sc_avg)/avg_cnt
    reedsq_sc_avg = reedsq_sc_avg + (reedsq_sc - reedsq_sc_avg)/avg_cnt
    asph_sc_avg = asph_sc_avg + (asph_sc - asph_sc_avg)/avg_cnt
    prol_sc_avg = prol_sc_avg + (prol_sc - prol_sc_avg)/avg_cnt
    gt_ev_avg = gt_ev_avg + (gt_ev - gt_ev_avg)/avg_cnt
    gt_ev_sc_avg = gt_ev_sc_avg + (gt_ev_sc - gt_ev_sc_avg)/avg_cnt


    !Double sums
    !Looping over the strictly upper triangular form
!   do j = 2, num_atoms
!       rj = coordinates(:,j)
!       rjm2 = dot_product(rj,rj)
!       do i = 1, (j-1)
!           ri = coordinates(:,i)
!           rim2 = dot_product(ri,ri)
!           ri_dot_rj = dot_product(ri,rj)
!           rij = rj - ri
!           rijm = norm2(rij)
!           irijm = 1.0_rp/rijm
!           irijm2 = irijm*irijm

!           if (rijm >= 2.0_rp) then
!             C1 =  1.0_rp + (2.0_rp/3.0_rp)*irijm2
!             C2 =  1.0_rp - 2.0_rp*irijm2
!             consij = 0.75_rp*irijm
!           else
!             C1 = 1.0_rp - 9.0_rp*rijm/(32.0_rp)
!             C2 = 3.0_rp*rijm/(32.0_rp)
!             consij = 1.0_rp
!           end if

!           tr = consij*( C1 + C2*rij(1)*rij(1)*irijm2 ) &
!              + consij*( C1 + C2*rij(2)*rij(2)*irijm2 ) &
!              + consij*( C1 + C2*rij(3)*rij(3)*irijm2 )

!           dkirk = dkirk + 2*tr !Twice since symmetric

!           tr = ri_dot_rj*irijm &
!               + (irijm*irijm2/10)*(4*ri_dot_rj*(rim2+rjm2) &
!                   - rim2*rjm2 - 7*ri_dot_rj**2)
!           iv = iv + 2*tr !Symmetric

!       end do
!   end do

!   dkirk = dkirk + 3*num_atoms !Diagonal components
!   dkirk = dkirk*rnb*rnb/3
!   rh = 1.0_rp/dkirk

!   radsq_tot = rgsq*num_atoms
!   iv = (10.0_rp/3.0_rp)*num_atoms + radsq_tot/(1.0_rp + 0.75_rp*iv/radsq_tot)
    

    end subroutine

!******************************************************************************

subroutine calc_shape(ev1, ev2, ev3, asph, prol)
    !! Given three eigen values of the gyration tensor, calculates asphericity
    !! and prolateness. Note that ev1 >= ev2 >= ev3.

    real(rp), intent(in)  :: ev1
    real(rp), intent(in)  :: ev2
    real(rp), intent(in)  :: ev3
    real(rp), intent(out) :: asph
    real(rp), intent(out) :: prol
    real(rp) :: rgsq
    real(rp) :: ev_av, ev1mav, ev2mav, ev3mav, evmav_sumsq

    rgsq = ev1 + ev2 + ev3
    ev_av = rgsq*math_third
    ev1mav = ev1 - ev_av
    ev2mav = ev2 - ev_av
    ev3mav = ev3 - ev_av
    evmav_sumsq = ev1mav**2 + ev2mav**2 + ev3mav**2

    asph = 1.5_rp*evmav_sumsq/(rgsq*rgsq)
    prol = 4*(ev1mav*ev2mav*ev3mav)/((2*math_third)*evmav_sumsq)**1.5

    end subroutine

!******************************************************************************

subroutine dsyevc3(a, w)

    !! Calculates the eigenvalues of a symmetric 3x3 matrix A using Cardano's
    !! analytical algorithm.
    !! Only the diagonal and upper triangular parts of A are accessed. The access
    !! is read-only.
    !!
    !! Copyright (C) 2006  Joachim Kopp
    ! ----------------------------------------------------------------------------
    ! Parameters:
    !   A: The symmetric input matrix
    !   W: Storage buffer for eigenvalues
    ! ----------------------------------------------------------------------------
    ! .. Arguments ..
      REAL(RP), DIMENSION(3,3), INTENT(IN) :: A
      REAL(RP), DIMENSION(3), INTENT(OUT) ::  W(3)

     !.. Local Variables ..
      REAL(RP) ::  M, C1, C0
      REAL(RP) ::  DE, DD, EE, FF
      REAL(RP) ::  P, SQRTP, Q, C, S, PHI
  
     !Determine coefficients of characteristic poynomial. We write
     !      | A   D   F  |
     ! A =  | D*  B   E  |
     !      | F*  E*  C  |

      DE    = A(1,2) * A(2,3)
      DD    = A(1,2)**2
      EE    = A(2,3)**2
      FF    = A(1,3)**2
      M     = A(1,1) + A(2,2) + A(3,3)
      C1    = ( A(1,1)*A(2,2) + A(1,1)*A(3,3) + A(2,2)*A(3,3) ) &
               - (DD + EE + FF)
      C0    = A(3,3)*DD + A(1,1)*EE + A(2,2)*FF - A(1,1)*A(2,2)*A(3,3) &
               - 2.0_RP * A(1,3)*DE

      P     = M**2 - 3.0_RP * C1
      Q     = M*(P - (3.0_RP/2.0_RP)*C1) - (27.0_RP/2.0_RP)*C0
      SQRTP = SQRT(ABS(P))

      PHI   = 27.0_RP * ( 0.25_RP * C1**2 * (P - C1) &
                + C0 * (Q + (27.0_RP/4.0_RP)*C0) )
      PHI   = (1.0_RP/3.0_RP) * ATAN2(SQRT(ABS(PHI)), Q)

      C     = SQRTP * COS(PHI)
      S     = (1.0_RP/MATH_SQRT3) * SQRTP * SIN(PHI)

      W(2) = (1.0_RP/3.0_RP) * (M - C)
      W(3) = W(2) + S
      W(1) = W(2) + C
      W(2) = W(2) - S

      END SUBROUTINE

!******************************************************************************

end module m_stats_io
