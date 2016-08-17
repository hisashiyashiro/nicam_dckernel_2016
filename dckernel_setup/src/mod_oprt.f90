!-------------------------------------------------------------------------------
!>
!! Operator module
!!
!! @par Description
!!         This module contains the subroutines for differential oeprators.
!!
!! @author  H.Tomita
!!
!! @par History
!! @li      2004-02-17 (H.Tomita)    Imported from igdc-4.33
!! @li      2006-04-17 (H.Tomita)    Add sub[OPRT_divergence]
!! @li      2006-08-11 (H.Tomita)    Implementatio of miura scheme(2004) with Thurbuen(1996)'s limiter.
!!                                   Add sub[OPRT_divergence2]
!!                                       sub[OPRT_divergence2_prep]
!!                                       sub[OPRT_divergence2_all].
!! @li      2006-08-22 (Y.Niwa)      divide the rows for the calc. of clap and clap_pl due to the rule of XLF.
!! @li      2007-01-26 (H.Tomita)    Optimization of sub[oprt_diffusion].
!! @li      2007-11-28 (T.Mitsui)    bugfix in oprt_divergence2, _all
!! @li      2008-01-24 (Y.Niwa)      add OPRT_divergence2{_prep,,_all}_rev
!! @li      2008-04-28 (T.Mitsui)    bug fix in OPRT_divergence2{_all}_rev
!! @li      2009-09-04 (H.Taniguchi) bug fix in OPRT_divergence2{_all}_rev Zero clear of wrk[_pl] is needed.
!! @li      2010-06-08 (S.Iga)       new grid is implemented (see, string XTMS)
!! @li      2011-09-27 (T.Seiki)     merge optimization by RIST and M.Terai
!! @li      2012-01-17 (M.Terai)     update optimization(case6) in div2rev
!! @li      2012-05-01 (T.Yamaura)   bug fix in div2rev
!! @li      2012-06-28 (M.Terai)     Removed wrapper subroutine to invoked the operators directly macro var.
!!
!<
module mod_oprt
  !-----------------------------------------------------------------------------
  !
  !++ used modules
  !
  use mod_precision
  use mod_misc
!ESC!  use mod_stdio
!ESC!  use mod_prof
!ESC!  use mod_adm, only: &
!ESC!     ADM_nxyz,           &
!ESC!     TI    => ADM_TI,    &
!ESC!     TJ    => ADM_TJ,    &
!ESC!     AI    => ADM_AI,    &
!ESC!     AIJ   => ADM_AIJ,   &
!ESC!     AJ    => ADM_AJ,    &
!ESC!     K0    => ADM_KNONE, &
!ESC!     vlink => ADM_vlink, &
!ESC!     ADM_lall,           &
!ESC!     ADM_lall_pl,        &
!ESC!     ADM_kall,           &
!ESC!     ADM_jall,           &
!ESC!     ADM_iall,           &
!ESC!     ADM_gall_pl
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  public :: OPRT_divergence_setup
  public :: OPRT_rotation_setup
  public :: OPRT_gradient_setup
  public :: OPRT_laplacian_setup
  public :: OPRT_diffusion_setup

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine OPRT_divergence_setup( &
       GMTR_p,   GMTR_p_pl,  &
       GMTR_t,   GMTR_t_pl,  &
       GMTR_a,   GMTR_a_pl,  &
       coef_div, coef_div_pl )
!ESC!    use mod_adm, only: &
!ESC!       ADM_have_pl,  &
!ESC!       ADM_have_sgp, &
!ESC!       ADM_jmin,     &
!ESC!       ADM_jmax,     &
!ESC!       ADM_imin,     &
!ESC!       ADM_imax,     &
!ESC!       ADM_gslf_pl,  &
!ESC!       ADM_gmin_pl,  &
!ESC!       ADM_gmax_pl
!ESC!    use mod_gmtr, only: &
!ESC!       P_RAREA => GMTR_p_RAREA, &
!ESC!       W1      => GMTR_t_W1,    &
!ESC!       W2      => GMTR_t_W2,    &
!ESC!       W3      => GMTR_t_W3,    &
!ESC!       HNX     => GMTR_a_HNX,   &
!ESC!       GMTR_p_nmax,             &
!ESC!       GMTR_t_nmax,             &
!ESC!       GMTR_a_nmax,             &
!ESC!       GMTR_a_nmax_pl
    implicit none

    real(RP), intent(in)  :: GMTR_p     (ADM_iall,ADM_jall,K0,ADM_lall   ,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_p_pl  (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_t     (ADM_iall,ADM_jall,K0,ADM_lall   ,TI:TJ,GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_t_pl  (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_a     (ADM_iall,ADM_jall,K0,ADM_lall   ,AI:AJ,GMTR_a_nmax   )
    real(RP), intent(in)  :: GMTR_a_pl  (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_a_nmax_pl)
!    real(RP), intent(out) :: coef_div   (ADM_nxyz,ADM_gall,0:6        ,ADM_lall   )
!    real(RP), intent(out) :: coef_div_pl(ADM_nxyz,         0:ADM_vlink,ADM_lall_pl)
    real(RP), intent(out) :: coef_div   (ADM_iall,ADM_jall,0:6        ,ADM_nxyz,ADM_lall   )
    real(RP), intent(out) :: coef_div_pl(                  0:ADM_vlink,ADM_nxyz,ADM_lall_pl)

    integer  :: ij, ijp1, ijm1

    real(RP) :: coef
    integer  :: i, j, l, d, n, v, hn
    !---------------------------------------------------------------------------

    write(ADM_LOG_FID,*) '*** setup coefficient of divergence operator'

    coef_div   (:,:,:,:,:) = 0.0_RP
    coef_div_pl(:,:,:)     = 0.0_RP

    do l = 1, ADM_lall
    do d = 1, ADM_nxyz
       hn = d + HNX - 1

       do j = ADM_jmin, ADM_jmax
       do i = ADM_imin, ADM_imax

          ! ij
!          coef_div(d,ij,0,l) = &
          coef_div(i,j,0,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q1 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q1 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q2 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q2 * b2
                                  + GMTR_t(i-1,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q3 * b2
                                  - GMTR_t(i-1,j  ,k0,l,TI,W2) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q3 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q4 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q4 * b4
                                  - GMTR_t(i-1,j-1,k0,l,TI,W3) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q5 * b4
                                  - GMTR_t(i-1,j-1,k0,l,TI,W3) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) & ! Q5 * b5
                                  - GMTR_t(i  ,j-1,k0,l,TJ,W3) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) & ! Q6 * b5
                                  + GMTR_t(i  ,j-1,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q6 * b6
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ip1j
!          coef_div(d,ij,1,l) = &
          coef_div(i,j,1,d,l) = ( - GMTR_t(i  ,j-1,k0,l,TJ,W2) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) & ! Q6 * b5
                                  + GMTR_t(i  ,j-1,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q6 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q1 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q1 * b1
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ip1jp1
!          coef_div(d,ij,2,l) = &
          coef_div(i,j,2,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q1 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q1 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q2 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q2 * b2
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ijp1
!          coef_div(d,ij,3,l) = &
          coef_div(i,j,3,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q2 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q2 * b2
                                  + GMTR_t(i-1,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q3 * b2
                                  - GMTR_t(i-1,j  ,k0,l,TI,W3) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q3 * b3
                                ) * 0.5_RP*GMTR_p(i,j,k0,l,P_RAREA)
          ! im1j
!          coef_div(d,ij,4,l) = &
          coef_div(i,j,4,d,l) = ( + GMTR_t(i-1,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q3 * b2
                                  - GMTR_t(i-1,j  ,k0,l,TI,W1) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q3 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q4 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q4 * b4
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! im1jm1
!          coef_div(d,ij,5,l) = &
          coef_div(i,j,5,d,l) = ( - GMTR_t(i-1,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q4 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q4 * b4
                                  - GMTR_t(i-1,j-1,k0,l,TI,W1) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q5 * b4
                                  - GMTR_t(i-1,j-1,k0,l,TI,W1) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) & ! Q5 * b5
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ijm1
!          coef_div(d,ij,6,l) = &
          coef_div(i,j,6,d,l) = ( - GMTR_t(i-1,j-1,k0,l,TI,W2) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q5 * b4
                                  - GMTR_t(i-1,j-1,k0,l,TI,W2) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) & ! Q5 * b5
                                  - GMTR_t(i  ,j-1,k0,l,TJ,W1) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) & ! Q6 * b5
                                  + GMTR_t(i  ,j-1,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q6 * b6
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
       enddo
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          j = ADM_jmin
          i = ADM_imin

          ! ij
!          coef_div(d,ij,0,l) = &
          coef_div(i,j,0,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q1 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q1 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q2 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q2 * b2
                                  + GMTR_t(i-1,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q3 * b2
                                  - GMTR_t(i-1,j  ,k0,l,TI,W2) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q3 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q4 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q4 * b4
                                  - GMTR_t(i  ,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q6 * b4
                                  + GMTR_t(i  ,j-1,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q6 * b6
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ip1j
!          coef_div(d,ij,1,l) = &
          coef_div(i,j,1,d,l) = ( - GMTR_t(i  ,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q6 * b4
                                  + GMTR_t(i  ,j-1,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q6 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q1 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q1 * b1
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ip1jp1
!          coef_div(d,ij,2,l) = &
          coef_div(i,j,2,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q1 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q1 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q2 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q2 * b2
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ijp1
!          coef_div(d,ij,3,l) = &
          coef_div(i,j,3,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q2 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q2 * b2
                                  + GMTR_t(i-1,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q3 * b2
                                  - GMTR_t(i-1,j  ,k0,l,TI,W3) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q3 * b3
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! im1j
!          coef_div(d,ij,4,l) = &
          coef_div(i,j,4,d,l) = ( + GMTR_t(i-1,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q3 * b2
                                  - GMTR_t(i-1,j  ,k0,l,TI,W1) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q3 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q4 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q4 * b4
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! im1jm1
!          coef_div(d,ij,5,l) = &
          coef_div(i,j,5,d,l) = ( - GMTR_t(i-1,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q4 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q4 * b4
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ijm1
!          coef_div(d,ij,6,l) = &
          coef_div(i,j,6,d,l) = ( - GMTR_t(i  ,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q6 * b4
                                  + GMTR_t(i  ,j-1,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q6 * b6
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
       endif

    enddo ! loop d
    enddo ! loop l

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl
       do l = 1, ADM_lall_pl
       do d = 1, ADM_nxyz
          hn = d + HNX - 1

          coef = 0.0_RP
          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             if( ijp1 == ADM_gmax_pl+1 ) ijp1 = ADM_gmin_pl

             coef = coef + ( GMTR_t_pl(ij,k0,l,W1) * GMTR_a_pl(ij  ,k0,l,hn) &
                           + GMTR_t_pl(ij,k0,l,W1) * GMTR_a_pl(ijp1,k0,l,hn) )
          enddo
!          coef_div_pl(d,0,l) &
          coef_div_pl(0,d,l) = coef * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             ijm1 = v - 1
             if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl
             if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

!          coef_div_pl(d,v-1,l) &
             coef_div_pl(v-1,d,l) = ( + GMTR_t_pl(ijm1,k0,l,W3) * GMTR_a_pl(ijm1,k0,l,hn) &
                                      + GMTR_t_pl(ijm1,k0,l,W3) * GMTR_a_pl(ij  ,k0,l,hn) &
                                      + GMTR_t_pl(ij  ,k0,l,W2) * GMTR_a_pl(ij  ,k0,l,hn) &
                                      + GMTR_t_pl(ij  ,k0,l,W2) * GMTR_a_pl(ijp1,k0,l,hn) &
                                    ) * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)
          enddo
       enddo ! loop d
       enddo ! loop l
    endif

    return
  end subroutine OPRT_divergence_setup

  !-----------------------------------------------------------------------------
  subroutine OPRT_rotation_setup( &
       GMTR_p,   GMTR_p_pl,  &
       GMTR_t,   GMTR_t_pl,  &
       GMTR_a,   GMTR_a_pl,  &
       coef_rot, coef_rot_pl )
!ESC!    use mod_adm, only: &
!ESC!       ADM_have_pl,  &
!ESC!       ADM_have_sgp, &
!ESC!       ADM_jmin,     &
!ESC!       ADM_jmax,     &
!ESC!       ADM_imin,     &
!ESC!       ADM_imax,     &
!ESC!       ADM_gslf_pl,  &
!ESC!       ADM_gmin_pl,  &
!ESC!       ADM_gmax_pl
!ESC!    use mod_gmtr, only: &
!ESC!       P_RAREA => GMTR_p_RAREA, &
!ESC!       W1      => GMTR_t_W1,    &
!ESC!       W2      => GMTR_t_W2,    &
!ESC!       W3      => GMTR_t_W3,    &
!ESC!       HTX     => GMTR_a_HTX,   &
!ESC!       GMTR_p_nmax,             &
!ESC!       GMTR_t_nmax,             &
!ESC!       GMTR_a_nmax,             &
!ESC!       GMTR_a_nmax_pl
    implicit none

    real(RP), intent(in)  :: GMTR_p     (ADM_iall,ADM_jall,K0,ADM_lall   ,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_p_pl  (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_t     (ADM_iall,ADM_jall,K0,ADM_lall   ,TI:TJ,GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_t_pl  (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_a     (ADM_iall,ADM_jall,K0,ADM_lall   ,AI:AJ,GMTR_a_nmax   )
    real(RP), intent(in)  :: GMTR_a_pl  (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_a_nmax_pl)
!    real(RP), intent(out) :: coef_rot   (ADM_nxyz,ADM_gall,0:6        ,ADM_lall   )
!    real(RP), intent(out) :: coef_rot_pl(ADM_nxyz,         0:ADM_vlink,ADM_lall_pl)
    real(RP), intent(out) :: coef_rot   (ADM_iall,ADM_jall,0:6        ,ADM_nxyz,ADM_lall   )
    real(RP), intent(out) :: coef_rot_pl(                  0:ADM_vlink,ADM_nxyz,ADM_lall_pl)

    integer  :: ij, ijp1, ijm1

    real(RP) :: coef
    integer  :: i, j, l, d, n, v, ht
    !---------------------------------------------------------------------------

    write(ADM_LOG_FID,*) '*** setup coefficient of rotation operator'

    coef_rot    (:,:,:,:,:) = 0.0_RP
    coef_rot_pl (:,:,:)     = 0.0_RP

    do l = 1, ADM_lall
    do d = 1, ADM_nxyz
       ht = d + HTX - 1

       do j = ADM_jmin, ADM_jmax
       do i = ADM_imin, ADM_imax

          ! ij
!          coef_rot(d,ij,0,l) &
          coef_rot(i,j,0,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AI ,ht) & ! Q1 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AIJ,ht) & ! Q1 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AIJ,ht) & ! Q2 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AJ ,ht) & ! Q2 * b2
                                  + GMTR_t(i-1,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AJ ,ht) & ! Q3 * b2
                                  - GMTR_t(i-1,j  ,k0,l,TI,W2) * GMTR_a(i-1,j  ,k0,l,AI ,ht) & ! Q3 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j  ,k0,l,AI ,ht) & ! Q4 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j-1,k0,l,AIJ,ht) & ! Q4 * b4
                                  - GMTR_t(i-1,j-1,k0,l,TI,W3) * GMTR_a(i-1,j-1,k0,l,AIJ,ht) & ! Q5 * b4
                                  - GMTR_t(i-1,j-1,k0,l,TI,W3) * GMTR_a(i  ,j-1,k0,l,AJ ,ht) & ! Q5 * b5
                                  - GMTR_t(i  ,j-1,k0,l,TJ,W3) * GMTR_a(i  ,j-1,k0,l,AJ ,ht) & ! Q6 * b5
                                  + GMTR_t(i  ,j-1,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AI ,ht) & ! Q6 * b6
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ip1j
!          coef_rot(d,ij,1,l) &
          coef_rot(i,j,1,d,l) = ( - GMTR_t(i  ,j-1,k0,l,TJ,W2) * GMTR_a(i  ,j-1,k0,l,AJ ,ht) & ! Q6 * b5
                                  + GMTR_t(i  ,j-1,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AI ,ht) & ! Q6 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AI ,ht) & ! Q1 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AIJ,ht) & ! Q1 * b1
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ip1jp1
!          coef_rot(d,ij,2,l) &
          coef_rot(i,j,2,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AI ,ht) & ! Q1 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AIJ,ht) & ! Q1 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AIJ,ht) & ! Q2 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AJ ,ht) & ! Q2 * b2
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ijp1
!          coef_rot(d,ij,3,l) &
          coef_rot(i,j,3,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AIJ,ht) & ! Q2 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AJ ,ht) & ! Q2 * b2
                                  + GMTR_t(i-1,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AJ ,ht) & ! Q3 * b2
                                  - GMTR_t(i-1,j  ,k0,l,TI,W3) * GMTR_a(i-1,j  ,k0,l,AI ,ht) & ! Q3 * b3
                                ) * 0.5_RP*GMTR_p(i,j,k0,l,P_RAREA)
          ! im1j
!          coef_rot(d,ij,4,l) &
          coef_rot(i,j,4,d,l) = ( + GMTR_t(i-1,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AJ ,ht) & ! Q3 * b2
                                  - GMTR_t(i-1,j  ,k0,l,TI,W1) * GMTR_a(i-1,j  ,k0,l,AI ,ht) & ! Q3 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j  ,k0,l,AI ,ht) & ! Q4 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j-1,k0,l,AIJ,ht) & ! Q4 * b4
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! im1jm1
!          coef_rot(d,ij,5,l) &
          coef_rot(i,j,5,d,l) = ( - GMTR_t(i-1,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j  ,k0,l,AI ,ht) & ! Q4 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j-1,k0,l,AIJ,ht) & ! Q4 * b4
                                  - GMTR_t(i-1,j-1,k0,l,TI,W1) * GMTR_a(i-1,j-1,k0,l,AIJ,ht) & ! Q5 * b4
                                  - GMTR_t(i-1,j-1,k0,l,TI,W1) * GMTR_a(i  ,j-1,k0,l,AJ ,ht) & ! Q5 * b5
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ijm1
!          coef_rot(d,ij,6,l) &
          coef_rot(i,j,6,d,l) = ( - GMTR_t(i-1,j-1,k0,l,TI,W2) * GMTR_a(i-1,j-1,k0,l,AIJ,ht) & ! Q5 * b4
                                  - GMTR_t(i-1,j-1,k0,l,TI,W2) * GMTR_a(i  ,j-1,k0,l,AJ ,ht) & ! Q5 * b5
                                  - GMTR_t(i  ,j-1,k0,l,TJ,W1) * GMTR_a(i  ,j-1,k0,l,AJ ,ht) & ! Q6 * b5
                                  + GMTR_t(i  ,j-1,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AI ,ht) & ! Q6 * b6
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
       enddo
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          j = ADM_jmin
          i = ADM_imin

          ! ij
!          coef_rot(d,ij,0,l) &
          coef_rot(i,j,0,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AI ,ht) & ! Q1 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AIJ,ht) & ! Q1 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AIJ,ht) & ! Q2 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AJ ,ht) & ! Q2 * b2
                                  + GMTR_t(i-1,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AJ ,ht) & ! Q3 * b2
                                  - GMTR_t(i-1,j  ,k0,l,TI,W2) * GMTR_a(i-1,j  ,k0,l,AI ,ht) & ! Q3 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j  ,k0,l,AI ,ht) & ! Q4 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j-1,k0,l,AIJ,ht) & ! Q4 * b4
                                  - GMTR_t(i  ,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j-1,k0,l,AIJ,ht) & ! Q6 * b4
                                  + GMTR_t(i  ,j-1,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AI ,ht) & ! Q6 * b6
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ip1j
!          coef_rot(d,ij,1,l) &
          coef_rot(i,j,1,d,l) = ( - GMTR_t(i  ,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j-1,k0,l,AIJ,ht) & ! Q6 * b4
                                  + GMTR_t(i  ,j-1,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AI ,ht) & ! Q6 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AI ,ht) & ! Q1 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AIJ,ht) & ! Q1 * b1
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ip1jp1
!          coef_rot(d,ij,2,l) &
          coef_rot(i,j,2,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AI ,ht) & ! Q1 * b6
                                  + GMTR_t(i  ,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AIJ,ht) & ! Q1 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AIJ,ht) & ! Q2 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AJ ,ht) & ! Q2 * b2
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ijp1
!          coef_rot(d,ij,3,l) &
          coef_rot(i,j,3,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AIJ,ht) & ! Q2 * b1
                                  + GMTR_t(i  ,j  ,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AJ ,ht) & ! Q2 * b2
                                  + GMTR_t(i-1,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AJ ,ht) & ! Q3 * b2
                                  - GMTR_t(i-1,j  ,k0,l,TI,W3) * GMTR_a(i-1,j  ,k0,l,AI ,ht) & ! Q3 * b3
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! im1j
!          coef_rot(d,ij,4,l) &
          coef_rot(i,j,4,d,l) = ( + GMTR_t(i-1,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AJ ,ht) & ! Q3 * b2
                                  - GMTR_t(i-1,j  ,k0,l,TI,W1) * GMTR_a(i-1,j  ,k0,l,AI ,ht) & ! Q3 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j  ,k0,l,AI ,ht) & ! Q4 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j-1,k0,l,AIJ,ht) & ! Q4 * b4
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! im1jm1
!          coef_rot(d,ij,5,l) &
          coef_rot(i,j,5,d,l) = ( - GMTR_t(i-1,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j  ,k0,l,AI ,ht) & ! Q4 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j-1,k0,l,AIJ,ht) & ! Q4 * b4
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ijm1
!          coef_rot(d,ij,6,l) &
          coef_rot(i,j,6,d,l) = ( - GMTR_t(i  ,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j-1,k0,l,AIJ,ht) & ! Q6 * b4
                                  + GMTR_t(i  ,j-1,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AI ,ht) & ! Q6 * b6
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
       endif

    enddo ! loop d
    enddo ! loop l

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl
       do l = 1, ADM_lall_pl
       do d = 1, ADM_nxyz
          ht = d + HTX - 1

          coef = 0.0_RP
          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             if( ijp1 == ADM_gmax_pl+1 ) ijp1 = ADM_gmin_pl

             coef = coef + ( GMTR_t_pl(ij,k0,l,W1) * GMTR_a_pl(ij  ,k0,l,ht) &
                           + GMTR_t_pl(ij,k0,l,W1) * GMTR_a_pl(ijp1,k0,l,ht) )
          enddo
!          coef_rot_pl(d,0,l) &
          coef_rot_pl(0,d,l) = coef * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             ijm1 = v - 1
             if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl
             if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

!             coef_rot_pl(d,v-1,l) &
             coef_rot_pl(v-1,d,l) = ( + GMTR_t_pl(ijm1,k0,l,W3) * GMTR_a_pl(ijm1,k0,l,ht) &
                                      + GMTR_t_pl(ijm1,k0,l,W3) * GMTR_a_pl(ij  ,k0,l,ht) &
                                      + GMTR_t_pl(ij  ,k0,l,W2) * GMTR_a_pl(ij  ,k0,l,ht) &
                                      + GMTR_t_pl(ij  ,k0,l,W2) * GMTR_a_pl(ijp1,k0,l,ht) &
                                    ) * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)
          enddo
       enddo ! loop d
       enddo ! loop l
    endif

    return
  end subroutine OPRT_rotation_setup

  !-----------------------------------------------------------------------------
  subroutine OPRT_gradient_setup( &
       GMTR_p,    GMTR_p_pl,   &
       GMTR_t,    GMTR_t_pl,   &
       GMTR_a,    GMTR_a_pl,   &
       coef_grad, coef_grad_pl )
!ESC!    use mod_adm, only: &
!ESC!       ADM_have_pl,  &
!ESC!       ADM_have_sgp, &
!ESC!       ADM_jmin,     &
!ESC!       ADM_jmax,     &
!ESC!       ADM_imin,     &
!ESC!       ADM_imax,     &
!ESC!       ADM_gslf_pl,  &
!ESC!       ADM_gmin_pl,  &
!ESC!       ADM_gmax_pl
!ESC!    use mod_gmtr, only: &
!ESC!       P_RAREA => GMTR_p_RAREA, &
!ESC!       W1      => GMTR_t_W1,    &
!ESC!       W2      => GMTR_t_W2,    &
!ESC!       W3      => GMTR_t_W3,    &
!ESC!       HNX     => GMTR_a_HNX,   &
!ESC!       GMTR_p_nmax,             &
!ESC!       GMTR_t_nmax,             &
!ESC!       GMTR_a_nmax,             &
!ESC!       GMTR_a_nmax_pl
    implicit none

    real(RP), intent(in)  :: GMTR_p      (ADM_iall,ADM_jall,K0,ADM_lall   ,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_p_pl   (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_t      (ADM_iall,ADM_jall,K0,ADM_lall   ,TI:TJ,GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_t_pl   (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_a      (ADM_iall,ADM_jall,K0,ADM_lall   ,AI:AJ,GMTR_a_nmax   )
    real(RP), intent(in)  :: GMTR_a_pl   (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_a_nmax_pl)
!    real(RP), intent(out) :: coef_grad   (ADM_nxyz,ADM_gall,0:6        ,ADM_lall   )
!    real(RP), intent(out) :: coef_grad_pl(ADM_nxyz,         0:ADM_vlink,ADM_lall_pl)
    real(RP), intent(out) :: coef_grad   (ADM_iall,ADM_jall,0:6        ,ADM_nxyz,ADM_lall   )
    real(RP), intent(out) :: coef_grad_pl(                  0:ADM_vlink,ADM_nxyz,ADM_lall_pl)

    integer  :: ij, ijp1, ijm1

    real(RP) :: coef
    integer  :: i, j, l, d, n, v, hn
    !---------------------------------------------------------------------------

    write(ADM_LOG_FID,*) '*** setup coefficient of gradient operator'

    coef_grad   (:,:,:,:,:) = 0.0_RP
    coef_grad_pl(:,:,:)     = 0.0_RP

    do l = 1, ADM_lall
    do d = 1, ADM_nxyz
       hn = d + HNX - 1

       do j = ADM_jmin, ADM_jmax
       do i = ADM_imin, ADM_imax

          ! ij
!          coef_grad(d,ij,0,l) &
          coef_grad(i,j,0,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q1 * b6
                                   + GMTR_t(i  ,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q1 * b1
                                   + GMTR_t(i  ,j  ,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q2 * b1
                                   + GMTR_t(i  ,j  ,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q2 * b2
                                   + GMTR_t(i-1,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q3 * b2
                                   - GMTR_t(i-1,j  ,k0,l,TI,W2) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q3 * b3
                                   - GMTR_t(i-1,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q4 * b3
                                   - GMTR_t(i-1,j-1,k0,l,TI,W3) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q4 * b4
                                   - GMTR_t(i-1,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q5 * b4
                                   - GMTR_t(i-1,j-1,k0,l,TI,W3) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) & ! Q5 * b5
                                   - GMTR_t(i  ,j-1,k0,l,TJ,W3) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) & ! Q6 * b5
                                   + GMTR_t(i  ,j-1,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q6 * b6
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,hn)                     & ! P0 * b1
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,hn)                     & ! P0 * b2
                                   + 2.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,hn)                     & ! P0 * b3
                                   + 2.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,hn)                     & ! P0 * b4
                                   + 2.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,hn)                     & ! P0 * b5
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,hn)                     & ! P0 * b6
                                 ) * 0.5_RP * GMTR_p(i  ,j  ,k0,l,P_RAREA)
          ! ip1j
!          coef_grad(d,ij,1,l) &
          coef_grad(i,j,1,d,l) = ( - GMTR_t(i  ,j-1,k0,l,TJ,W2) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) & ! Q6 * b5
                                   + GMTR_t(i  ,j-1,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q6 * b6
                                   + GMTR_t(i  ,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q1 * b6
                                   + GMTR_t(i  ,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q1 * b1
                                 ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ip1jp1
!          coef_grad(d,ij,2,l) &
          coef_grad(i,j,2,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q1 * b6
                                   + GMTR_t(i  ,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q1 * b1
                                   + GMTR_t(i  ,j  ,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q2 * b1
                                   + GMTR_t(i  ,j  ,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q2 * b2
                                 ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ijp1
!          coef_grad(d,ij,3,l) &
          coef_grad(i,j,3,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q2 * b1
                                   + GMTR_t(i  ,j  ,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q2 * b2
                                   + GMTR_t(i-1,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q3 * b2
                                   - GMTR_t(i-1,j  ,k0,l,TI,W3) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q3 * b3
                                 ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! im1j
!          coef_grad(d,ij,4,l) &
          coef_grad(i,j,4,d,l) = ( + GMTR_t(i-1,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q3 * b2
                                   - GMTR_t(i-1,j  ,k0,l,TI,W1) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q3 * b3
                                   - GMTR_t(i-1,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q4 * b3
                                   - GMTR_t(i-1,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q4 * b4
                                 ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! im1jm1
!          coef_grad(d,ij,5,l) &
          coef_grad(i,j,5,d,l) = ( - GMTR_t(i-1,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q4 * b3
                                  - GMTR_t(i-1,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q4 * b4
                                  - GMTR_t(i-1,j-1,k0,l,TI,W1) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q5 * b4
                                  - GMTR_t(i-1,j-1,k0,l,TI,W1) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) & ! Q5 * b5
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ijm1
!          coef_grad(d,ij,6,l) &
          coef_grad(i,j,6,d,l) = ( - GMTR_t(i-1,j-1,k0,l,TI,W2) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q5 * b4
                                  - GMTR_t(i-1,j-1,k0,l,TI,W2) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) & ! Q5 * b5
                                  - GMTR_t(i  ,j-1,k0,l,TJ,W1) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) & ! Q6 * b5
                                  + GMTR_t(i  ,j-1,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q6 * b6
                                ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
       enddo
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          j = ADM_jmin
          i = ADM_imin

          ! ij
!          coef_grad(d,ij,0,l) &
          coef_grad(i,j,0,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q1 * b6
                                   + GMTR_t(i  ,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q1 * b1
                                   + GMTR_t(i  ,j  ,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q2 * b1
                                   + GMTR_t(i  ,j  ,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q2 * b2
                                   + GMTR_t(i-1,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q3 * b2
                                   - GMTR_t(i-1,j  ,k0,l,TI,W2) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q3 * b3
                                   - GMTR_t(i-1,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q4 * b3
                                   - GMTR_t(i-1,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q4 * b4
                                   - GMTR_t(i  ,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q6 * b4
                                   + GMTR_t(i  ,j-1,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q6 * b6
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,hn)                     & ! P0 * b1
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,hn)                     & ! P0 * b2
                                   + 2.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,hn)                     & ! P0 * b3
                                   + 2.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,hn)                     & ! P0 * b4
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,hn)                     & ! P0 * b6
                                 ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ip1j
!          coef_grad(d,ij,1,l) &
          coef_grad(i,j,1,d,l) = ( - GMTR_t(i  ,j-1,k0,l,TJ,W2) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q6 * b4
                                   + GMTR_t(i  ,j-1,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q6 * b6
                                   + GMTR_t(i  ,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q1 * b6
                                   + GMTR_t(i  ,j  ,k0,l,TI,W2) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q1 * b1
                                 ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ip1jp1
!          coef_grad(d,ij,2,l) &
          coef_grad(i,j,2,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q1 * b6
                                   + GMTR_t(i  ,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q1 * b1
                                   + GMTR_t(i  ,j  ,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q2 * b1
                                   + GMTR_t(i  ,j  ,k0,l,TJ,W2) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q2 * b2
                                 ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ijp1
!          coef_grad(d,ij,3,l) &
          coef_grad(i,j,3,d,l) = ( + GMTR_t(i  ,j  ,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) & ! Q2 * b1
                                   + GMTR_t(i  ,j  ,k0,l,TJ,W3) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q2 * b2
                                   + GMTR_t(i-1,j  ,k0,l,TI,W3) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q3 * b2
                                   - GMTR_t(i-1,j  ,k0,l,TI,W3) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q3 * b3
                                 ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! im1j
!          coef_grad(d,ij,4,l) &
          coef_grad(i,j,4,d,l) = ( + GMTR_t(i-1,j  ,k0,l,TI,W1) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) & ! Q3 * b2
                                   - GMTR_t(i-1,j  ,k0,l,TI,W1) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q3 * b3
                                   - GMTR_t(i-1,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q4 * b3
                                   - GMTR_t(i-1,j-1,k0,l,TJ,W3) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q4 * b4
                                 ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! im1jm1
!          coef_grad(d,ij,5,l) &
          coef_grad(i,j,5,d,l) = ( - GMTR_t(i-1,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j  ,k0,l,AI ,hn) & ! Q4 * b3
                                   - GMTR_t(i-1,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q4 * b4
                                 ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          ! ijm1
!          coef_grad(d,ij,6,l) &
          coef_grad(i,j,6,d,l) = ( - GMTR_t(i  ,j-1,k0,l,TJ,W1) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) & ! Q6 * b4
                                   + GMTR_t(i  ,j-1,k0,l,TJ,W1) * GMTR_a(i  ,j  ,k0,l,AI ,hn) & ! Q6 * b6
                                 ) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
       endif
    enddo
    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do d = 1, ADM_nxyz
          hn = d + HNX - 1

          coef = 0.0_RP
          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl

             coef = coef + 2.0_RP * ( GMTR_t_pl(ij,k0,l,W1) - 1.0_RP ) * GMTR_a_pl(ijp1,k0,l,hn)
          enddo
!          coef_grad_pl(d,0,l) &
          coef_grad_pl(0,d,l) = coef * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             ijm1 = v - 1
             if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl
             if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

!             coef_grad_pl(d,v-1,l) &
             coef_grad_pl(v-1,d,l) = ( + GMTR_t_pl(ijm1,k0,l,W3) * GMTR_a_pl(ijm1,k0,l,hn) &
                                       + GMTR_t_pl(ijm1,k0,l,W3) * GMTR_a_pl(ij  ,k0,l,hn) &
                                       + GMTR_t_pl(ij  ,k0,l,W2) * GMTR_a_pl(ij  ,k0,l,hn) &
                                       + GMTR_t_pl(ij  ,k0,l,W2) * GMTR_a_pl(ijp1,k0,l,hn) &
                                     ) * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)
          enddo
       enddo ! loop d
       enddo ! loop l
    endif

    return
  end subroutine OPRT_gradient_setup

  !-----------------------------------------------------------------------------
  subroutine OPRT_laplacian_setup( &
       GMTR_p,   GMTR_p_pl,  &
       GMTR_t,   GMTR_t_pl,  &
       GMTR_a,   GMTR_a_pl,  &
       coef_lap, coef_lap_pl )
!ESC!    use mod_adm, only: &
!ESC!       ADM_have_pl,  &
!ESC!       ADM_have_sgp, &
!ESC!       ADM_jmin,     &
!ESC!       ADM_jmax,     &
!ESC!       ADM_imin,     &
!ESC!       ADM_imax,     &
!ESC!       ADM_gslf_pl,  &
!ESC!       ADM_gmin_pl,  &
!ESC!       ADM_gmax_pl
!ESC!    use mod_gmtr, only: &
!ESC!       P_RAREA => GMTR_p_RAREA, &
!ESC!       T_RAREA => GMTR_t_RAREA, &
!ESC!       HNX     => GMTR_a_HNX,   &
!ESC!       TNX     => GMTR_a_TNX,   &
!ESC!       TN2X    => GMTR_a_TN2X,  &
!ESC!       GMTR_p_nmax,             &
!ESC!       GMTR_t_nmax,             &
!ESC!       GMTR_a_nmax,             &
!ESC!       GMTR_a_nmax_pl
    implicit none

    real(RP), intent(in)  :: GMTR_p     (ADM_iall,ADM_jall,K0,ADM_lall   ,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_p_pl  (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_t     (ADM_iall,ADM_jall,K0,ADM_lall   ,TI:TJ,GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_t_pl  (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_a     (ADM_iall,ADM_jall,K0,ADM_lall   ,AI:AJ,GMTR_a_nmax   )
    real(RP), intent(in)  :: GMTR_a_pl  (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_a_nmax_pl)
    real(RP), intent(out) :: coef_lap   (ADM_iall,ADM_jall,0:6        ,ADM_lall   )
    real(RP), intent(out) :: coef_lap_pl(                  0:ADM_vlink,ADM_lall_pl)

    integer  :: ij, ijp1, ijm1

    integer  :: i, j, l, d, n, v, hn, tn, tn2
    !---------------------------------------------------------------------------

    write(ADM_LOG_FID,*) '*** setup coefficient of laplacian operator'

    coef_lap    (:,:,:,:) = 0.0_RP
    coef_lap_pl (:,:)     = 0.0_RP

    do l = 1, ADM_lall

       do j = ADM_jmin, ADM_jmax
       do i = ADM_imin, ADM_imax

          coef_lap(i,j,:,l) = 0.0_RP

          do d = 1, ADM_nxyz
             hn = d + HNX - 1
             tn = d + TNX - 1

             ! ij
             coef_lap(i,j,0,l) = coef_lap(i,j,0,l) &
                               + GMTR_t(i  ,j  ,k0,l,TI,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 2.0_RP * GMTR_a(i+1,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 2.0_RP * GMTR_a(i+1,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) )

             coef_lap(i,j,0,l) = coef_lap(i,j,0,l) &
                               + GMTR_t(i  ,j  ,k0,l,TJ,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j+1,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j+1,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) )

             coef_lap(i,j,0,l) = coef_lap(i,j,0,l) &
                               + GMTR_t(i-1,j  ,k0,l,TI,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 2.0_RP * GMTR_a(i-1,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 2.0_RP * GMTR_a(i-1,j  ,k0,l,AIJ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) )

             coef_lap(i,j,0,l) = coef_lap(i,j,0,l) &
                               + GMTR_t(i-1,j-1,k0,l,TJ,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 2.0_RP * GMTR_a(i-1,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 2.0_RP * GMTR_a(i-1,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) )

             coef_lap(i,j,0,l) = coef_lap(i,j,0,l) &
                               + GMTR_t(i-1,j-1,k0,l,TI,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 2.0_RP * GMTR_a(i-1,j-1,k0,l,AI ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   - 2.0_RP * GMTR_a(i-1,j-1,k0,l,AI ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) )

             coef_lap(i,j,0,l) = coef_lap(i,j,0,l) &
                               + GMTR_t(i  ,j-1,k0,l,TJ,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j-1,k0,l,AIJ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j-1,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) )

             ! ip1j
             coef_lap(i,j,1,l) = coef_lap(i,j,1,l) &
                               + GMTR_t(i  ,j-1,k0,l,TJ,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j-1,k0,l,AIJ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j-1,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) )

             coef_lap(i,j,1,l) = coef_lap(i,j,1,l) &
                               + GMTR_t(i  ,j  ,k0,l,TI,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i+1,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i+1,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) )

             ! ip1jp1
             coef_lap(i,j,2,l) = coef_lap(i,j,2,l) &
                               + GMTR_t(i  ,j  ,k0,l,TI,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i+1,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i+1,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) )

             coef_lap(i,j,2,l) = coef_lap(i,j,2,l) &
                               + GMTR_t(i  ,j  ,k0,l,TJ,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i  ,j+1,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j+1,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) )

             ! ijp1
             coef_lap(i,j,3,l) = coef_lap(i,j,3,l) &
                               + GMTR_t(i  ,j  ,k0,l,TJ,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j+1,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j+1,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) )

             coef_lap(i,j,3,l) = coef_lap(i,j,3,l) &
                               + GMTR_t(i-1,j  ,k0,l,TI,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i-1,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   + 2.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AIJ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 2.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) )

             ! im1j
             coef_lap(i,j,4,l) = coef_lap(i,j,4,l) &
                               + GMTR_t(i-1,j  ,k0,l,TI,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AIJ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) )

             coef_lap(i,j,4,l) = coef_lap(i,j,4,l) &
                               + GMTR_t(i-1,j-1,k0,l,TJ,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i-1,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 2.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 2.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) )

             ! im1jm1
             coef_lap(i,j,5,l) = coef_lap(i,j,5,l) &
                               + GMTR_t(i-1,j-1,k0,l,TJ,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 2.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 2.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) )

             coef_lap(i,j,5,l) = coef_lap(i,j,5,l) &
                               + GMTR_t(i-1,j-1,k0,l,TI,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i-1,j-1,k0,l,AI ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j-1,k0,l,AI ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) )

             ! ijm1
             coef_lap(i,j,6,l) = coef_lap(i,j,6,l) &
                               + GMTR_t(i-1,j-1,k0,l,TI,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 2.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j-1,k0,l,AI ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   + 2.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j-1,k0,l,AI ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) )

             coef_lap(i,j,6,l) = coef_lap(i,j,6,l) &
                               + GMTR_t(i  ,j-1,k0,l,TJ,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i  ,j-1,k0,l,AIJ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i  ,j-1,k0,l,AJ ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j-1,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) )
          enddo
       enddo
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          j = ADM_jmin
          i = ADM_imin

          coef_lap(i,j,:,l) = 0.0_RP

          do d = 1, ADM_nxyz
             hn = d + HNX - 1
             tn = d + TNX - 1

             ! ij
             coef_lap(i,j,0,l) = coef_lap(i,j,0,l) &
                               + GMTR_t(i  ,j  ,k0,l,TI,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 2.0_RP * GMTR_a(i+1,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 2.0_RP * GMTR_a(i+1,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) )

             coef_lap(i,j,0,l) = coef_lap(i,j,0,l) &
                               + GMTR_t(i  ,j  ,k0,l,TJ,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j+1,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j+1,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) )

             coef_lap(i,j,0,l) = coef_lap(i,j,0,l) &
                               + GMTR_t(i-1,j  ,k0,l,TI,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 2.0_RP * GMTR_a(i-1,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 2.0_RP * GMTR_a(i-1,j  ,k0,l,AIJ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) )

             coef_lap(i,j,0,l) = coef_lap(i,j,0,l) &
                               + GMTR_t(i-1,j-1,k0,l,TJ,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 2.0_RP * GMTR_a(i-1,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 2.0_RP * GMTR_a(i-1,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) )

             coef_lap(i,j,0,l) = coef_lap(i,j,0,l) &
                               + GMTR_t(i  ,j-1,k0,l,TJ,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j-1,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) )

             ! ip1j
             coef_lap(i,j,1,l) = coef_lap(i,j,1,l) &
                               + GMTR_t(i  ,j-1,k0,l,TJ,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i  ,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j-1,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) )

             coef_lap(i,j,1,l) = coef_lap(i,j,1,l) &
                               + GMTR_t(i  ,j  ,k0,l,TI,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i+1,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i+1,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) )

             ! ip1jp1
             coef_lap(i,j,2,l) = coef_lap(i,j,2,l) &
                               + GMTR_t(i  ,j  ,k0,l,TI,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i+1,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i+1,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) )

             coef_lap(i,j,2,l) = coef_lap(i,j,2,l) &
                               + GMTR_t(i  ,j  ,k0,l,TJ,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i  ,j+1,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j+1,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) )

             ! ijp1
             coef_lap(i,j,3,l) = coef_lap(i,j,3,l) &
                               + GMTR_t(i  ,j  ,k0,l,TJ,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i  ,j+1,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j+1,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) )

             coef_lap(i,j,3,l) = coef_lap(i,j,3,l) &
                               + GMTR_t(i-1,j  ,k0,l,TI,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i-1,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   + 2.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AIJ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 2.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) )

             ! im1j
             coef_lap(i,j,4,l) = coef_lap(i,j,4,l) &
                               + GMTR_t(i-1,j  ,k0,l,TI,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i-1,j  ,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AJ ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AIJ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AJ ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) )

             coef_lap(i,j,4,l) = coef_lap(i,j,4,l) &
                               + GMTR_t(i-1,j-1,k0,l,TJ,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i-1,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 2.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 2.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) )

             ! im1jm1
             coef_lap(i,j,5,l) = coef_lap(i,j,5,l) &
                               + GMTR_t(i-1,j-1,k0,l,TJ,T_RAREA) &
                               * ( - 1.0_RP * GMTR_a(i-1,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 2.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j  ,k0,l,AI ,hn) &
                                   - 1.0_RP * GMTR_a(i-1,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 2.0_RP * GMTR_a(i-1,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 1.0_RP * GMTR_a(i-1,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) )

             ! ijm1
             coef_lap(i,j,6,l) = coef_lap(i,j,6,l) &
                               + GMTR_t(i  ,j-1,k0,l,TJ,T_RAREA) &
                               * ( + 1.0_RP * GMTR_a(i  ,j-1,k0,l,AIJ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   + 2.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i-1,j-1,k0,l,AIJ,hn) &
                                   - 1.0_RP * GMTR_a(i  ,j-1,k0,l,AIJ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   - 2.0_RP * GMTR_a(i  ,j  ,k0,l,AI ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) &
                                   + 1.0_RP * GMTR_a(i  ,j-1,k0,l,AJ ,tn) * GMTR_a(i  ,j  ,k0,l,AI ,hn) )
          enddo
       endif

       coef_lap(:,:,0,l) = coef_lap(:,:,0,l) * GMTR_p(:,:,k0,l,P_RAREA) / 12.0_RP
       coef_lap(:,:,1,l) = coef_lap(:,:,1,l) * GMTR_p(:,:,k0,l,P_RAREA) / 12.0_RP
       coef_lap(:,:,2,l) = coef_lap(:,:,2,l) * GMTR_p(:,:,k0,l,P_RAREA) / 12.0_RP
       coef_lap(:,:,3,l) = coef_lap(:,:,3,l) * GMTR_p(:,:,k0,l,P_RAREA) / 12.0_RP
       coef_lap(:,:,4,l) = coef_lap(:,:,4,l) * GMTR_p(:,:,k0,l,P_RAREA) / 12.0_RP
       coef_lap(:,:,5,l) = coef_lap(:,:,5,l) * GMTR_p(:,:,k0,l,P_RAREA) / 12.0_RP
       coef_lap(:,:,6,l) = coef_lap(:,:,6,l) * GMTR_p(:,:,k0,l,P_RAREA) / 12.0_RP
    enddo

    if ( ADM_have_pl ) then
       n  = ADM_gslf_pl

       do l = 1,ADM_lall_pl

          coef_lap_pl(:,l) = 0.0_RP

          do d = 1, ADM_nxyz
             hn  = d + HNX  - 1
             tn  = d + TNX  - 1
             tn2 = d + TN2X - 1

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                ijm1 = v - 1
                if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl
                if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

                coef_lap_pl(0,l) = coef_lap_pl(0,l) &
                                 + GMTR_t_pl(ijm1,k0,l,T_RAREA) &
                                 * ( + 1.0_RP * GMTR_a_pl(ijm1,k0,l,tn ) * GMTR_a_pl(ij,k0,l,hn) &
                                     - 2.0_RP * GMTR_a_pl(ijm1,k0,l,tn2) * GMTR_a_pl(ij,k0,l,hn) &
                                     - 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn ) * GMTR_a_pl(ij,k0,l,hn) )

                coef_lap_pl(0,l) = coef_lap_pl(0,l) &
                                 + GMTR_t_pl(ij  ,k0,l,T_RAREA) &
                                 * ( + 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn ) * GMTR_a_pl(ij,k0,l,hn) &
                                     - 2.0_RP * GMTR_a_pl(ij  ,k0,l,tn2) * GMTR_a_pl(ij,k0,l,hn) &
                                     - 1.0_RP * GMTR_a_pl(ijp1,k0,l,tn ) * GMTR_a_pl(ij,k0,l,hn) )
             enddo

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                ijm1 = v - 1
                if( ijp1 == ADM_gmax_pl + 1 ) ijp1 = ADM_gmin_pl
                if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

                coef_lap_pl(v-1,l) = coef_lap_pl(v-1,l) &
                                   + GMTR_t_pl(ijm1,k0,l,T_RAREA) &
                                   * ( - 2.0_RP * GMTR_a_pl(ijm1,k0,l,tn ) * GMTR_a_pl(ijm1,k0,l,hn) &
                                       + 1.0_RP * GMTR_a_pl(ijm1,k0,l,tn2) * GMTR_a_pl(ijm1,k0,l,hn) &
                                       - 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn ) * GMTR_a_pl(ijm1,k0,l,hn) &
                                       - 2.0_RP * GMTR_a_pl(ijm1,k0,l,tn ) * GMTR_a_pl(ij  ,k0,l,hn) &
                                       + 1.0_RP * GMTR_a_pl(ijm1,k0,l,tn2) * GMTR_a_pl(ij  ,k0,l,hn) &
                                       - 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn ) * GMTR_a_pl(ij  ,k0,l,hn) )

                coef_lap_pl(v-1,l) = coef_lap_pl(v-1,l) &
                                   + GMTR_t_pl(ij  ,k0,l,T_RAREA) &
                                   * ( + 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn ) * GMTR_a_pl(ij  ,k0,l,hn) &
                                       + 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn2) * GMTR_a_pl(ij  ,k0,l,hn) &
                                       + 2.0_RP * GMTR_a_pl(ijp1,k0,l,tn ) * GMTR_a_pl(ij  ,k0,l,hn) &
                                       + 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn ) * GMTR_a_pl(ijp1,k0,l,hn) &
                                       + 1.0_RP * GMTR_a_pl(ij  ,k0,l,tn2) * GMTR_a_pl(ijp1,k0,l,hn) &
                                       + 2.0_RP * GMTR_a_pl(ijp1,k0,l,tn ) * GMTR_a_pl(ijp1,k0,l,hn) )
             enddo
          enddo ! d loop

          do v = ADM_gslf_pl, ADM_gmax_pl
             coef_lap_pl(v-1,l) = coef_lap_pl(v-1,l) * GMTR_p_pl(n,k0,l,P_RAREA) / 12.0_RP
          enddo

       enddo ! l loop
    endif

    return
  end subroutine OPRT_laplacian_setup

  !-----------------------------------------------------------------------------
  subroutine OPRT_diffusion_setup( &
       GMTR_p,    GMTR_p_pl,    &
       GMTR_t,    GMTR_t_pl,    &
       GMTR_a,    GMTR_a_pl,    &
       coef_intp, coef_intp_pl, &
       coef_diff, coef_diff_pl  )
!ESC!    use mod_adm, only: &
!ESC!       ADM_have_pl,  &
!ESC!       ADM_have_sgp, &
!ESC!       ADM_jmin,     &
!ESC!       ADM_jmax,     &
!ESC!       ADM_imin,     &
!ESC!       ADM_imax,     &
!ESC!       ADM_gslf_pl,  &
!ESC!       ADM_gmin_pl,  &
!ESC!       ADM_gmax_pl
!ESC!    use mod_gmtr, only: &
!ESC!       P_RAREA => GMTR_p_RAREA, &
!ESC!       T_RAREA => GMTR_t_RAREA, &
!ESC!       HNX     => GMTR_a_HNX,   &
!ESC!       TNX     => GMTR_a_TNX,   &
!ESC!       TN2X    => GMTR_a_TN2X,  &
!ESC!       GMTR_p_nmax,             &
!ESC!       GMTR_t_nmax,             &
!ESC!       GMTR_a_nmax,             &
!ESC!       GMTR_a_nmax_pl
    implicit none

    real(RP), intent(in)  :: GMTR_p      (ADM_iall,ADM_jall,K0,ADM_lall   ,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_p_pl   (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_p_nmax   )
    real(RP), intent(in)  :: GMTR_t      (ADM_iall,ADM_jall,K0,ADM_lall   ,TI:TJ,GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_t_pl   (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_t_nmax   )
    real(RP), intent(in)  :: GMTR_a      (ADM_iall,ADM_jall,K0,ADM_lall   ,AI:AJ,GMTR_a_nmax   )
    real(RP), intent(in)  :: GMTR_a_pl   (ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_a_nmax_pl)
!    real(RP), intent(out) :: coef_intp   (ADM_nxyz,ADM_gall   ,1:3,TI:TJ,ADM_lall   )
!    real(RP), intent(out) :: coef_intp_pl(ADM_nxyz,ADM_gall_pl,1:3,      ADM_lall_pl)
!    real(RP), intent(out) :: coef_diff   (ADM_nxyz,ADM_gall,1:6        ,ADM_lall   )
!    real(RP), intent(out) :: coef_diff_pl(ADM_nxyz,         1:ADM_vlink,ADM_lall_pl)
    real(RP), intent(out) :: coef_intp   (ADM_iall,ADM_jall,1:3        ,ADM_nxyz,TI:TJ,ADM_lall   )
    real(RP), intent(out) :: coef_intp_pl(ADM_gall_pl      ,1:3        ,ADM_nxyz,      ADM_lall_pl)
    real(RP), intent(out) :: coef_diff   (ADM_iall,ADM_jall,1:6        ,ADM_nxyz,      ADM_lall   )
    real(RP), intent(out) :: coef_diff_pl(                  1:ADM_vlink,ADM_nxyz,      ADM_lall_pl)

    integer  :: ij, ijp1

    integer :: i, j, l, d, n, v, hn, tn, tn2
    !---------------------------------------------------------------------------

    write(ADM_LOG_FID,*) '*** setup coefficient of diffusion operator'

    coef_intp   (:,:,:,:,:,:) = 0.0_RP
    coef_intp_pl(:,:,:,:)     = 0.0_RP
    coef_diff   (:,:,:,:,:)   = 0.0_RP
    coef_diff_pl(:,:,:)       = 0.0_RP

    do l = 1, ADM_lall
       do d = 1, ADM_nxyz
          hn = d + HNX - 1
          tn = d + TNX - 1

          do j = ADM_jmin-1, ADM_jmax
          do i = ADM_imin-1, ADM_imax

!             coef_intp(d,ij,1,TI,l) = + GMTR_a(ij  ,k0,l,AIJ,tn) - GMTR_a(ij  ,k0,l,AI ,tn)
!             coef_intp(d,ij,2,TI,l) = - GMTR_a(ij  ,k0,l,AI ,tn) - GMTR_a(ip1j,k0,l,AJ ,tn)
!             coef_intp(d,ij,3,TI,l) = - GMTR_a(ip1j,k0,l,AJ ,tn) + GMTR_a(ij  ,k0,l,AIJ,tn)
             coef_intp(i,j,1,d,TI,l) = + GMTR_a(i  ,j  ,k0,l,AIJ,tn) - GMTR_a(i  ,j  ,k0,l,AI ,tn)
             coef_intp(i,j,2,d,TI,l) = - GMTR_a(i  ,j  ,k0,l,AI ,tn) - GMTR_a(i+1,j  ,k0,l,AJ ,tn)
             coef_intp(i,j,3,d,TI,l) = - GMTR_a(i+1,j  ,k0,l,AJ ,tn) + GMTR_a(i  ,j  ,k0,l,AIJ,tn)

!             coef_intp(d,ij,1,TJ,l) = + GMTR_a(ij  ,k0,l,AJ ,tn) - GMTR_a(ij  ,k0,l,AIJ,tn)
!             coef_intp(d,ij,2,TJ,l) = - GMTR_a(ij  ,k0,l,AIJ,tn) + GMTR_a(ijp1,k0,l,AI ,tn)
!             coef_intp(d,ij,3,TJ,l) = + GMTR_a(ijp1,k0,l,AI ,tn) + GMTR_a(ij  ,k0,l,AJ ,tn)
             coef_intp(i,j,1,d,TJ,l) = + GMTR_a(i  ,j  ,k0,l,AJ ,tn) - GMTR_a(i  ,j  ,k0,l,AIJ,tn)
             coef_intp(i,j,2,d,TJ,l) = - GMTR_a(i  ,j  ,k0,l,AIJ,tn) + GMTR_a(i  ,j+1,k0,l,AI ,tn)
             coef_intp(i,j,3,d,TJ,l) = + GMTR_a(i  ,j+1,k0,l,AI ,tn) + GMTR_a(i  ,j  ,k0,l,AJ ,tn)

!             coef_intp(d,ij,:,TI,l) = coef_intp(d,ij,:,TI,l) * 0.5_RP * GMTR_t(ij,k0,l,TI,T_RAREA)
!             coef_intp(d,ij,:,TJ,l) = coef_intp(d,ij,:,TJ,l) * 0.5_RP * GMTR_t(ij,k0,l,TJ,T_RAREA)
             coef_intp(i,j,:,d,TI,l) = coef_intp(i,j,:,d,TI,l) * 0.5_RP * GMTR_t(i,j,k0,l,TI,T_RAREA)
             coef_intp(i,j,:,d,TJ,l) = coef_intp(i,j,:,d,TJ,l) * 0.5_RP * GMTR_t(i,j,k0,l,TJ,T_RAREA)
          enddo
          enddo

          do j = ADM_jmin, ADM_jmax
          do i = ADM_imin, ADM_imax

!             coef_diff(d,ij,1,l) = + GMTR_a(i  ,j  ,k0,l,AIJ,hn) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
!             coef_diff(d,ij,2,l) = + GMTR_a(i  ,j  ,k0,l,AJ ,hn) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
!             coef_diff(d,ij,3,l) = - GMTR_a(i-1,j  ,k0,l,AI ,hn) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
!             coef_diff(d,ij,4,l) = - GMTR_a(i-1,j-1,k0,l,AIJ,hn) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
!             coef_diff(d,ij,5,l) = - GMTR_a(i  ,j-1,k0,l,AJ ,hn) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
!             coef_diff(d,ij,6,l) = + GMTR_a(i  ,j  ,k0,l,AI ,hn) * 0.5_RP * GMTR_p(ij,k0,l,P_RAREA)
             coef_diff(i,j,1,d,l) = + GMTR_a(i  ,j  ,k0,l,AIJ,hn) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
             coef_diff(i,j,2,d,l) = + GMTR_a(i  ,j  ,k0,l,AJ ,hn) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
             coef_diff(i,j,3,d,l) = - GMTR_a(i-1,j  ,k0,l,AI ,hn) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
             coef_diff(i,j,4,d,l) = - GMTR_a(i-1,j-1,k0,l,AIJ,hn) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
             coef_diff(i,j,5,d,l) = - GMTR_a(i  ,j-1,k0,l,AJ ,hn) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
             coef_diff(i,j,6,d,l) = + GMTR_a(i  ,j  ,k0,l,AI ,hn) * 0.5_RP * GMTR_p(i,j,k0,l,P_RAREA)
          enddo
          enddo
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          coef_diff(ADM_imin,ADM_jmin,5,:,l) = 0.0_RP
       endif

    enddo ! l loop

    if ( ADM_have_pl ) then
       n  = ADM_gslf_pl

       do l = 1, ADM_lall_pl

          coef_intp_pl(:,:,:,l) = 0.0_RP

          do d = 1, ADM_nxyz
             hn  = d + HNX  - 1
             tn  = d + TNX  - 1
             tn2 = d + TN2X - 1

             do v = ADM_gmin_pl, ADM_gmax_pl
                ij   = v
                ijp1 = v + 1
                if( ijp1 == ADM_gmax_pl+1 ) ijp1 = ADM_gmin_pl

!                coef_intp_pl(d,v,1,l) = - GMTR_a_pl(ijp1,k0,l,tn ) + GMTR_a_pl(ij  ,k0,l,tn )
!                coef_intp_pl(d,v,2,l) = + GMTR_a_pl(ij  ,k0,l,tn ) + GMTR_a_pl(ij  ,k0,l,tn2)
!                coef_intp_pl(d,v,3,l) = + GMTR_a_pl(ij  ,k0,l,tn2) - GMTR_a_pl(ijp1,k0,l,tn )
                coef_intp_pl(v,1,d,l) = - GMTR_a_pl(ijp1,k0,l,tn ) + GMTR_a_pl(ij  ,k0,l,tn )
                coef_intp_pl(v,2,d,l) = + GMTR_a_pl(ij  ,k0,l,tn ) + GMTR_a_pl(ij  ,k0,l,tn2)
                coef_intp_pl(v,3,d,l) = + GMTR_a_pl(ij  ,k0,l,tn2) - GMTR_a_pl(ijp1,k0,l,tn )

!                coef_intp_pl(d,v,:,l) = coef_intp_pl(d,v,:,l) * 0.5_RP * GMTR_t_pl(v,k0,l,T_RAREA)
                coef_intp_pl(v,:,d,l) = coef_intp_pl(v,:,d,l) * 0.5_RP * GMTR_t_pl(v,k0,l,T_RAREA)

!                coef_diff_pl(d,v-1,l) = GMTR_a_pl(v,k0,l,hn) * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)
                coef_diff_pl(v-1,d,l) = GMTR_a_pl(v,k0,l,hn) * 0.5_RP * GMTR_p_pl(n,k0,l,P_RAREA)
             enddo
          enddo
       enddo ! l loop
    endif

    return
  end subroutine OPRT_diffusion_setup

end module mod_oprt
!-------------------------------------------------------------------------------
