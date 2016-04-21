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
!ESC!  use mod_adm, only: &
!ESC!     ADM_LOG_FID,      &
!ESC!     ADM_NSYS,         &
!ESC!     ADM_MAXFNAME,     &
!ESC!     TI  => ADM_TI,    &
!ESC!     TJ  => ADM_TJ,    &
!ESC!     AI  => ADM_AI,    &
!ESC!     AIJ => ADM_AIJ,   &
!ESC!     AJ  => ADM_AJ,    &
!ESC!     K0  => ADM_KNONE, &
!ESC!     ADM_nxyz,         &
!ESC!     ADM_lall,         &
!ESC!     ADM_lall_pl,      &
!ESC!     ADM_gall,         &
!ESC!     ADM_gall_pl,      &
!ESC!     ADM_kall
!ESC!  use mod_gmtr, only: &
!ESC!     P_RAREA => GMTR_P_RAREA, &
!ESC!     T_RAREA => GMTR_T_RAREA, &
!ESC!     W1      => GMTR_T_W1,    &
!ESC!     W2      => GMTR_T_W2,    &
!ESC!     W3      => GMTR_T_W3,    &
!ESC!     HNX     => GMTR_A_HNX,   &
!ESC!     HNY     => GMTR_A_HNY,   &
!ESC!     HNZ     => GMTR_A_HNZ,   &
!ESC!     HTX     => GMTR_A_HTX,   &
!ESC!     HTY     => GMTR_A_HTY,   &
!ESC!     HTZ     => GMTR_A_HTZ,   &
!ESC!     TNX     => GMTR_A_TNX,   &
!ESC!     TNY     => GMTR_A_TNY,   &
!ESC!     TNZ     => GMTR_A_TNZ,   &
!ESC!     TN2X    => GMTR_A_TN2X,  &
!ESC!     TN2Y    => GMTR_A_TN2Y,  &
!ESC!     TN2Z    => GMTR_A_TN2Z
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: OPRT_diffusion_putcoef
  public :: OPRT_diffusion

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !

#ifdef _FIXEDINDEX_
  real(RP), public              :: cinterp_TN    (ADM_gall   ,ADM_lall   ,AI:AJ,ADM_nxyz)
  real(RP), public              :: cinterp_TN_pl (ADM_gall_pl,ADM_lall_pl,2    ,ADM_nxyz)
  real(RP), public              :: cinterp_HN    (ADM_gall   ,ADM_lall   ,AI:AJ,ADM_nxyz)
  real(RP), public              :: cinterp_HN_pl (ADM_gall_pl,ADM_lall_pl,      ADM_nxyz)
  real(RP), public              :: cinterp_TRA   (ADM_gall   ,ADM_lall   ,TI:TJ         )
  real(RP), public              :: cinterp_TRA_pl(ADM_gall_pl,ADM_lall_pl               )
  real(RP), public              :: cinterp_PRA   (ADM_gall   ,ADM_lall                  )
  real(RP), public              :: cinterp_PRA_pl(ADM_gall_pl,ADM_lall_pl               )
#else
  real(RP), public, allocatable :: cinterp_TN    (:,:,:,:) ! coefficient for diffusion operator
  real(RP), public, allocatable :: cinterp_TN_pl (:,:,:,:)
  real(RP), public, allocatable :: cinterp_HN    (:,:,:,:)
  real(RP), public, allocatable :: cinterp_HN_pl (:,:,:)
  real(RP), public, allocatable :: cinterp_TRA   (:,:,:)
  real(RP), public, allocatable :: cinterp_TRA_pl(:,:)
  real(RP), public, allocatable :: cinterp_PRA   (:,:)
  real(RP), public, allocatable :: cinterp_PRA_pl(:,:)
#endif

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine OPRT_diffusion_putcoef( &
       cinterp_TN_in , cinterp_TN_pl_in , &
       cinterp_HN_in , cinterp_HN_pl_in , &
       cinterp_TRA_in, cinterp_TRA_pl_in, &
       cinterp_PRA_in, cinterp_PRA_pl_in  )
    implicit none

    real(RP), intent(in)  :: cinterp_TN_in    (ADM_gall   ,ADM_lall   ,AI:AJ,ADM_nxyz)
    real(RP), intent(in)  :: cinterp_TN_pl_in (ADM_gall_pl,ADM_lall_pl,2    ,ADM_nxyz)
    real(RP), intent(in)  :: cinterp_HN_in    (ADM_gall   ,ADM_lall   ,AI:AJ,ADM_nxyz)
    real(RP), intent(in)  :: cinterp_HN_pl_in (ADM_gall_pl,ADM_lall_pl,      ADM_nxyz)
    real(RP), intent(in)  :: cinterp_TRA_in   (ADM_gall   ,ADM_lall   ,TI:TJ         )
    real(RP), intent(in)  :: cinterp_TRA_pl_in(ADM_gall_pl,ADM_lall_pl               )
    real(RP), intent(in)  :: cinterp_PRA_in   (ADM_gall   ,ADM_lall                  )
    real(RP), intent(in)  :: cinterp_PRA_pl_in(ADM_gall_pl,ADM_lall_pl               )
    !---------------------------------------------------------------------------

#ifndef _FIXEDINDEX_
    allocate( cinterp_TN    (ADM_gall   ,ADM_lall   ,AI:AJ,ADM_nxyz) )
    allocate( cinterp_TN_pl (ADM_gall_pl,ADM_lall_pl,2    ,ADM_nxyz) )
    allocate( cinterp_HN    (ADM_gall   ,ADM_lall   ,AI:AJ,ADM_nxyz) )
    allocate( cinterp_HN_pl (ADM_gall_pl,ADM_lall_pl,      ADM_nxyz) )
    allocate( cinterp_TRA   (ADM_gall   ,ADM_lall   ,TI:TJ         ) )
    allocate( cinterp_TRA_pl(ADM_gall_pl,ADM_lall_pl               ) )
    allocate( cinterp_PRA   (ADM_gall   ,ADM_lall                  ) )
    allocate( cinterp_PRA_pl(ADM_gall_pl,ADM_lall_pl               ) )
#endif

    cinterp_TN    (:,:,:,:) = cinterp_TN_in    (:,:,:,:)
    cinterp_TN_pl (:,:,:,:) = cinterp_TN_pl_in (:,:,:,:)
    cinterp_HN    (:,:,:,:) = cinterp_HN_in    (:,:,:,:)
    cinterp_HN_pl (:,:,  :) = cinterp_HN_pl_in (:,:,  :)
    cinterp_TRA   (:,:,:  ) = cinterp_TRA_in   (:,:,:  )
    cinterp_TRA_pl(:,:    ) = cinterp_TRA_pl_in(:,:    )
    cinterp_PRA   (:,:    ) = cinterp_PRA_in   (:,:    )
    cinterp_PRA_pl(:,:    ) = cinterp_PRA_pl_in(:,:    )

    return
  end subroutine OPRT_diffusion_putcoef

  !-----------------------------------------------------------------------------
  subroutine OPRT_diffusion( &
       dscl, dscl_pl, &
       scl,  scl_pl,  &
       kh,   kh_pl    )
!ESC!    use mod_adm, only: &
!ESC!       ADM_have_pl,  &
!ESC!       ADM_have_sgp, &
!ESC!       ADM_gall_1d,  &
!ESC!       ADM_gmin,     &
!ESC!       ADM_gmax,     &
!ESC!       ADM_gslf_pl,  &
!ESC!       ADM_gmin_pl,  &
!ESC!       ADM_gmax_pl
    implicit none

    real(RP), intent(out) :: dscl   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: dscl_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: scl    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: scl_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: kh     (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: kh_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP)  :: vxt    (ADM_gall   ,TI:TJ)
    real(RP)  :: vyt    (ADM_gall   ,TI:TJ)
    real(RP)  :: vzt    (ADM_gall   ,TI:TJ)
    real(RP)  :: flux   (ADM_gall   ,AI:AJ)
    real(RP)  :: vxt_pl (ADM_gall_pl)
    real(RP)  :: vyt_pl (ADM_gall_pl)
    real(RP)  :: vzt_pl (ADM_gall_pl)
    real(RP)  :: flux_pl(ADM_gall_pl)

    real(RP) :: u1, u2, u3, smean

    integer :: nstart1, nstart2, nstart3, nstart4, nend
    integer :: ij
    integer :: ip1j, ijp1, ip1jp1
    integer :: im1j, ijm1, im1jm1

    integer :: gall, gall_1d, gmin, kall
    integer :: n, k, l, v
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('OPRT_diffusion')

    gall    = ADM_gall
    gall_1d = ADM_gall_1d
    gmin    = ADM_gmin
    kall    = ADM_kall

    nstart1 = suf(ADM_gmin-1,ADM_gmin-1)
    nstart2 = suf(ADM_gmin-1,ADM_gmin  )
    nstart3 = suf(ADM_gmin  ,ADM_gmin-1)
    nstart4 = suf(ADM_gmin  ,ADM_gmin  )
    nend    = suf(ADM_gmax  ,ADM_gmax  )

    do l = 1, ADM_lall
       !$omp parallel default(none),private(n,k,ij,ip1j,ip1jp1,ijp1,im1j,ijm1,im1jm1,smean,u1,u2,u3), &
       !$omp shared(gall,gall_1d,gmin,kall,nstart1,nstart2,nstart3,nstart4,nend,l,ADM_have_sgp, &
       !$omp dscl,scl,kh,vxt,vyt,vzt,flux,cinterp_TN,cinterp_HN,cinterp_PRA,cinterp_TRA)
       do k = 1, kall

          !$omp do
          do n = nstart1, nend
             ij     = n
             ip1j   = n + 1
             ip1jp1 = n + 1 + gall_1d

             smean = ( scl(ij,k,l) + scl(ip1j,k,l) + scl(ip1jp1,k,l) ) / 3.0_RP

             u1 = 0.5_RP * (scl(ij    ,k,l)+scl(ip1j  ,k,l)) - smean
             u2 = 0.5_RP * (scl(ip1j  ,k,l)+scl(ip1jp1,k,l)) - smean
             u3 = 0.5_RP * (scl(ip1jp1,k,l)+scl(ij    ,k,l)) - smean

             vxt(n,TI) = ( - u1 * cinterp_TN(ij  ,l,AI ,1) &
                           - u2 * cinterp_TN(ip1j,l,AJ ,1) &
                           + u3 * cinterp_TN(ij  ,l,AIJ,1) ) * cinterp_TRA(ij,l,TI)
             vyt(n,TI) = ( - u1 * cinterp_TN(ij  ,l,AI ,2) &
                           - u2 * cinterp_TN(ip1j,l,AJ ,2) &
                           + u3 * cinterp_TN(ij  ,l,AIJ,2) ) * cinterp_TRA(ij,l,TI)
             vzt(n,TI) = ( - u1 * cinterp_TN(ij  ,l,AI ,3) &
                           - u2 * cinterp_TN(ip1j,l,AJ ,3) &
                           + u3 * cinterp_TN(ij  ,l,AIJ,3) ) * cinterp_TRA(ij,l,TI)
          enddo
          !$omp end do nowait

          !$omp do
          do n = nstart1, nend
             ij     = n
             ijp1   = n     + gall_1d
             ip1jp1 = n + 1 + gall_1d

             smean = ( scl(ij,k,l) + scl(ip1jp1,k,l) + scl(ijp1,k,l) ) / 3.0_RP

             u1 = 0.5_RP * (scl(ij    ,k,l)+scl(ip1jp1,k,l)) - smean
             u2 = 0.5_RP * (scl(ip1jp1,k,l)+scl(ijp1  ,k,l)) - smean
             u3 = 0.5_RP * (scl(ijp1  ,k,l)+scl(ij    ,k,l)) - smean

             vxt(n,TJ) = ( - u1 * cinterp_TN(ij  ,l,AIJ,1) &
                           + u2 * cinterp_TN(ijp1,l,AI ,1) &
                           + u3 * cinterp_TN(ij  ,l,AJ ,1) ) * cinterp_TRA(ij,l,TJ)
             vyt(n,TJ) = ( - u1 * cinterp_TN(ij  ,l,AIJ,2) &
                           + u2 * cinterp_TN(ijp1,l,AI ,2) &
                           + u3 * cinterp_TN(ij  ,l,AJ ,2) ) * cinterp_TRA(ij,l,TJ)
             vzt(n,TJ) = ( - u1 * cinterp_TN(ij  ,l,AIJ,3) &
                           + u2 * cinterp_TN(ijp1,l,AI ,3) &
                           + u3 * cinterp_TN(ij  ,l,AJ ,3) ) * cinterp_TRA(ij,l,TJ)
          enddo
          !$omp end do

          if ( ADM_have_sgp(l) ) then ! pentagon
             !$omp master
             vxt(suf(gmin-1,gmin-1),TI) = vxt(suf(gmin,gmin-1),TJ)
             vyt(suf(gmin-1,gmin-1),TI) = vyt(suf(gmin,gmin-1),TJ)
             vzt(suf(gmin-1,gmin-1),TI) = vzt(suf(gmin,gmin-1),TJ)
             !$omp end master
             !$omp barrier
          endif

          !$omp do
          do n = nstart2, nend
             ij     = n
             ip1j   = n + 1
             ijp1   = n     + gall_1d
             ip1jp1 = n + 1 + gall_1d
             im1j   = n - 1
             ijm1   = n     - gall_1d

             flux(n,AI ) = 0.25_RP * ( (vxt(ijm1,TJ)+vxt(ij,TI)) * cinterp_HN(ij,l,AI ,1) &
                                     + (vyt(ijm1,TJ)+vyt(ij,TI)) * cinterp_HN(ij,l,AI ,2) &
                                     + (vzt(ijm1,TJ)+vzt(ij,TI)) * cinterp_HN(ij,l,AI ,3) &
                                     ) * (kh(ij,k,l)+kh(ip1j,k,l))
          enddo
          !$omp end do nowait

          !$omp do
          do n = nstart1, nend
             ij     = n
             ip1jp1 = n + 1 + gall_1d

             flux(n,AIJ) = 0.25_RP * ( (vxt(ij ,TI)+vxt(ij ,TJ)) * cinterp_HN(ij,l,AIJ,1) &
                                     + (vyt(ij ,TI)+vyt(ij ,TJ)) * cinterp_HN(ij,l,AIJ,2) &
                                     + (vzt(ij ,TI)+vzt(ij ,TJ)) * cinterp_HN(ij,l,AIJ,3) &
                                     ) * (kh(ij,k,l)+kh(ip1jp1,k,l))
          enddo
          !$omp end do nowait

          !$omp do
          do n = nstart3, nend
             ij     = n
             ijp1   = n     + gall_1d
             im1j   = n - 1

             flux(n,AJ ) = 0.25_RP * ( (vxt(ij,TJ)+vxt(im1j,TI)) * cinterp_HN(ij,l,AJ ,1) &
                                     + (vyt(ij,TJ)+vyt(im1j,TI)) * cinterp_HN(ij,l,AJ ,2) &
                                     + (vzt(ij,TJ)+vzt(im1j,TI)) * cinterp_HN(ij,l,AJ ,3) &
                                     ) * (kh(ij,k,l)+kh(ijp1,k,l))
          enddo
          !$omp end do

          if ( ADM_have_sgp(l) ) then ! pentagon
             !$omp master
             flux(suf(gmin,gmin-1),AJ) = 0.0_RP
             !$omp end master
             !$omp barrier
          endif

          !$omp do
          do n = nstart4, nend
             ij     = n
             im1j   = n - 1
             im1jm1 = n - 1 - gall_1d
             ijm1   = n     - gall_1d

             dscl(n,k,l) = ( flux(ij,AI ) - flux(im1j  ,AI ) &
                           + flux(ij,AIJ) - flux(im1jm1,AIJ) &
                           + flux(ij,AJ ) - flux(ijm1  ,AJ ) ) * cinterp_PRA(ij,l)
          enddo
          !$omp end do nowait

          !$omp do
          do n = 1, nstart4-1
             dscl(n,k,l) = 0.0_RP
          enddo
          !$omp end do nowait

          !$omp do
          do n = nend+1, gall
             dscl(n,k,l) = 0.0_RP
          enddo
          !$omp end do

       enddo
       !$omp end parallel
    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             if( ijp1 > ADM_gmax_pl ) ijp1 = ADM_gmin_pl

             smean = ( scl_pl(n,k,l) + scl_pl(ij,k,l) + scl_pl(ijp1,k,l) ) / 3.0_RP

             u1 = 0.5_RP * (scl_pl(n   ,k,l)+scl_pl(ij  ,k,l)) - smean
             u2 = 0.5_RP * (scl_pl(ij  ,k,l)+scl_pl(ijp1,k,l)) - smean
             u3 = 0.5_RP * (scl_pl(ijp1,k,l)+scl_pl(n   ,k,l)) - smean

             vxt_pl(v) = ( + u1 * cinterp_TN_pl(ij  ,l,1,1) &
                           + u2 * cinterp_TN_pl(ij  ,l,2,1) &
                           - u3 * cinterp_TN_pl(ijp1,l,1,1) ) * cinterp_TRA_pl(ij,l)
             vyt_pl(v) = ( + u1 * cinterp_TN_pl(ij  ,l,1,2) &
                           + u2 * cinterp_TN_pl(ij  ,l,2,2) &
                           - u3 * cinterp_TN_pl(ijp1,l,1,2) ) * cinterp_TRA_pl(ij,l)
             vzt_pl(v) = ( + u1 * cinterp_TN_pl(ij  ,l,1,3) &
                           + u2 * cinterp_TN_pl(ij  ,l,2,3) &
                           - u3 * cinterp_TN_pl(ijp1,l,1,3) ) * cinterp_TRA_pl(ij,l)
          enddo

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijm1 = v - 1
             if( ijm1 == ADM_gmin_pl - 1 ) ijm1 = ADM_gmax_pl

             flux_pl(v) = 0.25_RP * ( (vxt_pl(ij)+vxt_pl(ijm1)) * cinterp_HN_pl(ij,l,1) &
                                    + (vyt_pl(ij)+vyt_pl(ijm1)) * cinterp_HN_pl(ij,l,2) &
                                    + (vzt_pl(ij)+vzt_pl(ijm1)) * cinterp_HN_pl(ij,l,3) &
                                    ) * (kh_pl(n,k,l)+kh_pl(ij,k,l))
          enddo

          dscl_pl(:,k,l) = 0.0_RP
          do v = ADM_gmin_pl, ADM_gmax_pl
             dscl_pl(n,k,l) = dscl_pl(n,k,l) + flux_pl(v)
          enddo
          dscl_pl(n,k,l) = dscl_pl(n,k,l) * cinterp_PRA_pl(n,l)

       enddo
       enddo

    else
       dscl_pl(:,:,:) = 0.0_RP
    endif

    call DEBUG_rapend('OPRT_diffusion')

    return
  end subroutine OPRT_diffusion

  !-----------------------------------------------------------------------------
  integer function suf(i,j)
!ESC!    use mod_adm, only: &
!ESC!       ADM_gall_1d
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    suf = ADM_gall_1d * (j-1) + i

  end function suf

end module mod_oprt
!-------------------------------------------------------------------------------
