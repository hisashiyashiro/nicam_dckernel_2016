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
  public :: OPRT_diffusion

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
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
  subroutine OPRT_diffusion( &
       dscl,      dscl_pl,      &
       scl,       scl_pl,       &
       kh,        kh_pl,        &
       coef_intp, coef_intp_pl, &
       coef_diff, coef_diff_pl  )
!ESC!    use mod_adm, only: &
!ESC!       ADM_have_pl,  &
!ESC!       ADM_have_sgp, &
!ESC!       ADM_imin,     &
!ESC!       ADM_imax,     &
!ESC!       ADM_jmin,     &
!ESC!       ADM_jmax,     &
!ESC!       ADM_gslf_pl,  &
!ESC!       ADM_gmin_pl,  &
!ESC!       ADM_gmax_pl
!ESC!    use mod_grd, only: &
!ESC!       XDIR => GRD_XDIR, &
!ESC!       YDIR => GRD_YDIR, &
!ESC!       ZDIR => GRD_ZDIR
    implicit none

    real(RP), intent(out) :: dscl        (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(out) :: dscl_pl     (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: scl         (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: scl_pl      (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: kh          (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: kh_pl       (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: coef_intp   (ADM_iall,ADM_jall,1:3,    ADM_nxyz,TI:TJ,ADM_lall   )
    real(RP), intent(in)  :: coef_intp_pl(ADM_gall_pl      ,1:3,    ADM_nxyz,      ADM_lall_pl)
    real(RP), intent(in)  :: coef_diff   (ADM_iall,ADM_jall,1:6,    ADM_nxyz,      ADM_lall   )
    real(RP), intent(in)  :: coef_diff_pl(              1:ADM_vlink,ADM_nxyz,      ADM_lall_pl)

    real(RP) :: vt   (ADM_iall,ADM_jall,ADM_nxyz,TI:TJ)
    real(RP) :: vt_pl(ADM_gall_pl      ,ADM_nxyz)

    integer  :: imin, imax, jmin, jmax, kall, lall
    integer  :: ij, ijp1, ijm1

    integer  :: i, j, k, l, n, v, d
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('OPRT_diffusion')

    imin = ADM_imin
    imax = ADM_imax
    jmin = ADM_jmin
    jmax = ADM_jmax
    kall = ADM_kall
    lall = ADM_lall

    !$omp parallel default(none),private(i,j,k,l,d), &
    !$omp shared(imin,imax,jmin,jmax,kall,lall,ADM_have_sgp,dscl,scl,kh,vt,coef_intp,coef_diff)
    do l = 1, lall
    do k = 1, kall
       !$omp do schedule(static) collapse(2)
       do d = XDIR, ZDIR
       do j = jmin-1, jmax
       do i = imin-1, imax
          vt(i,j,d,TI) = ( ( + 2.0_RP * coef_intp(i,j,1,d,TI,l) &
                             - 1.0_RP * coef_intp(i,j,2,d,TI,l) &
                             - 1.0_RP * coef_intp(i,j,3,d,TI,l) ) * scl(i  ,j  ,k,l) &
                         + ( - 1.0_RP * coef_intp(i,j,1,d,TI,l) &
                             + 2.0_RP * coef_intp(i,j,2,d,TI,l) &
                             - 1.0_RP * coef_intp(i,j,3,d,TI,l) ) * scl(i+1,j  ,k,l) &
                         + ( - 1.0_RP * coef_intp(i,j,1,d,TI,l) &
                             - 1.0_RP * coef_intp(i,j,2,d,TI,l) &
                             + 2.0_RP * coef_intp(i,j,3,d,TI,l) ) * scl(i+1,j+1,k,l) &
                         ) / 3.0_RP
       enddo
       enddo
       enddo
       !$omp end do nowait

       !$omp do schedule(static) collapse(2)
       do d = XDIR, ZDIR
       do j = jmin-1, jmax
       do i = imin-1, imax
          vt(i,j,d,TJ) = ( ( + 2.0_RP * coef_intp(i,j,1,d,TJ,l) &
                             - 1.0_RP * coef_intp(i,j,2,d,TJ,l) &
                             - 1.0_RP * coef_intp(i,j,3,d,TJ,l) ) * scl(i  ,j  ,k,l) &
                         + ( - 1.0_RP * coef_intp(i,j,1,d,TJ,l) &
                             + 2.0_RP * coef_intp(i,j,2,d,TJ,l) &
                             - 1.0_RP * coef_intp(i,j,3,d,TJ,l) ) * scl(i+1,j+1,k,l) &
                         + ( - 1.0_RP * coef_intp(i,j,1,d,TJ,l) &
                             - 1.0_RP * coef_intp(i,j,2,d,TJ,l) &
                             + 2.0_RP * coef_intp(i,j,3,d,TJ,l) ) * scl(i  ,j+1,k,l) &
                         ) / 3.0_RP
       enddo
       enddo
       enddo
       !$omp end do

       if ( ADM_have_sgp(l) ) then ! pentagon
          !$omp master
          vt(imin-1,jmin-1,:,TI) = vt(imin,jmin-1,:,TJ)
          !$omp end master
       endif

       !$omp do schedule(static)
       do j = jmin, jmax
       do i = imin, imax
          dscl(i,j,k,l) = ( coef_diff(i,j,1,XDIR,l) * ( vt(i  ,j  ,XDIR,TI) + vt(i  ,j  ,XDIR,TJ) ) &
                          + coef_diff(i,j,1,YDIR,l) * ( vt(i  ,j  ,YDIR,TI) + vt(i  ,j  ,YDIR,TJ) ) &
                          + coef_diff(i,j,1,ZDIR,l) * ( vt(i  ,j  ,ZDIR,TI) + vt(i  ,j  ,ZDIR,TJ) ) &
                          ) * 0.5_RP * ( kh(i  ,j  ,k,l) + kh(i+1,j+1,k,l) )
       enddo
       enddo
       !$omp end do

       !$omp do schedule(static)
       do j = jmin, jmax
       do i = imin, imax
          dscl(i,j,k,l) = dscl(i,j,k,l) &
                        + ( coef_diff(i,j,2,XDIR,l) * ( vt(i  ,j  ,XDIR,TJ) + vt(i-1,j  ,XDIR,TI) ) &
                          + coef_diff(i,j,2,YDIR,l) * ( vt(i  ,j  ,YDIR,TJ) + vt(i-1,j  ,YDIR,TI) ) &
                          + coef_diff(i,j,2,ZDIR,l) * ( vt(i  ,j  ,ZDIR,TJ) + vt(i-1,j  ,ZDIR,TI) ) &
                          ) * 0.5_RP * ( kh(i  ,j  ,k,l) + kh(i  ,j+1,k,l) )
       enddo
       enddo
       !$omp end do

       !$omp do schedule(static)
       do j = jmin, jmax
       do i = imin, imax
          dscl(i,j,k,l) = dscl(i,j,k,l) &
                        + ( coef_diff(i,j,3,XDIR,l) * ( vt(i-1,j  ,XDIR,TI) + vt(i-1,j-1,XDIR,TJ) ) &
                          + coef_diff(i,j,3,YDIR,l) * ( vt(i-1,j  ,YDIR,TI) + vt(i-1,j-1,YDIR,TJ) ) &
                          + coef_diff(i,j,3,ZDIR,l) * ( vt(i-1,j  ,ZDIR,TI) + vt(i-1,j-1,ZDIR,TJ) ) &
                          ) * 0.5_RP * ( kh(i-1,j  ,k,l) + kh(i  ,j  ,k,l) )
       enddo
       enddo
       !$omp end do

       !$omp do schedule(static)
       do j = jmin, jmax
       do i = imin, imax
          dscl(i,j,k,l) = dscl(i,j,k,l) &
                        + ( coef_diff(i,j,4,XDIR,l) * ( vt(i-1,j-1,XDIR,TJ) + vt(i-1,j-1,XDIR,TI) ) &
                          + coef_diff(i,j,4,YDIR,l) * ( vt(i-1,j-1,YDIR,TJ) + vt(i-1,j-1,YDIR,TI) ) &
                          + coef_diff(i,j,4,ZDIR,l) * ( vt(i-1,j-1,ZDIR,TJ) + vt(i-1,j-1,ZDIR,TI) ) &
                          ) * 0.5_RP * ( kh(i-1,j-1,k,l) + kh(i  ,j  ,k,l) )
       enddo
       enddo
       !$omp end do

       !$omp do schedule(static)
       do j = jmin, jmax
       do i = imin, imax
          dscl(i,j,k,l) = dscl(i,j,k,l) &
                        + ( coef_diff(i,j,5,XDIR,l) * ( vt(i-1,j-1,XDIR,TI) + vt(i  ,j-1,XDIR,TJ) ) &
                          + coef_diff(i,j,5,YDIR,l) * ( vt(i-1,j-1,YDIR,TI) + vt(i  ,j-1,YDIR,TJ) ) &
                          + coef_diff(i,j,5,ZDIR,l) * ( vt(i-1,j-1,ZDIR,TI) + vt(i  ,j-1,ZDIR,TJ) ) &
                          ) * 0.5_RP * ( kh(i  ,j-1,k,l) + kh(i  ,j  ,k,l) )
       enddo
       enddo
       !$omp end do

       !$omp do schedule(static)
       do j = jmin, jmax
       do i = imin, imax
          dscl(i,j,k,l) = dscl(i,j,k,l) &
                        + ( coef_diff(i,j,6,XDIR,l) * ( vt(i  ,j-1,XDIR,TJ) + vt(i  ,j  ,XDIR,TI) ) &
                          + coef_diff(i,j,6,YDIR,l) * ( vt(i  ,j-1,YDIR,TJ) + vt(i  ,j  ,YDIR,TI) ) &
                          + coef_diff(i,j,6,ZDIR,l) * ( vt(i  ,j-1,ZDIR,TJ) + vt(i  ,j  ,ZDIR,TI) ) &
                          ) * 0.5_RP * ( kh(i  ,j  ,k,l) + kh(i+1,j  ,k,l) )
       enddo
       enddo
       !$omp end do

       !$omp workshare
       dscl(:,jmin-1,k,l) = 0.0_RP
       dscl(:,jmax+1,k,l) = 0.0_RP
       dscl(imin-1,:,k,l) = 0.0_RP
       dscl(imax+1,:,k,l) = 0.0_RP
       !$omp end workshare
    enddo ! k loop
    enddo ! l loop
    !$omp end parallel

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
       do k = 1, ADM_kall

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             if( ijp1 == ADM_gmax_pl+1 ) ijp1 = ADM_gmin_pl

             do d = 1, ADM_nxyz
                vt_pl(ij,d) = ( ( + 2.0_RP * coef_intp_pl(v,1,d,l) &
                                  - 1.0_RP * coef_intp_pl(v,2,d,l) &
                                  - 1.0_RP * coef_intp_pl(v,3,d,l) ) * scl_pl(n   ,k,l) &
                              + ( - 1.0_RP * coef_intp_pl(v,1,d,l) &
                                  + 2.0_RP * coef_intp_pl(v,2,d,l) &
                                  - 1.0_RP * coef_intp_pl(v,3,d,l) ) * scl_pl(ij  ,k,l) &
                              + ( - 1.0_RP * coef_intp_pl(v,1,d,l) &
                                  - 1.0_RP * coef_intp_pl(v,2,d,l) &
                                  + 2.0_RP * coef_intp_pl(v,3,d,l) ) * scl_pl(ijp1,k,l) &
                              ) / 3.0_RP
             enddo
          enddo

          dscl_pl(:,k,l) = 0.0_RP

          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijm1 = v - 1
             if( ijm1 == ADM_gmin_pl-1 ) ijm1 = ADM_gmax_pl ! cyclic condition

             dscl_pl(n,k,l) = dscl_pl(n,k,l) &
                            + ( coef_diff_pl(v-1,XDIR,l) * ( vt_pl(ijm1,XDIR) + vt_pl(ij,XDIR) ) &
                              + coef_diff_pl(v-1,YDIR,l) * ( vt_pl(ijm1,YDIR) + vt_pl(ij,YDIR) ) &
                              + coef_diff_pl(v-1,ZDIR,l) * ( vt_pl(ijm1,ZDIR) + vt_pl(ij,ZDIR) ) &
                              ) * 0.5_RP * ( kh_pl(n,k,l) + kh_pl(ij,k,l) )
          enddo

       enddo
       enddo
    else
       dscl_pl(:,:,:) = 0.0_RP
    endif

    call DEBUG_rapend('OPRT_diffusion')

    return
  end subroutine OPRT_diffusion

end module mod_oprt
!-------------------------------------------------------------------------------
