!-------------------------------------------------------------------------------
!
!+  Vertical Implicit module
!
!-------------------------------------------------------------------------------
module mod_vi
  !-----------------------------------------------------------------------------
  !
  !++ Description:
  !       This module is for the vertical implicit scheme of non-hydorostatic
  !       model.
  !
  !
  !++ Current Corresponding Author : H.Tomita
  !
  !++ History:
  !      Version   Date       Comment
  !      -----------------------------------------------------------------------
  !      0.00      04-02-17   Imported from igdc-4.34
  !                06-08-11   Add averaged rhog for tracer advection.
  !                11-05-07   Y.Yamada: Implementation of ES tuning cord by NEC.
  !                             Modified line: (20110405 NEC)
  !                                                or (!ftr< vi_small_step.r???)
  !                11-11-28   Y.Yamada: Merge Terai-san timer code
  !                                                    into the original code.
  !                11-12-29   Y.Yamada: Delete ES tuning and merge fjtimer
  !                12-3-9    S.Iga: tuned (phase4-1)
  !
  !      -----------------------------------------------------------------------
  !
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_misc
!ESC!  use mod_adm, only: &
!ESC!     ADM_LOG_FID
!ESC!  use mod_adm, only: &
!ESC!     ADM_lall,    &
!ESC!     ADM_lall_pl, &
!ESC!     ADM_gall,    &
!ESC!     ADM_gall_pl, &
!ESC!     ADM_kall
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: vi_rhow_solver_putcoef
  public :: vi_rhow_solver

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
#ifdef _FIXEDINDEX_
  real(RP), public              :: Mc              (ADM_gall   ,ADM_kall,ADM_lall   )
  real(RP), public              :: Mc_pl           (ADM_gall_pl,ADM_kall,ADM_lall_pl)
  real(RP), public              :: Ml              (ADM_gall   ,ADM_kall,ADM_lall   )
  real(RP), public              :: Ml_pl           (ADM_gall_pl,ADM_kall,ADM_lall_pl)
  real(RP), public              :: Mu              (ADM_gall   ,ADM_kall,ADM_lall   )
  real(RP), public              :: Mu_pl           (ADM_gall_pl,ADM_kall,ADM_lall_pl)
  real(RP), public              :: VMTR_RGSGAM2    (ADM_gall   ,ADM_kall,ADM_lall   )
  real(RP), public              :: VMTR_RGSGAM2_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
  real(RP), public              :: VMTR_RGSGAM2H   (ADM_gall   ,ADM_kall,ADM_lall   )
  real(RP), public              :: VMTR_RGSGAM2H_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
  real(RP), public              :: VMTR_RGAMH      (ADM_gall   ,ADM_kall,ADM_lall   )
  real(RP), public              :: VMTR_RGAMH_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
  real(RP), public              :: VMTR_RGAM       (ADM_gall   ,ADM_kall,ADM_lall   )
  real(RP), public              :: VMTR_RGAM_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
  real(RP), public              :: VMTR_GSGAM2H    (ADM_gall   ,ADM_kall,ADM_lall   )
  real(RP), public              :: VMTR_GSGAM2H_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
#else
  real(RP), public, allocatable :: Mc              (:,:,:)
  real(RP), public, allocatable :: Mc_pl           (:,:,:)
  real(RP), public, allocatable :: Ml              (:,:,:)
  real(RP), public, allocatable :: Ml_pl           (:,:,:)
  real(RP), public, allocatable :: Mu              (:,:,:)
  real(RP), public, allocatable :: Mu_pl           (:,:,:)
  real(RP), public, allocatable :: VMTR_RGSGAM2    (:,:,:)
  real(RP), public, allocatable :: VMTR_RGSGAM2_pl (:,:,:)
  real(RP), public, allocatable :: VMTR_RGSGAM2H   (:,:,:)
  real(RP), public, allocatable :: VMTR_RGSGAM2H_pl(:,:,:)
  real(RP), public, allocatable :: VMTR_RGAMH      (:,:,:)
  real(RP), public, allocatable :: VMTR_RGAMH_pl   (:,:,:)
  real(RP), public, allocatable :: VMTR_RGAM       (:,:,:)
  real(RP), public, allocatable :: VMTR_RGAM_pl    (:,:,:)
  real(RP), public, allocatable :: VMTR_GSGAM2H    (:,:,:)
  real(RP), public, allocatable :: VMTR_GSGAM2H_pl (:,:,:)
#endif

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine vi_rhow_solver_putcoef( &
       Mc_in           , Mc_pl_in           , &
       Ml_in           , Ml_pl_in           , &
       Mu_in           , Mu_pl_in           , &
       VMTR_RGSGAM2_in , VMTR_RGSGAM2_pl_in , &
       VMTR_RGSGAM2H_in, VMTR_RGSGAM2H_pl_in, &
       VMTR_RGAMH_in   , VMTR_RGAMH_pl_in   , &
       VMTR_RGAM_in    , VMTR_RGAM_pl_in    , &
       VMTR_GSGAM2H_in , VMTR_GSGAM2H_pl_in   )
    implicit none

    real(RP), intent(in)  :: Mc_in              (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: Mc_pl_in           (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: Ml_in              (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: Ml_pl_in           (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: Mu_in              (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: Mu_pl_in           (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: VMTR_RGSGAM2_in    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: VMTR_RGSGAM2_pl_in (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: VMTR_RGSGAM2H_in   (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: VMTR_RGSGAM2H_pl_in(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: VMTR_RGAMH_in      (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: VMTR_RGAMH_pl_in   (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: VMTR_RGAM_in       (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: VMTR_RGAM_pl_in    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)  :: VMTR_GSGAM2H_in    (ADM_gall   ,ADM_kall,ADM_lall   )
    real(RP), intent(in)  :: VMTR_GSGAM2H_pl_in (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    !---------------------------------------------------------------------------

#ifndef _FIXEDINDEX_
    allocate( Mc              (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( Mc_pl           (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( Ml              (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( Ml_pl           (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( Mu              (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( Mu_pl           (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( VMTR_RGSGAM2    (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( VMTR_RGSGAM2_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( VMTR_RGSGAM2H   (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( VMTR_RGSGAM2H_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( VMTR_RGAMH      (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( VMTR_RGAMH_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( VMTR_RGAM       (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( VMTR_RGAM_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
    allocate( VMTR_GSGAM2H    (ADM_gall   ,ADM_kall,ADM_lall   ) )
    allocate( VMTR_GSGAM2H_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
#endif

    Mc              (:,:,:) = Mc_in              (:,:,:)
    Mc_pl           (:,:,:) = Mc_pl_in           (:,:,:)
    Ml              (:,:,:) = Ml_in              (:,:,:)
    Ml_pl           (:,:,:) = Ml_pl_in           (:,:,:)
    Mu              (:,:,:) = Mu_in              (:,:,:)
    Mu_pl           (:,:,:) = Mu_pl_in           (:,:,:)
    VMTR_RGSGAM2    (:,:,:) = VMTR_RGSGAM2_in    (:,:,:)
    VMTR_RGSGAM2_pl (:,:,:) = VMTR_RGSGAM2_pl_in (:,:,:)
    VMTR_RGSGAM2H   (:,:,:) = VMTR_RGSGAM2H_in   (:,:,:)
    VMTR_RGSGAM2H_pl(:,:,:) = VMTR_RGSGAM2H_pl_in(:,:,:)
    VMTR_RGAMH      (:,:,:) = VMTR_RGAMH_in      (:,:,:)
    VMTR_RGAMH_pl   (:,:,:) = VMTR_RGAMH_pl_in   (:,:,:)
    VMTR_RGAM       (:,:,:) = VMTR_RGAM_in       (:,:,:)
    VMTR_RGAM_pl    (:,:,:) = VMTR_RGAM_pl_in    (:,:,:)
    VMTR_GSGAM2H    (:,:,:) = VMTR_GSGAM2H_in    (:,:,:)
    VMTR_GSGAM2H_pl (:,:,:) = VMTR_GSGAM2H_pl_in (:,:,:)

    return
  end subroutine vi_rhow_solver_putcoef

  !-----------------------------------------------------------------------------
  subroutine vi_rhow_solver( &
       rhogw,  rhogw_pl,  &
       rhogw0, rhogw0_pl, &
       preg0,  preg0_pl,  &
       rhog0,  rhog0_pl,  &
       Sr,     Sr_pl,     &
       Sw,     Sw_pl,     &
       Sp,     Sp_pl,     &
       dt                 )
!ESC!    use mod_adm, only: &
!ESC!       ADM_have_pl, &
!ESC!       ADM_gall,    &
!ESC!       ADM_gall_pl, &
!ESC!       ADM_lall,    &
!ESC!       ADM_lall_pl, &
!ESC!       ADM_kall,    &
!ESC!       ADM_kmin,    &
!ESC!       ADM_kmax
!ESC!    use mod_cnst, only: &
!ESC!       GRAV  => CNST_EGRAV, &
!ESC!       Rdry  => CNST_RAIR,  &
!ESC!       CVdry => CNST_CV
!ESC!    use mod_grd, only: &
!ESC!       GRD_rdgzh, &
!ESC!       GRD_afac,  &
!ESC!       GRD_bfac
!ESC!    use mod_vmtr, only: &
!ESC!       VMTR_RGSGAM2,     &
!ESC!       VMTR_RGSGAM2_pl,  &
!ESC!       VMTR_RGSGAM2H,    &
!ESC!       VMTR_RGSGAM2H_pl, &
!ESC!       VMTR_RGAMH,       &
!ESC!       VMTR_RGAMH_pl,    &
!ESC!       VMTR_RGAM,        &
!ESC!       VMTR_RGAM_pl,     &
!ESC!       VMTR_GSGAM2H,     &
!ESC!       VMTR_GSGAM2H_pl
!ESC!    use mod_runconf, only: &
!ESC!       NON_HYDRO_ALPHA
    implicit none

    real(RP), intent(inout) :: rhogw    (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*w          ( G^1/2 x gam2 ), n+1
    real(RP), intent(inout) :: rhogw_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)

    real(RP), intent(in)    :: rhogw0   (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho*w          ( G^1/2 x gam2 )
    real(RP), intent(in)    :: rhogw0_pl(ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: preg0    (ADM_gall   ,ADM_kall,ADM_lall   ) ! pressure prime ( G^1/2 x gam2 )
    real(RP), intent(in)    :: preg0_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: rhog0    (ADM_gall   ,ADM_kall,ADM_lall   ) ! rho            ( G^1/2 x gam2 )
    real(RP), intent(in)    :: rhog0_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: Sr       (ADM_gall   ,ADM_kall,ADM_lall   ) ! source term for rho  at the full level
    real(RP), intent(in)    :: Sr_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: Sw       (ADM_gall   ,ADM_kall,ADM_lall   ) ! source term for rhow at the half level
    real(RP), intent(in)    :: Sw_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: Sp       (ADM_gall   ,ADM_kall,ADM_lall   ) ! source term for pres at the full level
    real(RP), intent(in)    :: Sp_pl    (ADM_gall_pl,ADM_kall,ADM_lall_pl)
    real(RP), intent(in)    :: dt

    real(RP) :: Sall    (ADM_gall,   ADM_kall)
    real(RP) :: Sall_pl (ADM_gall_pl,ADM_kall)
    real(RP) :: beta    (ADM_gall   )
    real(RP) :: beta_pl (ADM_gall_pl)
    real(RP) :: gamma   (ADM_gall,   ADM_kall)
    real(RP) :: gamma_pl(ADM_gall_pl,ADM_kall)

    real(RP) :: alfa
    real(RP) :: CVovRt2 ! Cv / R / dt**2

    integer :: g, k, l
    !---------------------------------------------------------------------------

    call DEBUG_rapstart('____vi_rhow_solver')

    alfa = real(NON_HYDRO_ALPHA,kind=RP)
    CVovRt2 = CVdry / Rdry / (dt*dt)

    do l = 1, ADM_lall
       ! calc Sall
       do k  = ADM_kmin+1, ADM_kmax
       do g = 1, ADM_gall
          Sall(g,k) = (   ( rhogw0(g,k,  l)*alfa + dt * Sw(g,k,  l) ) * VMTR_RGAMH  (g,k,  l)**2            &
                      - ( ( preg0 (g,k,  l)      + dt * Sp(g,k,  l) ) * VMTR_RGSGAM2(g,k,  l)               &
                        - ( preg0 (g,k-1,l)      + dt * Sp(g,k-1,l) ) * VMTR_RGSGAM2(g,k-1,l)               &
                        ) * dt * GRD_rdgzh(k)                                                               &
                      - ( ( rhog0 (g,k,  l)      + dt * Sr(g,k,  l) ) * VMTR_RGAM(g,k,  l)**2 * GRD_afac(k) &
                        + ( rhog0 (g,k-1,l)      + dt * Sr(g,k-1,l) ) * VMTR_RGAM(g,k-1,l)**2 * GRD_bfac(k) &
                        ) * dt * 0.5_RP * GRAV                                                               &
                      ) * CVovRt2
       enddo
       enddo

       ! boundary conditions
       do g = 1, ADM_gall
          rhogw(g,ADM_kmin,  l) = rhogw(g,ADM_kmin,  l) * VMTR_RGSGAM2H(g,ADM_kmin,  l)
          rhogw(g,ADM_kmax+1,l) = rhogw(g,ADM_kmax+1,l) * VMTR_RGSGAM2H(g,ADM_kmax+1,l)
          Sall (g,ADM_kmin+1)   = Sall (g,ADM_kmin+1) - Ml(g,ADM_kmin+1,l) * rhogw(g,ADM_kmin,  l)
          Sall (g,ADM_kmax  )   = Sall (g,ADM_kmax  ) - Mu(g,ADM_kmax,  l) * rhogw(g,ADM_kmax+1,l)
       enddo

       !---< solve tri-daigonal matrix >

       ! condition at ADM_kmin+1
       k = ADM_kmin+1
       do g = 1, ADM_gall
          beta (g)     = Mc(g,k,l)
          rhogw(g,k,l) = Sall(g,k) / beta(g)
       enddo

       ! forward
       do k = ADM_kmin+2, ADM_kmax
       do g = 1, ADM_gall
          gamma(g,k)   = Mu(g,k-1,l) / beta(g)
          beta (g)     = Mc(g,k,l) - Ml(g,k,l) * gamma(g,k) ! update beta
          rhogw(g,k,l) = ( Sall(g,k) - Ml(g,k,l) * rhogw(g,k-1,l) ) / beta(g)
       enddo
       enddo

       ! backward
       do k = ADM_kmax-1, ADM_kmin+1, -1
       do g = 1, ADM_gall
          rhogw(g,k  ,l) = rhogw(g,k  ,l) - gamma(g,k+1) * rhogw(g,k+1,l)
          rhogw(g,k+1,l) = rhogw(g,k+1,l) * VMTR_GSGAM2H(g,k+1,l) ! return value ( G^1/2 x gam2 )
       enddo
       enddo

       ! boundary treatment
       do g = 1, ADM_gall
          rhogw(g,ADM_kmin  ,l) = rhogw(g,ADM_kmin  ,l) * VMTR_GSGAM2H(g,ADM_kmin  ,l)
          rhogw(g,ADM_kmin+1,l) = rhogw(g,ADM_kmin+1,l) * VMTR_GSGAM2H(g,ADM_kmin+1,l)
          rhogw(g,ADM_kmax+1,l) = rhogw(g,ADM_kmax+1,l) * VMTR_GSGAM2H(g,ADM_kmax+1,l)
       enddo
    enddo

    if ( ADM_have_pl ) then
       do l = 1, ADM_lall_pl
          do k  = ADM_kmin+1, ADM_kmax
          do g = 1, ADM_gall_pl
             Sall_pl(g,k) = (   ( rhogw0_pl(g,k,  l)*alfa + dt * Sw_pl(g,k,  l) ) * VMTR_RGAMH_pl  (g,k,  l)**2            &
                            - ( ( preg0_pl (g,k,  l)      + dt * Sp_pl(g,k,  l) ) * VMTR_RGSGAM2_pl(g,k,  l)               &
                              - ( preg0_pl (g,k-1,l)      + dt * Sp_pl(g,k-1,l) ) * VMTR_RGSGAM2_pl(g,k-1,l)               &
                              ) * dt * GRD_rdgzh(k)                                                                        &
                            - ( ( rhog0_pl (g,k,  l)      + dt * Sr_pl(g,k,  l) ) * VMTR_RGAM_pl(g,k,  l)**2 * GRD_afac(k) &
                              + ( rhog0_pl (g,k-1,l)      + dt * Sr_pl(g,k-1,l) ) * VMTR_RGAM_pl(g,k-1,l)**2 * GRD_bfac(k) &
                              ) * dt * 0.5_RP * GRAV                                                                        &
                            ) * CVovRt2
          enddo
          enddo

          do g = 1, ADM_gall_pl
             rhogw_pl(g,ADM_kmin,  l) = rhogw_pl(g,ADM_kmin,  l) * VMTR_RGSGAM2H_pl(g,ADM_kmin,  l)
             rhogw_pl(g,ADM_kmax+1,l) = rhogw_pl(g,ADM_kmax+1,l) * VMTR_RGSGAM2H_pl(g,ADM_kmax+1,l)
             Sall_pl (g,ADM_kmin+1)   = Sall_pl (g,ADM_kmin+1) - Ml_pl(g,ADM_kmin+1,l) * rhogw_pl(g,ADM_kmin,  l)
             Sall_pl (g,ADM_kmax  )   = Sall_pl (g,ADM_kmax  ) - Mu_pl(g,ADM_kmax,  l) * rhogw_pl(g,ADM_kmax+1,l)
          enddo

          k = ADM_kmin+1
          do g = 1, ADM_gall_pl
             beta_pl (g)     = Mc_pl(g,k,l)
             rhogw_pl(g,k,l) = Sall_pl(g,k) / beta_pl(g)
          enddo

          do k = ADM_kmin+2, ADM_kmax
          do g = 1, ADM_gall_pl
             gamma_pl(g,k)   = Mu_pl(g,k-1,l) / beta_pl(g)
             beta_pl (g)     = Mc_pl(g,k,l) - Ml_pl(g,k,l) * gamma_pl(g,k) ! update beta
             rhogw_pl(g,k,l) = ( Sall_pl(g,k) - Ml_pl(g,k,l) * rhogw_pl(g,k-1,l) ) / beta_pl(g)
          enddo
          enddo

          do k = ADM_kmax-1, ADM_kmin+1, -1
          do g = 1, ADM_gall_pl
             rhogw_pl(g,k,l) = rhogw_pl(g,k,l) - gamma_pl(g,k+1) * rhogw_pl(g,k+1,l)
          enddo
          enddo

          do k = ADM_kmin, ADM_kmax+1
          do g = 1, ADM_gall_pl
             rhogw_pl(g,k,l) = rhogw_pl(g,k,l) * VMTR_GSGAM2H_pl(g,k,l)
          enddo
          enddo
       enddo
    endif

    call DEBUG_rapend('____vi_rhow_solver')

    return
  end subroutine vi_rhow_solver

end module mod_vi
!-------------------------------------------------------------------------------

