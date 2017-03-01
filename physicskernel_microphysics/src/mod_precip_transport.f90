!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Precipitation Transport
!!
!! @par Description
!!          Module for the precipitation transport
!!
!! @author NICAM developers
!!
!<
!-------------------------------------------------------------------------------
module mod_precip_transport
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_debug
!ESC!  use mod_stdio
!ESC!  use mod_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: precip_transport_new

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine precip_transport_new( &
       ijdim,              &
       rhog,               &
       rhogvx,             &
       rhogvy,             &
       rhogvz,             &
       rhogw,              &
       rhoge,              &
       rhogq,              &
       rho,                &
       tem,                &
       pre,                &
       vx,                 &
       vy,                 &
       vz,                 &
       w,                  &
       q,                  &
       qd,                 &
       z,                  &
       Vterm,              &
       precipitating_flag, &
       precip,             &
       precip_rhoe,        &
       precip_lh_heat,     &
       precip_rhophi,      &
       precip_rhokin,      &
       frain,              &
       gsgam2,             &
       gsgam2h,            &
       rgs,                &
       rgsh,               &
       ix,                 &
       iy,                 &
       iz,                 &
       jx,                 &
       jy,                 &
       jz,                 &
       dt,                 &
       precip_trc          )
!ESC!    use mod_adm, only: &
!ESC!       kdim => ADM_kall, &
!ESC!       kmin => ADM_kmin, &
!ESC!       kmax => ADM_kmax
!ESC!    use mod_const, only: &
!ESC!       CONST_GRAV
!ESC!    use mod_grd, only: &
!ESC!       GRD_gz ,   &
!ESC!       GRD_gzh,   &
!ESC!       GRD_dgz,   &
!ESC!       GRD_dgzh,  &
!ESC!       GRD_rdgz,  &
!ESC!       GRD_rdgzh, &
!ESC!       GRD_afact, &
!ESC!       GRD_bfact, &
!ESC!       GRD_cfact, &
!ESC!       GRD_dfact
!ESC!    use mod_runconf, only: &
!ESC!       nqmax => TRC_VMAX, &
!ESC!       NQW_STR,           &
!ESC!       NQW_END,           &
!ESC!       I_QC,              &
!ESC!       I_QR,              &
!ESC!       I_QI,              &
!ESC!       I_QS,              &
!ESC!       I_QG,              &
!ESC!       CVW,               &
!ESC!       LHF,               &
!ESC!       PRCIP_TRN_ECORRECT
!ESC!    use mod_cnvvar, only: &
!ESC!       cnvvar_rhogkin_in
    use mod_thrmdyn, only: &
       thrmdyn_qd, &
       thrmdyn_tempre
    use mod_vadv1d, only: &
       vadv1d_prep,       &
       vadv1d_getflux_new
    implicit none

    integer,  intent(in)    :: ijdim
    real(RP), intent(inout) :: rhog              (ijdim,kdim)
    real(RP), intent(inout) :: rhogvx            (ijdim,kdim)
    real(RP), intent(inout) :: rhogvy            (ijdim,kdim)
    real(RP), intent(inout) :: rhogvz            (ijdim,kdim)
    real(RP), intent(inout) :: rhogw             (ijdim,kdim)
    real(RP), intent(inout) :: rhoge             (ijdim,kdim)
    real(RP), intent(inout) :: rhogq             (ijdim,kdim,nqmax)
    real(RP), intent(inout) :: rho               (ijdim,kdim)
    real(RP), intent(inout) :: tem               (ijdim,kdim)
    real(RP), intent(inout) :: pre               (ijdim,kdim)
    real(RP), intent(in)    :: vx                (ijdim,kdim)
    real(RP), intent(in)    :: vy                (ijdim,kdim)
    real(RP), intent(in)    :: vz                (ijdim,kdim)
    real(RP), intent(in)    :: w                 (ijdim,kdim)
    real(RP), intent(inout) :: q                 (ijdim,kdim,nqmax)
    real(RP), intent(out)   :: qd                (ijdim,kdim)
    real(RP), intent(in)    :: z                 (ijdim,kdim)
    real(RP), intent(in)    :: Vterm             (ijdim,kdim,nqmax)
    logical,  intent(in)    :: precipitating_flag(nqmax)
    real(RP), intent(out)   :: precip            (ijdim,2)
    real(RP), intent(out)   :: precip_rhoe       (ijdim)
    real(RP), intent(out)   :: precip_lh_heat    (ijdim)
    real(RP), intent(out)   :: precip_rhophi     (ijdim)
    real(RP), intent(out)   :: precip_rhokin     (ijdim)
    real(RP), intent(out)   :: frain             (ijdim,kdim)
    real(RP), intent(in)    :: gsgam2            (ijdim,kdim)
    real(RP), intent(in)    :: gsgam2h           (ijdim,kdim)
    real(RP), intent(in)    :: rgs               (ijdim,kdim)
    real(RP), intent(in)    :: rgsh              (ijdim,kdim)
    real(RP), intent(in)    :: ix                (ijdim)
    real(RP), intent(in)    :: iy                (ijdim)
    real(RP), intent(in)    :: iz                (ijdim)
    real(RP), intent(in)    :: jx                (ijdim)
    real(RP), intent(in)    :: jy                (ijdim)
    real(RP), intent(in)    :: jz                (ijdim)
    real(RP), intent(in)    :: dt
    real(RP), intent(inout), optional :: precip_trc(ijdim,nqmax) ! [Add] 2012/02/01 T.Seiki

    real(RP) :: rhogkin       (ijdim,kdim)
    real(RP) :: rhogkin_h     (ijdim,kdim)
    real(RP) :: rhogkin_v     (ijdim,kdim)

    real(RP) :: rhoq          (ijdim,kdim)
    real(RP) :: rhoeq         (ijdim,kdim)
    real(RP) :: rhophiq       (ijdim,kdim)
    real(RP) :: rhokin_h      (ijdim,kdim)
    real(RP) :: rhokin_v      (ijdim,kdim)
    real(RP) :: rhouq         (ijdim,kdim)
    real(RP) :: rhovq         (ijdim,kdim)
    real(RP) :: rhowq         (ijdim,kdim)

    real(RP) :: fprec_q       (ijdim,kdim)
    real(RP) :: fprec_rhoe    (ijdim,kdim)
    real(RP) :: fprec_rhophi  (ijdim,kdim)
    real(RP) :: fprec_rhokin_h(ijdim,kdim)
    real(RP) :: fprec_rhokin_v(ijdim,kdim)
    real(RP) :: fprec_rhou    (ijdim,kdim)
    real(RP) :: fprec_rhov    (ijdim,kdim)
    real(RP) :: fprec_rhow    (ijdim,kdim)

    real(RP) :: drhoq         (ijdim,kdim,nqmax)
    real(RP) :: drhoe         (ijdim,kdim)
    real(RP) :: drhophi       (ijdim,kdim)
    real(RP) :: drhokin_h     (ijdim,kdim)
    real(RP) :: drhokin_v     (ijdim,kdim)
    real(RP) :: drhogu        (ijdim,kdim)
    real(RP) :: drhogv        (ijdim,kdim)
    real(RP) :: drhogw        (ijdim,kdim)

    real(RP) :: kin_h0        (ijdim,kdim)
    real(RP) :: kin_h         (ijdim,kdim)
    real(RP) :: vx_t          (ijdim,kdim)
    real(RP) :: vy_t          (ijdim,kdim)
    real(RP) :: vz_t          (ijdim,kdim)
    real(RP) :: kin_v0        (ijdim,kdim)
    real(RP) :: kin_v         (ijdim,kdim)
    real(RP) :: w_t           (ijdim,kdim)

    real(RP) :: zdis0         (ijdim,kdim)
    integer  :: kcell         (ijdim,kdim)
    integer  :: kcell_max     (kdim)
    integer  :: kcell_min     (kdim)
    real(RP) :: zdis0h        (ijdim,kdim)
    integer  :: kcellh        (ijdim,kdim)
    integer  :: kcellh_max    (kdim)
    integer  :: kcellh_min    (kdim)

    real(RP) :: Vtermh        (ijdim,kdim)
    real(RP) :: qh            (ijdim,kdim)
    real(RP) :: rhogh         (ijdim,kdim)
    real(RP) :: ein           (ijdim,kdim)
    real(RP) :: tmp           (ijdim,kdim)
    real(RP) :: tmp_h         (ijdim,kdim)
    real(RP) :: tmp_v         (ijdim,kdim)
    real(RP) :: tmp2          (ijdim,kdim)

    real(RP) :: C2Wfact       (ijdim,kdim,2)
    real(RP) :: W2Cfact       (ijdim,kdim,2)

    real(RP) :: GRD_gz_shift  (kdim)
    real(RP) :: GRAV

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    call PROF_rapstart('____Precipitation')

    GRAV = CONST_GRAV

    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kmin,kmax,C2Wfact,GRD_afact,GRD_bfact,GSGAM2,GSGAM2H)
    do k  = kmin, kmax+1
    do ij = 1, ijdim
       C2Wfact(ij,k,1) = GRD_afact(k) / GSGAM2(ij,k  ) * GSGAM2H(ij,k)
       C2Wfact(ij,k,2) = GRD_bfact(k) / GSGAM2(ij,k-1) * GSGAM2H(ij,k)
    enddo
    enddo
    !$omp end parallel do

    !$omp parallel do default(none),private(ij), &
    !$omp shared(ijdim,kmin,C2Wfact)
    do ij = 1, ijdim
       C2Wfact(ij,kmin-1,1) = 0.0_RP
       C2Wfact(ij,kmin-1,2) = 0.0_RP
    enddo
    !$omp end parallel do

    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kmin,kmax,W2Cfact,GRD_cfact,GRD_dfact,GSGAM2,GSGAM2H)
    do k  = kmin-1, kmax
    do ij = 1, ijdim
       W2Cfact(ij,k,1) = GRD_cfact(k) * GSGAM2(ij,k) / GSGAM2H(ij,k+1)
       W2Cfact(ij,k,2) = GRD_dfact(k) * GSGAM2(ij,k) / GSGAM2H(ij,k  )
    enddo
    enddo
    !$omp end parallel do

    !$omp parallel do default(none),private(ij), &
    !$omp shared(ijdim,kmax,W2Cfact)
    do ij = 1, ijdim
       W2Cfact(ij,kmax+1,1) = 0.0_RP
       W2Cfact(ij,kmax+1,2) = 0.0_RP
    enddo
    !$omp end parallel do

    call cnvvar_rhogkin_in( ijdim,            & ! [IN]
                            kdim,             & ! [IN]
                            rhog     (:,:),   & ! [IN]
                            rhogvx   (:,:),   & ! [IN]
                            rhogvy   (:,:),   & ! [IN]
                            rhogvz   (:,:),   & ! [IN]
                            rhogw    (:,:),   & ! [IN]
                            C2Wfact  (:,:,:), & ! [IN]
                            W2Cfact  (:,:,:), & ! [IN]
                            rhogkin  (:,:),   & ! [OUT]
                            rhogkin_h(:,:),   & ! [OUT]
                            rhogkin_v(:,:)    ) ! [OUT]

    do k  = kmin, kmax+1
       GRD_gz_shift(k) = GRD_gz(k-1)
    enddo

    !$omp parallel do default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,nqmax,drhoq)
    do nq = 1, nqmax
    do k  = 1, kdim
    do ij = 1, ijdim
       drhoq(ij,k,nq) = 0.0_RP
    enddo
    enddo
    enddo
    !$omp end parallel do

    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,drhoe,drhophi,drhokin_h,drhokin_v,drhogu,drhogv,drhogw)
    do k  = 1, kdim
    do ij = 1, ijdim
       drhoe    (ij,k) = 0.0_RP
       drhophi  (ij,k) = 0.0_RP
       drhokin_h(ij,k) = 0.0_RP
       drhokin_v(ij,k) = 0.0_RP
       drhogu   (ij,k) = 0.0_RP
       drhogv   (ij,k) = 0.0_RP
       drhogw   (ij,k) = 0.0_RP
    enddo
    enddo
    !$omp end parallel do

    !$omp parallel do default(none),private(ij), &
    !$omp shared(ijdim,precip,precip_rhoe,precip_lh_heat,precip_rhophi,precip_rhokin)
    do ij = 1, ijdim
       precip        (ij,1) = 0.0_RP
       precip        (ij,2) = 0.0_RP
       precip_rhoe   (ij)   = 0.0_RP
       precip_lh_heat(ij)   = 0.0_RP
       precip_rhophi (ij)   = 0.0_RP
       precip_rhokin (ij)   = 0.0_RP
    enddo
    !$omp end parallel do

    do nq = 1, nqmax

       if( .NOT. precipitating_flag(nq) ) cycle

       call vadv1d_prep( ijdim, kdim, kmin, kmax, & !--- [IN]
                         GRD_dgz  (:),            & !--- [IN]
                         GRD_gzh  (:),            & !--- [IN]
                         Vterm    (:,:,nq),       & !--- [IN]
                         zdis0    (:,:),          & !--- [OUT] [bugfix] H.Yashiro 20120606
                         kcell    (:,:),          & !--- [OUT]
                         kcell_max(:),            & !--- [OUT]
                         kcell_min(:),            & !--- [OUT]
                         dt                       ) !--- [IN]

       !----- mass
       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(nq,ijdim,kdim,rhoq,rhogq,rgs)
       do k  = 1, kdim
       do ij = 1, ijdim
          rhoq(ij,k) = rhogq(ij,k,nq) * rgs(ij,k)
       enddo
       enddo
       !$omp end parallel do

       call vadv1d_getflux_new( ijdim, kdim, kmin, kmax, & !--- [IN]
                                GRD_dgz  (:),            & !--- [IN]
                                rhoq     (:,:),          & !--- [IN]
                                zdis0    (:,:),          & !--- [IN] [bugfix] H.Yashiro 20120606
                                kcell    (:,:),          & !--- [IN]
                                kcell_max(:),            & !--- [IN]
                                kcell_min(:),            & !--- [IN]
                                fprec_q  (:,:)           ) !--- [OUT]

       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(nq,ijdim,kmin,kmax,drhoq,fprec_q,GRD_rdgz)
       do k  = kmin, kmax
       do ij = 1, ijdim
          drhoq(ij,k,nq) = -( fprec_q(ij,k+1)-fprec_q(ij,k) ) * GRD_rdgz(k)
       enddo
       enddo
       !$omp end parallel do

       ! Only for mass tracer
       if (       nq >= NQW_STR &
            .AND. nq <= NQW_END ) then

          !--- internal energy
          !--- potential energy
          !--- horizontal kinetic energy
          !--- momentum u
          !--- momentum v

          !$omp parallel do default(none),private(ij,k), &
          !$omp shared(nq,ijdim,kdim,rhoeq,rhophiq,rhokin_h,rhouq,rhovq,rhogq,q,rgs, &
          !$omp        tem,z,rhogkin_h,vx,vy,vz,ix,iy,iz,jx,jy,jz,CVW,GRAV)
          do k  = 1, kdim
          do ij = 1, ijdim
             rhoeq   (ij,k) = rhogq(ij,k,nq) * rgs(ij,k) * CVW(nq) * tem(ij,k)
             rhophiq (ij,k) = rhogq(ij,k,nq) * rgs(ij,k) * GRAV * z(ij,k)
             rhokin_h(ij,k) =     q(ij,k,nq) * rgs(ij,k) * rhogkin_h(ij,k)
             rhouq   (ij,k) = rhogq(ij,k,nq) * rgs(ij,k) * ( vx(ij,k) * ix(ij) &
                                                           + vy(ij,k) * iy(ij) &
                                                           + vz(ij,k) * iz(ij) )
             rhovq   (ij,k) = rhogq(ij,k,nq) * rgs(ij,k) * ( vx(ij,k) * jx(ij) &
                                                           + vy(ij,k) * jy(ij) &
                                                           + vz(ij,k) * jz(ij) )
          enddo
          enddo
          !$omp end parallel do

          call vadv1d_getflux_new( ijdim, kdim, kmin, kmax, &
                                   GRD_dgz  (:),            &
                                   rhoeq    (:,:),          &
                                   zdis0    (:,:),          &
                                   kcell    (:,:),          &
                                   kcell_max(:),            &
                                   kcell_min(:),            &
                                   fprec_rhoe(:,:)          )

          call vadv1d_getflux_new( ijdim, kdim, kmin, kmax, &
                                   GRD_dgz  (:),            &
                                   rhophiq  (:,:),          &
                                   zdis0    (:,:),          &
                                   kcell    (:,:),          &
                                   kcell_max(:),            &
                                   kcell_min(:),            &
                                   fprec_rhophi(:,:)        )

          call vadv1d_getflux_new( ijdim, kdim, kmin, kmax, &
                                   GRD_dgz  (:),            &
                                   rhokin_h (:,:),          &
                                   zdis0    (:,:),          &
                                   kcell    (:,:),          &
                                   kcell_max(:),            &
                                   kcell_min(:),            &
                                   fprec_rhokin_h(:,:)      )

          call vadv1d_getflux_new( ijdim, kdim, kmin, kmax, &
                                   GRD_dgz  (:),            &
                                   rhouq    (:,:),          &
                                   zdis0    (:,:),          &
                                   kcell    (:,:),          &
                                   kcell_max(:),            &
                                   kcell_min(:),            &
                                   fprec_rhou(:,:)          )

          call vadv1d_getflux_new( ijdim, kdim, kmin, kmax, &
                                   GRD_dgz  (:),            &
                                   rhovq    (:,:),          &
                                   zdis0    (:,:),          &
                                   kcell    (:,:),          &
                                   kcell_max(:),            &
                                   kcell_min(:),            &
                                   fprec_rhov(:,:)          )

          !$omp parallel default(none),private(ij,k), &
          !$omp shared(ijdim,kmin,kmax,drhoe,drhophi,drhokin_h,drhogu,drhogv,GRD_rdgz, &
          !$omp        fprec_rhoe,fprec_rhophi,fprec_rhokin_h,fprec_rhou,fprec_rhov)

          !$omp do
          do k  = kmin, kmax
          do ij = 1, ijdim
             drhoe    (ij,k) = drhoe    (ij,k) -( fprec_rhoe    (ij,k+1)-fprec_rhoe    (ij,k) ) * GRD_rdgz(k)
          enddo
          enddo
          !$omp end do nowait

          !$omp do
          do k  = kmin, kmax
          do ij = 1, ijdim
             drhophi  (ij,k) = drhophi  (ij,k) -( fprec_rhophi  (ij,k+1)-fprec_rhophi  (ij,k) ) * GRD_rdgz(k)
          enddo
          enddo
          !$omp end do nowait

          !$omp do
          do k  = kmin, kmax
          do ij = 1, ijdim
             drhokin_h(ij,k) = drhokin_h(ij,k) -( fprec_rhokin_h(ij,k+1)-fprec_rhokin_h(ij,k) ) * GRD_rdgz(k)
          enddo
          enddo
          !$omp end do nowait

          !$omp do
          do k  = kmin, kmax
          do ij = 1, ijdim
             drhogu   (ij,k) = drhogu   (ij,k) -( fprec_rhou    (ij,k+1)-fprec_rhou    (ij,k) ) * GRD_rdgz(k)
          enddo
          enddo
          !$omp end do nowait

          !$omp do
          do k  = kmin, kmax
          do ij = 1, ijdim
             drhogv   (ij,k) = drhogv   (ij,k) -( fprec_rhov    (ij,k+1)-fprec_rhov    (ij,k) ) * GRD_rdgz(k)
          enddo
          enddo
          !$omp end do

          !$omp end parallel

          !--- half level

          !$omp parallel default(none),private(ij,k), &
          !$omp shared(ijdim,kmin,kmax,Vtermh,Vterm,nq)

          !$omp do
          do k  = kmin+1, kmax-1
          do ij = 1,      ijdim
             Vtermh(ij,k) = 0.5_RP * ( Vterm(ij,k,nq) + Vterm(ij,k-1,nq) )
          enddo
          enddo
          !$omp end do nowait

          !$omp do
          do ij = 1, ijdim
             Vtermh(ij,kmin-1) = 0.0_RP
             Vtermh(ij,kmin  ) = 0.0_RP
             Vtermh(ij,kmax  ) = 0.0_RP
             Vtermh(ij,kmax+1) = 0.0_RP
          enddo
          !$omp end do

          !$omp end parallel

          call vadv1d_prep( ijdim, kdim, kmin+1, kmax, &
                            GRD_dgzh    (:),           &
                            GRD_gz_shift(:),           &
                            Vtermh      (:,:),         &
                            zdis0h      (:,:),         & ! [bugfix] H.Yashiro 20120606
                            kcellh      (:,:),         &
                            kcellh_max  (:),           &
                            kcellh_min  (:),           &
                            dt                         )

          !$omp parallel default(none),private(ij,k), &
          !$omp shared(nq,ijdim,kmin,kmax,qh,q,rgsh)

          !$omp do
          do k  = kmin+1, kmax
          do ij = 1, ijdim
             qh(ij,k) = 0.5_RP * ( q(ij,k,nq) + q(ij,k-1,nq) ) * rgsh(ij,k)
          enddo
          enddo
          !$omp end do nowait

          !$omp do
          do ij = 1, ijdim
             qh(ij,kmin-1) = 0.0_RP
             qh(ij,kmin  ) = 0.0_RP
             qh(ij,kmax+1) = 0.0_RP
          enddo
          !$omp end do

          !$omp end parallel

          !--- vertical kinetic energy
          !--- moment w
          !--- half level

          !$omp parallel do default(none),private(ij,k), &
          !$omp shared(nq,ijdim,kmin,kmax,rhokin_v,rhowq,qh,rhogkin_v,rhogw)
          do k  = kmin+1, kmax
          do ij = 1, ijdim
             rhokin_v(ij,k) = qh(ij,k) * rhogkin_v(ij,k)
             rhowq   (ij,k) = qh(ij,k) * rhogw    (ij,k)
          enddo
          enddo
          !$omp end parallel do

          call vadv1d_getflux_new( ijdim, kdim, kmin+1, kmax, &
                                   GRD_dgzh  (:),             &
                                   rhokin_v  (:,:),           &
                                   zdis0h    (:,:),           & ! [bugfix] H.Yashiro 20120606
                                   kcellh    (:,:),           &
                                   kcellh_max(:),             &
                                   kcellh_min(:),             &
                                   fprec_rhokin_v(:,:)        )

          call vadv1d_getflux_new( ijdim, kdim, kmin+1, kmax, &
                                   GRD_dgzh  (:),             &
                                   rhowq     (:,:),           &
                                   zdis0h    (:,:),           &
                                   kcellh    (:,:),           &
                                   kcellh_max(:),             &
                                   kcellh_min(:),             &
                                   fprec_rhow(:,:)            )

          !$omp parallel do default(none),private(ij,k), &
          !$omp shared(ijdim,kmin,kmax,drhokin_v,drhogw,fprec_rhokin_v,fprec_rhow,GRD_rdgzh)
          do k  = kmin+1, kmax
          do ij = 1, ijdim
             drhokin_v(ij,k) = drhokin_v(ij,k) -( fprec_rhokin_v(ij,k)-fprec_rhokin_v(ij,k-1) ) * GRD_rdgzh(k)
             drhogw   (ij,k) = drhogw   (ij,k) -( fprec_rhow    (ij,k)-fprec_rhow    (ij,k-1) ) * GRD_rdgzh(k)
          enddo
          enddo
          !$omp end parallel do

          ! precipitation on the ground
          if ( nq == I_QC ) then

             !$omp parallel do default(none),private(ij), &
             !$omp shared(ijdim,kmin,precip,fprec_q,dt)
             do ij = 1, ijdim
                precip(ij,1) = precip(ij,1) - fprec_q(ij,kmin) / dt
             enddo
             !$omp end parallel do

          elseif( nq == I_QR ) then

             !$omp parallel do default(none),private(ij), &
             !$omp shared(ijdim,kmin,precip,fprec_q,dt)
             do ij = 1, ijdim
                precip(ij,1) = precip(ij,1) - fprec_q(ij,kmin) / dt
             enddo
             !$omp end parallel do

             !$omp parallel default(none),private(ij,k), &
             !$omp shared(ijdim,kmin,kmax,frain,fprec_q,dt)

             !$omp do
             do k  = kmin, kmax
             do ij = 1, ijdim
                frain(ij,k) = -fprec_q(ij,k) / dt
             enddo
             enddo
             !$omp end do nowait

             !$omp do
             do ij = 1, ijdim
                frain(ij,kmin-1) = 0.0_RP
                frain(ij,kmax+1) = 0.0_RP
             enddo
             !$omp end do

             !$omp end parallel

          elseif( nq == I_QI .OR. nq == I_QS .OR. nq == I_QG ) then

             !$omp parallel do default(none),private(ij), &
             !$omp shared(ijdim,kmin,precip,fprec_q,dt,precip_lh_heat,LHF)
             do ij = 1, ijdim
                precip(ij,2) = precip(ij,2) - fprec_q(ij,kmin) / dt
                precip_lh_heat(ij) = precip_lh_heat(ij) + fprec_q(ij,kmin) * LHF / dt
             enddo
             !$omp end parallel do

          endif

          if ( present(precip_trc) ) then ! [Add] 2012/02/01 T.Seiki

             !$omp parallel do default(none),private(ij), &
             !$omp shared(nq,ijdim,kmin,precip_trc,fprec_q,dt)
             do ij = 1, ijdim
                precip_trc(ij,nq) = precip_trc(ij,nq) - fprec_q(ij,kmin) / dt
             enddo
             !$omp end parallel do

          endif

          !$omp parallel do default(none),private(ij), &
          !$omp shared(ijdim,kmin,precip_rhoe,precip_rhophi,precip_rhokin,     &
          !$omp        fprec_rhoe,fprec_rhophi,fprec_rhokin_h,fprec_rhokin_v,dt)
          do ij = 1, ijdim
             precip_rhoe  (ij) = precip_rhoe  (ij) - fprec_rhoe    (ij,kmin) / dt
             precip_rhophi(ij) = precip_rhophi(ij) - fprec_rhophi  (ij,kmin) / dt
             precip_rhokin(ij) = precip_rhokin(ij) - fprec_rhokin_h(ij,kmin) / dt &
                                                   - fprec_rhokin_v(ij,kmin) / dt
          enddo
          !$omp end parallel do

       endif

    enddo ! tracer LOOP

    ! Change in internal energy comes from precipitation and dissipation of kinetic energy due to drag force.
    ! See Ooyama(2001) (3.13)

    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,rhoge,rhogkin_h,rhogkin_v,rhogvx,rhogvy,rhogvz,rhogw, &
    !$omp        drhoe,drhophi,drhokin_h,drhokin_v,drhogu,drhogv,drhogw,ix,iy,iz,jx,jy,jz)
    do k  = 1, kdim
    do ij = 1, ijdim
       rhoge    (ij,k) = rhoge    (ij,k) + drhoe    (ij,k) + drhophi(ij,k)
       rhogkin_h(ij,k) = rhogkin_h(ij,k) + drhokin_h(ij,k)
       rhogkin_v(ij,k) = rhogkin_v(ij,k) + drhokin_v(ij,k)
       rhogvx   (ij,k) = rhogvx   (ij,k) + drhogu   (ij,k) * ix(ij) + drhogv(ij,k) * jx(ij)
       rhogvy   (ij,k) = rhogvy   (ij,k) + drhogu   (ij,k) * iy(ij) + drhogv(ij,k) * jy(ij)
       rhogvz   (ij,k) = rhogvz   (ij,k) + drhogu   (ij,k) * iz(ij) + drhogv(ij,k) * jz(ij)
       rhogw    (ij,k) = rhogw    (ij,k) + drhogw   (ij,k)
    enddo
    enddo
    !$omp end parallel do

    do nq = 1, nqmax
       if ( nq >= NQW_STR .AND. nq <= NQW_END ) then
          !$omp parallel do default(none),private(ij,k), &
          !$omp shared(nq,ijdim,kdim,rhogq,rhog,rhoge,drhoq,z,GRAV)
          do k  = 1, kdim
          do ij = 1, ijdim
             rhogq(ij,k,nq) = rhogq(ij,k,nq) + drhoq(ij,k,nq)
             rhog (ij,k)    = rhog (ij,k)    + drhoq(ij,k,nq)
             rhoge(ij,k)    = rhoge(ij,k)    - drhoq(ij,k,nq) * GRAV * z(ij,k)
          enddo
          enddo
          !$omp end parallel do
       else
          !$omp parallel do default(none),private(ij,k), &
          !$omp shared(nq,ijdim,kdim,rhogq,drhoq)
          do k  = 1, kdim
          do ij = 1, ijdim
             rhogq(ij,k,nq) = rhogq(ij,k,nq) + drhoq(ij,k,nq)
          enddo
          enddo
          !$omp end parallel do
       endif
    enddo

    if ( PRCIP_TRN_ECORRECT == 'KIN2EIN' ) then

       !$omp parallel default(none),private(ij,k), &
       !$omp shared(ijdim,kmin,kmax,tmp2,rhogkin_h,rhogkin_v,W2Cfact)

       !$omp do
       do k  = kmin, kmax
       do ij = 1, ijdim
          tmp2(ij,k) = rhogkin_h(ij,k) + ( W2Cfact(ij,k,1) * rhogkin_v(ij,k+1) &
                                         + W2Cfact(ij,k,2) * rhogkin_v(ij,k  ) )
       enddo
       enddo
       !$omp end do

       !$omp do
       do ij = 1, ijdim
          tmp2(ij,kmin-1) = 0.0_RP
          tmp2(ij,kmax+1) = 0.0_RP
       enddo
       !$omp end do

       !$omp end parallel

       call cnvvar_rhogkin_in( ijdim,          & ! [IN]
                               kdim,           & ! [IN]
                               rhog   (:,:),   & ! [IN]
                               rhogvx (:,:),   & ! [IN]
                               rhogvy (:,:),   & ! [IN]
                               rhogvz (:,:),   & ! [IN]
                               rhogw  (:,:),   & ! [IN]
                               C2Wfact(:,:,:), & ! [IN]
                               W2Cfact(:,:,:), & ! [IN]
                               tmp    (:,:),   & ! [OUT]
                               tmp_h  (:,:),   & ! [OUT]
                               tmp_v  (:,:)    ) ! [OUT]

       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(nq,ijdim,kmin,kmax,rhoge,tmp2,tmp)
       do k  = kmin, kmax
       do ij = 1, ijdim
          rhoge(ij,k) = rhoge(ij,k) + ( tmp2(ij,k) - tmp(ij,k) )
       enddo
       enddo
       !$omp end parallel do

    elseif( PRCIP_TRN_ECORRECT == 'KIN2KIN' ) then

       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kmin,kmax,kin_h0,vx_t,vy_t,vz_t,kin_h, &
       !$omp        rhogkin_h,rhogvx,rhogvy,rhogvz,rhog)
       do k  = kmin, kmax
       do ij = 1, ijdim
          kin_h0(ij,k) = rhogkin_h(ij,k) / rhog(ij,k)
          vx_t  (ij,k) = rhogvx   (ij,k) / rhog(ij,k)
          vy_t  (ij,k) = rhogvy   (ij,k) / rhog(ij,k)
          vz_t  (ij,k) = rhogvz   (ij,k) / rhog(ij,k)

          kin_h (ij,k) = 0.5_RP * ( vx_t(ij,k)*vx_t(ij,k) &
                                  + vy_t(ij,k)*vy_t(ij,k) &
                                  + vz_t(ij,k)*vz_t(ij,k) )
       enddo
       enddo
       !$omp end parallel do

       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kmin,kmax,vx_t,vy_t,vz_t,kin_h0,kin_h)
       do k  = kmin, kmax
       do ij = 1, ijdim
          if ( kin_h(ij,k) > 1.E-20_RP ) then
             vx_t(ij,k) = vx_t(ij,k) * sqrt( abs( kin_h0(ij,k) / kin_h(ij,k) ) )
             vy_t(ij,k) = vy_t(ij,k) * sqrt( abs( kin_h0(ij,k) / kin_h(ij,k) ) )
             vz_t(ij,k) = vz_t(ij,k) * sqrt( abs( kin_h0(ij,k) / kin_h(ij,k) ) )
          endif
       enddo
       enddo
       !$omp end parallel do

       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kmin,kmax,rhogvx,rhogvy,rhogvz,vx_t,vy_t,vz_t,rhog)
       do k  = kmin, kmax
       do ij = 1, ijdim
          rhogvx(ij,k) = vx_t(ij,k) * rhog(ij,k)
          rhogvy(ij,k) = vy_t(ij,k) * rhog(ij,k)
          rhogvz(ij,k) = vz_t(ij,k) * rhog(ij,k)
       enddo
       enddo
       !$omp end parallel do

       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kmin,kmax,rhogh,rhog,C2Wfact)
       do k  = kmin, kmax+1
       do ij = 1, ijdim
          rhogh(ij,k) = ( C2Wfact(ij,k,1) * rhog(ij,k  ) &
                        + C2Wfact(ij,k,2) * rhog(ij,k-1) )
       enddo
       enddo
       !$omp end parallel do

       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kmin,kmax,kin_v0,w_t,kin_v,rhogkin_v,rhogw,rhogh)
       do k  = kmin, kmax+1
       do ij = 1, ijdim
          kin_v0(ij,k) = rhogkin_v(ij,k) / rhogh(ij,k)
          w_t   (ij,k) = rhogw    (ij,k) / rhogh(ij,k)
          kin_v (ij,k) = 0.5_RP * w_t(ij,k)**2
       enddo
       enddo
       !$omp end parallel do

       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kmin,kmax,w_t,kin_v0,kin_v)
       do k  = kmin, kmax+1
       do ij = 1, ijdim
          if ( kin_v(ij,k) > 1.E-20_RP ) then
             w_t(ij,k) = w_t(ij,k) * sqrt( abs( kin_v0(ij,k) / kin_v(ij,k) ) )
          endif
       enddo
       enddo
       !$omp end parallel do

       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kmin,kmax,rhogw,w_t,rhogh)
       do k  = kmin, kmax+1
       do ij = 1, ijdim
          rhogw(ij,k) = w_t(ij,k) * rhogh(ij,k)
       enddo
       enddo
       !$omp end parallel do

    else
       write(*,*) 'Error in PRCIP_TRN_ECORRECT: ', trim(PRCIP_TRN_ECORRECT)
    endif

    !$omp parallel do default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kdim,nqmax,q,rhogq,rhog)
    do nq = 1, nqmax
    do k  = 1, kdim
    do ij = 1, ijdim
       q(ij,k,nq) = rhogq(ij,k,nq) / rhog(ij,k)
    enddo
    enddo
    enddo
    !$omp end parallel do

    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,rho,ein,rhog,rhoge,gsgam2)
    do k  = 1, kdim
    do ij = 1, ijdim
       rho(ij,k) = rhog (ij,k) / gsgam2(ij,k)
       ein(ij,k) = rhoge(ij,k) / rhog(ij,k)
    enddo
    enddo
    !$omp end parallel do

    call thrmdyn_qd( ijdim,     & ! [IN]
                     kdim,      & ! [IN]
                     q (:,:,:), & ! [IN]
                     qd(:,:)    ) ! [OUT]

    call thrmdyn_tempre( ijdim,      & ! [IN]
                         kdim,       & ! [IN]
                         ein(:,:),   & ! [IN]
                         rho(:,:),   & ! [IN]
                         q  (:,:,:), & ! [IN]
                         tem(:,:),   & ! [OUT]
                         pre(:,:)    ) ! [OUT]

    call PROF_rapend  ('____Precipitation')

    return
  end subroutine precip_transport_new

end module mod_precip_transport
