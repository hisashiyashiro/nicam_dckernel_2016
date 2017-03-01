!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics driver
!!
!! @par Description
!!          driver of cloud microphysics schemes
!!
!! @author NICAM developers
!!
!<
!-------------------------------------------------------------------------------
module mod_mp_driver
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_debug
  use mod_precision
!ESC!  use mod_stdio
!ESC!  use mod_prof
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: mp_init
  public :: mp_driver
!ESC!  public :: mp_terminal_velocity
!ESC!  public :: mp_diag_volume
!ESC!  public :: mp_effective_radius

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
!ESC!  character(len=H_SHORT), private :: MP_TYPE = 'NONE'

  logical, private :: opt_radius_explicit = .false.
  logical, private :: opt_volume_explicit = .false.

  real(RP), private :: TSICE = 273.15_RP ! temperature of ice/liq = 0 (for tentative use)
  real(RP), private :: TWICE = 258.15_RP ! temperature of ice/liq = 1 (for tentative use)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> setup
  subroutine mp_init( MP_TYPE_in )
!ESC!    use mod_adm, only: &
!ESC!       ADM_proc_stop
!ESC!    use mod_mp_kessler, only: &
!ESC!       mp_kessler_init
!ESC!    use mod_mp_g98, only: &
!ESC!       mp_g98_init
!ESC!    use mod_mp_wsm6, only: &
!ESC!       mp_wsm6_init
!ESC!    use mod_mp_nsw5, only: &
!ESC!       mp_nsw5_init
    use mod_mp_nsw6, only: &
       mp_nsw6_init
!ESC!    use mod_mp_ndw6, only: &
!ESC!       mp_ndw6_init
!ESC!    use mod_mp_lin, only: &
!ESC!       mp_lin_init
!ESC!    use mod_mp_lsc, only: &
!ESC!       mp_lsc_init
    implicit none

    character(len=*), intent(in) :: MP_TYPE_in

    namelist / nm_mp_driver_init / &
       opt_radius_explicit, &
       opt_volume_explicit, &
       TSICE,               &
       TWICE

    integer  :: ierr
    !---------------------------------------------------------------------------

!ESC!    MP_TYPE = trim(MP_TYPE_in)

!ESC!    !--- read parameters
!ESC!    write(IO_FID_LOG,*)
!ESC!    write(IO_FID_LOG,*) '+++ Module[Cloud Microphysics Driver]/Category[nhm physics]'
!ESC!    rewind(IO_FID_CONF)
!ESC!    read(IO_FID_CONF,nml=nm_mp_driver_init,iostat=ierr)
!ESC!    if ( ierr < 0 ) then
!ESC!       write(IO_FID_LOG,*) '*** nm_mp_driver_init is not specified. use default.'
!ESC!    elseif( ierr > 0 ) then
!ESC!       write(*,          *) 'xxx Not appropriate names in namelist nm_mp_driver_init. STOP.'
!ESC!       write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist nm_mp_driver_init. STOP.'
!ESC!       call ADM_proc_stop
!ESC!    endif
!ESC!    write(IO_FID_LOG,nml=nm_mp_driver_init)

    select case(MP_TYPE)
!ESC!    case('NONE')
!ESC!       write(IO_FID_LOG,*) '*** no microphysics'
!ESC!    case('SAT_ONLY_REV')
!ESC!       write(IO_FID_LOG,*) '*** microphysics type = SAT_ONLY_REV'
!ESC!    case('SAT_ONLY_PSEUDO')
!ESC!       write(IO_FID_LOG,*) '*** microphysics type = SAT_ONLY_PSEUDO'
!ESC!    case('KESSLER')
!ESC!       write(IO_FID_LOG,*) '*** microphysics type = KESSLER'
!ESC!       call mp_kessler_init
!ESC!    case('G98')
!ESC!       write(IO_FID_LOG,*) '*** microphysics type = G98'
!ESC!       call mp_g98_init
!ESC!    case('WSM6')
!ESC!       write(IO_FID_LOG,*) '*** microphysics type = WSM6'
!ESC!       call mp_wsm6_init
!ESC!    case('NSW5')
!ESC!       write(IO_FID_LOG,*) '*** microphysics type = NSW5'
!ESC!       call mp_nsw5_init
    case('NSW6')
       write(IO_FID_LOG,*) '*** microphysics type = NSW6'
       call mp_nsw6_init
!ESC!    case('NDW6')
!ESC!       write(IO_FID_LOG,*) '*** microphysics type = NDW6'
!ESC!       call mp_ndw6_init
!ESC!    case('LIN')
!ESC!       write(IO_FID_LOG,*) '*** microphysics type = LIN'
!ESC!       call mp_lin_init
!ESC!    case('LSC')
!ESC!       write(IO_FID_LOG,*) '*** microphysics type = Large scale condensation'
!ESC!       call mp_lsc_init
    case default
       write(*,*) 'xxx Not appropriate type. MP_TYPE = ', trim(MP_TYPE)
       call ADM_proc_stop
    end select

    return
  end subroutine mp_init

  !-----------------------------------------------------------------------------
  subroutine mp_driver( &
       ijdim,          &
       l_region,       &
       rhog,           &
       rhogvx,         &
       rhogvy,         &
       rhogvz,         &
       rhogw,          &
       rhoge,          &
       rhogq,          &
       vx,             &
       vy,             &
       vz,             &
       w,              &
       unccn,          &
       rho,            &
       tem,            &
       pre,            &
       q,              &
       qd,             &
       precip,         &
       ISO1_precip,    &
       ISO2_precip,    &
       precip_rhoe,    &
       precip_lh_heat, &
       precip_rhophi,  &
       precip_rhokin,  &
       re_liquid,      &
       re_solid,       &
       re_cld,         &
       rctop,          &
       rwtop,          &
       tctop,          &
       frhoge_af,      &
       frhogqv_af,     &
       frhoge_rad,     &
       rhogqke,        &
       gsgam2,         &
       gsgam2h,        &
       gam2,           &
       gam2h,          &
       ix,             &
       iy,             &
       iz,             &
       jx,             &
       jy,             &
       jz,             &
       z,              &
       zh,             &
       dt,             &
       ct,             &
       GDCLW,          &
       GDCFRC,         &
       GPREC,          &
       CBMFX           )
!ESC!    use mod_adm, only: &
!ESC!       kdim => ADM_kall
!ESC!    use mod_const, only: &
!ESC!       CONST_UNDEF
!ESC!    use mod_runconf, only: &
!ESC!       nqmax => TRC_VMAX, &
!ESC!       I_QC, I_QR,        &
!ESC!       opt_2moment_water, &
!ESC!       MP_DIV_NUM,        &
!ESC!       ISOTOPE               ! [add] K.Yoshimura 20120414
!ESC!    use mod_satadjust, only: &
!ESC!       SATURATION_setrange, &
!ESC!       SATURATION_adjustment
!ESC!    use mod_mp_lsc, only: &
!ESC!       mp_lsc
!ESC!    use mod_mp_kessler, only: &
!ESC!       mp_kessler
!ESC!    use mod_mp_g98, only: &
!ESC!       mp_g98
!ESC!    use mod_mp_lin, only: &
!ESC!       mp_lin
!ESC!    use mod_mp_wsm6, only: &
!ESC!       mp_wsm6
!ESC!    use mod_mp_nsw5, only: &
!ESC!       mp_nsw5
    use mod_mp_nsw6, only: &
       mp_nsw6
!ESC!    use mod_mp_ndw6, only: &
!ESC!       mp_ndw6
!ESC!    use mod_history, only: &
!ESC!       history_in
    implicit none

    integer,  intent(in)    :: ijdim
    integer,  intent(in)    :: l_region
    real(RP), intent(inout) :: rhog          (ijdim,kdim)
    real(RP), intent(inout) :: rhogvx        (ijdim,kdim)
    real(RP), intent(inout) :: rhogvy        (ijdim,kdim)
    real(RP), intent(inout) :: rhogvz        (ijdim,kdim)
    real(RP), intent(inout) :: rhogw         (ijdim,kdim)
    real(RP), intent(inout) :: rhoge         (ijdim,kdim)
    real(RP), intent(inout) :: rhogq         (ijdim,kdim,nqmax)
    real(RP), intent(in)    :: vx            (ijdim,kdim)
    real(RP), intent(in)    :: vy            (ijdim,kdim)
    real(RP), intent(in)    :: vz            (ijdim,kdim)
    real(RP), intent(in)    :: w             (ijdim,kdim)
    real(RP), intent(in)    :: unccn         (ijdim,kdim) ! CCN
    real(RP), intent(inout) :: rho           (ijdim,kdim)
    real(RP), intent(inout) :: tem           (ijdim,kdim)
    real(RP), intent(inout) :: pre           (ijdim,kdim)
    real(RP), intent(inout) :: q             (ijdim,kdim,nqmax)
    real(RP), intent(out)   :: qd            (ijdim,kdim)
    real(RP), intent(inout) :: precip        (ijdim,2)
    real(RP), intent(inout) :: ISO1_precip   (ijdim,2)    ! [add] K.Yoshimura 20110414
    real(RP), intent(inout) :: ISO2_precip   (ijdim,2)    ! [add] K.Yoshimura 20110414
    real(RP), intent(inout) :: precip_rhoe   (ijdim)
    real(RP), intent(inout) :: precip_lh_heat(ijdim)
    real(RP), intent(inout) :: precip_rhophi (ijdim)
    real(RP), intent(inout) :: precip_rhokin (ijdim)
    real(RP), intent(out)   :: re_liquid     (ijdim,kdim) ! Effective Radius
    real(RP), intent(out)   :: re_solid      (ijdim,kdim) ! Effective Radius of solid
    real(RP), intent(out)   :: re_cld        (ijdim,kdim) ! Effective Radius
    real(RP), intent(out)   :: rctop         (ijdim,1)    ! Effective Radius of Cloud Top
    real(RP), intent(out)   :: rwtop         (ijdim,1)    ! Effective Radius of Warm-Cloud Top
    real(RP), intent(out)   :: tctop         (ijdim,1)    ! Cloud Top Temperature
    real(RP), intent(in)    :: frhoge_af     (ijdim,kdim)
    real(RP), intent(in)    :: frhogqv_af    (ijdim,kdim)
    real(RP), intent(in)    :: frhoge_rad    (ijdim,kdim) ! energy tendency by radiation
    real(RP), intent(in)    :: rhogqke       (ijdim,kdim) ! rhog*2*TKE
    real(RP), intent(in)    :: gsgam2        (ijdim,kdim)
    real(RP), intent(in)    :: gsgam2h       (ijdim,kdim)
    real(RP), intent(in)    :: gam2          (ijdim,kdim)
    real(RP), intent(in)    :: gam2h         (ijdim,kdim)
    real(RP), intent(in)    :: ix            (ijdim)
    real(RP), intent(in)    :: iy            (ijdim)
    real(RP), intent(in)    :: iz            (ijdim)
    real(RP), intent(in)    :: jx            (ijdim)
    real(RP), intent(in)    :: jy            (ijdim)
    real(RP), intent(in)    :: jz            (ijdim)
    real(RP), intent(in)    :: z             (ijdim,kdim)
    real(RP), intent(in)    :: zh            (ijdim,kdim)
    real(RP), intent(in)    :: dt
    real(RP), intent(in)    :: ct
    real(RP), intent(out)   :: GDCLW         (ijdim,kdim)
    real(RP), intent(out)   :: GDCFRC        (ijdim,kdim)
    real(RP), intent(inout) :: GPREC         (ijdim,kdim) ! rain flux
    real(RP), intent(in)    :: CBMFX         (ijdim,kdim)

    real(RP) :: precip_trc        (ijdim,nqmax)

    real(RP) :: precip_sum        (ijdim,2)
    real(RP) :: ISO1_precip_sum   (ijdim,2)
    real(RP) :: ISO2_precip_sum   (ijdim,2)
    real(RP) :: precip_rhoe_sum   (ijdim)
    real(RP) :: precip_lh_heat_sum(ijdim)
    real(RP) :: precip_rhophi_sum (ijdim)
    real(RP) :: precip_rhokin_sum (ijdim)
    real(RP) :: precip_trc_sum    (ijdim,nqmax)
    real(RP) :: GDCLW_sum         (ijdim,kdim)
    real(RP) :: GDCFRC_sum        (ijdim,kdim)
    real(RP) :: GPREC_sum         (ijdim,kdim)

    real(RP) :: re_rain           (ijdim,kdim) ! Effective Radius
    real(RP) :: re_ice            (ijdim,kdim) ! Effective Radius
    real(RP) :: re_snow           (ijdim,kdim) ! Effective Radius
    real(RP) :: re_graupel        (ijdim,kdim) ! Effective Radius
    real(RP) :: rctop_cld         (ijdim,1)    ! Effective Radius of Cloud Top
    real(RP) :: rwtop_cld         (ijdim,1)    ! Effective Radius of Warm-Cloud Top
    real(RP) :: tctop_cld         (ijdim,1)    ! Cloud Top Temperature

    real(RP) :: qke               (ijdim,kdim)

    real(RP) :: fraction_mp ! 1 / MP_DIV_NUM
    real(RP) :: dt_mp

    integer  :: k ,ij, m
    !---------------------------------------------------------------------------

    fraction_mp = 1.0_RP / real(MP_DIV_NUM,kind=8)
    dt_mp       = dt * fraction_mp

    !$omp parallel workshare
    precip_sum        (:,:) = 0.0_RP
    ISO1_precip_sum   (:,:) = 0.0_RP
    ISO2_precip_sum   (:,:) = 0.0_RP
    precip_rhoe_sum   (:)   = 0.0_RP
    precip_lh_heat_sum(:)   = 0.0_RP
    precip_rhophi_sum (:)   = 0.0_RP
    precip_rhokin_sum (:)   = 0.0_RP
    precip_trc_sum    (:,:) = 0.0_RP
    GPREC_sum         (:,:) = 0.0_RP
    GDCLW_sum         (:,:) = 0.0_RP
    GDCFRC_sum        (:,:) = 0.0_RP
    !$omp end parallel workshare

    if ( MP_TYPE /= 'NDW6') then
       !$omp parallel workshare
       re_rain   (:,:) = CONST_UNDEF
       re_ice    (:,:) = CONST_UNDEF
       re_snow   (:,:) = CONST_UNDEF
       re_graupel(:,:) = CONST_UNDEF
       !$omp end parallel workshare
    endif

    do m = 1, MP_DIV_NUM
       !$omp parallel workshare
       precip        (:,:) = 0.0_RP
       ISO1_precip   (:,:) = 0.0_RP
       ISO2_precip   (:,:) = 0.0_RP
       precip_rhoe   (:)   = 0.0_RP
       precip_lh_heat(:)   = 0.0_RP
       precip_rhophi (:)   = 0.0_RP
       precip_rhokin (:)   = 0.0_RP
       precip_trc    (:,:) = 0.0_RP
       GPREC         (:,:) = 0.0_RP
       GDCLW         (:,:) = 0.0_RP
       GDCFRC        (:,:) = 0.0_RP
       !$omp end parallel workshare

       if ( MP_TYPE == 'NONE' ) then

          !$omp parallel workshare
          re_liquid(:,:) = 10.E-6_RP ! not supported
          re_solid (:,:) = 20.E-6_RP ! not supported
          !$omp end parallel workshare

!ESC!       elseif( MP_TYPE == 'NDW6' ) then
!ESC!
!ESC!          !$omp parallel workshare
!ESC!          qke(:,:) = rhogqke(:,:) / rhog(:,:)
!ESC!          !$omp end parallel workshare
!ESC!
!ESC!          call mp_ndw6( ijdim,          & ! [IN]
!ESC!                        l_region,       & ! [IN]
!ESC!                        m, MP_DIV_NUM,  & ! [IN] [Add] 09/04/14 T.Mitsui
!ESC!                        rhog,           & ! [INOUT]
!ESC!                        rhogvx,         & ! [INOUT]
!ESC!                        rhogvy,         & ! [INOUT]
!ESC!                        rhogvz,         & ! [INOUT]
!ESC!                        rhogw,          & ! [INOUT]
!ESC!                        rhoge,          & ! [INOUT]
!ESC!                        rhogq,          & ! [INOUT]
!ESC!                        vx,             & ! [IN]
!ESC!                        vy,             & ! [IN]
!ESC!                        vz,             & ! [IN]
!ESC!                        w,              & ! [IN]
!ESC!                        rho,            & ! [INOUT]
!ESC!                        tem,            & ! [INOUT]
!ESC!                        pre,            & ! [INOUT]
!ESC!                        q,              & ! [INOUT]
!ESC!                        qd,             & ! [OUT]
!ESC!                        precip,         & ! [OUT]
!ESC!                        precip_rhoe,    & !---- OUT: [Add] 09/09/03 T.Mitsui
!ESC!                        precip_lh_heat, & ! [OUT] : [Add] 09/09/03 T.Mitsui
!ESC!                        precip_rhophi,  & ! [OUT] : [Add] 09/09/03 T.Mitsui
!ESC!                        precip_rhokin,  & ! [OUT] : [Add] 09/09/03 T.Mitsui
!ESC!                        precip_trc,     & ! [OUT] : [Add] 12/02/01 T.Seiki
!ESC!                        GPREC,          & ! [OUT] :
!ESC!                        re_liquid,      & ! [OUT] :
!ESC!                        re_solid,       & ! [OUT] : ! 08/05/30 [Add] T.Mitsui
!ESC!                        rctop,          & ! [OUT] :
!ESC!                        rwtop,          & ! [OUT] :
!ESC!                        tctop,          & ! [OUT] :
!ESC!                        re_cld,         & ! [OUT] :
!ESC!                        re_rain,        & ! [OUT] : ! 08/05/30 [Add] T.Mitsui
!ESC!                        re_ice,         & ! [OUT] : ! 08/05/30 [Add] T.Mitsui
!ESC!                        re_snow,        & ! [OUT] : ! 08/05/30 [Add] T.Mitsui
!ESC!                        re_graupel,     & ! [OUT] : ! 08/05/30 [Add] T.Mitsui
!ESC!                        rctop_cld,      & ! [OUT] :
!ESC!                        rwtop_cld,      & ! [OUT] :
!ESC!                        tctop_cld,      & ! [OUT] :
!ESC!                        frhoge_af,      & ! [IN]  : energy tendency by additional forcing [Add] 10/08/03 T.Mitsui
!ESC!                        frhogqv_af,     & ! [IN]  : vapor  tendency by additional forcing [Add] 10/08/03 T.Mitsui
!ESC!                        frhoge_rad,     & ! [IN]  : energy tendency by radiation [Add] 09/08/18 T.Mitsui
!ESC!                        qke,            & ! [IN]  : 2*TKE [Add] 09/08/18 T.Mitsui
!ESC!                        gsgam2,         & ! [IN]
!ESC!                        gsgam2h,        & ! [IN]
!ESC!                        gam2,           & ! [IN]
!ESC!                        gam2h,          & ! [IN]
!ESC!                        ix,             & ! [IN]
!ESC!                        iy,             & ! [IN]
!ESC!                        iz,             & ! [IN]
!ESC!                        jx,             & ! [IN]
!ESC!                        jy,             & ! [IN]
!ESC!                        jz,             & ! [IN]
!ESC!                        z,              & ! [IN]
!ESC!                        dt_mp,          & ! [IN]
!ESC!                        ct+dt_mp*m      ) ! [IN] 09/04/14 [Add] T.Mitsui
!ESC!
!ESC!       elseif( MP_TYPE == 'KESSLER' ) then
!ESC!
!ESC!          call mp_kessler( ijdim,                 & ! [IN]
!ESC!                           rhog          (:,:),   & ! [INOUT]
!ESC!                           rhogvx        (:,:),   & ! [INOUT]
!ESC!                           rhogvy        (:,:),   & ! [INOUT]
!ESC!                           rhogvz        (:,:),   & ! [INOUT]
!ESC!                           rhogw         (:,:),   & ! [INOUT]
!ESC!                           rhoge         (:,:),   & ! [INOUT]
!ESC!                           rhogq         (:,:,:), & ! [INOUT]
!ESC!                           vx            (:,:),   & ! [IN]
!ESC!                           vy            (:,:),   & ! [IN]
!ESC!                           vz            (:,:),   & ! [IN]
!ESC!                           w             (:,:),   & ! [IN]
!ESC!                           rho           (:,:),   & ! [INOUT]
!ESC!                           tem           (:,:),   & ! [INOUT]
!ESC!                           pre           (:,:),   & ! [INOUT]
!ESC!                           q             (:,:,:), & ! [INOUT]
!ESC!                           qd            (:,:),   & ! [OUT]
!ESC!                           precip        (:,1),   & ! [OUT]
!ESC!                           precip_rhoe   (:),     & ! [OUT]
!ESC!                           precip_lh_heat(:),     & ! [OUT]
!ESC!                           precip_rhophi (:),     & ! [OUT]
!ESC!                           precip_rhokin (:),     & ! [OUT]
!ESC!                           GPREC         (:,:),   & ! [OUT]
!ESC!                           gsgam2        (:,:),   & ! [IN]
!ESC!                           gsgam2h       (:,:),   & ! [IN]
!ESC!                           gam2          (:,:),   & ! [IN]
!ESC!                           gam2h         (:,:),   & ! [IN]
!ESC!                           ix            (:),     & ! [IN]
!ESC!                           iy            (:),     & ! [IN]
!ESC!                           iz            (:),     & ! [IN]
!ESC!                           jx            (:),     & ! [IN]
!ESC!                           jy            (:),     & ! [IN]
!ESC!                           jz            (:),     & ! [IN]
!ESC!                           z             (:,:),   & ! [IN]
!ESC!                           dt_mp                  ) ! [IN]
!ESC!
!ESC!          !$omp parallel workshare
!ESC!          precip(:,2) = 0.0_RP ! snow
!ESC!          !$omp end parallel workshare
!ESC!
!ESC!          call SATURATION_setrange( 10.0_RP, 5.0_RP ) ! [IN]
!ESC!
!ESC!          call SATURATION_adjustment( ijdim,             & ! [IN]
!ESC!                                      kdim,              & ! [IN]
!ESC!                                      rhog  (:,:),       & ! [IN]
!ESC!                                      rhoge (:,:),       & ! [INOUT]
!ESC!                                      rhogq (:,:,:),     & ! [INOUT]
!ESC!                                      tem   (:,:),       & ! [INOUT]
!ESC!                                      q     (:,:,:),     & ! [INOUT]
!ESC!                                      qd    (:,:),       & ! [IN]
!ESC!                                      gsgam2(:,:),       & ! [IN]
!ESC!                                      ice_adjust=.false. ) ! [IN]
!ESC!
!ESC!          !$omp parallel workshare
!ESC!          re_liquid(:,:) = 10.E-6_RP ! not supported
!ESC!          re_solid (:,:) = 20.E-6_RP ! not supported
!ESC!          re_cld   (:,:) = 10.E-6_RP ! not supported
!ESC!          rctop    (:,:) = 10.E-6_RP ! not supported
!ESC!          rwtop    (:,:) = 10.E-6_RP ! not supported
!ESC!          rctop_cld(:,:) = 10.E-6_RP ! not supported
!ESC!          rwtop_cld(:,:) = 10.E-6_RP ! not supported
!ESC!          !$omp end parallel workshare
!ESC!
!ESC!          do k  = 1, kdim
!ESC!          do ij = 1, ijdim
!ESC!             if ( q(ij,k,I_QC) + q(ij,k,I_QR) > 1.E-5_RP ) then
!ESC!                tctop    (ij,1) = tem(ij,k)
!ESC!                tctop_cld(ij,1) = tem(ij,k)
!ESC!             endif
!ESC!          enddo
!ESC!          enddo
!ESC!
!ESC!       elseif( MP_TYPE == 'G98' ) then
!ESC!
!ESC!          call mp_g98( ijdim,          & ! [IN]
!ESC!                       rhog,           & ! [INOUT]
!ESC!                       rhogvx,         & ! [INOUT]
!ESC!                       rhogvy,         & ! [INOUT]
!ESC!                       rhogvz,         & ! [INOUT]
!ESC!                       rhogw,          & ! [INOUT]
!ESC!                       rhoge,          & ! [INOUT]
!ESC!                       rhogq,          & ! [INOUT]
!ESC!                       vx,             & ! [IN]
!ESC!                       vy,             & ! [IN]
!ESC!                       vz,             & ! [IN]
!ESC!                       w,              & ! [IN]
!ESC!                       unccn,          & ! [IN]
!ESC!                       rho,            & ! [INOUT]
!ESC!                       tem,            & ! [INOUT]
!ESC!                       pre,            & ! [INOUT]
!ESC!                       q,              & ! [INOUT]
!ESC!                       qd,             & ! [OUT]
!ESC!                       precip(:,1),    & ! [OUT]
!ESC!                       precip_rhoe,    & ! [OUT]
!ESC!                       precip_lh_heat, & ! [OUT]
!ESC!                       precip_rhophi,  & ! [OUT]
!ESC!                       precip_rhokin,  & ! [OUT]
!ESC!                       GPREC,          & ! [OUT]
!ESC!                       re_liquid,      & ! [OUT]
!ESC!                       rctop,          & ! [OUT]
!ESC!                       rwtop,          & ! [OUT]
!ESC!                       tctop,          & ! [OUT]
!ESC!                       re_cld,         & ! [OUT]
!ESC!                       rctop_cld,      & ! [OUT]
!ESC!                       rwtop_cld,      & ! [OUT]
!ESC!                       tctop_cld,      & ! [OUT]
!ESC!                       gsgam2,         & ! [IN]
!ESC!                       gsgam2h,        & ! [IN]
!ESC!                       gam2,           & ! [IN]
!ESC!                       gam2h,          & ! [IN]
!ESC!                       ix,             & ! [IN]
!ESC!                       iy,             & ! [IN]
!ESC!                       iz,             & ! [IN]
!ESC!                       jx,             & ! [IN]
!ESC!                       jy,             & ! [IN]
!ESC!                       jz,             & ! [IN]
!ESC!                       z,              & ! [IN]
!ESC!                       dt_mp           ) ! [IN]
!ESC!
!ESC!          !$omp parallel workshare
!ESC!          precip  (:,2) = 0.0_RP    ! snow
!ESC!          re_solid(:,:) = 20.E-6_RP ! not supported
!ESC!          !$omp end parallel workshare
!ESC!
!ESC!          call SATURATION_setrange( 273.15_RP, 258.15_RP ) ! [IN]
!ESC!
!ESC!          call SATURATION_adjustment( ijdim,             & ! [IN]
!ESC!                                      kdim,              & ! [IN]
!ESC!                                      rhog  (:,:),       & ! [IN]
!ESC!                                      rhoge (:,:),       & ! [INOUT]
!ESC!                                      rhogq (:,:,:),     & ! [INOUT]
!ESC!                                      tem   (:,:),       & ! [INOUT]
!ESC!                                      q     (:,:,:),     & ! [INOUT]
!ESC!                                      qd    (:,:),       & ! [IN]
!ESC!                                      gsgam2(:,:),       & ! [IN]
!ESC!                                      ice_adjust=.false. ) ! [IN]
!ESC!
!ESC!       elseif( MP_TYPE == 'WSM6' ) then
!ESC!
!ESC!          call mp_wsm6( ijdim,          & ! [IN]
!ESC!                        rhog,           & ! [INOUT]
!ESC!                        rhogvx,         & ! [INOUT]
!ESC!                        rhogvy,         & ! [INOUT]
!ESC!                        rhogvz,         & ! [INOUT]
!ESC!                        rhogw,          & ! [INOUT]
!ESC!                        rhoge,          & ! [INOUT]
!ESC!                        rhogq,          & ! [INOUT]
!ESC!                        vx,             & ! [IN]
!ESC!                        vy,             & ! [IN]
!ESC!                        vz,             & ! [IN]
!ESC!                        w,              & ! [IN]
!ESC!                        unccn,          & ! [IN]
!ESC!                        rho,            & ! [INOUT]
!ESC!                        tem,            & ! [INOUT]
!ESC!                        pre,            & ! [INOUT]
!ESC!                        q,              & ! [INOUT]
!ESC!                        qd,             & ! [OUT]
!ESC!                        precip,         & ! [OUT]
!ESC!                        precip_rhoe,    & ! [OUT]
!ESC!                        precip_lh_heat, & ! [OUT]
!ESC!                        precip_rhophi,  & ! [OUT]
!ESC!                        precip_rhokin,  & ! [OUT]
!ESC!                        GPREC,          & ! [OUT]
!ESC!                        re_liquid,      & ! [OUT]
!ESC!                        rctop,          & ! [OUT]
!ESC!                        rwtop,          & ! [OUT]
!ESC!                        tctop,          & ! [OUT]
!ESC!                        re_cld,         & ! [OUT]
!ESC!                        rctop_cld,      & ! [OUT]
!ESC!                        rwtop_cld,      & ! [OUT]
!ESC!                        tctop_cld,      & ! [OUT]
!ESC!                        gsgam2,         & ! [IN]
!ESC!                        gsgam2h,        & ! [IN]
!ESC!                        gam2,           & ! [IN]
!ESC!                        gam2h,          & ! [IN]
!ESC!                        ix,             & ! [IN]
!ESC!                        iy,             & ! [IN]
!ESC!                        iz,             & ! [IN]
!ESC!                        jx,             & ! [IN]
!ESC!                        jy,             & ! [IN]
!ESC!                        jz,             & ! [IN]
!ESC!                        z,              & ! [IN]
!ESC!                        dt_mp           ) ! [IN]
!ESC!
!ESC!          !$omp parallel workshare
!ESC!          re_solid(:,:) = 20.E-6_RP ! not supported
!ESC!          !$omp end parallel workshare
!ESC!
!ESC!          call SATURATION_setrange( 100.0_RP, 90.0_RP ) ! [IN]
!ESC!
!ESC!          call SATURATION_adjustment( ijdim,             & ! [IN]
!ESC!                                      kdim,              & ! [IN]
!ESC!                                      rhog  (:,:),       & ! [IN]
!ESC!                                      rhoge (:,:),       & ! [INOUT]
!ESC!                                      rhogq (:,:,:),     & ! [INOUT]
!ESC!                                      tem   (:,:),       & ! [INOUT]
!ESC!                                      q     (:,:,:),     & ! [INOUT]
!ESC!                                      qd    (:,:),       & ! [IN]
!ESC!                                      gsgam2(:,:),       & ! [IN]
!ESC!                                      ice_adjust=.false. ) ! [IN]
!ESC!
!ESC!       elseif( MP_TYPE == 'NSW5' ) then
!ESC!
!ESC!          call mp_nsw5( ijdim,          & ! [IN]
!ESC!                        rhog,           & ! [INOUT]
!ESC!                        rhogvx,         & ! [INOUT]
!ESC!                        rhogvy,         & ! [INOUT]
!ESC!                        rhogvz,         & ! [INOUT]
!ESC!                        rhogw,          & ! [INOUT]
!ESC!                        rhoge,          & ! [INOUT]
!ESC!                        rhogq,          & ! [INOUT]
!ESC!                        vx,             & ! [IN]
!ESC!                        vy,             & ! [IN]
!ESC!                        vz,             & ! [IN]
!ESC!                        w,              & ! [IN]
!ESC!                        unccn,          & ! [IN]
!ESC!                        rho,            & ! [INOUT]
!ESC!                        tem,            & ! [INOUT]
!ESC!                        pre,            & ! [INOUT]
!ESC!                        q,              & ! [INOUT]
!ESC!                        qd,             & ! [OUT]
!ESC!                        precip,         & ! [OUT]
!ESC!                        precip_rhoe,    & ! [OUT]
!ESC!                        precip_lh_heat, & ! [OUT]
!ESC!                        precip_rhophi,  & ! [OUT]
!ESC!                        precip_rhokin,  & ! [OUT]
!ESC!                        GPREC,          & ! [OUT]
!ESC!                        re_liquid,      & ! [OUT]
!ESC!                        rctop,          & ! [OUT]
!ESC!                        rwtop,          & ! [OUT]
!ESC!                        tctop,          & ! [OUT]
!ESC!                        re_cld,         & ! [OUT]
!ESC!                        rctop_cld,      & ! [OUT]
!ESC!                        rwtop_cld,      & ! [OUT]
!ESC!                        tctop_cld,      & ! [OUT]
!ESC!                        gsgam2,         & ! [IN]
!ESC!                        gsgam2h,        & ! [IN]
!ESC!                        gam2,           & ! [IN]
!ESC!                        gam2h,          & ! [IN]
!ESC!                        ix,             & ! [IN]
!ESC!                        iy,             & ! [IN]
!ESC!                        iz,             & ! [IN]
!ESC!                        jx,             & ! [IN]
!ESC!                        jy,             & ! [IN]
!ESC!                        jz,             & ! [IN]
!ESC!                        z,              & ! [IN]
!ESC!                        dt_mp           ) ! [IN]
!ESC!
!ESC!          !$omp parallel workshare
!ESC!          re_solid(:,:) = 20.E-6_RP ! not supported
!ESC!          !$omp end parallel workshare
!ESC!
!ESC!          call SATURATION_setrange( 273.16D0, 233.16D0 ) ! [IN]
!ESC!
!ESC!          call SATURATION_adjustment( ijdim,             & ! [IN]
!ESC!                                      kdim,              & ! [IN]
!ESC!                                      rhog  (:,:),       & ! [IN]
!ESC!                                      rhoge (:,:),       & ! [INOUT]
!ESC!                                      rhogq (:,:,:),     & ! [INOUT]
!ESC!                                      tem   (:,:),       & ! [INOUT]
!ESC!                                      q     (:,:,:),     & ! [INOUT]
!ESC!                                      qd    (:,:),       & ! [IN]
!ESC!                                      gsgam2(:,:),       & ! [IN]
!ESC!                                      ice_adjust=.true.  ) ! [IN]
!ESC!
       elseif( MP_TYPE == 'NSW6' ) then

          call mp_nsw6( ijdim,                 & ! [IN]
                        l_region,              & ! [IN]
                        rhog          (:,:),   & ! [INOUT]
                        rhogvx        (:,:),   & ! [INOUT]
                        rhogvy        (:,:),   & ! [INOUT]
                        rhogvz        (:,:),   & ! [INOUT]
                        rhogw         (:,:),   & ! [INOUT]
                        rhoge         (:,:),   & ! [INOUT]
                        rhogq         (:,:,:), & ! [INOUT]
                        vx            (:,:),   & ! [IN]
                        vy            (:,:),   & ! [IN]
                        vz            (:,:),   & ! [IN]
                        w             (:,:),   & ! [IN]
                        unccn         (:,:),   & ! [IN]
                        rho           (:,:),   & ! [INOUT]
                        tem           (:,:),   & ! [INOUT]
                        pre           (:,:),   & ! [INOUT]
                        q             (:,:,:), & ! [INOUT]
                        qd            (:,:),   & ! [OUT]
                        precip        (:,:),   & ! [OUT]
                        precip_rhoe   (:),     & ! [OUT]
                        precip_lh_heat(:),     & ! [OUT]
                        precip_rhophi (:),     & ! [OUT]
                        precip_rhokin (:),     & ! [OUT]
                        GPREC         (:,:),   & ! [OUT]
                        re_liquid     (:,:),   & ! [OUT]
                        rctop         (:,:),   & ! [OUT]
                        rwtop         (:,:),   & ! [OUT]
                        tctop         (:,:),   & ! [OUT]
                        re_cld        (:,:),   & ! [OUT]
                        rctop_cld     (:,:),   & ! [OUT]
                        rwtop_cld     (:,:),   & ! [OUT]
                        tctop_cld     (:,:),   & ! [OUT]
                        gsgam2        (:,:),   & ! [IN]
                        gsgam2h       (:,:),   & ! [IN]
                        gam2          (:,:),   & ! [IN]
                        gam2h         (:,:),   & ! [IN]
                        ix            (:),     & ! [IN]
                        iy            (:),     & ! [IN]
                        iz            (:),     & ! [IN]
                        jx            (:),     & ! [IN]
                        jy            (:),     & ! [IN]
                        jz            (:),     & ! [IN]
                        z             (:,:),   & ! [IN]
                        dt_mp                  ) ! [IN]

          !$omp parallel workshare
          re_solid(:,:) = 20.E-6_RP ! not supported
          !$omp end parallel workshare

!ESC!       elseif( MP_TYPE == 'LIN' ) then
!ESC!
!ESC!          call mp_lin( ijdim,          & ! [IN]
!ESC!                       rhog,           & ! [INOUT]
!ESC!                       rhogvx,         & ! [INOUT]
!ESC!                       rhogvy,         & ! [INOUT]
!ESC!                       rhogvz,         & ! [INOUT]
!ESC!                       rhogw,          & ! [INOUT]
!ESC!                       rhoge,          & ! [INOUT]
!ESC!                       rhogq,          & ! [INOUT]
!ESC!                       vx,             & ! [IN]
!ESC!                       vy,             & ! [IN]
!ESC!                       vz,             & ! [IN]
!ESC!                       w,              & ! [IN]
!ESC!                       unccn,          & ! [IN]
!ESC!                       rho,            & ! [INOUT]
!ESC!                       tem,            & ! [INOUT]
!ESC!                       pre,            & ! [INOUT]
!ESC!                       q,              & ! [INOUT]
!ESC!                       qd,             & ! [OUT]
!ESC!                       precip,         & ! [OUT]
!ESC!                       precip_rhoe,    & ! [OUT]
!ESC!                       precip_lh_heat, & ! [OUT]
!ESC!                       precip_rhophi,  & ! [OUT]
!ESC!                       precip_rhokin,  & ! [OUT]
!ESC!                       GPREC,          & ! [OUT]
!ESC!                       re_liquid,      & ! [OUT]
!ESC!                       rctop,          & ! [OUT]
!ESC!                       rwtop,          & ! [OUT]
!ESC!                       tctop,          & ! [OUT]
!ESC!                       re_cld,         & ! [OUT]
!ESC!                       rctop_cld,      & ! [OUT]
!ESC!                       rwtop_cld,      & ! [OUT]
!ESC!                       tctop_cld,      & ! [OUT]
!ESC!                       gsgam2,         & ! [IN]
!ESC!                       gsgam2h,        & ! [IN]
!ESC!                       gam2,           & ! [IN]
!ESC!                       gam2h,          & ! [IN]
!ESC!                       ix,             & ! [IN]
!ESC!                       iy,             & ! [IN]
!ESC!                       iz,             & ! [IN]
!ESC!                       jx,             & ! [IN]
!ESC!                       jy,             & ! [IN]
!ESC!                       jz,             & ! [IN]
!ESC!                       z,              & ! [IN]
!ESC!                       dt_mp           ) ! [IN]
!ESC!
!ESC!          !$omp parallel workshare
!ESC!          re_solid(:,:) = 20.E-6_RP ! not supported
!ESC!          !$omp end parallel workshare
!ESC!
!ESC!       elseif( MP_TYPE == 'LSC' ) then
!ESC!
!ESC!          call mp_lsc( ijdim,              & ! [IN]
!ESC!                       rhog       (:,:),   & ! [INOUT]
!ESC!                       rhogvx     (:,:),   & ! [INOUT]
!ESC!                       rhogvy     (:,:),   & ! [INOUT]
!ESC!                       rhogvz     (:,:),   & ! [INOUT]
!ESC!                       rhogw      (:,:),   & ! [INOUT]
!ESC!                       rhoge      (:,:),   & ! [INOUT]
!ESC!                       rhogq      (:,:,:), & ! [INOUT]
!ESC!                       precip     (:,:),   & ! [OUT]
!ESC!                       ISO1_precip(:,:),   & ! [OUT]
!ESC!                       ISO2_precip(:,:),   & ! [OUT]
!ESC!                       precip_rhoe(:),     & ! [OUT]
!ESC!                       GDCLW      (:,:),   & ! [OUT]
!ESC!                       GDCFRC     (:,:),   & ! [OUT]
!ESC!                       GPREC      (:,:),   & ! [OUT]
!ESC!                       re_liquid  (:,:),   & ! [OUT]
!ESC!                       rctop      (:,:),   & ! [OUT]
!ESC!                       rwtop      (:,:),   & ! [OUT]
!ESC!                       tctop      (:,:),   & ! [OUT]
!ESC!                       CBMFX      (:,:),   & ! [IN]
!ESC!                       rho        (:,:),   & ! [IN]
!ESC!                       tem        (:,:),   & ! [IN]
!ESC!                       pre        (:,:),   & ! [IN]
!ESC!                       gsgam2     (:,:),   & ! [IN]
!ESC!                       gsgam2h    (:,:),   & ! [IN]
!ESC!                       q          (:,:,:), & ! [IN]
!ESC!                       unccn      (:,:),   & ! [IN]
!ESC!                       z          (:,:),   & ! [IN]
!ESC!                       zh         (:,:),   & ! [IN]
!ESC!                       dt_mp               ) ! [IN]
!ESC!
!ESC!          !$omp parallel workshare
!ESC!          re_solid (:,:) = 20.E-6_RP ! not supported
!ESC!          re_cld   (:,:) = re_liquid(:,:) ! only qc is prognostic
!ESC!          rctop_cld(:,:) = rctop(:,:)
!ESC!          rwtop_cld(:,:) = rwtop(:,:)
!ESC!          tctop_cld(:,:) = tctop(:,:)
!ESC!          !$omp end parallel workshare
       endif

       !$omp parallel workshare
       precip_sum        (:,:) = precip_sum        (:,:) + precip        (:,:)
       ISO1_precip_sum   (:,:) = ISO1_precip_sum   (:,:) + ISO1_precip   (:,:)
       ISO2_precip_sum   (:,:) = ISO2_precip_sum   (:,:) + ISO2_precip   (:,:)
       precip_rhoe_sum   (:)   = precip_rhoe_sum   (:)   + precip_rhoe   (:)
       precip_lh_heat_sum(:)   = precip_lh_heat_sum(:)   + precip_lh_heat(:)
       precip_rhophi_sum (:)   = precip_rhophi_sum (:)   + precip_rhophi (:)
       precip_rhokin_sum (:)   = precip_rhokin_sum (:)   + precip_rhokin (:)
       precip_trc_sum    (:,:) = precip_trc_sum    (:,:) + precip_trc    (:,:)
       GDCLW_sum         (:,:) = GDCLW_sum         (:,:) + GDCLW         (:,:)
       GDCFRC_sum        (:,:) = GDCFRC_sum        (:,:) + GDCFRC        (:,:)
       GPREC_sum         (:,:) = GPREC_sum         (:,:) + GPREC         (:,:)
       !$omp end parallel workshare
    enddo

    !$omp parallel workshare
    precip        (:,:) = precip_sum        (:,:) * fraction_mp
    ISO1_precip   (:,:) = ISO1_precip_sum   (:,:) * fraction_mp
    ISO2_precip   (:,:) = ISO2_precip_sum   (:,:) * fraction_mp
    precip_rhoe   (:)   = precip_rhoe_sum   (:)   * fraction_mp
    precip_lh_heat(:)   = precip_lh_heat_sum(:)   * fraction_mp
    precip_rhophi (:)   = precip_rhophi_sum (:)   * fraction_mp
    precip_rhokin (:)   = precip_rhokin_sum (:)   * fraction_mp
    precip_trc    (:,:) = precip_trc_sum    (:,:) * fraction_mp
    GDCLW         (:,:) = GDCLW_sum         (:,:) * fraction_mp
    GDCFRC        (:,:) = GDCFRC_sum        (:,:) * fraction_mp
    GPREC         (:,:) = GPREC_sum         (:,:) * fraction_mp
    !$omp end parallel workshare

    ! [Move] 2016/03/29 T.Seiki
!!$    call history_in( 'ml_re_cld'    , re_cld     )
!!$    call history_in( 'ml_re_rain'   , re_rain    )
!!$    call history_in( 'ml_re_ice'    , re_ice     )
!!$    call history_in( 'ml_re_snow'   , re_snow    )
!!$    call history_in( 'ml_re_graupel', re_graupel )
!ESC!    call history_in( 'ml_re_liq'    , re_liquid  )
!ESC!    call history_in( 'ml_re_sol'    , re_solid   )
!ESC!
!ESC!    call history_in( 'sl_rctop_cld' , rctop_cld  )
!ESC!    call history_in( 'sl_rwtop_cld' , rwtop_cld  )
!ESC!    call history_in( 'sl_tctop_cld' , tctop_cld  )

    return
  end subroutine mp_driver

!ESC!  !-----------------------------------------------------------------------------
!ESC!  subroutine mp_terminal_velocity( &
!ESC!       ijdim, kdim,      & ! in
!ESC!       rho, tem, pre, w, & ! in
!ESC!       q,                & ! in
!ESC!       vt_qc, vt_qr, vt_qi, vt_qs, vt_qg, vt_qh,& ! out
!ESC!       vt_nc, vt_nr, vt_ni, vt_ns, vt_ng, vt_nh ) ! out
!ESC!    use mod_adm, only: &
!ESC!         kmin => adm_kmin, &
!ESC!         kmax => adm_kmax
!ESC!    use mod_mp_g98, only: &
!ESC!         mp_g98_terminal_velocity
!ESC!    use mod_mp_nsw5, only: &
!ESC!         mp_nsw5_terminal_velocity
!ESC!    use mod_mp_nsw6, only: &
!ESC!         mp_nsw6_terminal_velocity
!ESC!    use mod_mp_ndw6, only: &
!ESC!         mp_ndw6_terminal_velocity, &
!ESC!         mp_ndw6_diag_number ! [Add] 09/04/14 T.Mitsui
!ESC!    use mod_runconf, only: &
!ESC!         nqmax=>TRC_VMAX,  &
!ESC!         MP_TYPE,          &
!ESC!         opt_2moment_water, & ! [Add] 09/04/14 T.Mitsui
!ESC!         I_QV,              & ! [Add] 09/04/14 T.Mitsui
!ESC!         I_QC, I_QR, I_QI, I_QS, I_QG, &
!ESC!         I_NC, I_NR, I_NI, I_NS, I_NG
!ESC!    implicit none
!ESC!
!ESC!    integer, intent(in)  :: ijdim
!ESC!    integer, intent(in)  :: kdim
!ESC!
!ESC!    real(RP), intent(in)  :: rho(ijdim,kdim)   ! density
!ESC!    real(RP), intent(in)  :: tem(ijdim,kdim)   ! temperature
!ESC!    real(RP), intent(in)  :: pre(ijdim,kdim)   ! pressure
!ESC!    real(RP), intent(in)  :: w  (ijdim,kdim)   ! vertical velocity
!ESC!    real(RP), intent(in)  :: q  (ijdim,kdim,nqmax) ! mixing ratio of hydrometeors
!ESC!
!ESC!    real(RP), intent(out) :: vt_qc(ijdim,kdim) ! terminal velocity of cloud mass
!ESC!    real(RP), intent(out) :: vt_qr(ijdim,kdim) ! terminal velocity of rain mass
!ESC!    real(RP), intent(out) :: vt_qi(ijdim,kdim) ! terminal velocity of ice mass
!ESC!    real(RP), intent(out) :: vt_qs(ijdim,kdim) ! terminal velocity of snow mass
!ESC!    real(RP), intent(out) :: vt_qg(ijdim,kdim) ! terminal velocity of graupel mass
!ESC!    real(RP), intent(out) :: vt_qh(ijdim,kdim) ! terminal velocity of hail mass
!ESC!
!ESC!    real(RP), intent(out) :: vt_nc(ijdim,kdim) ! terminal velocity of cloud mass
!ESC!    real(RP), intent(out) :: vt_nr(ijdim,kdim) ! terminal velocity of rain mass
!ESC!    real(RP), intent(out) :: vt_ni(ijdim,kdim) ! terminal velocity of ice mass
!ESC!    real(RP), intent(out) :: vt_ns(ijdim,kdim) ! terminal velocity of snow mass
!ESC!    real(RP), intent(out) :: vt_ng(ijdim,kdim) ! terminal velocity of graupel mass
!ESC!    real(RP), intent(out) :: vt_nh(ijdim,kdim) ! terminal velocity of hail mass
!ESC!    ! [Add] 09/04/14 T.Mitsui
!ESC!    real(RP) :: qnc(ijdim,kdim)
!ESC!    real(RP) :: qnr(ijdim,kdim)
!ESC!    real(RP) :: qni(ijdim,kdim)
!ESC!    real(RP) :: qns(ijdim,kdim)
!ESC!    real(RP) :: qng(ijdim,kdim)
!ESC!    real(RP) :: rrho(ijdim,kdim)
!ESC!
!ESC!    vt_qc(:,:)=0.0_RP
!ESC!    vt_qr(:,:)=0.0_RP
!ESC!    vt_qi(:,:)=0.0_RP
!ESC!    vt_qs(:,:)=0.0_RP
!ESC!    vt_qg(:,:)=0.0_RP
!ESC!    vt_qh(:,:)=0.0_RP
!ESC!    if( trim(MP_TYPE) /= 'NDW6' ) then
!ESC!       vt_nc(:,:)=0.0_RP
!ESC!       vt_nr(:,:)=0.0_RP
!ESC!       vt_ni(:,:)=0.0_RP
!ESC!       vt_ns(:,:)=0.0_RP
!ESC!       vt_ng(:,:)=0.0_RP
!ESC!    endif
!ESC!    vt_nh(:,:)=0.0_RP
!ESC!
!ESC!    if( MP_TYPE == 'NONE' ) then
!ESC!       ! Notice, T.Mitsui never use these mp_schemes.
!ESC!       !         Please calculate for yourself who want to.
!ESC!    elseif( MP_TYPE == 'KESSLER') then
!ESC!    elseif( MP_TYPE == 'LIN') then
!ESC!    elseif( MP_TYPE == 'WSM6') then
!ESC!    elseif( MP_TYPE == 'G98') then
!ESC!       call mp_g98_terminal_velocity( &
!ESC!            ijdim, kdim, & ! in
!ESC!            rho, tem,    & ! in
!ESC!            q(:,:,I_QC), q(:,:,I_QR),  & ! in
!ESC!            vt_qc, vt_qr, vt_qi, vt_qs ) ! out
!ESC!    elseif( MP_TYPE == 'NSW5') then
!ESC!       call mp_nsw5_terminal_velocity( &
!ESC!            ijdim, kdim, & ! in
!ESC!            rho, tem,    & ! in
!ESC!            q(:,:,I_QC), q(:,:,I_QR),  & ! in
!ESC!            q(:,:,I_QI), q(:,:,I_QS),  & ! in
!ESC!            vt_qc, vt_qr, vt_qi, vt_qs ) ! out
!ESC!    elseif( MP_TYPE == 'NSW6') then
!ESC!       call mp_nsw6_terminal_velocity( &
!ESC!            ijdim, kdim, & ! in
!ESC!            rho, tem,    & ! in
!ESC!            w,    & ! [Add] 2016/06/04 WS.ROH
!ESC!            q(:,:,I_QC), q(:,:,I_QR),  & ! in
!ESC!            q(:,:,I_QI), q(:,:,I_QS), q(:,:,I_QG),& ! in
!ESC!            vt_qc, vt_qr, vt_qi, vt_qs, vt_qg     ) ! out
!ESC!    elseif( MP_TYPE == 'NDW6' ) then
!ESC!       ! 09/04/14 T.Mitsui [Add] conditioning
!ESC!       if(opt_2moment_water) then
!ESC!          qnc(:,:) = q(:,:,I_NC)
!ESC!          qnr(:,:) = q(:,:,I_NR)
!ESC!          qni(:,:) = q(:,:,I_NI)
!ESC!          qns(:,:) = q(:,:,I_NS)
!ESC!          qng(:,:) = q(:,:,I_NG)
!ESC!       else
!ESC!          call mp_ndw6_diag_number( &
!ESC!               ijdim, kdim,   & ! in
!ESC!               kmin, kmax,    & ! in
!ESC!               rho, tem,      & ! in
!ESC!               q(:,:,I_QV),   & ! in
!ESC!               q(:,:,I_QC),   & ! in
!ESC!               q(:,:,I_QR),   & ! in
!ESC!               q(:,:,I_QI),   & ! in
!ESC!               q(:,:,I_QS),   & ! in
!ESC!               q(:,:,I_QG),   & ! in
!ESC!               qnc, qnr, qni, qns, qng ) ! out
!ESC!          rrho(:,:) = 1.0_RP/rho(:,:)
!ESC!          qnc(:,:) = qnc(:,:)*rrho(:,:)
!ESC!          qnr(:,:) = qnr(:,:)*rrho(:,:)
!ESC!          qni(:,:) = qni(:,:)*rrho(:,:)
!ESC!          qns(:,:) = qns(:,:)*rrho(:,:)
!ESC!          qng(:,:) = qng(:,:)*rrho(:,:)
!ESC!       endif
!ESC!
!ESC!       call mp_ndw6_terminal_velocity( &
!ESC!            ijdim, kdim,   & ! in
!ESC!            ! [Mod] 10/08/03 T.Mitsui, add argument "pre"
!ESC!!!$         rho, tem,      & ! in
!ESC!            rho, tem, pre, & ! in
!ESC!            q(:,:,I_QC), q(:,:,I_QR), q(:,:,I_QI), q(:,:,I_QS), q(:,:,I_QG), & ! in
!ESC!            ! [Mod] 09/04/14 T.Mitsui
!ESC!!!$         q(:,:,I_NC), q(:,:,I_NR), q(:,:,I_NI), q(:,:,I_NS), q(:,:,I_NG), & ! in
!ESC!            qnc(:,:), qnr(:,:), qni(:,:), qns(:,:), qng(:,:), & ! in
!ESC!            vt_qc, vt_qr, vt_qi, vt_qs, vt_qg, & ! out
!ESC!            vt_nc, vt_nr, vt_ni, vt_ns, vt_ng  ) ! out
!ESC!    elseif( MP_TYPE == 'LSC') then
!ESC!       vt_qr(:,:)=-5.d0 ! given in namelist variable as VTERM
!ESC!       vt_qs(:,:)=-5.d0 ! given in namelist variable as VTERM
!ESC!    else
!ESC!       write(*,*) &
!ESC!            'Msg : Sub[mp_terminal_velocity]/Mod[mp_driver]'
!ESC!       write(*,*) &
!ESC!            ' *** WARNING : Not appropriate type!! CHECK!!'
!ESC!    endif
!ESC!
!ESC!    return
!ESC!  end subroutine mp_terminal_velocity
!ESC!  !-----------------------------------------------------------------------------
!ESC!  ! 09/04/14 [Add] T.Mitsui, used for radiation
!ESC!  ! Temporarily values are imported from each subroutine and set density as constant.
!ESC!  ! Please prepare subroutines for yourself to synchronize parameters with cloud microphysics.
!ESC!  subroutine mp_diag_volume(      &
!ESC!       ijdim, kmin, kmax, kdim,   &
!ESC!       l_region,                  &
!ESC!       rho, tem, pre,             &
!ESC!       q, qc_cum, cfrac, cfrac_cum, &
!ESC!       cloud_volume, cumulus_volume)
!ESC!    use mod_const, only: &
!ESC!       DWATR => CONST_DWATR, &
!ESC!       DICE  => CONST_DICE
!ESC!    use mod_runconf, only: &
!ESC!       nqmax=>TRC_VMAX,   &
!ESC!       MP_TYPE,           &
!ESC!       RAIN_TYPE,         &
!ESC!       opt_2moment_water, &
!ESC!       HYDRO_MAX,         &
!ESC!       NQW_STR, NQW_END,  &
!ESC!       I_QV,              &
!ESC!       I_QC, I_QR, I_QI, I_QS, I_QG, I_QH, &
!ESC!       I_NC, I_NR, I_NI, I_NS, I_NG, I_NH
!ESC!    use mod_mp_ndw6, only: &
!ESC!       mp_ndw6_diag_volume, &
!ESC!       mp_ndw6_diag_number
!ESC!    use mod_mp_g98, only: &
!ESC!       mp_g98_diag_volume
!ESC!    use mod_mp_lsc, only: &
!ESC!       mp_lsc_diag_volume
!ESC!    implicit none
!ESC!
!ESC!    integer, intent(in)  :: ijdim
!ESC!    integer, intent(in)  :: kmin
!ESC!    integer, intent(in)  :: kmax
!ESC!    integer, intent(in)  :: kdim
!ESC!    integer, intent(in)  :: l_region
!ESC!    real(RP), intent(in)  :: q             (ijdim,kdim,nqmax)     ! [kg/kg] mixing ratio
!ESC!    real(RP), intent(in)  :: qc_cum        (ijdim,kdim)           ! [kg/kg] cumulus mixing ratio
!ESC!    real(RP), intent(in)  :: cfrac         (ijdim,kdim)           ! [0:1] cloud fraction
!ESC!    real(RP), intent(in)  :: cfrac_cum     (ijdim,1)              ! [0:1] cumulus cloud fraction
!ESC!    real(RP), intent(in)  :: rho           (ijdim,kdim)           ! [kg/m3] air density
!ESC!    real(RP), intent(in)  :: tem           (ijdim,kdim)           ! [K] temperature
!ESC!    real(RP), intent(in)  :: pre           (ijdim,kdim)           ! [Pa] Pressure
!ESC!    real(RP), intent(out) :: cloud_volume  (ijdim,kdim,HYDRO_MAX) ! [m3/m3]
!ESC!    real(RP), intent(out) :: cumulus_volume(ijdim,kdim,2)         ! [m3/m3]
!ESC!
!ESC!    real(RP) :: nc(ijdim,kdim)
!ESC!    real(RP) :: nr(ijdim,kdim)
!ESC!    real(RP) :: ni(ijdim,kdim)
!ESC!    real(RP) :: ns(ijdim,kdim)
!ESC!    real(RP) :: ng(ijdim,kdim)
!ESC!    real(RP) :: FLIQ
!ESC!
!ESC!    integer, parameter :: I_rd_QC = 2
!ESC!    integer, parameter :: I_rd_QR = 3
!ESC!    integer, parameter :: I_rd_QI = 4
!ESC!    integer, parameter :: I_rd_QS = 5
!ESC!    integer, parameter :: I_rd_QG = 6
!ESC!    integer, parameter :: I_rd_QH = 7
!ESC!
!ESC!    integer :: ij, k, nq
!ESC!    !---------------------------------------------------------------------------
!ESC!
!ESC!    cloud_volume(:,:,1) = 0.0_RP
!ESC!
!ESC!    if ( opt_volume_explicit ) then
!ESC!
!ESC!       if ( MP_TYPE == 'NDW6' ) then
!ESC!
!ESC!          if ( opt_2moment_water ) then
!ESC!             nc(:,:) = rho(:,:) * q(:,:,I_NC)
!ESC!             nr(:,:) = rho(:,:) * q(:,:,I_NR)
!ESC!             ni(:,:) = rho(:,:) * q(:,:,I_NI)
!ESC!             ns(:,:) = rho(:,:) * q(:,:,I_NS)
!ESC!             ng(:,:) = rho(:,:) * q(:,:,I_NG)
!ESC!          else
!ESC!             call mp_ndw6_diag_number( ijdim,       & ! [IN]
!ESC!                                       kdim,        & ! [IN]
!ESC!                                       kmin,        & ! [IN]
!ESC!                                       kmax,        & ! [IN]
!ESC!                                       rho(:,:),    & ! [IN]
!ESC!                                       tem(:,:),    & ! [IN]
!ESC!                                       q(:,:,I_QV), & ! [IN]
!ESC!                                       q(:,:,I_QC), & ! [IN]
!ESC!                                       q(:,:,I_QR), & ! [IN]
!ESC!                                       q(:,:,I_QI), & ! [IN]
!ESC!                                       q(:,:,I_QS), & ! [IN]
!ESC!                                       q(:,:,I_QG), & ! [IN]
!ESC!                                       nc(:,:),     & ! [OUT]
!ESC!                                       nr(:,:),     & ! [OUT]
!ESC!                                       ni(:,:),     & ! [OUT]
!ESC!                                       ns(:,:),     & ! [OUT]
!ESC!                                       ng(:,:)      ) ! [OUT]
!ESC!          endif
!ESC!
!ESC!          call mp_ndw6_diag_volume( ijdim,                  & ! [IN]
!ESC!                                    kdim,                   & ! [IN]
!ESC!                                    kmin,                   & ! [IN]
!ESC!                                    kmax,                   & ! [IN]
!ESC!                                    rho         (:,:),      & ! [IN]
!ESC!                                    nc          (:,:),      & ! [IN]
!ESC!                                    nr          (:,:),      & ! [IN]
!ESC!                                    ni          (:,:),      & ! [IN]
!ESC!                                    ns          (:,:),      & ! [IN]
!ESC!                                    ng          (:,:),      & ! [IN]
!ESC!                                    q           (:,:,I_QC), & ! [IN]
!ESC!                                    q           (:,:,I_QR), & ! [IN]
!ESC!                                    q           (:,:,I_QI), & ! [IN]
!ESC!                                    q           (:,:,I_QS), & ! [IN]
!ESC!                                    q           (:,:,I_QG), & ! [IN]
!ESC!                                    cloud_volume(:,:,I_rd_QC), & ! [OUT]
!ESC!                                    cloud_volume(:,:,I_rd_QR), & ! [OUT]
!ESC!                                    cloud_volume(:,:,I_rd_QI), & ! [OUT]
!ESC!                                    cloud_volume(:,:,I_rd_QS), & ! [OUT]
!ESC!                                    cloud_volume(:,:,I_rd_QG)  ) ! [OUT]
!ESC!
!ESC!          cloud_volume(:,:,I_rd_QH) = 0.0_RP
!ESC!
!ESC!       elseif( MP_TYPE == 'NSW6' ) then
!ESC!
!ESC!          cloud_volume(:,:,I_rd_QC) = rho(:,:) * q(:,:,I_QC) / DWATR
!ESC!          cloud_volume(:,:,I_rd_QR) = rho(:,:) * q(:,:,I_QR) / DWATR
!ESC!          ! [Fix] 2016/05/13 T.Seiki: here we should use optically effective volume
!ESC!          !                           see Fu etal.(1996),JC or Seiki etal.(2014),JGR
!ESC!!          cloud_volume(:,:,I_rd_QI) = rho(:,:) * q(:,:,I_QI) / 400.0_RP ! rho_g
!ESC!!          cloud_volume(:,:,I_rd_QS) = rho(:,:) * q(:,:,I_QS) / 100.0_RP ! rho_s
!ESC!!          cloud_volume(:,:,I_rd_QG) = rho(:,:) * q(:,:,I_QG) / 400.0_RP ! rho_g
!ESC!          cloud_volume(:,:,I_rd_QI) = rho(:,:) * q(:,:,I_QI) / DICE
!ESC!          cloud_volume(:,:,I_rd_QS) = rho(:,:) * q(:,:,I_QS) / DICE
!ESC!          cloud_volume(:,:,I_rd_QG) = rho(:,:) * q(:,:,I_QG) / DICE
!ESC!          cloud_volume(:,:,I_rd_QH) = 0.0_RP
!ESC!
!ESC!       elseif( MP_TYPE == 'NSW5' ) then
!ESC!
!ESC!          cloud_volume(:,:,I_rd_QC) = rho(:,:) * q(:,:,I_QC) / DWATR
!ESC!          cloud_volume(:,:,I_rd_QR) = rho(:,:) * q(:,:,I_QR) / DWATR
!ESC!          cloud_volume(:,:,I_rd_QI) = rho(:,:) * q(:,:,I_QI) / 400.0_RP ! rho_g
!ESC!          cloud_volume(:,:,I_rd_QS) = rho(:,:) * q(:,:,I_QS) / 100.0_RP ! rho_s
!ESC!          cloud_volume(:,:,I_rd_QG) = 0.0_RP
!ESC!          cloud_volume(:,:,I_rd_QH) = 0.0_RP
!ESC!
!ESC!       elseif( MP_TYPE == 'WSM6' ) then
!ESC!
!ESC!          cloud_volume(:,:,I_rd_QC) = rho(:,:) * q(:,:,I_QC) / DWATR
!ESC!          cloud_volume(:,:,I_rd_QR) = rho(:,:) * q(:,:,I_QR) / DWATR
!ESC!          cloud_volume(:,:,I_rd_QI) = rho(:,:) * q(:,:,I_QI) / 500.0_RP ! rho_g
!ESC!          cloud_volume(:,:,I_rd_QS) = rho(:,:) * q(:,:,I_QS) / 100.0_RP ! rho_s
!ESC!          cloud_volume(:,:,I_rd_QG) = rho(:,:) * q(:,:,I_QG) / 500.0_RP ! rho_g
!ESC!          cloud_volume(:,:,I_rd_QH) = 0.0_RP
!ESC!
!ESC!       elseif( MP_TYPE == 'LIN' ) then
!ESC!
!ESC!          cloud_volume(:,:,I_rd_QC) = rho(:,:) * q(:,:,I_QC) / DWATR
!ESC!          cloud_volume(:,:,I_rd_QR) = rho(:,:) * q(:,:,I_QR) / DWATR
!ESC!          cloud_volume(:,:,I_rd_QI) = rho(:,:) * q(:,:,I_QI) / 400.0_RP ! rho_g
!ESC!          cloud_volume(:,:,I_rd_QS) = rho(:,:) * q(:,:,I_QS) / 100.0_RP ! rho_s
!ESC!          cloud_volume(:,:,I_rd_QG) = rho(:,:) * q(:,:,I_QG) / 400.0_RP ! rho_g
!ESC!          cloud_volume(:,:,I_rd_QH) = 0.0_RP
!ESC!
!ESC!       elseif( MP_TYPE == 'KESSLER' ) then
!ESC!
!ESC!          cloud_volume(:,:,I_rd_QC) = rho(:,:) * q(:,:,I_QC) / DWATR
!ESC!          cloud_volume(:,:,I_rd_QR) = 0.0_RP
!ESC!          cloud_volume(:,:,I_rd_QI) = 0.0_RP
!ESC!          cloud_volume(:,:,I_rd_QS) = 0.0_RP
!ESC!          cloud_volume(:,:,I_rd_QG) = 0.0_RP
!ESC!          cloud_volume(:,:,I_rd_QH) = 0.0_RP
!ESC!
!ESC!       elseif( MP_TYPE == 'G98' ) then
!ESC!
!ESC!          call mp_g98_diag_volume( ijdim,                  & ! [IN]
!ESC!                                   kdim,                   & ! [IN]
!ESC!                                   kmin,                   & ! [IN]
!ESC!                                   kmax,                   & ! [IN]
!ESC!                                   rho         (:,:),      & ! [IN]
!ESC!                                   tem         (:,:),      & ! [IN]
!ESC!                                   pre         (:,:),      & ! [IN]
!ESC!                                   q           (:,:,I_QC), & ! [IN]
!ESC!                                   q           (:,:,I_QR), & ! [IN]
!ESC!                                   cloud_volume(:,:,I_rd_QC), & ! [OUT]
!ESC!                                   cloud_volume(:,:,I_rd_QR), & ! [OUT]
!ESC!                                   cloud_volume(:,:,I_rd_QI), & ! [OUT]
!ESC!                                   cloud_volume(:,:,I_rd_QS)  ) ! [OUT]
!ESC!
!ESC!          cloud_volume(:,:,I_rd_QG) = 0.0_RP
!ESC!          cloud_volume(:,:,I_rd_QH) = 0.0_RP
!ESC!
!ESC!       elseif( MP_TYPE == 'LSC' ) then
!ESC!
!ESC!          call mp_lsc_diag_volume( ijdim,                  & ! [IN]
!ESC!                                   kdim,                   & ! [IN]
!ESC!                                   kmin,                   & ! [IN]
!ESC!                                   kmax,                   & ! [IN]
!ESC!                                   rho         (:,:),      & ! [IN]
!ESC!                                   tem         (:,:),      & ! [IN]
!ESC!                                   q           (:,:,I_QC), & ! [IN]
!ESC!                                   cfrac       (:,:),      & ! [IN]
!ESC!                                   cloud_volume(:,:,I_rd_QC), & ! [OUT]
!ESC!                                   cloud_volume(:,:,I_rd_QI)  ) ! [OUT]
!ESC!
!ESC!          cloud_volume(:,:,I_rd_QR) = 0.0_RP
!ESC!          cloud_volume(:,:,I_rd_QS) = 0.0_RP
!ESC!          cloud_volume(:,:,I_rd_QG) = 0.0_RP
!ESC!          cloud_volume(:,:,I_rd_QH) = 0.0_RP
!ESC!
!ESC!       elseif( MP_TYPE == 'NONE' ) then
!ESC!
!ESC!          cloud_volume(:,:,I_rd_QC) = 0.0_RP
!ESC!          cloud_volume(:,:,I_rd_QR) = 0.0_RP
!ESC!          cloud_volume(:,:,I_rd_QI) = 0.0_RP
!ESC!          cloud_volume(:,:,I_rd_QS) = 0.0_RP
!ESC!          cloud_volume(:,:,I_rd_QG) = 0.0_RP
!ESC!          cloud_volume(:,:,I_rd_QH) = 0.0_RP
!ESC!
!ESC!       endif
!ESC!
!ESC!    else ! set default
!ESC!
!ESC!       do nq = NQW_STR+1, NQW_END
!ESC!          if ( nq == I_QC ) then
!ESC!             cloud_volume(:,:,nq) = rho(:,:) * q(:,:,nq) / DWATR
!ESC!          elseif( nq == I_QR ) then
!ESC!             cloud_volume(:,:,nq) = rho(:,:) * q(:,:,nq) / DWATR
!ESC!          elseif( nq == I_QI ) then
!ESC!             cloud_volume(:,:,nq) = rho(:,:) * q(:,:,nq) / 400.0_RP
!ESC!          elseif( nq == I_QS ) then
!ESC!             cloud_volume(:,:,nq) = rho(:,:) * q(:,:,nq) / 100.0_RP
!ESC!          elseif( nq == I_QG ) then
!ESC!             cloud_volume(:,:,nq) = rho(:,:) * q(:,:,nq) / 400.0_RP
!ESC!          elseif( nq == I_QH ) then
!ESC!             cloud_volume(:,:,nq) = rho(:,:) * q(:,:,nq) / DICE
!ESC!          else
!ESC!             cloud_volume(:,:,nq) = 0.0_RP
!ESC!          endif
!ESC!       enddo
!ESC!
!ESC!    endif
!ESC!
!ESC!    ! cumulus_volume with liquid/ice separation
!ESC!    do k  = kmin, kmax
!ESC!    do ij = 1, ijdim
!ESC!       FLIQ = ( tem(ij,k)-TWICE ) / ( TSICE-TWICE )
!ESC!       FLIQ = min( max( FLIQ, 0.0_RP ), 1.0_RP )
!ESC!
!ESC!       if ( cfrac_cum(ij,1) >= 1.D-12 ) then
!ESC!          cumulus_volume(ij,k,1) = rho(ij,k) * (      FLIQ ) * qc_cum(ij,k) / ( cfrac_cum(ij,1) * DWATR )
!ESC!          cumulus_volume(ij,k,2) = rho(ij,k) * ( 1.0_RP-FLIQ ) * qc_cum(ij,k) / ( cfrac_cum(ij,1) * DICE  )
!ESC!       else
!ESC!          cumulus_volume(ij,k,1) = 0.0_RP
!ESC!          cumulus_volume(ij,k,2) = 0.0_RP
!ESC!       endif
!ESC!    enddo
!ESC!    enddo
!ESC!
!ESC!    return
!ESC!  end subroutine mp_diag_volume
!ESC!
!ESC!  !-----------------------------------------------------------------------------
!ESC!  subroutine mp_effective_radius( &
!ESC!       ijdim,    &
!ESC!       kmin,     &
!ESC!       kmax,     &
!ESC!       kdim,     &
!ESC!       l_region, &
!ESC!       q,        &
!ESC!       GDCFRC,   &
!ESC!       rho,      &
!ESC!       tem,      &
!ESC!       pre,      &
!ESC!       ccn,      &
!ESC!       w,        &
!ESC!       re        )
!ESC!    use mod_runconf, only: &
!ESC!       nqmax => TRC_VMAX,            &
!ESC!       MP_TYPE,                      &
!ESC!       opt_2moment_water,            &
!ESC!       HYDRO_MAX,                    &
!ESC!       I_QV,                         &
!ESC!       I_QC, I_QR, I_QI, I_QS, I_QG, &
!ESC!       I_NC, I_NR, I_NI, I_NS, I_NG
!ESC!    use mod_mp_ndw6, only: &
!ESC!       mp_ndw6_effective_radius, &
!ESC!       mp_ndw6_diag_number
!ESC!    use mod_mp_nsw6, only: &
!ESC!       mp_nsw6_effective_radius
!ESC!    use mod_mp_nsw5, only: &
!ESC!       mp_nsw5_effective_radius
!ESC!    use mod_mp_g98, only: &
!ESC!       mp_g98_effective_radius
!ESC!    use mod_mp_lsc, only: &
!ESC!       mp_lsc_effective_radius
!ESC!    use mod_history, only: &
!ESC!       history_in
!ESC!    implicit none
!ESC!
!ESC!    integer,  intent(in)  :: ijdim
!ESC!    integer,  intent(in)  :: kmin
!ESC!    integer,  intent(in)  :: kmax
!ESC!    integer,  intent(in)  :: kdim
!ESC!    integer,  intent(in)  :: l_region
!ESC!    real(RP), intent(in)  :: q     (ijdim,kdim,nqmax)     ! mixing ratio [kg/kg]
!ESC!    real(RP), intent(in)  :: GDCFRC(ijdim,kdim)           ! cloud fraction[0:1]
!ESC!    real(RP), intent(in)  :: rho   (ijdim,kdim)           ! air density [kg/m3]
!ESC!    real(RP), intent(in)  :: tem   (ijdim,kdim)           ! temperature [K]
!ESC!    real(RP), intent(in)  :: pre   (ijdim,kdim)           ! pressure [Pa]
!ESC!    real(RP), intent(in)  :: ccn   (ijdim,kdim)           ! cloud condensation nuclei[1/m3]
!ESC!    real(RP), intent(in)  :: w     (ijdim,kdim)           ! vertical velocity [m/s] 2016/03/29 T.Seiki
!ESC!    real(RP), intent(out) :: re    (ijdim,kdim,HYDRO_MAX) ! effective radius[m]
!ESC!
!ESC!    ! number concentration[1/m3]
!ESC!    real(RP) :: nc(ijdim,kdim)
!ESC!    real(RP) :: nr(ijdim,kdim)
!ESC!    real(RP) :: ni(ijdim,kdim)
!ESC!    real(RP) :: ns(ijdim,kdim)
!ESC!    real(RP) :: ng(ijdim,kdim)
!ESC!
!ESC!    integer, parameter :: I_rd_QC = 2
!ESC!    integer, parameter :: I_rd_QR = 3
!ESC!    integer, parameter :: I_rd_QI = 4
!ESC!    integer, parameter :: I_rd_QS = 5
!ESC!    integer, parameter :: I_rd_QG = 6
!ESC!    integer, parameter :: I_rd_QH = 7
!ESC!    !---------------------------------------------------------------------------
!ESC!
!ESC!
!ESC!    re(:,:,:)       =  0.0_RP
!ESC!    re(:,:,I_rd_QC) =  8.E-6_RP
!ESC!    re(:,:,I_rd_QI) = 20.E-6_RP
!ESC!
!ESC!    if ( .NOT. opt_radius_explicit ) then
!ESC!       re(:,:,I_rd_QR) = 500.E-6_RP ! 0.5 mm
!ESC!       re(:,:,I_rd_QS) = 500.E-6_RP ! 0.5 mm
!ESC!       re(:,:,I_rd_QG) = 500.E-6_RP ! 0.5 mm
!ESC!       re(:,:,I_rd_QH) = 500.E-6_RP ! 0.5 mm
!ESC!       return
!ESC!    endif
!ESC!
!ESC!    if ( MP_TYPE == 'NDW6' ) then
!ESC!       if ( opt_2moment_water ) then
!ESC!          nc(:,:) = rho(:,:) * q(:,:,I_NC)
!ESC!          nr(:,:) = rho(:,:) * q(:,:,I_NR)
!ESC!          ni(:,:) = rho(:,:) * q(:,:,I_NI)
!ESC!          ns(:,:) = rho(:,:) * q(:,:,I_NS)
!ESC!          ng(:,:) = rho(:,:) * q(:,:,I_NG)
!ESC!       else
!ESC!          call mp_ndw6_diag_number( ijdim,             & ! [IN]
!ESC!                                    kdim, kmin, kmax,  & ! [IN]
!ESC!                                    rho, tem,          & ! [IN]
!ESC!                                    q(:,:,I_QV),       & ! [IN]
!ESC!                                    q(:,:,I_QC),       & ! [IN]
!ESC!                                    q(:,:,I_QR),       & ! [IN]
!ESC!                                    q(:,:,I_QI),       & ! [IN]
!ESC!                                    q(:,:,I_QS),       & ! [IN]
!ESC!                                    q(:,:,I_QG),       & ! [IN]
!ESC!                                    nc, nr, ni, ns, ng ) ! [OUT]
!ESC!       endif
!ESC!
!ESC!       call mp_ndw6_effective_radius( ijdim,              & ! [IN]
!ESC!                                      kdim, kmin, kmax,   & ! [IN]
!ESC!                                      rho, tem, pre,      & ! [IN]
!ESC!                                      nc, nr, ni, ns, ng, & ! [IN]
!ESC!                                      q (:,:,I_QC),       & ! [IN]
!ESC!                                      q (:,:,I_QR),       & ! [IN]
!ESC!                                      q (:,:,I_QI),       & ! [IN]
!ESC!                                      q (:,:,I_QS),       & ! [IN]
!ESC!                                      q (:,:,I_QG),       & ! [IN]
!ESC!                                      re(:,:,I_rd_QC),    & ! [OUT]
!ESC!                                      re(:,:,I_rd_QR),    & ! [OUT]
!ESC!                                      re(:,:,I_rd_QI),    & ! [OUT]
!ESC!                                      re(:,:,I_rd_QS),    & ! [OUT]
!ESC!                                      re(:,:,I_rd_QG)     ) ! [OUT]
!ESC!    elseif( MP_TYPE == 'KESSLER' ) then
!ESC!       ! Tentative, constant value
!ESC!       re(:,:,I_rd_QC) = 8.E-6_RP
!ESC!    elseif( MP_TYPE == 'G98'  ) then
!ESC!       call mp_g98_effective_radius( ijdim,              & ! [IN]
!ESC!                                     kdim, kmin, kmax,   & ! [IN]
!ESC!                                     rho, tem, pre,      & ! [IN]
!ESC!                                     ccn,                & ! [IN]
!ESC!                                     q (:,:,I_QC),       & ! [IN]
!ESC!                                     q (:,:,I_QR),       & ! [IN]
!ESC!                                     re(:,:,I_rd_QC),    & ! [OUT]
!ESC!                                     re(:,:,I_rd_QR),    & ! [OUT]
!ESC!                                     re(:,:,I_rd_QI),    & ! [OUT]
!ESC!                                     re(:,:,I_rd_QS)     ) ! [OUT]
!ESC!    elseif( MP_TYPE == 'WSM6' ) then ! TODO
!ESC!       ! Tentative, constant value
!ESC!       re(:,:,I_rd_QC) =   8.E-6_RP
!ESC!       re(:,:,I_rd_QR) = 100.E-6_RP
!ESC!       re(:,:,I_rd_QI) =  20.E-6_RP
!ESC!       re(:,:,I_rd_QS) = 100.E-6_RP
!ESC!       re(:,:,I_rd_QG) = 100.E-6_RP
!ESC!    elseif( MP_TYPE == 'NSW5' ) then
!ESC!       call mp_nsw5_effective_radius( ijdim,              & ! [IN]
!ESC!                                      kdim, kmin, kmax,   & ! [IN]
!ESC!                                      rho, tem, pre,      & ! [IN]
!ESC!                                      ccn,                & ! [IN]
!ESC!                                      q (:,:,I_QC),       & ! [IN]
!ESC!                                      q (:,:,I_QR),       & ! [IN]
!ESC!                                      q (:,:,I_QI),       & ! [IN]
!ESC!                                      q (:,:,I_QS),       & ! [IN]
!ESC!                                      re(:,:,I_rd_QC),    & ! [OUT]
!ESC!                                      re(:,:,I_rd_QR),    & ! [OUT]
!ESC!                                      re(:,:,I_rd_QI),    & ! [OUT]
!ESC!                                      re(:,:,I_rd_QS)     ) ! [OUT]
!ESC!    elseif( MP_TYPE == 'NSW6' ) then
!ESC!       call mp_nsw6_effective_radius( ijdim,              & ! [IN]
!ESC!                                      kdim, kmin, kmax,   & ! [IN]
!ESC!                                      rho, tem, pre,      & ! [IN]
!ESC!                                      ccn,                & ! [IN]
!ESC!                                      w,                  & ! [IN] [Add] 2016/03/29 T.Seiki
!ESC!                                      q (:,:,I_QC),       & ! [IN]
!ESC!                                      q (:,:,I_QR),       & ! [IN]
!ESC!                                      q (:,:,I_QI),       & ! [IN]
!ESC!                                      q (:,:,I_QS),       & ! [IN]
!ESC!                                      q (:,:,I_QG),       & ! [IN]
!ESC!                                      re(:,:,I_rd_QC),    & ! [OUT]
!ESC!                                      re(:,:,I_rd_QR),    & ! [OUT]
!ESC!                                      re(:,:,I_rd_QI),    & ! [OUT]
!ESC!                                      re(:,:,I_rd_QS),    & ! [OUT]
!ESC!                                      re(:,:,I_rd_QG)     ) ! [OUT]
!ESC!    elseif( MP_TYPE == 'LIN'  ) then  ! TODO
!ESC!       ! Tentative, constant value
!ESC!       re(:,:,I_rd_QC) =   8.E-6_RP
!ESC!       re(:,:,I_rd_QR) = 100.E-6_RP
!ESC!       re(:,:,I_rd_QI) =  20.E-6_RP
!ESC!       re(:,:,I_rd_QS) = 100.E-6_RP
!ESC!       re(:,:,I_rd_QG) = 100.E-6_RP
!ESC!    elseif( MP_TYPE == 'LSC'  ) then  ! TODO
!ESC!       call mp_lsc_effective_radius( ijdim,              & ! [IN]
!ESC!                                     kdim, kmin, kmax,   & ! [IN]
!ESC!                                     rho, tem, pre,      & ! [IN]
!ESC!                                     q (:,:,I_QC),       & ! [IN]
!ESC!                                     ccn,                & ! [IN]
!ESC!                                     GDCFRC,             & ! [IN]
!ESC!                                     re(:,:,I_rd_QC),    & ! [OUT]
!ESC!                                     re(:,:,I_rd_QI)     ) ! [OUT]
!ESC!    endif
!ESC!
!ESC!    ! [Add] 2016/03/29 T.Seiki
!ESC!    call history_in( 'ml_re_cld'    , re(:,:,I_rd_QC) )
!ESC!    call history_in( 'ml_re_rain'   , re(:,:,I_rd_QR) )
!ESC!    call history_in( 'ml_re_ice'    , re(:,:,I_rd_QI) )
!ESC!    call history_in( 'ml_re_snow'   , re(:,:,I_rd_QS) )
!ESC!    call history_in( 'ml_re_graupel', re(:,:,I_rd_QG) )
!ESC!
!ESC!    return
!ESC!  end subroutine mp_effective_radius

end module mod_mp_driver
