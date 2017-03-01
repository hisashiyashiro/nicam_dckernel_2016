!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Radiative Transfer driver
!!
!! @par Description
!!          driver of radiative transfer schemes
!!
!! @author NICAM developers
!!
!<
!-------------------------------------------------------------------------------
module mod_rd_driver
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
  !++ Public procedures
  !
  public :: rd_init
  public :: rd_driver
!ESC!  public :: rd_surface
!ESC!  public :: rd_tendency
!ESC!  public :: rd_heatingrate
  public :: rd_solarincidence
!ESC!  public :: rd_ccover

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
!ESC!  character(len=H_SHORT), private :: RD_TYPE

  ! NAMELIST NM_RD_TIME
  real(DP),         private :: TINTV = 5.0_DP
  character(len=4), private :: TUNIT = 'MIN'

  ! NAMELIST NM_RD_SOLARCONSTANT
  logical,  private :: OSOLCON_HIST = .false.   ! historical Solar Constant data
  real(RP), private :: SOLCON       = 1365.0_RP ! climatorogical Solar Constant value

  real(DP), private :: RD_TIME_PREV
  real(DP), private :: RD_TIME_ORIGIN

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> setup
  subroutine rd_init( RD_TYPE_in )
!ESC!    use mod_adm, only: &
!ESC!       ADM_proc_stop,        &
!ESC!       ijdim => ADM_gall_in, &
!ESC!       kmin  => ADM_kmin,    &
!ESC!       kmax  => ADM_kmax
!ESC!    use mod_time, only: &
!ESC!       TIME_CTIME
!ESC!    use mod_runconf, only: &
!ESC!       NCRF,  &
!ESC!       NRBND, &
!ESC!       NRDIR
    use mod_rd_mstrnx, only: &
       rd_mstrnx_init
!ESC!    use mod_rd_mstrnx_ar5, only: &
!ESC!       rd_mstrnx_ar5_init
!ESC!    use mod_rd_mstrnx_ar5_crm, only: &
!ESC!       rd_mstrnx_ar5_crm_init
!ESC!    use mod_rd_mstrnx_cmip6, only: &
!ESC!       rd_mstrnx_cmip6_init
    implicit none

    NAMELIST / NM_RD_TIME / &
       TINTV, &
       TUNIT

    NAMELIST / NM_RD_SOLARCONSTANT / &
       OSOLCON_HIST, &
       SOLCON

    character(len=*), intent(in) :: RD_TYPE_in

    integer :: ierr
    !---------------------------------------------------------------------------

!ESC!    RD_TYPE = trim(RD_TYPE_in)

!ESC!    rewind(IO_FID_CONF)
!ESC!    read(IO_FID_CONF,nml=NM_RD_TIME,iostat=ierr)
!ESC!    if ( ierr < 0 ) then
!ESC!       write(IO_FID_LOG,*) '*** NM_RD_TIME is not specified. use default.'
!ESC!    elseif( ierr > 0 ) then
!ESC!       write(*,         *) 'xxx Not appropriate names in namelist NM_RD_TIME. STOP.'
!ESC!       write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist NM_RD_TIME. STOP.'
!ESC!       call ADM_proc_stop
!ESC!    endif
!ESC!    write(IO_FID_LOG,nml=NM_RD_TIME)
!ESC!
!ESC!    rewind(IO_FID_CONF)
!ESC!    read(IO_FID_CONF,nml=NM_RD_SOLARCONSTANT,iostat=ierr)
!ESC!    if ( ierr < 0 ) then
!ESC!       write(IO_FID_LOG,*) '*** NM_RD_SOLARCONSTANT is not specified. use default.'
!ESC!    elseif( ierr > 0 ) then
!ESC!       write(*,         *) 'xxx Not appropriate names in namelist NM_RD_SOLARCONSTANT. STOP.'
!ESC!       write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist NM_RD_SOLARCONSTANT. STOP.'
!ESC!       call ADM_proc_stop
!ESC!    endif
!ESC!    write(IO_FID_LOG,nml=NM_RD_SOLARCONSTANT)

    RD_TIME_PREV   = TIME_CTIME
    RD_TIME_ORIGIN = TIME_CTIME



    select case(RD_TYPE)
    case('NONE')
       !--- nothing
    case('MSTRNX')
       call rd_mstrnx_init( ADM_gall_in )
!ESC!    case('MSTRNX_AR5')
!ESC!       call rd_mstrnx_ar5_init    ( kmin, kmax, NCRF, NRDIR, NRBND )
!ESC!    case('MSTRNX_AR5_CRM')
!ESC!       call rd_mstrnx_ar5_crm_init( kmin, kmax, NCRF, NRDIR, NRBND )
!ESC!    case('MSTRNX_CMIP6')
!ESC!       call rd_mstrnx_cmip6_init( ijdim )
    case default
       write(*,*) 'xxx [rd_init/rd_driver] Not appropriate RD_TYPE. CHECK! ', trim(RD_TYPE)
       call ADM_proc_stop
    end select

    return
  end subroutine rd_init

  !-----------------------------------------------------------------------------
  subroutine rd_driver( &
       ijdim,          &
       rho,            &
       pre,            &
       tem,            &
       qv,             &
       q_clw,          &
       q_cli,          &
       qr,             &
       cloud_volume,   &
       cumulus_volume, &
       re_all,         &
       cfrac,          &
       q_cumclw,       &
       q_cumcli,       &
       cumfrac,        &
       qo3,            &
       tem_sfc,        &
       pre_sfc,        &
       zsfc,           &
       lat,            &
       lon,            &
       gralb,          &
       outqld,         &
       unccn,          &
       z,              &
       zh,             &
       RFLXSD,         &
       RFLXSU,         &
       RFLXLD,         &
       RFLXLU,         &
       DRFLXS,         &
       DRFLXL,         &
       RFSFCD,         &
       dfq_isccp,      &
       TAUC_ALL,       &
       TAUCL,          &
       TAUCI,          &
       TAUCLK,         &
       TAUCIK,         &
       RCEFF,          &
       RCEFF_SOLID,    &
       RCEFF_CLD,      &
       AERDFS,         &
       AERDFL,         &
       AERDFS_TRP,     &
       AERDFL_TRP,     &
       tbb_11um,       &
       aot_ext_vis,    & ! [Add] 2016/05/18 T.Seiki
       aot_abs_vis,    & ! [Add] 2016/05/18 T.Seiki
       update_flag     )
!ESC!    use mod_const, only: &
!ESC!       CONST_GRAV
    use mod_calendar, only: &
       CALENDAR_ointvl, &
       CALENDAR_ssaft,  &
       CALENDAR_ss2cc,  &
       CALENDAR_ss2yd,  &
       CALENDAR_dayyr
!ESC!    use mod_adm, only: &
!ESC!       kdim => ADM_kall, &
!ESC!       ADM_kmin,         &
!ESC!       ADM_kmax
!ESC!    use mod_grd, only: &
!ESC!       GRD_afact, &
!ESC!       GRD_bfact
!ESC!    use mod_time, only: &
!ESC!       TIME_CTIME
!ESC!    use mod_runconf, only: &
!ESC!       NCRF,                 &
!ESC!       NRBND,                &
!ESC!       NRDIR,                &
!ESC!       KAPCL,                &
!ESC!       NTAU  => NTAU_ISCCP,  &
!ESC!       NPRES => NPRES_ISCCP, &
!ESC!       HYDRO_MAX
    use mod_rd_mstrnx, only: &
       rd_mstrnx
!ESC!    use mod_rd_mstrnx_ar5, only: &
!ESC!       rd_mstrnx_ar5
!ESC!    use mod_rd_mstrnx_ar5_crm, only: &
!ESC!       rd_mstrnx_ar5_crm
!ESC!    use mod_rd_mstrnx_cmip6, only: &
!ESC!       rd_mstrnx_cmip6
    implicit none

    integer,  intent(in)    :: ijdim
    real(RP), intent(in)    :: rho           (ijdim,kdim)
    real(RP), intent(in)    :: pre           (ijdim,kdim)
    real(RP), intent(in)    :: tem           (ijdim,kdim)
    real(RP), intent(in)    :: qv            (ijdim,kdim)
    real(RP), intent(in)    :: q_clw         (ijdim,kdim)
    real(RP), intent(in)    :: q_cli         (ijdim,kdim)
    real(RP), intent(in)    :: qr            (ijdim,kdim)
    real(RP), intent(in)    :: cloud_volume  (ijdim,kdim,HYDRO_MAX)   ! [m3/m3]
    real(RP), intent(in)    :: cumulus_volume(ijdim,kdim,2)           ! [m3/m3]
    real(RP), intent(in)    :: re_all        (ijdim,kdim,HYDRO_MAX)   ! [m3/m3]
    real(RP), intent(in)    :: cfrac         (ijdim,kdim)             ! ratio of cloudy area
    real(RP), intent(in)    :: q_cumclw      (ijdim,kdim)             ! cloud water in cumulus
    real(RP), intent(in)    :: q_cumcli      (ijdim,kdim)
    real(RP), intent(in)    :: cumfrac       (ijdim)                  ! area rate of cumulus
    real(RP), intent(in)    :: qo3           (ijdim,kdim)
    real(RP), intent(in)    :: tem_sfc       (ijdim)
    real(RP), intent(in)    :: pre_sfc       (ijdim)
    real(RP), intent(in)    :: zsfc          (ijdim)
    real(RP), intent(in)    :: lat           (ijdim)
    real(RP), intent(in)    :: lon           (ijdim)
    real(RP), intent(in)    :: gralb         (ijdim,NRDIR,NRBND)
    real(RP), intent(in)    :: outqld        (ijdim,kdim,KAPCL)       ! aerosol mass conc
    real(RP), intent(in)    :: unccn         (ijdim,kdim)             ! CCN number conc.
    real(RP), intent(in)    :: z             (ijdim,kdim)
    real(RP), intent(in)    :: zh            (ijdim,kdim)
    real(RP), intent(inout) :: RFLXSD        (ijdim,kdim,NCRF)        ! downward short wave rad.(half lev)
    real(RP), intent(inout) :: RFLXSU        (ijdim,kdim,NCRF)        ! upward short wave(half lev)
    real(RP), intent(inout) :: RFLXLD        (ijdim,kdim,NCRF)        ! upward long wave rad.(half lev)
    real(RP), intent(inout) :: RFLXLU        (ijdim,kdim,NCRF)        ! upward long wave(half lev)
    real(RP), intent(inout) :: DRFLXL        (ijdim,kdim,NCRF)        ! long wave deriv.(half lev)
    real(RP), intent(inout) :: DRFLXS        (ijdim,kdim,NCRF)        ! short wave deriv.(half lev)
    real(RP), intent(inout) :: RFSFCD        (ijdim,NRDIR,NRBND,NCRF) ! surface downrad
    real(RP), intent(inout) :: dfq_isccp     (ijdim,NTAU,NPRES)       ! daytime fq_isccp
    real(RP), intent(inout) :: TAUC_ALL      (ijdim,HYDRO_MAX)
    real(RP), intent(inout) :: TAUCL         (ijdim)                  ! optical thickness of cloud (liquid)
    real(RP), intent(inout) :: TAUCI         (ijdim)                  ! optical thickness of cloud (ice)
    real(RP), intent(inout) :: TAUCLK        (ijdim,kdim)
    real(RP), intent(inout) :: TAUCIK        (ijdim,kdim)
    real(RP), intent(inout) :: RCEFF         (ijdim,kdim)             ! cloud particle effective radius
    real(RP), intent(inout) :: RCEFF_SOLID   (ijdim,kdim)             ! cloud particle effective radius for solid particle(qi,qs,qg)
    real(RP), intent(inout) :: RCEFF_CLD     (ijdim,kdim)             ! cloud particle effective radius for qc
    real(RP), intent(inout) :: AERDFS        (ijdim,kdim,NCRF)        ! aerosol radiative forcing (shortwave)
    real(RP), intent(inout) :: AERDFL        (ijdim,kdim,NCRF)        ! aerosol radiative forcing (longwave)
    real(RP), intent(inout) :: AERDFS_TRP    (ijdim,NCRF)             ! aerosol radiative forcing (shortwave, tropopause)
    real(RP), intent(inout) :: AERDFL_TRP    (ijdim,NCRF)             ! aerosol radiative forcing (longwave,  tropopause)
    real(RP), intent(inout) :: tbb_11um      (ijdim)                  ! black body temperature @ 11um
    real(RP), intent(inout) :: aot_ext_vis   (ijdim,KAPCL)            ! aerosol optical thickness (absorption+scatter) [Add] 2016/05/18 T.Seiki
    real(RP), intent(inout) :: aot_abs_vis   (ijdim,KAPCL)            ! aerosol optical thickness (only absorption)    [Add] 2016/05/18 T.Seiki
    logical,  intent(out)   :: update_flag                            ! radiation is updated?

    real(RP) :: SINS (ijdim)      ! solar incidence
    real(RP) :: COSZ (ijdim)      ! COS(solar zenith angle)
    real(RP) :: tem_h(ijdim,kdim) ! temperature (half lev)
    real(RP) :: pre_h(ijdim,kdim) ! pressure    (half lev)

    real(DP)          :: TIME_NEXT
    character(len=20) :: HTIME0
    character(len=20) :: HTIME1

    integer  :: IYEAR, IDOY, NDOY
    real(DP) :: year_frac

    integer  :: kmin, kmax
    real(RP) :: GRAV

    integer  :: ij, k, ip
    !---------------------------------------------------------------------------

    kmin = ADM_kmin
    kmax = ADM_kmax

    GRAV = CONST_GRAV

    update_flag = Calendar_ointvl( TIME_CTIME,     & ! [IN]
                                   RD_TIME_PREV ,  & ! [IN]
                                   RD_TIME_ORIGIN, & ! [IN]
                                   TINTV,          & ! [IN]
                                   TUNIT           ) ! [IN]

    update_flag = .true.

    if ( update_flag ) then
       RD_TIME_PREV = TIME_CTIME

       call CALENDAR_ssaft( TIME_NEXT, TIME_CTIME, TINTV, TUNIT )
       call CALENDAR_ss2cc( HTIME0, TIME_CTIME )
       call CALENDAR_ss2cc( HTIME1, TIME_NEXT  )

       write(IO_FID_LOG,*) '*** Radiation calculation : ', trim(HTIME0), ' - ', trim(HTIME1)

       ! solar incidence
       call rd_solarincidence( ijdim,      & ! [IN]
                               SINS(:),    & ! [OUT]
                               COSZ(:),    & ! [OUT]
                               TIME_CTIME, & ! [IN]
                               TIME_NEXT,  & ! [IN]
                               LON (:),    & ! [IN]
                               LAT (:)     ) ! [INOUT]

       call CALENDAR_ss2yd( IYEAR, IDOY, TIME_CTIME )
       call CALENDAR_dayyr( NDOY, IYEAR )
       year_frac = real(IYEAR,kind=DP) + real(IDOY-1,kind=DP) / real(NDOY,kind=DP)

       ! vartical interpolation
       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kmin,kmax,tem_h,tem,GRD_afact,GRD_bfact)
       do k  = kmin+1, kmax
       do ij = 1, ijdim
          tem_h(ij,k) = ( GRD_afact(k) * tem(ij,k  ) &
                        + GRD_bfact(k) * tem(ij,k-1) ) ! [mod] 20161026 H.Yashiro afac is halved
       enddo
       enddo
       !$omp end parallel do

       !$omp parallel do default(none),private(ij), &
       !$omp shared(ijdim,kmin,kmax,tem_h,tem_sfc,tem)
       do ij = 1, ijdim
          tem_h(ij,kmin-1) = 0.0_RP ! not used
          tem_h(ij,kmin)   = tem_sfc(ij)
          tem_h(ij,kmax+1) = tem(ij,kmax)
       enddo
       !$omp end parallel do

       !--- pressure interpolation
       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kmin,kmax,pre_h,pre,GRD_afact,GRD_bfact)
       do k  = kmin+1, kmax
       do ij = 1, ijdim
          pre_h(ij,k) = pre(ij,k  )**GRD_afact(k) &
                      * pre(ij,k-1)**GRD_bfact(k) ! [mod] 20161026 H.Yashiro afac is halved
       enddo
       enddo
       !$omp end parallel do

       !$omp parallel do default(none),private(ij), &
       !$omp shared(ijdim,kmin,kmax,pre_h,pre_sfc,pre,rho,zh,GRAV)
       do ij = 1, ijdim
          pre_h(ij,kmin-1) = 0.0_RP ! not used
          pre_h(ij,kmin)   = pre_sfc(ij)
          pre_h(ij,kmax+1) = pre(ij,kmax) - rho(ij,kmax) * GRAV * 0.5_RP * ( zh(ij,kmax+1) -zh(ij,kmax) )
       enddo
       !$omp end parallel do

       !$omp parallel default(none),private(ij,ip), &
       !$omp shared(ijdim,tbb_11um,aot_ext_vis,aot_abs_vis)

       !$omp do
       do ij = 1, ijdim
          tbb_11um(ij) = 0.0_RP
       enddo
       !$omp end do nowait

       !$omp do
       do ip = 1, KAPCL
       do ij = 1, ijdim
          aot_ext_vis(ij,ip) = 0.0_RP ! [Add] 2016/05/18 T.Seiki
          aot_abs_vis(ij,ip) = 0.0_RP ! [Add] 2016/05/18 T.Seiki
       enddo
       enddo
       !$omp end do

       !$omp end parallel

       select case( RD_TYPE )
       case('NONE')

          RFLXSU     (:,:,:)   = 0.0_RP
          RFLXSD     (:,:,:)   = 0.0_RP
          RFLXLU     (:,:,:)   = 0.0_RP
          RFLXLD     (:,:,:)   = 0.0_RP
          DRFLXL     (:,:,:)   = 0.0_RP
          DRFLXS     (:,:,:)   = 0.0_RP
          RFSFCD     (:,:,:,:) = 0.0_RP
          dfq_isccp  (:,:,:)   = 0.0_RP
          TAUC_ALL   (:,:)     = 0.0_RP
          TAUCL      (:)       = 0.0_RP
          TAUCI      (:)       = 0.0_RP
          TAUCLK     (:,:)     = 0.0_RP
          TAUCIK     (:,:)     = 0.0_RP
          RCEFF      (:,:)     = 0.0_RP
          RCEFF_SOLID(:,:)     = 0.0_RP
          RCEFF_CLD  (:,:)     = 0.0_RP
          AERDFS     (:,:,:)   = 0.0_RP
          AERDFL     (:,:,:)   = 0.0_RP
          AERDFS_TRP (:,:)     = 0.0_RP
          AERDFL_TRP (:,:)     = 0.0_RP
          tbb_11um   (:)       = 0.0_RP

       case('MSTRNX')

          call rd_mstrnx        ( ijdim,          & ! [IN]
                                  RFLXSD,         & ! [OUT]
                                  RFLXSU,         & ! [OUT]
                                  RFLXLD,         & ! [OUT]
                                  RFLXLU,         & ! [OUT]
                                  DRFLXS,         & ! [OUT]
                                  DRFLXL,         & ! [OUT]
                                  RFSFCD,         & ! [OUT]
                                  dfq_isccp,      & ! [OUT]
                                  TAUC_ALL,       & ! [OUT]
                                  TAUCL,          & ! [OUT]
                                  TAUCI,          & ! [OUT]
                                  TAUCLK,         & ! [OUT]
                                  TAUCIK,         & ! [OUT]
                                  tem_h,          & ! [IN]
                                  pre_h,          & ! [IN]
                                  zh,             & ! [IN]
                                  tem,            & ! [IN]
                                  pre,            & ! [IN]
                                  qv,             & ! [IN]
                                  qo3,            & ! [IN]
                                  cloud_volume,   & ! [IN]
                                  cumulus_volume, & ! [IN]
                                  re_all,         & ! [IN]
                                  cfrac,          & ! [IN]
                                  cumfrac,        & ! [IN]
                                  SINS,           & ! [IN]
                                  COSZ,           & ! [IN]
                                  gralb,          & ! [IN]
                                  tem_sfc,        & ! [IN]
                                  OUTQLD,         & ! [IN]
                                  AERDFS,         & ! [OUT]
                                  AERDFL,         & ! [OUT]
                                  AERDFS_TRP,     & ! [OUT]
                                  AERDFL_TRP,     & ! [OUT]
                                  tbb_11um,       & ! [OUT]
                                  aot_ext_vis,    & ! [out] [add] 2016/05/18 T.Seiki
                                  aot_abs_vis     ) ! [out] [add] 2016/05/18 T.Seiki

!ESC!       case('MSTRNX_AR5_CRM')
!ESC!
!ESC!          call rd_mstrnx_ar5_crm( ijdim,          & ! [IN]
!ESC!                                  RFLXSU,         & ! [OUT]
!ESC!                                  RFLXSD,         & ! [OUT]
!ESC!                                  RFLXLU,         & ! [OUT]
!ESC!                                  RFLXLD,         & ! [OUT]
!ESC!                                  DRFLXS,         & ! [OUT]
!ESC!                                  DRFLXL,         & ! [OUT]
!ESC!                                  RFSFCD,         & ! [OUT]
!ESC!                                  dfq_isccp,      & ! [OUT]
!ESC!                                  TAUC_ALL,       & ! [OUT]
!ESC!                                  TAUCL,          & ! [OUT]
!ESC!                                  TAUCI,          & ! [OUT]
!ESC!                                  TAUCLK,         & ! [OUT]
!ESC!                                  TAUCIK,         & ! [OUT]
!ESC!                                  tem_h,          & ! [IN]
!ESC!                                  pre_h,          & ! [IN]
!ESC!                                  zh,             & ! [IN]
!ESC!                                  tem,            & ! [IN]
!ESC!                                  pre,            & ! [IN]
!ESC!                                  qv,             & ! [IN]
!ESC!                                  qo3,            & ! [IN]
!ESC!                                  cloud_volume,   & ! [IN]
!ESC!                                  cumulus_volume, & ! [IN]
!ESC!                                  re_all,         & ! [IN]
!ESC!                                  cfrac,          & ! [IN]
!ESC!                                  cumfrac,        & ! [IN]
!ESC!                                  SINS,           & ! [IN]
!ESC!                                  COSZ,           & ! [IN]
!ESC!                                  gralb,          & ! [IN]
!ESC!                                  tem_sfc,        & ! [IN]
!ESC!                                  OUTQLD,         & ! [IN]
!ESC!                                  AERDFS,         & ! [OUT]
!ESC!                                  AERDFL,         & ! [OUT]
!ESC!                                  AERDFS_TRP,     & ! [OUT]
!ESC!                                  AERDFL_TRP,     & ! [OUT]
!ESC!                                  tbb_11um        ) ! [OUT]
!ESC!
!ESC!       case('MSTRNX_AR5')
!ESC!
!ESC!          call rd_mstrnx_ar5    ( ijdim,          & ! [IN]
!ESC!                                  RFLXSU,         & ! [OUT]
!ESC!                                  RFLXSD,         & ! [OUT]
!ESC!                                  RFLXLU,         & ! [OUT]
!ESC!                                  RFLXLD,         & ! [OUT]
!ESC!                                  DRFLXS,         & ! [OUT]
!ESC!                                  DRFLXL,         & ! [OUT]
!ESC!                                  RFSFCD,         & ! [OUT]
!ESC!                                  dfq_isccp,      & ! [OUT]
!ESC!                                  TAUCL,          & ! [OUT]
!ESC!                                  TAUCI,          & ! [OUT]
!ESC!                                  TAUCLK,         & ! [OUT]
!ESC!                                  TAUCIK,         & ! [OUT]
!ESC!                                  RCEFF,          & ! [INOUT]
!ESC!                                  RCEFF_SOLID,    & ! [INOUT]
!ESC!                                  RCEFF_CLD,      & ! [INOUT]
!ESC!                                  tem_h,          & ! [IN]
!ESC!                                  pre_h,          & ! [IN]
!ESC!                                  zh,             & ! [IN]
!ESC!                                  tem,            & ! [IN]
!ESC!                                  pre,            & ! [IN]
!ESC!                                  qv,             & ! [IN]
!ESC!                                  qo3,            & ! [IN]
!ESC!                                  q_clw,          & ! [IN]
!ESC!                                  q_cli,          & ! [IN]
!ESC!                                  qr,             & ! [IN]
!ESC!                                  cfrac,          & ! [IN]
!ESC!                                  q_cumclw,       & ! [IN]
!ESC!                                  q_cumcli,       & ! [IN]
!ESC!                                  cumfrac,        & ! [IN]
!ESC!                                  SINS,           & ! [IN]
!ESC!                                  COSZ,           & ! [INOUT]
!ESC!                                  gralb,          & ! [IN]
!ESC!                                  tem_sfc,        & ! [IN]
!ESC!                                  OUTQLD,         & ! [IN]
!ESC!                                  UNCCN,          & ! [IN]
!ESC!                                  AERDFS,         & ! [INOUT]
!ESC!                                  AERDFL,         & ! [INOUT]
!ESC!                                  AERDFS_TRP,     & ! [INOUT]
!ESC!                                  AERDFL_TRP      ) ! [INOUT]
!ESC!
!ESC!       case('MSTRNX_CMIP6')
!ESC!
!ESC!          call rd_mstrnx_cmip6  ( ijdim,          & ! [IN]
!ESC!                                  RFLXSD,         & ! [OUT]
!ESC!                                  RFLXSU,         & ! [OUT]
!ESC!                                  RFLXLD,         & ! [OUT]
!ESC!                                  RFLXLU,         & ! [OUT]
!ESC!                                  DRFLXS,         & ! [OUT]
!ESC!                                  DRFLXL,         & ! [OUT]
!ESC!                                  RFSFCD,         & ! [OUT]
!ESC!                                  dfq_isccp,      & ! [OUT]
!ESC!                                  TAUC_ALL,       & ! [OUT]
!ESC!                                  TAUCL,          & ! [OUT]
!ESC!                                  TAUCI,          & ! [OUT]
!ESC!                                  TAUCLK,         & ! [OUT]
!ESC!                                  TAUCIK,         & ! [OUT]
!ESC!                                  tem_h,          & ! [IN]
!ESC!                                  pre_h,          & ! [IN]
!ESC!                                  zh,             & ! [IN]
!ESC!                                  tem,            & ! [IN]
!ESC!                                  pre,            & ! [IN]
!ESC!                                  z,              & ! [IN]
!ESC!                                  qv,             & ! [IN]
!ESC!                                  qo3,            & ! [IN]
!ESC!                                  cloud_volume,   & ! [IN]
!ESC!                                  cumulus_volume, & ! [IN]
!ESC!                                  re_all,         & ! [IN]
!ESC!                                  cfrac,          & ! [IN]
!ESC!                                  cumfrac,        & ! [IN]
!ESC!                                  SINS,           & ! [IN]
!ESC!                                  COSZ,           & ! [IN]
!ESC!                                  zsfc,           & ! [IN]
!ESC!                                  lat,            & ! [IN]
!ESC!                                  lon,            & ! [IN]
!ESC!                                  year_frac,      & ! [IN]
!ESC!                                  gralb,          & ! [IN]
!ESC!                                  tem_sfc,        & ! [IN]
!ESC!                                  OUTQLD,         & ! [IN]
!ESC!                                  AERDFS,         & ! [OUT]
!ESC!                                  AERDFL,         & ! [OUT]
!ESC!                                  AERDFS_TRP,     & ! [OUT]
!ESC!                                  AERDFL_TRP,     & ! [OUT]
!ESC!                                  tbb_11um        ) ! [OUT]

       end select

    endif ! update_flag?

    return
  end subroutine rd_driver

!ESC!  !-----------------------------------------------------------------------------
!ESC!  subroutine rd_tendency( &
!ESC!       ijdim,     &
!ESC!       RFLXSD,    &
!ESC!       RFLXSU,    &
!ESC!       RFLXLD,    &
!ESC!       RFLXLU,    &
!ESC!       DRFLXS,    &
!ESC!       DRFLXL,    &
!ESC!       RFLXSG,    &
!ESC!       RFLXLG,    &
!ESC!       rhog,      &
!ESC!       frhoge,    &
!ESC!       frhogetot, &
!ESC!       gam2h      )
!ESC!    use mod_const, only: &
!ESC!       CONST_CVdry
!ESC!    use mod_adm, only: &
!ESC!       kdim => ADM_kall, &
!ESC!       ADM_kmin,         &
!ESC!       ADM_kmax
!ESC!    use mod_grd, only: &
!ESC!       GRD_rdgz
!ESC!    use mod_runconf, only: &
!ESC!       NCRF
!ESC!    use mod_history, only: &
!ESC!       history_in
!ESC!    implicit none
!ESC!
!ESC!    integer,  intent(in)    :: ijdim
!ESC!    real(RP), intent(in)    :: RFLXSD   (ijdim,kdim,NCRF)
!ESC!    real(RP), intent(inout) :: RFLXSU   (ijdim,kdim,NCRF)
!ESC!    real(RP), intent(in)    :: RFLXLD   (ijdim,kdim,NCRF)
!ESC!    real(RP), intent(inout) :: RFLXLU   (ijdim,kdim,NCRF)
!ESC!    real(RP), intent(in)    :: DRFLXS   (ijdim,kdim,NCRF)
!ESC!    real(RP), intent(in)    :: DRFLXL   (ijdim,kdim,NCRF)
!ESC!    real(RP), intent(in)    :: RFLXSG   (ijdim)
!ESC!    real(RP), intent(in)    :: RFLXLG   (ijdim)
!ESC!    real(RP), intent(in)    :: rhog     (ijdim,kdim)
!ESC!    real(RP), intent(inout) :: frhoge   (ijdim,kdim)
!ESC!    real(RP), intent(inout) :: frhogetot(ijdim,kdim)
!ESC!    real(RP), intent(in)    :: gam2h    (ijdim,kdim)
!ESC!
!ESC!    real(RP) :: DRSG  (ijdim)
!ESC!    real(RP) :: DRLG  (ijdim)
!ESC!    real(RP) :: drhoge(ijdim,kdim)
!ESC!    real(RP) :: hr    (ijdim,kdim)
!ESC!    real(RP) :: rrhodz(ijdim,kdim)
!ESC!
!ESC!    integer  :: kmin, kmax, kdim_
!ESC!    real(RP) :: cvdry
!ESC!
!ESC!    integer  :: ij, k, ic
!ESC!    !---------------------------------------------------------------------------
!ESC!
!ESC!    kdim_ = kdim
!ESC!    kmin  = ADM_kmin
!ESC!    kmax  = ADM_kmax
!ESC!    cvdry = CONST_CVdry
!ESC!
!ESC!    if ( RD_TYPE /= 'NONE') then
!ESC!
!ESC!       !$omp parallel do default(none),private(ij), &
!ESC!       !$omp shared(ijdim,kmin,DRSG,RFLXSG,RFLXSU,DRFLXS)
!ESC!       do ij = 1, ijdim
!ESC!          DRSG(ij) = 0.0_RP
!ESC!          if ( DRFLXS(ij,kmin,1) /= 0.0_RP ) then
!ESC!             DRSG(ij) = ( RFLXSG(ij) - RFLXSU(ij,kmin,1) ) / DRFLXS(ij,kmin,1)
!ESC!          endif
!ESC!       enddo
!ESC!       !$omp end parallel do
!ESC!
!ESC!       !$omp parallel do default(none),private(ij), &
!ESC!       !$omp shared(ijdim,kmin,DRLG,RFLXLG,RFLXLU,DRFLXL)
!ESC!       do ij = 1, ijdim
!ESC!          DRLG(ij) = 0.0_RP
!ESC!          if ( DRFLXL(ij,kmin,1) /= 0.0_RP ) then
!ESC!             DRLG(ij) = ( RFLXLG(ij) - RFLXLU(ij,kmin,1) ) / DRFLXL(ij,kmin,1)
!ESC!          endif
!ESC!       enddo
!ESC!       !$omp end parallel do
!ESC!
!ESC!       !$omp parallel do default(none),private(ij,k,ic), &
!ESC!       !$omp shared(ijdim,kdim_,RFLXSU,RFLXLU,DRSG,DRLG,DRFLXS,DRFLXL)
!ESC!       do ic = 1, NCRF
!ESC!       do k  = 1, kdim_
!ESC!       do ij = 1, ijdim
!ESC!          RFLXSU(ij,k,ic) = RFLXSU(ij,k,ic) + DRSG(ij) * DRFLXS(ij,k,ic)
!ESC!          RFLXLU(ij,k,ic) = RFLXLU(ij,k,ic) + DRLG(ij) * DRFLXL(ij,k,ic)
!ESC!       enddo
!ESC!       enddo
!ESC!       enddo
!ESC!       !$omp end parallel do
!ESC!
!ESC!       call rd_heatingrate( ijdim,         & ! [IN]
!ESC!                            RFLXSD(:,:,1), & ! [IN]
!ESC!                            RFLXSU(:,:,1), & ! [IN]
!ESC!                            RFLXLD(:,:,1), & ! [IN]
!ESC!                            RFLXLU(:,:,1), & ! [IN]
!ESC!                            drhoge(:,:),   & ! [OUT]
!ESC!                            gam2h (:,:)    ) ! [IN]
!ESC!
!ESC!       !$omp parallel do default(none),private(ij,k), &
!ESC!       !$omp shared(ijdim,kmin,kmax,frhoge,frhogetot,drhoge)
!ESC!       do k  = kmin, kmax
!ESC!       do ij = 1, ijdim
!ESC!          frhoge   (ij,k) = frhoge   (ij,k) + drhoge(ij,k)
!ESC!          frhogetot(ij,k) = frhogetot(ij,k) + drhoge(ij,k)
!ESC!       enddo
!ESC!       enddo
!ESC!       !$omp end parallel do
!ESC!
!ESC!    endif
!ESC!
!ESC!    !$omp parallel do default(none),private(ij,k), &
!ESC!    !$omp shared(ijdim,kmin,kmax,rrhodz,GRD_rdgz,rhog,cvdry)
!ESC!    do k  = kmin, kmax
!ESC!    do ij = 1, ijdim
!ESC!       rrhodz(ij,k) = GRD_rdgz(k) / ( rhog(ij,k) * cvdry )
!ESC!    enddo
!ESC!    enddo
!ESC!    !$omp end parallel do
!ESC!
!ESC!    !$omp parallel workshare
!ESC!    hr(:,kmin-1) = 0.0_RP
!ESC!    hr(:,kmax+1) = 0.0_RP
!ESC!    !$omp end parallel workshare
!ESC!
!ESC!    !--- short wave heating rate (K/s)
!ESC!
!ESC!    !$omp parallel do default(none),private(ij,k), &
!ESC!    !$omp shared(ijdim,kmin,kmax,hr,RFLXSU,RFLXSD,rrhodz,gam2h)
!ESC!    do k  = kmin, kmax
!ESC!    do ij = 1, ijdim
!ESC!       hr(ij,k) = ( gam2h(ij,k) * RFLXSU(ij,k,1) - gam2h(ij,k+1) * RFLXSU(ij,k+1,1) &
!ESC!                  - gam2h(ij,k) * RFLXSD(ij,k,1) + gam2h(ij,k+1) * RFLXSD(ij,k+1,1) &
!ESC!                  ) * rrhodz(ij,k)
!ESC!    enddo
!ESC!    enddo
!ESC!    !$omp end parallel do
!ESC!    call history_in('ml_swhr',hr(:,:))
!ESC!
!ESC!
!ESC!    !$omp parallel do default(none),private(ij,k), &
!ESC!    !$omp shared(ijdim,kmin,kmax,hr,RFLXSU,RFLXSD,rrhodz,gam2h)
!ESC!    do k  = kmin, kmax
!ESC!    do ij = 1, ijdim
!ESC!       hr(ij,k) = ( gam2h(ij,k) * RFLXSU(ij,k,2) - gam2h(ij,k+1) * RFLXSU(ij,k+1,2) &
!ESC!                  - gam2h(ij,k) * RFLXSD(ij,k,2) + gam2h(ij,k+1) * RFLXSD(ij,k+1,2) &
!ESC!                  ) * rrhodz(ij,k)
!ESC!    enddo
!ESC!    enddo
!ESC!    !$omp end parallel do
!ESC!    call history_in('ml_swhr_c',hr(:,:))
!ESC!
!ESC!    !--- long wave heating rate [K/s]
!ESC!
!ESC!    !$omp parallel do default(none),private(ij,k), &
!ESC!    !$omp shared(ijdim,kmin,kmax,hr,RFLXLU,RFLXLD,rrhodz,gam2h)
!ESC!    do k  = kmin, kmax
!ESC!    do ij = 1, ijdim
!ESC!       hr(ij,k) = ( gam2h(ij,k) * RFLXLU(ij,k,1) - gam2h(ij,k+1) * RFLXLU(ij,k+1,1) &
!ESC!                  - gam2h(ij,k) * RFLXLD(ij,k,1) + gam2h(ij,k+1) * RFLXLD(ij,k+1,1) &
!ESC!                  ) * rrhodz(ij,k)
!ESC!    enddo
!ESC!    enddo
!ESC!    !$omp end parallel do
!ESC!    call history_in('ml_lwhr',hr(:,:))
!ESC!
!ESC!    !$omp parallel do default(none),private(ij,k), &
!ESC!    !$omp shared(ijdim,kmin,kmax,hr,RFLXLU,RFLXLD,rrhodz,gam2h)
!ESC!    do k  = kmin, kmax
!ESC!    do ij = 1, ijdim
!ESC!       hr(ij,k) = ( gam2h(ij,k) * RFLXLU(ij,k,2) - gam2h(ij,k+1) * RFLXLU(ij,k+1,2) &
!ESC!                  - gam2h(ij,k) * RFLXLD(ij,k,2) + gam2h(ij,k+1) * RFLXLD(ij,k+1,2) &
!ESC!                  ) * rrhodz(ij,k)
!ESC!    enddo
!ESC!    enddo
!ESC!    !$omp end parallel do
!ESC!    call history_in('ml_lwhr_c',hr(:,:))
!ESC!
!ESC!    return
!ESC!  end subroutine rd_tendency
!ESC!
!ESC!  !-----------------------------------------------------------------------------
!ESC!  subroutine rd_heatingrate( &
!ESC!       ijdim,    &
!ESC!       rflux_sd, &
!ESC!       rflux_su, &
!ESC!       rflux_ld, &
!ESC!       rflux_lu, &
!ESC!       drhoge,   &
!ESC!       gam2h     )
!ESC!    use mod_adm, only: &
!ESC!       kdim => ADM_kall, &
!ESC!       ADM_kmin,         &
!ESC!       ADM_kmax
!ESC!    use mod_grd, only: &
!ESC!       GRD_rdgz
!ESC!    implicit none
!ESC!
!ESC!    integer,  intent(in)  :: ijdim
!ESC!    real(RP), intent(in)  :: rflux_sd(ijdim,kdim)
!ESC!    real(RP), intent(in)  :: rflux_su(ijdim,kdim)
!ESC!    real(RP), intent(in)  :: rflux_ld(ijdim,kdim)
!ESC!    real(RP), intent(in)  :: rflux_lu(ijdim,kdim)
!ESC!    real(RP), intent(out) :: drhoge  (ijdim,kdim)
!ESC!    real(RP), intent(in)  :: gam2h   (ijdim,kdim)
!ESC!
!ESC!    integer  :: kmin, kmax
!ESC!
!ESC!    integer  :: ij, k
!ESC!    !---------------------------------------------------------------------------
!ESC!
!ESC!    kmin  = ADM_kmin
!ESC!    kmax  = ADM_kmax
!ESC!
!ESC!    !$omp parallel do default(none),private(ij,k), &
!ESC!    !$omp shared(ijdim,kmin,kmax,drhoge,rflux_su,rflux_sd,rflux_lu,rflux_ld,GRD_rdgz,gam2h)
!ESC!    do k  = kmin, kmax
!ESC!    do ij = 1, ijdim
!ESC!       drhoge(ij,k) = ( gam2h(ij,k) * rflux_su(ij,k) - gam2h(ij,k+1) * rflux_su(ij,k+1) &
!ESC!                      - gam2h(ij,k) * rflux_sd(ij,k) + gam2h(ij,k+1) * rflux_sd(ij,k+1) &
!ESC!                      + gam2h(ij,k) * rflux_lu(ij,k) - gam2h(ij,k+1) * rflux_lu(ij,k+1) &
!ESC!                      - gam2h(ij,k) * rflux_ld(ij,k) + gam2h(ij,k+1) * rflux_ld(ij,k+1) &
!ESC!                      ) * GRD_rdgz(k)
!ESC!    enddo
!ESC!    enddo
!ESC!    !$omp end parallel do
!ESC!
!ESC!    !$omp parallel workshare
!ESC!    drhoge(:,kmin-1) = 0.0_RP
!ESC!    drhoge(:,kmax+1) = 0.0_RP
!ESC!    !$omp end parallel workshare
!ESC!
!ESC!    return
!ESC!  end subroutine rd_heatingrate

  !-----------------------------------------------------------------------------
  subroutine rd_solarincidence( &
       ijdim, &
       SINS,  &
       COSZ,  &
       TIMES, &
       TIMEE, &
       ALON,  &
       ALAT   )
!ESC!    use mod_adm, only: &
!ESC!       ADM_proc_stop
!ESC!    use mod_const, only: &
!ESC!       CONST_PI
    use mod_calendar, only: &
       CALENDAR_secdy, &
       CALENDAR_ss2ds, &
       CALENDAR_ss2yd, &
       CALENDAR_ym2dd, &
       CALENDAR_ds2ss, &
       CALENDAR_dayyr, &
       CALENDAR_perpr
    implicit none

    integer,  intent(in)  :: ijdim       ! number of horizontal grid
    real(RP), intent(out) :: SINS(ijdim) ! flux of incidence
    real(RP), intent(out) :: COSZ(ijdim) ! cos(incidnt angle)
    real(DP), intent(in)  :: TIMES       ! start time
    real(DP), intent(in)  :: TIMEE       ! finish time
    real(RP), intent(in)  :: ALON(ijdim) ! longitude
    real(RP), intent(in)  :: ALAT(ijdim) ! latitude

    real(RP), save :: E         = 0.01672_RP ! eccentricity
    real(RP), save :: EPSD      =   23.45_RP ! orbital incline angle
    real(RP), save :: VPID      =  102.04_RP ! longitude of periherion
    integer,  save :: MDORIG(2) = (/3,21/)   ! date of origin
    logical,  save :: OYRAVR    = .false.    ! day-av solar incidence
    logical,  save :: ODYAVR    = .false.    ! day-av solar incidence
    real(RP), save :: AZET      = 0.41_RP    ! fact. ann-av incidence
    real(RP), save :: BZET      = 0.59_RP    ! fact. ann-av incidence
    integer,  save :: NHSUB     = 12         ! # of averaging angle
    logical,  save :: OSET      = .false.    ! special setting

    real(RP), save :: BETA, A1, A2, A3, B1, B2, B3
    real(RP), save :: EPS, VPI, ALM0, R0
    integer,  save :: IYEAR0, IMONTH0, IDAYS0
    logical,  save :: OPERP
    logical,  save :: OFIRST = .true. ! intial flag

    NAMELIST / NM_RD_SOLAR / &
       E     , &
       EPSD  , &
       VPID  , &
       MDORIG, &
       OYRAVR, &
       ODYAVR, &
       AZET  , &
       BZET  , &
       NHSUB , &
       OSET

    real(RP) :: DELTS, ANGHR, ANG0, COSZ0
    real(RP) :: RSECYR, ALM, EM, V, RR, SOLINS
    integer  :: IDAYS, NSECDY, NDAYYR, IDAYY, IDAY0, IYEAR, NHSUBX
    real(DP) :: TIMESX, TIMEEX, TIMEX, TIME0
    real(DP) :: RSEC

    integer  :: ij, IH
    integer  :: ierr
    !---------------------------------------------------------------------------

    if ( ofirst ) then
       ofirst = .false.

!ESC!       !--- read parameters
!ESC!       write(IO_FID_LOG,*)
!ESC!       write(IO_FID_LOG,*) '+++ Module[RD DRIVER]/Category[nhm physics]'
!ESC!       rewind(IO_FID_CONF)
!ESC!       read(IO_FID_CONF,nml=NM_RD_SOLAR,iostat=ierr)
!ESC!       if ( ierr < 0 ) then
!ESC!          write(IO_FID_LOG,*) '*** NM_RD_SOLAR is not specified. use default.'
!ESC!       elseif( ierr > 0 ) then
!ESC!          write(*,         *) 'xxx Not appropriate names in namelist NM_RD_SOLAR. STOP.'
!ESC!          write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist NM_RD_SOLAR. STOP.'
!ESC!          call ADM_proc_stop
!ESC!       endif
!ESC!       write(IO_FID_LOG,nml=NM_RD_SOLAR)

       !--- Ref. Berger(1978)
       BETA = sqrt( 1.0_RP - E**2 )
       A1   = -2.0_RP * ( 1.0_RP/2.0_RP * E + 1.0_RP/8.0_RP * E**3 )* ( 1.0_RP        + BETA )
       A2   = -2.0_RP * (-1.0_RP/4.0_RP * E**2 )                    * ( 1.0_RP/2.0_RP + BETA )
       A3   = -2.0_RP * ( 1.0_RP/8.0_RP * E**3 )                    * ( 1.0_RP/3.0_RP + BETA )
       B1   =  2.0_RP * E - 1.0_RP/4.0_RP * E**3
       B2   =  5.0_RP / 4.0_RP  * E**2
       B3   = 13.0_RP / 12.0_RP * E**3

       EPS  = EPSD / 180.0_RP * CONST_PI
       VPI  = ( VPID + 180.0_RP ) / 180.0_RP * CONST_PI
       ALM0 = A1*sin(-VPI) + A2*sin(-2.0_RP*VPI) + A3*sin(-3.0_RP*VPI)
       R0   = 1.0_RP - E**2

       call CALENDAR_perpr( IYEAR0, IMONTH0, IDAY0,   OPERP  )
       call CALENDAR_ym2dd( IDAYS0, IYEAR0,  IMONTH0, IDAY0  )
    endif

    !--- get seconds for 1 day
    call Calendar_SECDY( NSECDY )

    !--- Finish time
    if ( ODYAVR ) then
       TIMEEX = TIMES + real(NSECDY,kind=DP)
    else
       TIMEEX = TIMEE
    endif

    if ( TIMEEX /= TIMES ) then
       NHSUBX = NHSUB
    else
       NHSUBX = 1
    endif

    !$omp parallel workshare
    sins(:) = 0.0_RP
    COSZ(:) = 0.0_RP
    !$omp end parallel workshare

    if ( OSET ) then ! special set

       !$omp parallel workshare
       COSZ(:) = cos( ALAT(:) )
       SINS(:) = SOLCON
       !$omp end parallel workshare

    elseif( OYRAVR .and. ODYAVR ) then ! DYAVR : annual and daily avg.

       !$omp parallel workshare
       COSZ(:) = AZET + BZET*cos( ALAT(:) )**2
       SINS(:) = SOLCON / CONST_PI
       !$omp end parallel workshare

    elseif( OYRAVR ) then ! YRAVR : annual avg. with diurnal cycle

       do IH = 1, NHSUBX
          TIMEX = TIMES + ( TIMEEX - TIMES ) * (IH-1) / NHSUB
          call Calendar_SS2DS( IDAYS, RSEC, TIMEX )
          ANG0  = RSEC / real(NSECDY,kind=RP) * 2.0_RP * CONST_PI - CONST_PI

          !$omp parallel do default(none),private(ij,ANGHR,COSZ0), &
          !$omp shared(ijdim,ANG0,ALON,AZET,BZET,ALAT,COSZ,SINS)
          do ij = 1, ijdim
             ANGHR = ANG0 + ALON(ij)
             COSZ0 = ( AZET + BZET*cos( ALAT(ij) )**2 ) * cos(ANGHR)
             if ( COSZ0 > 0.0_RP ) then
                COSZ(ij) = COSZ(ij) + COSZ0
                SINS(ij) = SINS(ij) + 1.0_RP
             endif
          enddo
          !$omp end parallel do
       enddo

       !$omp parallel do default(none),private(ij), &
       !$omp shared(ijdim,COSZ,SINS,SOLCON,NHSUBX)
       do ij = 1, ijdim
          if ( SINS(ij) > 0.0_RP ) then
             COSZ(ij) = COSZ(ij) / SINS(ij)
             SINS(ij) = SOLCON * SINS(ij) / NHSUBX
          else
             COSZ(ij) = -1.0_RP
             SINS(ij) =  0.0_RP
          endif
       enddo
       !$omp end parallel do

    else ! MARCH : with seasonal cycle

       TIMESX = TIMES
       call CALENDAR_ss2yd( IYEAR, IDAYY, TIMES )
       call CALENDAR_ym2dd( IDAY0, IYEAR, MDORIG(1), MDORIG(2) )
       call CALENDAR_ds2ss( TIME0, IDAY0, 0.0_DP  )
       call CALENDAR_dayyr( NDAYYR, IYEAR )
       if ( OPERP ) then
          CALL CALENDAR_ym2dd( IDAY0,  IYEAR0, MDORIG(1), MDORIG(2) )
          CALL CALENDAR_ds2ss( TIME0,  IDAY0,  0.0_DP  )
          CALL CALENDAR_ss2ds( IDAYS,  RSEC,   TIMES )
          CALL CALENDAR_ds2ss( TIMESX, IDAYS0, RSEC  )
          CALL CALENDAR_dayyr( NDAYYR, IYEAR0 )
       endif

       RSECYR = real(NSECDY*NDAYYR,kind=RP)

       do IH = 1, NHSUBX
          TIMEX = TIMESX + ( TIMEEX - TIMES ) * (IH-1) / NHSUB
          call CALENDAR_SS2DS( IDAYS, RSEC, TIMEX )
          ANG0  = RSEC / real(NSECDY,kind=RP) * 2.0_RP * CONST_PI - CONST_PI

          ALM    = ( TIMEX - TIME0 ) / RSECYR * 2.0_RP * CONST_PI + ALM0
          EM     = ALM - VPI
          V      = EM + B1*sin(EM) + B2*sin(2.0_RP*EM) + B3*sin(3.0_RP*EM)
          DELTS  = asin( sin(EPS) * sin(V+VPI) )
          RR     = R0 / ( 1.0_RP + E*cos(V) )
          SOLINS = SOLCON / (RR**2)

          !$omp parallel do default(none),private(ij,ANGHR,COSZ0), &
          !$omp shared(ijdim,ANG0,ALON,ALAT,DELTS,COSZ,SINS)
          do ij = 1, ijdim
             ANGHR = ANG0 + ALON(ij)
             COSZ0 = sin( ALAT(ij) ) * sin( DELTS ) &
                   + cos( ALAT(ij) ) * cos( DELTS ) * cos( ANGHR )

             if ( COSZ0 > 0.0_RP ) then
                COSZ(ij) = COSZ(ij) + COSZ0
                SINS(ij) = SINS(ij) + 1.0_RP
             endif
          enddo
          !$omp end parallel do
       enddo

       !$omp parallel do default(none),private(ij), &
       !$omp shared(ijdim,COSZ,SINS,SOLINS,SOLCON,NHSUBX)
       do ij = 1, ijdim
          if ( SINS(ij) > 0.0_RP ) then
             COSZ(ij) = COSZ(ij) / SINS(ij)
             SINS(ij) = SOLINS * SINS(ij) / NHSUBX
          else
             COSZ(ij) = -1.0_RP
             SINS(ij) =  0.0_RP
          endif
       enddo
       !$omp end parallel do
    endif

    return
  end subroutine rd_solarincidence

!ESC!  !-----------------------------------------------------------------------------
!ESC!  subroutine rd_surface( &
!ESC!       ijdim,        &
!ESC!       rflux_sfc_sd, &
!ESC!       rflux_sfc_su, &
!ESC!       rflux_sfc_ld, &
!ESC!       rflux_sfc_lu, &
!ESC!       drflux_dT,    &
!ESC!       rflux_sfc_d,  &
!ESC!       albedo_sfc,   &
!ESC!       tem_sfc       )
!ESC!    use mod_const, only: &
!ESC!       CONST_STB
!ESC!    use mod_runconf, only: &
!ESC!       SF_TYPE,       &
!ESC!       NRBND,         &
!ESC!       NRDIR,         &
!ESC!       NRDIR_DIRECT,  &
!ESC!       NRDIR_DIFFUSE, &
!ESC!       NRBND_VIS,     &
!ESC!       NRBND_NIR,     &
!ESC!       NRBND_IR
!ESC!    implicit none
!ESC!
!ESC!    integer,  intent(in)  :: ijdim
!ESC!    real(RP), intent(out) :: rflux_sfc_sd(ijdim)             ! surface downward SW flux
!ESC!    real(RP), intent(out) :: rflux_sfc_su(ijdim)             ! surface upward   SW flux
!ESC!    real(RP), intent(out) :: rflux_sfc_ld(ijdim)             ! surface downward LW flux
!ESC!    real(RP), intent(out) :: rflux_sfc_lu(ijdim)             ! surface upward   LW flux
!ESC!    real(RP), intent(out) :: drflux_dT   (ijdim)             ! derivative of surface LW emission
!ESC!    real(RP), intent(in)  :: rflux_sfc_d (ijdim,NRDIR,NRBND) ! surface downward flux (DIR,DIF x VIS,NIR,IR)
!ESC!    real(RP), intent(in)  :: albedo_sfc  (ijdim,NRDIR,NRBND) ! surface albedo
!ESC!    real(RP), intent(in)  :: tem_sfc     (ijdim)             ! surface temperature
!ESC!
!ESC!    real(RP) :: emission ! surface LW emission
!ESC!
!ESC!    integer  :: ij
!ESC!    !---------------------------------------------------------------------------
!ESC!
!ESC!    if ( SF_TYPE == 'NO-FLUX' ) then
!ESC!       !$omp parallel do default(none),private(ij), &
!ESC!       !$omp shared(ijdim,rflux_sfc_ld,rflux_sfc_lu,rflux_sfc_sd,rflux_sfc_su)
!ESC!       do ij = 1, ijdim
!ESC!          rflux_sfc_ld(ij) = 0.0_RP
!ESC!          rflux_sfc_lu(ij) = 0.0_RP
!ESC!          rflux_sfc_sd(ij) = 0.0_RP
!ESC!          rflux_sfc_su(ij) = 0.0_RP
!ESC!       enddo
!ESC!       !$omp end parallel do
!ESC!    else
!ESC!
!ESC!       ! SW
!ESC!
!ESC!       !$omp parallel do default(none),private(ij), &
!ESC!       !$omp shared(ijdim,rflux_sfc_sd,rflux_sfc_su,rflux_sfc_d,albedo_sfc)
!ESC!       do ij = 1, ijdim
!ESC!          rflux_sfc_sd(ij) = rflux_sfc_d(ij,NRDIR_DIRECT ,NRBND_VIS) &
!ESC!                           + rflux_sfc_d(ij,NRDIR_DIFFUSE,NRBND_VIS) &
!ESC!                           + rflux_sfc_d(ij,NRDIR_DIRECT ,NRBND_NIR) &
!ESC!                           + rflux_sfc_d(ij,NRDIR_DIFFUSE,NRBND_NIR)
!ESC!
!ESC!          rflux_sfc_su(ij) = rflux_sfc_d(ij,NRDIR_DIRECT ,NRBND_VIS) * albedo_sfc(ij,NRDIR_DIRECT ,NRBND_VIS) &
!ESC!                           + rflux_sfc_d(ij,NRDIR_DIFFUSE,NRBND_VIS) * albedo_sfc(ij,NRDIR_DIFFUSE,NRBND_VIS) &
!ESC!                           + rflux_sfc_d(ij,NRDIR_DIRECT ,NRBND_NIR) * albedo_sfc(ij,NRDIR_DIRECT ,NRBND_NIR) &
!ESC!                           + rflux_sfc_d(ij,NRDIR_DIFFUSE,NRBND_NIR) * albedo_sfc(ij,NRDIR_DIFFUSE,NRBND_NIR)
!ESC!       enddo
!ESC!       !$omp end parallel do
!ESC!
!ESC!       ! LW
!ESC!
!ESC!       !$omp parallel do default(none),private(ij,emission), &
!ESC!       !$omp shared(ijdim,rflux_sfc_ld,rflux_sfc_lu,drflux_dT,rflux_sfc_d,albedo_sfc,tem_sfc)
!ESC!       do ij = 1, ijdim
!ESC!          emission = ( 1.0_RP - albedo_sfc(ij,NRDIR_DIFFUSE,NRBND_IR) ) * CONST_STB * tem_sfc(ij)**4
!ESC!
!ESC!          rflux_sfc_ld(ij) = rflux_sfc_d(ij,NRDIR_DIFFUSE,NRBND_IR )
!ESC!          rflux_sfc_lu(ij) = rflux_sfc_d(ij,NRDIR_DIFFUSE,NRBND_IR ) * albedo_sfc(ij,NRDIR_DIFFUSE,NRBND_IR ) &
!ESC!                           + emission
!ESC!
!ESC!          drflux_dT(ij) = 4.0_RP * emission / tem_sfc(ij)
!ESC!       enddo
!ESC!       !$omp end parallel do
!ESC!    endif
!ESC!
!ESC!    return
!ESC!  end subroutine rd_surface
!ESC!
!ESC!  !-----------------------------------------------------------------------------
!ESC!  ! cloud cover diagnosis
!ESC!  subroutine rd_ccover( &
!ESC!       ijdim,  &
!ESC!       CCOV,   &
!ESC!       GDCFRC, &
!ESC!       CUMFRC  )
!ESC!    use mod_adm, only: &
!ESC!       kdim => ADM_kall, &
!ESC!       kmin => ADM_kmin, &
!ESC!       kmax => ADM_kmax
!ESC!    implicit none
!ESC!
!ESC!    integer,  intent(in)  :: ijdim
!ESC!    real(RP), intent(out) :: CCOV  (ijdim)      ! cloud cover
!ESC!    real(RP), intent(in)  :: GDCFRC(ijdim,kdim) ! ratio of cloudy area
!ESC!    real(RP), intent(in)  :: CUMFRC(ijdim)      ! areal rate of cumulus
!ESC!
!ESC!    integer :: k
!ESC!    !---------------------------------------------------------------------------
!ESC!
!ESC!    ! 1 - total cloudiness
!ESC!    CCOV(:) = 1.0_RP
!ESC!    do k = kmin, kmax
!ESC!       CCOV(:) = CCOV(:) * ( 1.0_RP-GDCFRC(:,k) )
!ESC!    enddo
!ESC!
!ESC!    ! marge cumlus fraction
!ESC!    CCOV(:) = (        CUMFRC(:) ) &
!ESC!            + ( 1.0_RP-CUMFRC(:) ) * ( 1.0_RP-CCOV(:) )
!ESC!
!ESC!    return
!ESC!  end subroutine rd_ccover

end module mod_rd_driver
!-------------------------------------------------------------------------------
