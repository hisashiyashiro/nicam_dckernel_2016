!-------------------------------------------------------------------------------
!
!+  Program physics kernel driver (cloud radiation)
!
!-------------------------------------------------------------------------------
program physicskernel_radiation
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_debug
  use mod_satadjust, only: &
     SATURATION_setup
  use mod_rd_driver, only: &
     rd_init,  &
     rd_driver
  !-----------------------------------------------------------------------------
  implicit none

  real(DP), allocatable :: work1 (:,:)
  real(DP), allocatable :: work2 (:)
  real(DP), allocatable :: work3 (:,:,:)
  real(DP), allocatable :: work4 (:,:,:)
  real(DP), allocatable :: work5 (:,:,:)
  real(DP), allocatable :: work6 (:,:,:)
  real(DP), allocatable :: work7 (:,:,:)
  real(DP), allocatable :: work8 (:,:,:)
  real(DP), allocatable :: work9 (:,:)
  real(DP), allocatable :: work10(:,:)
  real(DP), allocatable :: work11(:,:,:,:)
  real(DP), allocatable :: work12(:,:)

  real(RP), allocatable :: rho           (:,:,:)
  real(RP), allocatable :: pre           (:,:,:)
  real(RP), allocatable :: tem           (:,:,:)
  real(RP), allocatable :: q_Lswp        (:,:,:,:)
  real(RP), allocatable :: q_clw         (:,:,:)
  real(RP), allocatable :: q_cli         (:,:,:)
  real(RP), allocatable :: qr            (:,:,:)
  real(RP), allocatable :: cloud_volume  (:,:,:)
  real(RP), allocatable :: cumulus_volume(:,:,:)
  real(RP), allocatable :: re_all        (:,:,:)
  real(RP), allocatable :: cfrac         (:,:,:)
  real(RP), allocatable :: q_cumclw      (:,:,:)
  real(RP), allocatable :: q_cumcli      (:,:,:)
  real(RP), allocatable :: cumfrac       (:,:,:)
  real(RP), allocatable :: qo3           (:,:,:)
  real(RP), allocatable :: tem_sfc       (:,:,:)
  real(RP), allocatable :: pre_sfc       (:,:,:)
  real(RP), allocatable :: zs            (:,:)
  real(RP), allocatable :: alat          (:,:)
  real(RP), allocatable :: alon          (:,:)
  real(RP), allocatable :: albedo_sfc    (:,:,:,:)
  real(RP), allocatable :: outqld        (:,:,:,:)
  real(RP), allocatable :: unccn         (:,:,:)
  real(RP), allocatable :: z             (:,:,:)
  real(RP), allocatable :: zh            (:,:,:)
  real(RP), allocatable :: RFLXSD        (:,:,:,:)
  real(RP), allocatable :: RFLXSU        (:,:,:,:)
  real(RP), allocatable :: RFLXLD        (:,:,:,:)
  real(RP), allocatable :: RFLXLU        (:,:,:,:)
  real(RP), allocatable :: DRFLXS        (:,:,:,:)
  real(RP), allocatable :: DRFLXL        (:,:,:,:)
  real(RP), allocatable :: RFSFCD        (:,:,:,:,:)
  real(RP), allocatable :: dfq_isccp     (:,:,:,:)
  real(RP), allocatable :: tauc_all      (:,:,:)
  real(RP), allocatable :: taucl         (:,:,:)
  real(RP), allocatable :: tauci         (:,:,:)
  real(RP), allocatable :: tauclk        (:,:,:)
  real(RP), allocatable :: taucik        (:,:,:)
  real(RP), allocatable :: rceff         (:,:,:)
  real(RP), allocatable :: rceff_solid   (:,:,:)
  real(RP), allocatable :: rceff_cld     (:,:,:)
  real(RP), allocatable :: aerdfs        (:,:,:,:)
  real(RP), allocatable :: aerdfl        (:,:,:,:)
  real(RP), allocatable :: aerdfs_trp    (:,:,:,:)
  real(RP), allocatable :: aerdfl_trp    (:,:,:,:)
  real(RP), allocatable :: tbb_11um      (:,:,:)
  real(RP), allocatable :: aot_ext_vis   (:,:,:)
  real(RP), allocatable :: aot_abs_vis   (:,:,:)

  real(RP), allocatable :: CHECK_RFLXSD     (:,:,:,:)
  real(RP), allocatable :: CHECK_RFLXSU     (:,:,:,:)
  real(RP), allocatable :: CHECK_RFLXLD     (:,:,:,:)
  real(RP), allocatable :: CHECK_RFLXLU     (:,:,:,:)
  real(RP), allocatable :: CHECK_DRFLXS     (:,:,:,:)
  real(RP), allocatable :: CHECK_DRFLXL     (:,:,:,:)
  real(RP), allocatable :: CHECK_RFSFCD     (:,:,:,:,:)
  real(RP), allocatable :: CHECK_dfq_isccp  (:,:,:,:)
  real(RP), allocatable :: CHECK_tauc_all   (:,:,:)
  real(RP), allocatable :: CHECK_taucl      (:,:,:)
  real(RP), allocatable :: CHECK_tauci      (:,:,:)
  real(RP), allocatable :: CHECK_tauclk     (:,:,:)
  real(RP), allocatable :: CHECK_taucik     (:,:,:)
  real(RP), allocatable :: CHECK_rceff      (:,:,:)
  real(RP), allocatable :: CHECK_rceff_solid(:,:,:)
  real(RP), allocatable :: CHECK_rceff_cld  (:,:,:)
  real(RP), allocatable :: CHECK_aerdfs     (:,:,:,:)
  real(RP), allocatable :: CHECK_aerdfl     (:,:,:,:)
  real(RP), allocatable :: CHECK_aerdfs_trp (:,:,:,:)
  real(RP), allocatable :: CHECK_aerdfl_trp (:,:,:,:)
  real(RP), allocatable :: CHECK_tbb_11um   (:,:,:)
  real(RP), allocatable :: CHECK_aot_ext_vis(:,:,:)
  real(RP), allocatable :: CHECK_aot_abs_vis(:,:,:)

  logical  :: rd_update_flag

  integer :: l, iteration
  !=============================================================================

  write(*,*) "[KERNEL] physicskernel_radiation"
  write(*,*) "*** Start  initialize"

  allocate( work1 (ADM_gall_in,ADM_kall)           )
  allocate( work2 (ADM_gall_in)                    )
  allocate( work3 (ADM_gall_in,ADM_kall,HYDRO_MAX) )
  allocate( work4 (ADM_gall_in,ADM_kall,2)         )
  allocate( work5 (ADM_gall_in,NRDIR,NRBND)        )
  allocate( work6 (ADM_gall_in,ADM_kall,KAPCL)     )
  allocate( work7 (ADM_gall_in,ADM_kall,NCRF)      )
  allocate( work8 (ADM_gall_in,NTAU,NPRES)         )
  allocate( work9 (ADM_gall_in,HYDRO_MAX)          )
  allocate( work10(ADM_gall_in,NCRF)               )
  allocate( work11(ADM_gall_in,NRDIR,NRBND,NCRF)   )
  allocate( work12(ADM_gall_in,KAPCL)              )

  allocate( rho              (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( pre              (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( tem              (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( q_Lswp           (ADM_gall_in,ADM_kall,TRC_VMAX,ADM_lall) )
  allocate( q_clw            (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( q_cli            (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( qr               (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( cloud_volume     (ADM_gall_in,ADM_kall,HYDRO_MAX)         )
  allocate( cumulus_volume   (ADM_gall_in,ADM_kall,2)                 )
  allocate( re_all           (ADM_gall_in,ADM_kall,HYDRO_MAX)         )
  allocate( cfrac            (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( q_cumclw         (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( q_cumcli         (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( cumfrac          (ADM_gall_in,ADM_KNONE,ADM_lall)         )
  allocate( qo3              (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( tem_sfc          (ADM_gall_in,ADM_KNONE,ADM_lall)         )
  allocate( pre_sfc          (ADM_gall_in,ADM_KNONE,ADM_lall)         )
  allocate( zs               (ADM_gall_in,ADM_lall)                   )
  allocate( alat             (ADM_gall_in,ADM_lall)                   )
  allocate( alon             (ADM_gall_in,ADM_lall)                   )
  allocate( albedo_sfc       (ADM_gall_in,NRDIR,NRBND,ADM_lall)       )
  allocate( outqld           (ADM_gall_in,ADM_kall,KAPCL,ADM_lall)    )
  allocate( unccn            (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( z                (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( zh               (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( RFLXSD           (ADM_gall_in,ADM_kall,NCRF,ADM_lall)     )
  allocate( RFLXSU           (ADM_gall_in,ADM_kall,NCRF,ADM_lall)     )
  allocate( RFLXLD           (ADM_gall_in,ADM_kall,NCRF,ADM_lall)     )
  allocate( RFLXLU           (ADM_gall_in,ADM_kall,NCRF,ADM_lall)     )
  allocate( DRFLXS           (ADM_gall_in,ADM_kall,NCRF,ADM_lall)     )
  allocate( DRFLXL           (ADM_gall_in,ADM_kall,NCRF,ADM_lall)     )
  allocate( RFSFCD           (ADM_gall_in,NRDIR,NRBND,NCRF,ADM_lall)  )
  allocate( dfq_isccp        (ADM_gall_in,NTAU,NPRES,ADM_lall)        )
  allocate( tauc_all         (ADM_gall_in,HYDRO_MAX,ADM_lall)         )
  allocate( taucl            (ADM_gall_in,ADM_KNONE,ADM_lall)         )
  allocate( tauci            (ADM_gall_in,ADM_KNONE,ADM_lall)         )
  allocate( tauclk           (ADM_gall_in,ADM_kall ,ADM_lall)         )
  allocate( taucik           (ADM_gall_in,ADM_kall ,ADM_lall)         )
  allocate( rceff            (ADM_gall_in,ADM_kall ,ADM_lall)         )
  allocate( rceff_solid      (ADM_gall_in,ADM_kall ,ADM_lall)         )
  allocate( rceff_cld        (ADM_gall_in,ADM_kall ,ADM_lall)         )
  allocate( aerdfs           (ADM_gall_in,ADM_kall ,ADM_lall,NCRF)    )
  allocate( aerdfl           (ADM_gall_in,ADM_kall ,ADM_lall,NCRF)    )
  allocate( aerdfs_trp       (ADM_gall_in,ADM_KNONE,ADM_lall,NCRF)    )
  allocate( aerdfl_trp       (ADM_gall_in,ADM_KNONE,ADM_lall,NCRF)    )
  allocate( tbb_11um         (ADM_gall_in,ADM_KNONE,ADM_lall)         )
  allocate( aot_ext_vis      (ADM_gall_in,KAPCL    ,ADM_lall)         )
  allocate( aot_abs_vis      (ADM_gall_in,KAPCL    ,ADM_lall)         )

  allocate( CHECK_RFLXSD     (ADM_gall_in,ADM_kall,NCRF,ADM_lall)     )
  allocate( CHECK_RFLXSU     (ADM_gall_in,ADM_kall,NCRF,ADM_lall)     )
  allocate( CHECK_RFLXLD     (ADM_gall_in,ADM_kall,NCRF,ADM_lall)     )
  allocate( CHECK_RFLXLU     (ADM_gall_in,ADM_kall,NCRF,ADM_lall)     )
  allocate( CHECK_DRFLXS     (ADM_gall_in,ADM_kall,NCRF,ADM_lall)     )
  allocate( CHECK_DRFLXL     (ADM_gall_in,ADM_kall,NCRF,ADM_lall)     )
  allocate( CHECK_RFSFCD     (ADM_gall_in,NRDIR,NRBND,NCRF,ADM_lall)  )
  allocate( CHECK_dfq_isccp  (ADM_gall_in,NTAU,NPRES,ADM_lall)        )
  allocate( CHECK_tauc_all   (ADM_gall_in,HYDRO_MAX,ADM_lall)         )
  allocate( CHECK_taucl      (ADM_gall_in,ADM_KNONE,ADM_lall)         )
  allocate( CHECK_tauci      (ADM_gall_in,ADM_KNONE,ADM_lall)         )
  allocate( CHECK_tauclk     (ADM_gall_in,ADM_kall ,ADM_lall)         )
  allocate( CHECK_taucik     (ADM_gall_in,ADM_kall ,ADM_lall)         )
  allocate( CHECK_rceff      (ADM_gall_in,ADM_kall ,ADM_lall)         )
  allocate( CHECK_rceff_solid(ADM_gall_in,ADM_kall ,ADM_lall)         )
  allocate( CHECK_rceff_cld  (ADM_gall_in,ADM_kall ,ADM_lall)         )
  allocate( CHECK_aerdfs     (ADM_gall_in,ADM_kall ,ADM_lall,NCRF)    )
  allocate( CHECK_aerdfl     (ADM_gall_in,ADM_kall ,ADM_lall,NCRF)    )
  allocate( CHECK_aerdfs_trp (ADM_gall_in,ADM_KNONE,ADM_lall,NCRF)    )
  allocate( CHECK_aerdfl_trp (ADM_gall_in,ADM_KNONE,ADM_lall,NCRF)    )
  allocate( CHECK_tbb_11um   (ADM_gall_in,ADM_KNONE,ADM_lall)         )
  allocate( CHECK_aot_ext_vis(ADM_gall_in,KAPCL    ,ADM_lall)         )
  allocate( CHECK_aot_abs_vis(ADM_gall_in,KAPCL    ,ADM_lall)         )

  l = SET_l

  !###############################################################################

   EX_rgnid = SET_rgnid

   call dumpio_syscheck
   call dumpio_mk_fname(EX_fname,'snapshot.radiation','pe',EX_rgnid-1,6)
   call dumpio_fopen(EX_fid,EX_fname,IO_FREAD)

   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall          , work1 (:,:)     ); rho           (1:ADM_gall_in,:,l)           = real(work1 (1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall          , work1 (:,:)     ); pre           (1:ADM_gall_in,:,l)           = real(work1 (1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall          , work1 (:,:)     ); tem           (1:ADM_gall_in,:,l)           = real(work1 (1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall          , work1 (:,:)     ); q_Lswp        (1:ADM_gall_in,:,I_QV,l)      = real(work1 (1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall          , work1 (:,:)     ); q_clw         (1:ADM_gall_in,:,l)           = real(work1 (1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall          , work1 (:,:)     ); q_cli         (1:ADM_gall_in,:,l)           = real(work1 (1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall          , work1 (:,:)     ); qr            (1:ADM_gall_in,:,l)           = real(work1 (1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall*HYDRO_MAX, work3 (:,:,:)   ); cloud_volume  (1:ADM_gall_in,:,:)           = real(work3 (1:ADM_gall_in,:,:)  ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall*2        , work4 (:,:,:)   ); cumulus_volume(1:ADM_gall_in,:,:)           = real(work4 (1:ADM_gall_in,:,:)  ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall*HYDRO_MAX, work3 (:,:,:)   ); re_all        (1:ADM_gall_in,:,:)           = real(work3 (1:ADM_gall_in,:,:)  ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall          , work1 (:,:)     ); cfrac         (1:ADM_gall_in,:,l)           = real(work1 (1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall          , work1 (:,:)     ); q_cumclw      (1:ADM_gall_in,:,l)           = real(work1 (1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall          , work1 (:,:)     ); q_cumcli      (1:ADM_gall_in,:,l)           = real(work1 (1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in                   , work2 (:)       ); cumfrac       (1:ADM_gall_in,ADM_KNONE,l)   = real(work2 (1:ADM_gall_in)      ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall          , work1 (:,:)     ); qo3           (1:ADM_gall_in,:,l)           = real(work1 (1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in                   , work2 (:)       ); tem_sfc       (1:ADM_gall_in,ADM_KNONE,l)   = real(work2 (1:ADM_gall_in)      ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in                   , work2 (:)       ); pre_sfc       (1:ADM_gall_in,ADM_KNONE,l)   = real(work2 (1:ADM_gall_in)      ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in                   , work2 (:)       ); zs            (1:ADM_gall_in,l)             = real(work2 (1:ADM_gall_in)      ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in                   , work2 (:)       ); alat          (1:ADM_gall_in,l)             = real(work2 (1:ADM_gall_in)      ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in                   , work2 (:)       ); alon          (1:ADM_gall_in,l)             = real(work2 (1:ADM_gall_in)      ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*NRDIR*NRBND       , work5 (:,:,:)   ); albedo_sfc    (1:ADM_gall_in,:,:,l)         = real(work5 (1:ADM_gall_in,:,:)  ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall*KAPCL    , work6 (:,:,:)   ); outqld        (1:ADM_gall_in,:,:,l)         = real(work6 (1:ADM_gall_in,:,:)  ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall          , work1 (:,:)     ); unccn         (1:ADM_gall_in,:,l)           = real(work1 (1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall          , work1 (:,:)     ); z             (1:ADM_gall_in,:,l)           = real(work1 (1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall          , work1 (:,:)     ); zh            (1:ADM_gall_in,:,l)           = real(work1 (1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall*NCRF     , work7 (:,:,:)   ); RFLXSD        (1:ADM_gall_in,:,:,l)         = real(work7 (1:ADM_gall_in,:,:)  ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall*NCRF     , work7 (:,:,:)   ); RFLXSU        (1:ADM_gall_in,:,:,l)         = real(work7 (1:ADM_gall_in,:,:)  ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall*NCRF     , work7 (:,:,:)   ); RFLXLD        (1:ADM_gall_in,:,:,l)         = real(work7 (1:ADM_gall_in,:,:)  ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall*NCRF     , work7 (:,:,:)   ); RFLXLU        (1:ADM_gall_in,:,:,l)         = real(work7 (1:ADM_gall_in,:,:)  ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall*NCRF     , work7 (:,:,:)   ); DRFLXS        (1:ADM_gall_in,:,:,l)         = real(work7 (1:ADM_gall_in,:,:)  ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall*NCRF     , work7 (:,:,:)   ); DRFLXL        (1:ADM_gall_in,:,:,l)         = real(work7 (1:ADM_gall_in,:,:)  ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*NRDIR*NRBND*NCRF  , work11(:,:,:,:) ); RFSFCD        (1:ADM_gall_in,:,:,:,l)       = real(work11(1:ADM_gall_in,:,:,:),kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*NTAU*NPRES        , work8 (:,:,:)   ); dfq_isccp     (1:ADM_gall_in,:,:,l)         = real(work8 (1:ADM_gall_in,:,:)  ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*HYDRO_MAX         , work9 (:,:)     ); tauc_all      (1:ADM_gall_in,:,l)           = real(work9 (1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in                   , work2 (:)       ); taucl         (1:ADM_gall_in,ADM_KNONE,l)   = real(work2 (1:ADM_gall_in)      ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in                   , work2 (:)       ); tauci         (1:ADM_gall_in,ADM_KNONE,l)   = real(work2 (1:ADM_gall_in)      ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall          , work1 (:,:)     ); tauclk        (1:ADM_gall_in,:,l)           = real(work1 (1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall          , work1 (:,:)     ); taucik        (1:ADM_gall_in,:,l)           = real(work1 (1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall          , work1 (:,:)     ); rceff         (1:ADM_gall_in,:,l)           = real(work1 (1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall          , work1 (:,:)     ); rceff_solid   (1:ADM_gall_in,:,l)           = real(work1 (1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall          , work1 (:,:)     ); rceff_cld     (1:ADM_gall_in,:,l)           = real(work1 (1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall*NCRF     , work7 (:,:,:)   ); aerdfs        (1:ADM_gall_in,:,l,:)         = real(work7 (1:ADM_gall_in,:,:)  ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall*NCRF     , work7 (:,:,:)   ); aerdfl        (1:ADM_gall_in,:,l,:)         = real(work7 (1:ADM_gall_in,:,:)  ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*NCRF              , work10(:,:)     ); aerdfs_trp    (1:ADM_gall_in,ADM_KNONE,l,:) = real(work10(1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*NCRF              , work10(:,:)     ); aerdfl_trp    (1:ADM_gall_in,ADM_KNONE,l,:) = real(work10(1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in                   , work2 (:)       ); tbb_11um      (1:ADM_gall_in,ADM_KNONE,l)   = real(work2 (1:ADM_gall_in)      ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*KAPCL             , work12(:,:)     ); aot_ext_vis   (1:ADM_gall_in,:,l)           = real(work12(1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*KAPCL             , work12(:,:)     ); aot_abs_vis   (1:ADM_gall_in,:,l)           = real(work12(1:ADM_gall_in,:)    ,kind=RP)

   call dumpio_fclose(EX_fid)

  !###############################################################################

  CONST_PI  = 4.E0_RP * atan( 1.0_RP )
  PI        = CONST_PI
  CONST_EPS = epsilon(0.0_RP)

  !--- vertical grid setup
  call GRD_setup

  call SATURATION_setup

  !---  radiation initialization
  call rd_init( RD_TYPE )

  !###############################################################################

  write(*,*) "*** Finish initialize"
  write(*,*) "*** Start  kernel"
  call PROF_rapstart('Radiation_kernel')

  do iteration = 1, SET_iteration
     write(*,*) '### Before ###'
     call PROF_valcheck( 'before', 'RFLXSD', RFLXSD(:,:,:,l) )
     call PROF_valcheck( 'before', 'RFLXSU', RFLXSU(:,:,:,l) )
     call PROF_valcheck( 'before', 'RFLXLD', RFLXLD(:,:,:,l) )
     call PROF_valcheck( 'before', 'RFLXLU', RFLXLU(:,:,:,l) )
     call PROF_valcheck( 'before', 'DRFLXS', DRFLXS(:,:,:,l) )
     call PROF_valcheck( 'before', 'DRFLXL', DRFLXL(:,:,:,l) )
     call PROF_valcheck( 'before', 'RFSFCD', RFSFCD(:,:,:,:,l) )

       call rd_driver( ADM_gall_in,                     & ! [IN]
                       rho           (:,:,l),           & ! [IN]
                       pre           (:,:,l),           & ! [IN]
                       tem           (:,:,l),           & ! [IN]
                       q_Lswp        (:,:,I_QV,l),      & ! [IN]
                       q_clw         (:,:,l),           & ! [IN]
                       q_cli         (:,:,l),           & ! [IN]
                       qr            (:,:,l),           & ! [IN]
                       cloud_volume  (:,:,:),           & ! [IN]
                       cumulus_volume(:,:,:),           & ! [IN]
                       re_all        (:,:,:),           & ! [IN]
                       cfrac         (:,:,l),           & ! [IN]
                       q_cumclw      (:,:,l),           & ! [IN]
                       q_cumcli      (:,:,l),           & ! [IN]
                       cumfrac       (:,ADM_KNONE,l),   & ! [IN]
                       qo3           (:,:,l),           & ! [IN]
                       tem_sfc       (:,ADM_KNONE,l),   & ! [IN]
                       pre_sfc       (:,ADM_KNONE,l),   & ! [IN]
                       zs            (:,l),             & ! [IN]
                       alat          (:,l),             & ! [IN]
                       alon          (:,l),             & ! [IN]
                       albedo_sfc    (:,:,:,l),         & ! [IN]
                       outqld        (:,:,:,l),         & ! [IN]
                       unccn         (:,:,l),           & ! [IN]
                       z             (:,:,l),           & ! [IN]
                       zh            (:,:,l),           & ! [IN]
                       RFLXSD        (:,:,:,l),         & ! [INOUT]
                       RFLXSU        (:,:,:,l),         & ! [INOUT]
                       RFLXLD        (:,:,:,l),         & ! [INOUT]
                       RFLXLU        (:,:,:,l),         & ! [INOUT]
                       DRFLXS        (:,:,:,l),         & ! [INOUT]
                       DRFLXL        (:,:,:,l),         & ! [INOUT]
                       RFSFCD        (:,:,:,:,l),       & ! [INOUT]
                       dfq_isccp     (:,:,:,l),         & ! [INOUT]
                       tauc_all      (:,:,l),           & ! [INOUT]
                       taucl         (:,ADM_KNONE,l),   & ! [INOUT]
                       tauci         (:,ADM_KNONE,l),   & ! [INOUT]
                       tauclk        (:,:,l),           & ! [INOUT]
                       taucik        (:,:,l),           & ! [INOUT]
                       rceff         (:,:,l),           & ! [INOUT]
                       rceff_solid   (:,:,l),           & ! [INOUT]
                       rceff_cld     (:,:,l),           & ! [INOUT]
                       aerdfs        (:,:,l,:),         & ! [INOUT]
                       aerdfl        (:,:,l,:),         & ! [INOUT]
                       aerdfs_trp    (:,ADM_KNONE,l,:), & ! [INOUT]
                       aerdfl_trp    (:,ADM_KNONE,l,:), & ! [INOUT]
                       tbb_11um      (:,ADM_KNONE,l),   & ! [INOUT]
                       aot_ext_vis   (:,:,l),           & ! [INOUT] [Add] 2016/05/18 T.Seiki
                       aot_abs_vis   (:,:,l),           & ! [INOUT] [Add] 2016/05/18 T.Seiki
                       rd_update_flag                   ) ! [OUT]

     write(*,*) '### After ###'
     call PROF_valcheck( 'after', 'RFLXSD', RFLXSD(:,:,:,l) )
     call PROF_valcheck( 'after', 'RFLXSU', RFLXSU(:,:,:,l) )
     call PROF_valcheck( 'after', 'RFLXLD', RFLXLD(:,:,:,l) )
     call PROF_valcheck( 'after', 'RFLXLU', RFLXLU(:,:,:,l) )
     call PROF_valcheck( 'after', 'DRFLXS', DRFLXS(:,:,:,l) )
     call PROF_valcheck( 'after', 'DRFLXL', DRFLXL(:,:,:,l) )
     call PROF_valcheck( 'after', 'RFSFCD', RFSFCD(:,:,:,:,l) )
  enddo

  call PROF_rapend('Radiation_kernel')
  write(*,*) "*** Finish kernel"

  !###############################################################################
  if ( SET_check .AND. SET_iteration == 1 ) then
   EX_rgnid = SET_rgnid

   call dumpio_mk_fname(EX_fname,'check.radiation','pe',EX_rgnid-1,6)
   call dumpio_fopen(EX_fid,EX_fname,IO_FREAD)

   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall*NCRF   , work7 (:,:,:)   ); CHECK_RFLXSD     (1:ADM_gall_in,:,:,l)         = real(work7 (1:ADM_gall_in,:,:)  ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall*NCRF   , work7 (:,:,:)   ); CHECK_RFLXSU     (1:ADM_gall_in,:,:,l)         = real(work7 (1:ADM_gall_in,:,:)  ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall*NCRF   , work7 (:,:,:)   ); CHECK_RFLXLD     (1:ADM_gall_in,:,:,l)         = real(work7 (1:ADM_gall_in,:,:)  ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall*NCRF   , work7 (:,:,:)   ); CHECK_RFLXLU     (1:ADM_gall_in,:,:,l)         = real(work7 (1:ADM_gall_in,:,:)  ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall*NCRF   , work7 (:,:,:)   ); CHECK_DRFLXS     (1:ADM_gall_in,:,:,l)         = real(work7 (1:ADM_gall_in,:,:)  ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall*NCRF   , work7 (:,:,:)   ); CHECK_DRFLXL     (1:ADM_gall_in,:,:,l)         = real(work7 (1:ADM_gall_in,:,:)  ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*NRDIR*NRBND*NCRF, work11(:,:,:,:) ); CHECK_RFSFCD     (1:ADM_gall_in,:,:,:,l)       = real(work11(1:ADM_gall_in,:,:,:),kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*NTAU*NPRES      , work8 (:,:,:)   ); CHECK_dfq_isccp  (1:ADM_gall_in,:,:,l)         = real(work8 (1:ADM_gall_in,:,:)  ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*HYDRO_MAX       , work9 (:,:)     ); CHECK_tauc_all   (1:ADM_gall_in,:,l)           = real(work9 (1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in                 , work2 (:)       ); CHECK_taucl      (1:ADM_gall_in,ADM_KNONE,l)   = real(work2 (1:ADM_gall_in)      ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in                 , work2 (:)       ); CHECK_tauci      (1:ADM_gall_in,ADM_KNONE,l)   = real(work2 (1:ADM_gall_in)      ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall        , work1 (:,:)     ); CHECK_tauclk     (1:ADM_gall_in,:,l)           = real(work1 (1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall        , work1 (:,:)     ); CHECK_taucik     (1:ADM_gall_in,:,l)           = real(work1 (1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall        , work1 (:,:)     ); CHECK_rceff      (1:ADM_gall_in,:,l)           = real(work1 (1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall        , work1 (:,:)     ); CHECK_rceff_solid(1:ADM_gall_in,:,l)           = real(work1 (1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall        , work1 (:,:)     ); CHECK_rceff_cld  (1:ADM_gall_in,:,l)           = real(work1 (1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall*NCRF   , work7 (:,:,:)   ); CHECK_aerdfs     (1:ADM_gall_in,:,l,:)         = real(work7 (1:ADM_gall_in,:,:)  ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall*NCRF   , work7 (:,:,:)   ); CHECK_aerdfl     (1:ADM_gall_in,:,l,:)         = real(work7 (1:ADM_gall_in,:,:)  ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*NCRF            , work10(:,:)     ); CHECK_aerdfs_trp (1:ADM_gall_in,ADM_KNONE,l,:) = real(work10(1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*NCRF            , work10(:,:)     ); CHECK_aerdfl_trp (1:ADM_gall_in,ADM_KNONE,l,:) = real(work10(1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in                 , work2 (:)       ); CHECK_tbb_11um   (1:ADM_gall_in,ADM_KNONE,l)   = real(work2 (1:ADM_gall_in)      ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*KAPCL           , work12(:,:)     ); CHECK_aot_ext_vis(1:ADM_gall_in,:,l)           = real(work12(1:ADM_gall_in,:)    ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*KAPCL           , work12(:,:)     ); CHECK_aot_abs_vis(1:ADM_gall_in,:,l)           = real(work12(1:ADM_gall_in,:)    ,kind=RP)

   call dumpio_fclose(EX_fid)

     CHECK_RFLXSD     (:,:,:,:)   = RFLXSD     (:,:,:,:)   - CHECK_RFLXSD     (:,:,:,:)
     CHECK_RFLXSU     (:,:,:,:)   = RFLXSU     (:,:,:,:)   - CHECK_RFLXSU     (:,:,:,:)
     CHECK_RFLXLD     (:,:,:,:)   = RFLXLD     (:,:,:,:)   - CHECK_RFLXLD     (:,:,:,:)
     CHECK_RFLXLU     (:,:,:,:)   = RFLXLU     (:,:,:,:)   - CHECK_RFLXLU     (:,:,:,:)
     CHECK_DRFLXS     (:,:,:,:)   = DRFLXS     (:,:,:,:)   - CHECK_DRFLXS     (:,:,:,:)
     CHECK_DRFLXL     (:,:,:,:)   = DRFLXL     (:,:,:,:)   - CHECK_DRFLXL     (:,:,:,:)
     CHECK_RFSFCD     (:,:,:,:,:) = RFSFCD     (:,:,:,:,:) - CHECK_RFSFCD     (:,:,:,:,:)
     CHECK_dfq_isccp  (:,:,:,:)   = dfq_isccp  (:,:,:,:)   - CHECK_dfq_isccp  (:,:,:,:)
     CHECK_tauc_all   (:,:,:)     = tauc_all   (:,:,:)     - CHECK_tauc_all   (:,:,:)
     CHECK_taucl      (:,:,:)     = taucl      (:,:,:)     - CHECK_taucl      (:,:,:)
     CHECK_tauci      (:,:,:)     = tauci      (:,:,:)     - CHECK_tauci      (:,:,:)
     CHECK_tauclk     (:,:,:)     = tauclk     (:,:,:)     - CHECK_tauclk     (:,:,:)
     CHECK_taucik     (:,:,:)     = taucik     (:,:,:)     - CHECK_taucik     (:,:,:)
     CHECK_rceff      (:,:,:)     = rceff      (:,:,:)     - CHECK_rceff      (:,:,:)
     CHECK_rceff_solid(:,:,:)     = rceff_solid(:,:,:)     - CHECK_rceff_solid(:,:,:)
     CHECK_rceff_cld  (:,:,:)     = rceff_cld  (:,:,:)     - CHECK_rceff_cld  (:,:,:)
     CHECK_aerdfs     (:,:,:,:)   = aerdfs     (:,:,:,:)   - CHECK_aerdfs     (:,:,:,:)
     CHECK_aerdfl     (:,:,:,:)   = aerdfl     (:,:,:,:)   - CHECK_aerdfl     (:,:,:,:)
     CHECK_aerdfs_trp (:,:,:,:)   = aerdfs_trp (:,:,:,:)   - CHECK_aerdfs_trp (:,:,:,:)
     CHECK_aerdfl_trp (:,:,:,:)   = aerdfl_trp (:,:,:,:)   - CHECK_aerdfl_trp (:,:,:,:)
     CHECK_tbb_11um   (:,:,:)     = tbb_11um   (:,:,:)     - CHECK_tbb_11um   (:,:,:)
     CHECK_aot_ext_vis(:,:,:)     = aot_ext_vis(:,:,:)     - CHECK_aot_ext_vis(:,:,:)
     CHECK_aot_abs_vis(:,:,:)     = aot_abs_vis(:,:,:)     - CHECK_aot_abs_vis(:,:,:)

     write(*,*) '### Check ###'
     call PROF_valcheck( 'check', 'check_RFLXSD     ', CHECK_RFLXSD     (:,:,:,:)   )
     call PROF_valcheck( 'check', 'check_RFLXSU     ', CHECK_RFLXSU     (:,:,:,:)   )
     call PROF_valcheck( 'check', 'check_RFLXLD     ', CHECK_RFLXLD     (:,:,:,:)   )
     call PROF_valcheck( 'check', 'check_RFLXLU     ', CHECK_RFLXLU     (:,:,:,:)   )
     call PROF_valcheck( 'check', 'check_DRFLXS     ', CHECK_DRFLXS     (:,:,:,:)   )
     call PROF_valcheck( 'check', 'check_DRFLXL     ', CHECK_DRFLXL     (:,:,:,:)   )
     call PROF_valcheck( 'check', 'check_RFSFCD     ', CHECK_RFSFCD     (:,:,:,:,:) )
     call PROF_valcheck( 'check', 'check_dfq_isccp  ', CHECK_dfq_isccp  (:,:,:,:)   )
     call PROF_valcheck( 'check', 'check_tauc_all   ', CHECK_tauc_all   (:,:,:)     )
     call PROF_valcheck( 'check', 'check_taucl      ', CHECK_taucl      (:,:,:)     )
     call PROF_valcheck( 'check', 'check_tauci      ', CHECK_tauci      (:,:,:)     )
     call PROF_valcheck( 'check', 'check_tauclk     ', CHECK_tauclk     (:,:,:)     )
     call PROF_valcheck( 'check', 'check_taucik     ', CHECK_taucik     (:,:,:)     )
     call PROF_valcheck( 'check', 'check_rceff      ', CHECK_rceff      (:,:,:)     )
     call PROF_valcheck( 'check', 'check_rceff_solid', CHECK_rceff_solid(:,:,:)     )
     call PROF_valcheck( 'check', 'check_rceff_cld  ', CHECK_rceff_cld  (:,:,:)     )
!     call PROF_valcheck( 'check', 'check_aerdfs     ', CHECK_aerdfs     (:,:,:,:)   ) not for use
!     call PROF_valcheck( 'check', 'check_aerdfl     ', CHECK_aerdfl     (:,:,:,:)   ) not for use
!     call PROF_valcheck( 'check', 'check_aerdfs_trp ', CHECK_aerdfs_trp (:,:,:,:)   ) not for use
!     call PROF_valcheck( 'check', 'check_aerdfl_trp ', CHECK_aerdfl_trp (:,:,:,:)   ) not for use
     call PROF_valcheck( 'check', 'check_tbb_11um   ', CHECK_tbb_11um   (:,:,:)     )
!     call PROF_valcheck( 'check', 'check_aot_ext_vis', CHECK_aot_ext_vis(:,:,:)     ) not for use
!     call PROF_valcheck( 'check', 'check_aot_abs_vis', CHECK_aot_abs_vis(:,:,:)     ) not for use

  endif
  !###############################################################################

  call PROF_rapreport

  stop
end program physicskernel_radiation
!-------------------------------------------------------------------------------
