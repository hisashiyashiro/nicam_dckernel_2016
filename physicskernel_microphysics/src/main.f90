!-------------------------------------------------------------------------------
!
!+  Program physics kernel driver (cloud microphysics)
!
!-------------------------------------------------------------------------------
program physicskernel_microphysics
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_debug
  use mod_satadjust, only: &
     SATURATION_setup
  use mod_mp_driver, only: &
     mp_init,  &
     mp_driver
  !-----------------------------------------------------------------------------
  implicit none

  real(DP), allocatable :: work1(:,:,:)
  real(DP), allocatable :: work2(:,:,:,:)
  real(DP), allocatable :: work3(:,:)
  real(DP), allocatable :: work4(:,:,:,:)

  real(RP), allocatable :: rhog             (:,:,:)
  real(RP), allocatable :: rhogvx           (:,:,:)
  real(RP), allocatable :: rhogvy           (:,:,:)
  real(RP), allocatable :: rhogvz           (:,:,:)
  real(RP), allocatable :: rhogw            (:,:,:)
  real(RP), allocatable :: rhoge            (:,:,:)
  real(RP), allocatable :: rhogq_Lswp       (:,:,:,:)
  real(RP), allocatable :: vx               (:,:,:)
  real(RP), allocatable :: vy               (:,:,:)
  real(RP), allocatable :: vz               (:,:,:)
  real(RP), allocatable :: w                (:,:,:)
  real(RP), allocatable :: unccn            (:,:,:)
  real(RP), allocatable :: rho              (:,:,:)
  real(RP), allocatable :: pre              (:,:,:)
  real(RP), allocatable :: tem              (:,:,:)
  real(RP), allocatable :: q_Lswp           (:,:,:,:)
  real(RP), allocatable :: precip_mp        (:,:,:,:)
  real(RP), allocatable :: precip1_mp       (:,:,:,:)
  real(RP), allocatable :: precip2_mp       (:,:,:,:)
  real(RP), allocatable :: rhoein_precip_mp (:,:,:)
  real(RP), allocatable :: lh_precip_mp     (:,:,:)
  real(RP), allocatable :: rhophi_precip_mp (:,:,:)
  real(RP), allocatable :: rhokin_precip_mp (:,:,:)
  real(RP), allocatable :: frhoge_af        (:,:,:)
  real(RP), allocatable :: frhogqv_af       (:,:,:)
  real(RP), allocatable :: frhoge_rad       (:,:,:)
  real(RP), allocatable :: qke              (:,:,:)
  real(RP), allocatable :: gsgam2           (:,:,:)
  real(RP), allocatable :: gsgam2h          (:,:,:)
  real(RP), allocatable :: gam2             (:,:,:)
  real(RP), allocatable :: gam2h            (:,:,:)
  real(RP), allocatable :: ix               (:,:)
  real(RP), allocatable :: iy               (:,:)
  real(RP), allocatable :: iz               (:,:)
  real(RP), allocatable :: jx               (:,:)
  real(RP), allocatable :: jy               (:,:)
  real(RP), allocatable :: jz               (:,:)
  real(RP), allocatable :: z                (:,:,:)
  real(RP), allocatable :: zh               (:,:,:)
  real(RP), allocatable :: GPREC            (:,:,:)
  real(RP), allocatable :: CBMFX            (:,:,:)
  real(RP), allocatable :: qd               (:,:,:)
  real(RP), allocatable :: rceff            (:,:,:)
  real(RP), allocatable :: rceff_solid      (:,:,:)
  real(RP), allocatable :: rceff_cld        (:,:,:)
  real(RP), allocatable :: rctop            (:,:,:)
  real(RP), allocatable :: rwtop            (:,:,:)
  real(RP), allocatable :: tctop            (:,:,:)
  real(RP), allocatable :: GDCLW            (:,:,:)
  real(RP), allocatable :: GDCFRC           (:,:,:)

  real(RP), allocatable :: CHECK_rhog             (:,:,:)
  real(RP), allocatable :: CHECK_rhogvx           (:,:,:)
  real(RP), allocatable :: CHECK_rhogvy           (:,:,:)
  real(RP), allocatable :: CHECK_rhogvz           (:,:,:)
  real(RP), allocatable :: CHECK_rhogw            (:,:,:)
  real(RP), allocatable :: CHECK_rhoge            (:,:,:)
  real(RP), allocatable :: CHECK_rhogq_Lswp       (:,:,:,:)
  real(RP), allocatable :: CHECK_rho              (:,:,:)
  real(RP), allocatable :: CHECK_pre              (:,:,:)
  real(RP), allocatable :: CHECK_tem              (:,:,:)
  real(RP), allocatable :: CHECK_q_Lswp           (:,:,:,:)
  real(RP), allocatable :: CHECK_qd               (:,:,:)
  real(RP), allocatable :: CHECK_precip_mp        (:,:,:,:)
  real(RP), allocatable :: CHECK_precip1_mp       (:,:,:,:)
  real(RP), allocatable :: CHECK_precip2_mp       (:,:,:,:)
  real(RP), allocatable :: CHECK_rhoein_precip_mp (:,:,:)
  real(RP), allocatable :: CHECK_lh_precip_mp     (:,:,:)
  real(RP), allocatable :: CHECK_rhophi_precip_mp (:,:,:)
  real(RP), allocatable :: CHECK_rhokin_precip_mp (:,:,:)
  real(RP), allocatable :: CHECK_rceff            (:,:,:)
  real(RP), allocatable :: CHECK_rceff_solid      (:,:,:)
  real(RP), allocatable :: CHECK_rceff_cld        (:,:,:)
  real(RP), allocatable :: CHECK_rctop            (:,:,:)
  real(RP), allocatable :: CHECK_rwtop            (:,:,:)
  real(RP), allocatable :: CHECK_tctop            (:,:,:)
  real(RP), allocatable :: CHECK_GDCLW            (:,:,:)
  real(RP), allocatable :: CHECK_GDCFRC           (:,:,:)
  real(RP), allocatable :: CHECK_GPREC            (:,:,:)

  integer :: l, nq, iteration
  !=============================================================================

  write(*,*) "[KERNEL] physicskernel_microphysics"
  write(*,*) "*** Start  initialize"

  allocate( work1(ADM_gall_in_orig,ADM_kall,ADM_lall)          )
  allocate( work2(ADM_gall_in_orig,ADM_kall,TRC_VMAX,ADM_lall) )
  allocate( work3(ADM_gall_in_orig,ADM_lall)                   )
  allocate( work4(ADM_gall_in_orig,ADM_KNONE,ADM_lall,2)       )

  allocate( rhog                   (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( rhogvx                 (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( rhogvy                 (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( rhogvz                 (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( rhogw                  (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( rhoge                  (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( rhogq_Lswp             (ADM_gall_in,ADM_kall,TRC_VMAX,ADM_lall) )
  allocate( vx                     (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( vy                     (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( vz                     (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( w                      (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( unccn                  (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( rho                    (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( pre                    (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( tem                    (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( q_Lswp                 (ADM_gall_in,ADM_kall,TRC_VMAX,ADM_lall) )
  allocate( precip_mp              (ADM_gall_in,ADM_KNONE,ADM_lall,2)       )
  allocate( precip1_mp             (ADM_gall_in,ADM_KNONE,ADM_lall,2)       )
  allocate( precip2_mp             (ADM_gall_in,ADM_KNONE,ADM_lall,2)       )
  allocate( rhoein_precip_mp       (ADM_gall_in,ADM_KNONE,ADM_lall)         )
  allocate( lh_precip_mp           (ADM_gall_in,ADM_KNONE,ADM_lall)         )
  allocate( rhophi_precip_mp       (ADM_gall_in,ADM_KNONE,ADM_lall)         )
  allocate( rhokin_precip_mp       (ADM_gall_in,ADM_KNONE,ADM_lall)         )
  allocate( frhoge_af              (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( frhogqv_af             (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( frhoge_rad             (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( qke                    (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( gsgam2                 (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( gsgam2h                (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( gam2                   (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( gam2h                  (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( ix                     (ADM_gall_in,ADM_lall)                   )
  allocate( iy                     (ADM_gall_in,ADM_lall)                   )
  allocate( iz                     (ADM_gall_in,ADM_lall)                   )
  allocate( jx                     (ADM_gall_in,ADM_lall)                   )
  allocate( jy                     (ADM_gall_in,ADM_lall)                   )
  allocate( jz                     (ADM_gall_in,ADM_lall)                   )
  allocate( z                      (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( zh                     (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( GPREC                  (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( CBMFX                  (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( qd                     (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( rceff                  (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( rceff_solid            (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( rceff_cld              (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( rctop                  (ADM_gall_in,ADM_KNONE,ADM_lall)         )
  allocate( rwtop                  (ADM_gall_in,ADM_KNONE,ADM_lall)         )
  allocate( tctop                  (ADM_gall_in,ADM_KNONE,ADM_lall)         )
  allocate( GDCLW                  (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( GDCFRC                 (ADM_gall_in,ADM_kall,ADM_lall)          )

  allocate( CHECK_rhog             (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( CHECK_rhogvx           (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( CHECK_rhogvy           (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( CHECK_rhogvz           (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( CHECK_rhogw            (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( CHECK_rhoge            (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( CHECK_rhogq_Lswp       (ADM_gall_in,ADM_kall,TRC_VMAX,ADM_lall) )
  allocate( CHECK_rho              (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( CHECK_pre              (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( CHECK_tem              (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( CHECK_q_Lswp           (ADM_gall_in,ADM_kall,TRC_VMAX,ADM_lall) )
  allocate( CHECK_qd               (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( CHECK_precip_mp        (ADM_gall_in,ADM_KNONE,ADM_lall,2)       )
  allocate( CHECK_precip1_mp       (ADM_gall_in,ADM_KNONE,ADM_lall,2)       )
  allocate( CHECK_precip2_mp       (ADM_gall_in,ADM_KNONE,ADM_lall,2)       )
  allocate( CHECK_rhoein_precip_mp (ADM_gall_in,ADM_KNONE,ADM_lall)         )
  allocate( CHECK_lh_precip_mp     (ADM_gall_in,ADM_KNONE,ADM_lall)         )
  allocate( CHECK_rhophi_precip_mp (ADM_gall_in,ADM_KNONE,ADM_lall)         )
  allocate( CHECK_rhokin_precip_mp (ADM_gall_in,ADM_KNONE,ADM_lall)         )
  allocate( CHECK_rceff            (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( CHECK_rceff_solid      (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( CHECK_rceff_cld        (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( CHECK_rctop            (ADM_gall_in,ADM_KNONE,ADM_lall)         )
  allocate( CHECK_rwtop            (ADM_gall_in,ADM_KNONE,ADM_lall)         )
  allocate( CHECK_tctop            (ADM_gall_in,ADM_KNONE,ADM_lall)         )
  allocate( CHECK_GDCLW            (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( CHECK_GDCFRC           (ADM_gall_in,ADM_kall,ADM_lall)          )
  allocate( CHECK_GPREC            (ADM_gall_in,ADM_kall,ADM_lall)          )

  l = SET_l

  !###############################################################################

   EX_rgnid = SET_rgnid

   call dumpio_syscheck
   call dumpio_mk_fname(EX_fname,'snapshot.microphysics','pe',EX_rgnid-1,6)
   call dumpio_fopen(EX_fid,EX_fname,IO_FREAD)

   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); rhog            (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); rhogvx          (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); rhogvy          (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); rhogvz          (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); rhogw           (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); rhoge           (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall*TRC_VMAX, work2(:,:,:,l)         ); rhogq_Lswp      (1:ADM_gall_in,:,:,l)         = real(work2(1:ADM_gall_in,:,:,l)        ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); vx              (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); vy              (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); vz              (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); w               (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); unccn           (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); rho             (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); tem             (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); pre             (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall*TRC_VMAX, work2(:,:,:,l)         ); q_Lswp          (1:ADM_gall_in,:,:,l)         = real(work2(1:ADM_gall_in,:,:,l)        ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*2                , work4(:,ADM_KNONE,l,:) ); precip_mp       (1:ADM_gall_in,ADM_KNONE,l,:) = real(work4(1:ADM_gall_in,ADM_KNONE,l,:),kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*2                , work4(:,ADM_KNONE,l,:) ); precip1_mp      (1:ADM_gall_in,ADM_KNONE,l,:) = real(work4(1:ADM_gall_in,ADM_KNONE,l,:),kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*2                , work4(:,ADM_KNONE,l,:) ); precip2_mp      (1:ADM_gall_in,ADM_KNONE,l,:) = real(work4(1:ADM_gall_in,ADM_KNONE,l,:),kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in                  , work3(:,l)             ); rhoein_precip_mp(1:ADM_gall_in,ADM_KNONE,l)   = real(work3(1:ADM_gall_in,l)            ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in                  , work3(:,l)             ); lh_precip_mp    (1:ADM_gall_in,ADM_KNONE,l)   = real(work3(1:ADM_gall_in,l)            ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in                  , work3(:,l)             ); rhophi_precip_mp(1:ADM_gall_in,ADM_KNONE,l)   = real(work3(1:ADM_gall_in,l)            ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in                  , work3(:,l)             ); rhokin_precip_mp(1:ADM_gall_in,ADM_KNONE,l)   = real(work3(1:ADM_gall_in,l)            ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); frhoge_af       (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); frhogqv_af      (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); frhoge_rad      (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); qke             (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); gsgam2          (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); gsgam2h         (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); gam2            (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); gam2h           (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in                  , work3(:,l)             ); ix              (1:ADM_gall_in,l)             = real(work3(1:ADM_gall_in,l)            ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in                  , work3(:,l)             ); iy              (1:ADM_gall_in,l)             = real(work3(1:ADM_gall_in,l)            ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in                  , work3(:,l)             ); iz              (1:ADM_gall_in,l)             = real(work3(1:ADM_gall_in,l)            ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in                  , work3(:,l)             ); jx              (1:ADM_gall_in,l)             = real(work3(1:ADM_gall_in,l)            ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in                  , work3(:,l)             ); jy              (1:ADM_gall_in,l)             = real(work3(1:ADM_gall_in,l)            ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in                  , work3(:,l)             ); jz              (1:ADM_gall_in,l)             = real(work3(1:ADM_gall_in,l)            ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); z               (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); zh              (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); GPREC           (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); CBMFX           (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)

   call dumpio_fclose(EX_fid)

  !###############################################################################

  CONST_PI  = 4.E0_RP * atan( 1.0_RP )
  PI        = CONST_PI
  CONST_EPS = epsilon(0.0_RP)
  EPS       = CONST_EPS

  !--- vertical grid setup
  call GRD_setup

  call SATURATION_setup

  !---  microphysics initialization
  call mp_init( MP_TYPE )

  !###############################################################################

  write(*,*) "*** Finish initialize"
  write(*,*) "*** Start  kernel"
  call PROF_rapstart('Cloud_Microphysics_kernel')

  do iteration = 1, SET_iteration
     write(*,*) '### Before ###'
     call PROF_valcheck( 'before', 'rhog  ', rhog  (:,:,:) )
     call PROF_valcheck( 'before', 'rhogvx', rhogvx(:,:,:) )
     call PROF_valcheck( 'before', 'rhogvy', rhogvy(:,:,:) )
     call PROF_valcheck( 'before', 'rhogvz', rhogvz(:,:,:) )
     call PROF_valcheck( 'before', 'rhogw ', rhogw (:,:,:) )
     call PROF_valcheck( 'before', 'rhoge ', rhoge (:,:,:) )
     do nq = 1, TRC_VMAX
        call PROF_valcheck( 'before', TRC_name(nq), rhogq_Lswp(:,:,nq,:) )
     enddo

       call mp_driver( ADM_gall_in,                       & ! [IN]
                       l,                                 & ! [IN]
                       rhog            (:,:,l),           & ! [INOUT]
                       rhogvx          (:,:,l),           & ! [INOUT]
                       rhogvy          (:,:,l),           & ! [INOUT]
                       rhogvz          (:,:,l),           & ! [INOUT]
                       rhogw           (:,:,l),           & ! [INOUT]
                       rhoge           (:,:,l),           & ! [INOUT]
                       rhogq_Lswp      (:,:,:,l),         & ! [INOUT]
                       vx              (:,:,l),           & ! [IN]
                       vy              (:,:,l),           & ! [IN]
                       vz              (:,:,l),           & ! [IN]
                       w               (:,:,l),           & ! [IN]
                       unccn           (:,:,l),           & ! [IN]
                       rho             (:,:,l),           & ! [INOUT]
                       tem             (:,:,l),           & ! [INOUT]
                       pre             (:,:,l),           & ! [INOUT]
                       q_Lswp          (:,:,:,l),         & ! [INOUT]
                       qd              (:,:,l),           & ! [OUT]
                       precip_mp       (:,ADM_KNONE,l,:), & ! [INOUT]
                       precip1_mp      (:,ADM_KNONE,l,:), & ! [INOUT]
                       precip2_mp      (:,ADM_KNONE,l,:), & ! [INOUT]
                       rhoein_precip_mp(:,ADM_KNONE,l),   & ! [INOUT]
                       lh_precip_mp    (:,ADM_KNONE,l),   & ! [INOUT]
                       rhophi_precip_mp(:,ADM_KNONE,l),   & ! [INOUT]
                       rhokin_precip_mp(:,ADM_KNONE,l),   & ! [INOUT]
                       rceff           (:,:,l),           & ! [OUT]
                       rceff_solid     (:,:,l),           & ! [OUT]
                       rceff_cld       (:,:,l),           & ! [OUT]
                       rctop           (:,:,l),           & ! [OUT]
                       rwtop           (:,:,l),           & ! [OUT]
                       tctop           (:,:,l),           & ! [OUT]
                       frhoge_af       (:,:,l),           & ! [IN]
                       frhogqv_af      (:,:,l),           & ! [IN]
                       frhoge_rad      (:,:,l),           & ! [IN]
                       qke             (:,:,l),           & ! [IN]
                       gsgam2          (:,:,l),           & ! [IN]
                       gsgam2h         (:,:,l),           & ! [IN]
                       gam2            (:,:,l),           & ! [IN]
                       gam2h           (:,:,l),           & ! [IN]
                       ix              (:,l),             & ! [IN]
                       iy              (:,l),             & ! [IN]
                       iz              (:,l),             & ! [IN]
                       jx              (:,l),             & ! [IN]
                       jy              (:,l),             & ! [IN]
                       jz              (:,l),             & ! [IN]
                       z               (:,:,l),           & ! [IN]
                       zh              (:,:,l),           & ! [IN]
                       TIME_DTL,                          & ! [IN]
                       TIME_CTIME,                        & ! [IN]
                       GDCLW           (:,:,l),           & ! [OUT]
                       GDCFRC          (:,:,l),           & ! [OUT]
                       GPREC           (:,:,l),           & ! [INOUT]
                       CBMFX           (:,:,l)            ) ! [IN]

     write(*,*) '### After ###'
     call PROF_valcheck( 'after', 'rhog  ', rhog  (:,:,:) )
     call PROF_valcheck( 'after', 'rhogvx', rhogvx(:,:,:) )
     call PROF_valcheck( 'after', 'rhogvy', rhogvy(:,:,:) )
     call PROF_valcheck( 'after', 'rhogvz', rhogvz(:,:,:) )
     call PROF_valcheck( 'after', 'rhogw ', rhogw (:,:,:) )
     call PROF_valcheck( 'after', 'rhoge ', rhoge (:,:,:) )
     do nq = 1, TRC_VMAX
        call PROF_valcheck( 'after', TRC_name(nq), rhogq_Lswp(:,:,nq,:) )
     enddo
  enddo

  call PROF_rapend('Cloud_Microphysics_kernel')
  write(*,*) "*** Finish kernel"

  !###############################################################################
  if ( SET_check .AND. SET_iteration == 1 ) then
   EX_rgnid = SET_rgnid

   call dumpio_mk_fname(EX_fname,'check.microphysics','pe',EX_rgnid-1,6)
   call dumpio_fopen(EX_fid,EX_fname,IO_FREAD)

   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); CHECK_rhog            (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); CHECK_rhogvx          (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); CHECK_rhogvy          (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); CHECK_rhogvz          (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); CHECK_rhogw           (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); CHECK_rhoge           (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall*TRC_VMAX, work2(:,:,:,l)         ); CHECK_rhogq_Lswp      (1:ADM_gall_in,:,:,l)         = real(work2(1:ADM_gall_in,:,:,l)        ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); CHECK_rho             (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); CHECK_tem             (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); CHECK_pre             (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall*TRC_VMAX, work2(:,:,:,l)         ); CHECK_q_Lswp          (1:ADM_gall_in,:,:,l)         = real(work2(1:ADM_gall_in,:,:,l)        ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); CHECK_qd              (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*2                , work4(:,ADM_KNONE,l,:) ); CHECK_precip_mp       (1:ADM_gall_in,ADM_KNONE,l,:) = real(work4(1:ADM_gall_in,ADM_KNONE,l,:),kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*2                , work4(:,ADM_KNONE,l,:) ); CHECK_precip1_mp      (1:ADM_gall_in,ADM_KNONE,l,:) = real(work4(1:ADM_gall_in,ADM_KNONE,l,:),kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*2                , work4(:,ADM_KNONE,l,:) ); CHECK_precip2_mp      (1:ADM_gall_in,ADM_KNONE,l,:) = real(work4(1:ADM_gall_in,ADM_KNONE,l,:),kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in                  , work3(:,l)             ); CHECK_rhoein_precip_mp(1:ADM_gall_in,ADM_KNONE,l)   = real(work3(1:ADM_gall_in,l)            ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in                  , work3(:,l)             ); CHECK_lh_precip_mp    (1:ADM_gall_in,ADM_KNONE,l)   = real(work3(1:ADM_gall_in,l)            ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in                  , work3(:,l)             ); CHECK_rhophi_precip_mp(1:ADM_gall_in,ADM_KNONE,l)   = real(work3(1:ADM_gall_in,l)            ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in                  , work3(:,l)             ); CHECK_rhokin_precip_mp(1:ADM_gall_in,ADM_KNONE,l)   = real(work3(1:ADM_gall_in,l)            ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); CHECK_rceff           (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); CHECK_rceff_solid     (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); CHECK_rceff_cld       (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in                  , work3(:,l)             ); CHECK_rctop           (1:ADM_gall_in,ADM_KNONE,l)   = real(work3(1:ADM_gall_in,l)            ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in                  , work3(:,l)             ); CHECK_rwtop           (1:ADM_gall_in,ADM_KNONE,l)   = real(work3(1:ADM_gall_in,l)            ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in                  , work3(:,l)             ); CHECK_tctop           (1:ADM_gall_in,ADM_KNONE,l)   = real(work3(1:ADM_gall_in,l)            ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); CHECK_GDCLW           (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); CHECK_GDCFRC          (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)
   call dumpio_read_data( EX_fid, ADM_gall_in*ADM_kall         , work1(:,:,l)           ); CHECK_GPREC           (1:ADM_gall_in,:,l)           = real(work1(1:ADM_gall_in,:,l)          ,kind=RP)

   call dumpio_fclose(EX_fid)

     CHECK_rhog  (:,:,:) = rhog  (:,:,:) - CHECK_rhog  (:,:,:)
     CHECK_rhogvx(:,:,:) = rhogvx(:,:,:) - CHECK_rhogvx(:,:,:)
     CHECK_rhogvy(:,:,:) = rhogvy(:,:,:) - CHECK_rhogvy(:,:,:)
     CHECK_rhogvz(:,:,:) = rhogvz(:,:,:) - CHECK_rhogvz(:,:,:)
     CHECK_rhogw (:,:,:) = rhogw (:,:,:) - CHECK_rhogw (:,:,:)
     CHECK_rhoge (:,:,:) = rhoge (:,:,:) - CHECK_rhoge (:,:,:)
     CHECK_rhogq_Lswp(:,:,:,:) = rhogq_Lswp(:,:,:,:) - CHECK_rhogq_Lswp(:,:,:,:)

     write(*,*) '### Check ###'
     call PROF_valcheck( 'check', 'check_rhog  ', CHECK_rhog  (:,:,:) )
     call PROF_valcheck( 'check', 'check_rhogvx', CHECK_rhogvx(:,:,:) )
     call PROF_valcheck( 'check', 'check_rhogvy', CHECK_rhogvy(:,:,:) )
     call PROF_valcheck( 'check', 'check_rhogvz', CHECK_rhogvz(:,:,:) )
     call PROF_valcheck( 'check', 'check_rhogw ', CHECK_rhogw (:,:,:) )
     call PROF_valcheck( 'check', 'check_rhoge ', CHECK_rhoge (:,:,:) )
     do nq = 1, TRC_VMAX
        call PROF_valcheck( 'check', TRC_name(nq), CHECK_rhogq_Lswp(:,:,nq,:) )
     enddo

  endif
  !###############################################################################

  call PROF_rapreport

  stop
end program physicskernel_microphysics
!-------------------------------------------------------------------------------
