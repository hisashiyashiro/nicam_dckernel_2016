!-------------------------------------------------------------------------------
!
!+  Program dynamics kernel driver (diffusion operator)
!
!-------------------------------------------------------------------------------
program dckernel_diffusion
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_misc
  use mod_oprt, only: &
     OPRT_diffusion
  !-----------------------------------------------------------------------------
  implicit none

  real(DP), allocatable :: ORG_dscl        (:,:,:,:)
  real(DP), allocatable :: ORG_dscl_pl     (:,:,:)
  real(DP), allocatable :: ORG_scl         (:,:,:,:)
  real(DP), allocatable :: ORG_scl_pl      (:,:,:)
  real(DP), allocatable :: ORG_kh          (:,:,:,:)
  real(DP), allocatable :: ORG_kh_pl       (:,:,:)
  real(DP), allocatable :: ORG_coef_intp   (:,:,:,:,:,:)
  real(DP), allocatable :: ORG_coef_intp_pl(:,:,:,:)
  real(DP), allocatable :: ORG_coef_diff   (:,:,:,:,:)
  real(DP), allocatable :: ORG_coef_diff_pl(:,:,:)

  real(RP), allocatable :: dscl        (:,:,:,:)
  real(RP), allocatable :: dscl_pl     (:,:,:)
  real(RP), allocatable :: scl         (:,:,:,:)
  real(RP), allocatable :: scl_pl      (:,:,:)
  real(RP), allocatable :: kh          (:,:,:,:)
  real(RP), allocatable :: kh_pl       (:,:,:)
  real(RP), allocatable :: coef_intp   (:,:,:,:,:,:)
  real(RP), allocatable :: coef_intp_pl(:,:,:,:)
  real(RP), allocatable :: coef_diff   (:,:,:,:,:)
  real(RP), allocatable :: coef_diff_pl(:,:,:)

  real(RP), allocatable :: check_dscl   (:,:,:,:)
  real(RP), allocatable :: check_dscl_pl(:,:,:)

  integer :: iteration
  !=============================================================================

  write(*,*) "[KERNEL] dckernel_diffusion"

  !###############################################################################

  write(*,*) "*** Start  initialize"

  allocate( dscl        (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )               )
  allocate( dscl_pl     (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)               )
  allocate( scl         (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )               )
  allocate( scl_pl      (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)               )
  allocate( kh          (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )               )
  allocate( kh_pl       (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)               )
  allocate( coef_intp   (ADM_iall,ADM_jall,1:3    ,ADM_nxyz,TI:TJ,ADM_lall   ) )
  allocate( coef_intp_pl(ADM_gall_pl      ,1:3    ,ADM_nxyz,      ADM_lall_pl) )
  allocate( coef_diff   (ADM_iall,ADM_jall,1:6    ,ADM_nxyz,ADM_lall   )       )
  allocate( coef_diff_pl(              1:ADM_vlink,ADM_nxyz,ADM_lall_pl)       )
  ! Todo : first touch with considering NUMA
  dscl        (:,:,:,:)     = 0.0_RP
  dscl_pl     (:,:,:)       = 0.0_RP
  scl         (:,:,:,:)     = 0.0_RP
  scl_pl      (:,:,:)       = 0.0_RP
  kh          (:,:,:,:)     = 0.0_RP
  kh_pl       (:,:,:)       = 0.0_RP
  coef_intp   (:,:,:,:,:,:) = 0.0_RP
  coef_intp_pl(:,:,:,:)     = 0.0_RP
  coef_diff   (:,:,:,:,:)   = 0.0_RP
  coef_diff_pl(:,:,:)       = 0.0_RP

  allocate( check_dscl    (ADM_iall,ADM_jall,ADM_kall,ADM_lall   ) )
  allocate( check_dscl_pl (ADM_gall_pl      ,ADM_kall,ADM_lall_pl) )
  check_dscl    (:,:,:,:)   = 0.0_RP
  check_dscl_pl (:,:,:)     = 0.0_RP

  !###############################################################################

  !---< read input data >---
  allocate( ORG_dscl        (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )               )
  allocate( ORG_dscl_pl     (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)               )
  allocate( ORG_scl         (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )               )
  allocate( ORG_scl_pl      (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)               )
  allocate( ORG_kh          (ADM_iall,ADM_jall,ADM_kall,ADM_lall   )               )
  allocate( ORG_kh_pl       (ADM_gall_pl      ,ADM_kall,ADM_lall_pl)               )
  allocate( ORG_coef_intp   (ADM_iall,ADM_jall,1:3    ,ADM_nxyz,TI:TJ,ADM_lall   ) )
  allocate( ORG_coef_intp_pl(ADM_gall_pl      ,1:3    ,ADM_nxyz,      ADM_lall_pl) )
  allocate( ORG_coef_diff   (ADM_iall,ADM_jall,1:6    ,ADM_nxyz,ADM_lall   )       )
  allocate( ORG_coef_diff_pl(              1:ADM_vlink,ADM_nxyz,ADM_lall_pl)       )

  call dumpio_syscheck
  call dumpio_mk_fname(EX_fname,'snapshot.dc_diffusion','pe',SET_prc_me-1,6)
  call dumpio_fopen(EX_fid,EX_fname,IO_FREAD)

  call dumpio_read_data( EX_fid, ADM_iall*ADM_jall*ADM_kall*ADM_lall   ,     ORG_dscl   (:,:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl      *ADM_kall*ADM_lall_pl,     ORG_dscl_pl(:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_iall*ADM_jall*ADM_kall*ADM_lall   ,     ORG_scl    (:,:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl      *ADM_kall*ADM_lall_pl,     ORG_scl_pl (:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_iall*ADM_jall*ADM_kall*ADM_lall   ,     ORG_kh     (:,:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl      *ADM_kall*ADM_lall_pl,     ORG_kh_pl  (:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_iall*ADM_jall*3*ADM_nxyz*2*ADM_lall   , ORG_coef_intp   (:,:,:,:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl      *3*ADM_nxyz*  ADM_lall_pl, ORG_coef_intp_pl(:,:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_iall*ADM_jall*6*ADM_nxyz*  ADM_lall   , ORG_coef_diff   (:,:,:,:,:)   )
  call dumpio_read_data( EX_fid,           ADM_vlink*ADM_nxyz*  ADM_lall_pl, ORG_coef_diff_pl(:,:,:)       )

  call dumpio_fclose(EX_fid)

  dscl        (:,:,:,:)     = real( ORG_dscl        (:,:,:,:)    , kind=RP )
  dscl_pl     (:,:,:)       = real( ORG_dscl_pl     (:,:,:)      , kind=RP )
  scl         (:,:,:,:)     = real( ORG_scl         (:,:,:,:)    , kind=RP )
  scl_pl      (:,:,:)       = real( ORG_scl_pl      (:,:,:)      , kind=RP )
  kh          (:,:,:,:)     = real( ORG_kh          (:,:,:,:)    , kind=RP )
  kh_pl       (:,:,:)       = real( ORG_kh_pl       (:,:,:)      , kind=RP )
  coef_intp   (:,:,:,:,:,:) = real( ORG_coef_intp   (:,:,:,:,:,:), kind=RP )
  coef_intp_pl(:,:,:,:)     = real( ORG_coef_intp_pl(:,:,:,:)    , kind=RP )
  coef_diff   (:,:,:,:,:)   = real( ORG_coef_diff   (:,:,:,:,:)  , kind=RP )
  coef_diff_pl(:,:,:)       = real( ORG_coef_diff_pl(:,:,:)      , kind=RP )

  deallocate( ORG_dscl         )
  deallocate( ORG_dscl_pl      )
  deallocate( ORG_scl          )
  deallocate( ORG_scl_pl       )
  deallocate( ORG_kh           )
  deallocate( ORG_kh_pl        )
  deallocate( ORG_coef_intp    )
  deallocate( ORG_coef_intp_pl )
  deallocate( ORG_coef_diff    )
  deallocate( ORG_coef_diff_pl )

  !###############################################################################

  check_dscl   (:,:,:,:) = dscl   (:,:,:,:)
  check_dscl_pl(:,:,:)   = dscl_pl(:,:,:)

  write(*,*) "*** Finish initialize"

  !###############################################################################

  write(*,*) "*** Start  kernel"
  call DEBUG_rapstart('DC_diffusion_kernel')

  do iteration = 1, SET_iteration



     call OPRT_diffusion( dscl     (:,:,:,:),     dscl_pl     (:,:,:),   & ! [OUT]
                          scl      (:,:,:,:),     scl_pl      (:,:,:),   & ! [IN]
                          kh       (:,:,:,:),     kh_pl       (:,:,:),   & ! [IN]
                          coef_intp(:,:,:,:,:,:), coef_intp_pl(:,:,:,:), & ! [IN]
                          coef_diff(:,:,:,:,:),   coef_diff_pl(:,:,:)    ) ! [IN]



     write(ADM_LOG_FID,*) '### Input ###'
     call DEBUG_valuecheck( 'check_dscl   ', check_dscl   (:,:,:,:) )
     call DEBUG_valuecheck( 'check_dscl_pl', check_dscl_pl(:,:,:)   )
     call DEBUG_valuecheck( 'scl          ', scl          (:,:,:,:) )
     call DEBUG_valuecheck( 'scl_pl       ', scl_pl       (:,:,:)   )
     call DEBUG_valuecheck( 'kh           ', kh           (:,:,:,:) )
     call DEBUG_valuecheck( 'kh_pl        ', kh_pl        (:,:,:)   )
     write(ADM_LOG_FID,*) '### Output ###'
     call DEBUG_valuecheck( 'dscl         ', dscl         (:,:,:,:) )
     call DEBUG_valuecheck( 'dscl_pl      ', dscl_pl      (:,:,:)   )
  enddo

  write(ADM_LOG_FID,*) '### Varidation by diff ###'
  check_dscl   (:,:,:,:) = check_dscl   (:,:,:,:) - dscl   (:,:,:,:)
  check_dscl_pl(:,:,:)   = check_dscl_pl(:,:,:)   - dscl_pl(:,:,:)
  call DEBUG_valuecheck( 'check_dscl   ', check_dscl   (:,:,:,:) )
  call DEBUG_valuecheck( 'check_dscl_pl', check_dscl_pl(:,:,:)   )

  call DEBUG_rapend('DC_diffusion_kernel')
  write(*,*) "*** Finish kernel"

  !###############################################################################

  call DEBUG_rapreport

  stop
end program dckernel_diffusion
!-------------------------------------------------------------------------------
