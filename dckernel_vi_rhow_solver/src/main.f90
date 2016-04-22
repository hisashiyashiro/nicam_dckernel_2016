!-------------------------------------------------------------------------------
!
!+  Program dynamics kernel driver (vi_rhow_solver operator)
!
!-------------------------------------------------------------------------------
program dckernel_vi_rhow_solver
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_misc
  use mod_vi, only: &
     vi_rhow_solver_putcoef, &
     vi_rhow_solver
  !-----------------------------------------------------------------------------
  implicit none

  real(RP), allocatable :: rhogw_prev      (:,:,:)
  real(RP), allocatable :: rhogw_prev_pl   (:,:,:)
  real(RP), allocatable :: rhogw           (:,:,:)
  real(RP), allocatable :: rhogw_pl        (:,:,:)
  real(RP), allocatable :: rhogw0          (:,:,:)
  real(RP), allocatable :: rhogw0_pl       (:,:,:)
  real(RP), allocatable :: preg0           (:,:,:)
  real(RP), allocatable :: preg0_pl        (:,:,:)
  real(RP), allocatable :: rhog0           (:,:,:)
  real(RP), allocatable :: rhog0_pl        (:,:,:)
  real(RP), allocatable :: Srho            (:,:,:)
  real(RP), allocatable :: Srho_pl         (:,:,:)
  real(RP), allocatable :: Sw              (:,:,:)
  real(RP), allocatable :: Sw_pl           (:,:,:)
  real(RP), allocatable :: Spre            (:,:,:)
  real(RP), allocatable :: Spre_pl         (:,:,:)
  real(RP), allocatable :: Mc              (:,:,:)
  real(RP), allocatable :: Mc_pl           (:,:,:)
  real(RP), allocatable :: Ml              (:,:,:)
  real(RP), allocatable :: Ml_pl           (:,:,:)
  real(RP), allocatable :: Mu              (:,:,:)
  real(RP), allocatable :: Mu_pl           (:,:,:)
  real(RP), allocatable :: VMTR_RGSGAM2    (:,:,:)
  real(RP), allocatable :: VMTR_RGSGAM2_pl (:,:,:)
  real(RP), allocatable :: VMTR_RGSGAM2H   (:,:,:)
  real(RP), allocatable :: VMTR_RGSGAM2H_pl(:,:,:)
  real(RP), allocatable :: VMTR_RGAMH      (:,:,:)
  real(RP), allocatable :: VMTR_RGAMH_pl   (:,:,:)
  real(RP), allocatable :: VMTR_RGAM       (:,:,:)
  real(RP), allocatable :: VMTR_RGAM_pl    (:,:,:)
  real(RP), allocatable :: VMTR_GSGAM2H    (:,:,:)
  real(RP), allocatable :: VMTR_GSGAM2H_pl (:,:,:)
  real(RP), allocatable :: check_rhogw     (:,:,:)
  real(RP), allocatable :: check_rhogw_pl  (:,:,:)

  integer :: iteration
  !=============================================================================

  write(*,*) "[KERNEL] dckernel_vi_rhow_solver"

  !###############################################################################

  write(*,*) "*** Start  initialize"

  allocate( rhogw_prev      (ADM_gall   ,ADM_kall,ADM_lall   ) )
  allocate( rhogw_prev_pl   (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
  allocate( rhogw           (ADM_gall   ,ADM_kall,ADM_lall   ) )
  allocate( rhogw_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
  allocate( rhogw0          (ADM_gall   ,ADM_kall,ADM_lall   ) )
  allocate( rhogw0_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
  allocate( preg0           (ADM_gall   ,ADM_kall,ADM_lall   ) )
  allocate( preg0_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
  allocate( rhog0           (ADM_gall   ,ADM_kall,ADM_lall   ) )
  allocate( rhog0_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
  allocate( Srho            (ADM_gall   ,ADM_kall,ADM_lall   ) )
  allocate( Srho_pl         (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
  allocate( Sw              (ADM_gall   ,ADM_kall,ADM_lall   ) )
  allocate( Sw_pl           (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
  allocate( Spre            (ADM_gall   ,ADM_kall,ADM_lall   ) )
  allocate( Spre_pl         (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
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
  ! Todo : first touch with considering NUMA
  rhogw_prev      (:,:,:)   = 0.0_RP
  rhogw_prev_pl   (:,:,:)   = 0.0_RP
  rhogw           (:,:,:)   = 0.0_RP
  rhogw_pl        (:,:,:)   = 0.0_RP
  rhogw0          (:,:,:)   = 0.0_RP
  rhogw0_pl       (:,:,:)   = 0.0_RP
  preg0           (:,:,:)   = 0.0_RP
  preg0_pl        (:,:,:)   = 0.0_RP
  rhog0           (:,:,:)   = 0.0_RP
  rhog0_pl        (:,:,:)   = 0.0_RP
  Srho            (:,:,:)   = 0.0_RP
  Srho_pl         (:,:,:)   = 0.0_RP
  Sw              (:,:,:)   = 0.0_RP
  Sw_pl           (:,:,:)   = 0.0_RP
  Spre            (:,:,:)   = 0.0_RP
  Spre_pl         (:,:,:)   = 0.0_RP
  Mc              (:,:,:)   = 0.0_RP
  Mc_pl           (:,:,:)   = 0.0_RP
  Ml              (:,:,:)   = 0.0_RP
  Ml_pl           (:,:,:)   = 0.0_RP
  Mu              (:,:,:)   = 0.0_RP
  Mu_pl           (:,:,:)   = 0.0_RP
  VMTR_RGSGAM2    (:,:,:)   = 0.0_RP
  VMTR_RGSGAM2_pl (:,:,:)   = 0.0_RP
  VMTR_RGSGAM2H   (:,:,:)   = 0.0_RP
  VMTR_RGSGAM2H_pl(:,:,:)   = 0.0_RP
  VMTR_RGAMH      (:,:,:)   = 0.0_RP
  VMTR_RGAMH_pl   (:,:,:)   = 0.0_RP
  VMTR_RGAM       (:,:,:)   = 0.0_RP
  VMTR_RGAM_pl    (:,:,:)   = 0.0_RP
  VMTR_GSGAM2H    (:,:,:)   = 0.0_RP
  VMTR_GSGAM2H_pl (:,:,:)   = 0.0_RP

  allocate( check_rhogw     (ADM_gall   ,ADM_kall,ADM_lall   ) )
  allocate( check_rhogw_pl  (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
  check_rhogw     (:,:,:)   = 0.0_RP
  check_rhogw_pl  (:,:,:)   = 0.0_RP

  call MISC_make_idstr(EX_fname,'snapshot.dc_vi_rhow_solver','pe',SET_prc_me)
  EX_fid = MISC_get_available_fid()
  open( unit   = EX_fid,         &
        file   = trim(EX_fname), &
        form   = 'unformatted',  &
        access = 'sequential',   &
        status = 'old'           )

     read(EX_fid) rhogw_prev      (:,:,:)
     read(EX_fid) rhogw_prev_pl   (:,:,:)
     read(EX_fid) rhogw           (:,:,:)
     read(EX_fid) rhogw_pl        (:,:,:)
     read(EX_fid) rhogw0          (:,:,:)
     read(EX_fid) rhogw0_pl       (:,:,:)
     read(EX_fid) preg0           (:,:,:)
     read(EX_fid) preg0_pl        (:,:,:)
     read(EX_fid) rhog0           (:,:,:)
     read(EX_fid) rhog0_pl        (:,:,:)
     read(EX_fid) Srho            (:,:,:)
     read(EX_fid) Srho_pl         (:,:,:)
     read(EX_fid) Sw              (:,:,:)
     read(EX_fid) Sw_pl           (:,:,:)
     read(EX_fid) Spre            (:,:,:)
     read(EX_fid) Spre_pl         (:,:,:)
     read(EX_fid) Mc              (:,:,:)
     read(EX_fid) Mc_pl           (:,:,:)
     read(EX_fid) Ml              (:,:,:)
     read(EX_fid) Ml_pl           (:,:,:)
     read(EX_fid) Mu              (:,:,:)
     read(EX_fid) Mu_pl           (:,:,:)
     read(EX_fid) VMTR_RGSGAM2    (:,:,:)
     read(EX_fid) VMTR_RGSGAM2_pl (:,:,:)
     read(EX_fid) VMTR_RGSGAM2H   (:,:,:)
     read(EX_fid) VMTR_RGSGAM2H_pl(:,:,:)
     read(EX_fid) VMTR_RGAMH      (:,:,:)
     read(EX_fid) VMTR_RGAMH_pl   (:,:,:)
     read(EX_fid) VMTR_RGAM       (:,:,:)
     read(EX_fid) VMTR_RGAM_pl    (:,:,:)
     read(EX_fid) VMTR_GSGAM2H    (:,:,:)
     read(EX_fid) VMTR_GSGAM2H_pl (:,:,:)
  close(EX_fid)

  call GRD_setup ! allocate GRD_rdgz

  !---< operator module setup >---
  call vi_rhow_solver_putcoef( Mc           (:,:,:), Mc_pl           (:,:,:), & ! [IN]
                               Ml           (:,:,:), Ml_pl           (:,:,:), & ! [IN]
                               Mu           (:,:,:), Mu_pl           (:,:,:), & ! [IN]
                               VMTR_RGSGAM2 (:,:,:), VMTR_RGSGAM2_pl (:,:,:), & ! [IN]
                               VMTR_RGSGAM2H(:,:,:), VMTR_RGSGAM2H_pl(:,:,:), & ! [IN]
                               VMTR_RGAMH   (:,:,:), VMTR_RGAMH_pl   (:,:,:), & ! [IN]
                               VMTR_RGAM    (:,:,:), VMTR_RGAM_pl    (:,:,:), & ! [IN]
                               VMTR_GSGAM2H (:,:,:), VMTR_GSGAM2H_pl (:,:,:)  ) ! [IN]

  check_rhogw   (:,:,:) = rhogw   (:,:,:)
  check_rhogw_pl(:,:,:) = rhogw_pl(:,:,:)

  write(*,*) "*** Finish initialize"

  !###############################################################################

  write(*,*) "*** Start  kernel"
  call DEBUG_rapstart('DC_vi_rhow_solver_kernel')

  do iteration = 1, SET_iteration



     ! restore previous value for multiple iteration
     rhogw   (:,:,:) = rhogw_prev   (:,:,:)
     rhogw_pl(:,:,:) = rhogw_prev_pl(:,:,:)

     call vi_rhow_solver( rhogw (:,:,:),  rhogw_pl (:,:,:), & ! [INOUT]
                          rhogw0(:,:,:),  rhogw0_pl(:,:,:), & ! [IN]
                          preg0 (:,:,:),  preg0_pl (:,:,:), & ! [IN]
                          rhog0 (:,:,:),  rhog0_pl (:,:,:), & ! [IN]
                          Srho  (:,:,:),  Srho_pl  (:,:,:), & ! [IN]
                          Sw    (:,:,:),  Sw_pl    (:,:,:), & ! [IN]
                          Spre  (:,:,:),  Spre_pl  (:,:,:), & ! [IN]
                          SET_dt                            ) ! [IN]



     write(ADM_LOG_FID,*) '### Input ###'
     EX_item =       'rhogw_prev   '
     EX_max  = maxval(rhogw_prev   (:,:,:))
     EX_min  = minval(rhogw_prev   (:,:,:))
     EX_sum  = sum   (rhogw_prev   (:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
     EX_item =       'rhogw_prev_pl'
     EX_max  = maxval(rhogw_prev_pl(:,:,:))
     EX_min  = minval(rhogw_prev_pl(:,:,:))
     EX_sum  = sum   (rhogw_prev_pl(:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
     EX_item =       'check_rhogw'
     EX_max  = maxval(check_rhogw(:,:,:))
     EX_min  = minval(check_rhogw(:,:,:))
     EX_sum  = sum   (check_rhogw(:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
     EX_item =       'check_rhogw_pl'
     EX_max  = maxval(check_rhogw_pl(:,:,:))
     EX_min  = minval(check_rhogw_pl(:,:,:))
     EX_sum  = sum   (check_rhogw_pl(:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
     EX_item =       'rhogw0   '
     EX_max  = maxval(rhogw0   (:,:,:))
     EX_min  = minval(rhogw0   (:,:,:))
     EX_sum  = sum   (rhogw0   (:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
     EX_item =       'rhogw0_pl'
     EX_max  = maxval(rhogw0_pl(:,:,:))
     EX_min  = minval(rhogw0_pl(:,:,:))
     EX_sum  = sum   (rhogw0_pl(:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
     EX_item =       'preg0    '
     EX_max  = maxval(preg0    (:,:,:))
     EX_min  = minval(preg0    (:,:,:))
     EX_sum  = sum   (preg0    (:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
     EX_item =       'preg0_pl '
     EX_max  = maxval(preg0_pl (:,:,:))
     EX_min  = minval(preg0_pl (:,:,:))
     EX_sum  = sum   (preg0_pl (:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
     EX_item =       'Srho     '
     EX_max  = maxval(Srho     (:,:,:))
     EX_min  = minval(Srho     (:,:,:))
     EX_sum  = sum   (Srho     (:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
     EX_item =       'Srho_pl  '
     EX_max  = maxval(Srho_pl  (:,:,:))
     EX_min  = minval(Srho_pl  (:,:,:))
     EX_sum  = sum   (Srho_pl  (:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
     EX_item =       'Sw       '
     EX_max  = maxval(Sw       (:,:,:))
     EX_min  = minval(Sw       (:,:,:))
     EX_sum  = sum   (Sw       (:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
     EX_item =       'Sw_pl    '
     EX_max  = maxval(Sw_pl    (:,:,:))
     EX_min  = minval(Sw_pl    (:,:,:))
     EX_sum  = sum   (Sw_pl    (:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
     EX_item =       'Spre     '
     EX_max  = maxval(Spre     (:,:,:))
     EX_min  = minval(Spre     (:,:,:))
     EX_sum  = sum   (Spre     (:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
     EX_item =       'Spre_pl  '
     EX_max  = maxval(Spre_pl  (:,:,:))
     EX_min  = minval(Spre_pl  (:,:,:))
     EX_sum  = sum   (Spre_pl  (:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
     EX_item =       'Mc       '
     EX_max  = maxval(Mc       (:,:,:))
     EX_min  = minval(Mc       (:,:,:))
     EX_sum  = sum   (Mc       (:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
     EX_item =       'Mc_pl    '
     EX_max  = maxval(Mc_pl    (:,:,:))
     EX_min  = minval(Mc_pl    (:,:,:))
     EX_sum  = sum   (Mc_pl    (:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
     EX_item =       'Ml       '
     EX_max  = maxval(Ml       (:,:,:))
     EX_min  = minval(Ml       (:,:,:))
     EX_sum  = sum   (Ml       (:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
     EX_item =       'Ml_pl    '
     EX_max  = maxval(Ml_pl    (:,:,:))
     EX_min  = minval(Ml_pl    (:,:,:))
     EX_sum  = sum   (Ml_pl    (:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
     EX_item =       'Mu       '
     EX_max  = maxval(Mu       (:,:,:))
     EX_min  = minval(Mu       (:,:,:))
     EX_sum  = sum   (Mu       (:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
     EX_item =       'Mu_pl    '
     EX_max  = maxval(Mu_pl    (:,:,:))
     EX_min  = minval(Mu_pl    (:,:,:))
     EX_sum  = sum   (Mu_pl    (:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum

     write(ADM_LOG_FID,*) '### Output ###'
     EX_item =       'rhogw    '
     EX_max  = maxval(rhogw    (:,:,:))
     EX_min  = minval(rhogw    (:,:,:))
     EX_sum  = sum   (rhogw    (:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
     EX_item =       'rhogw_pl '
     EX_max  = maxval(rhogw_pl (:,:,:))
     EX_min  = minval(rhogw_pl (:,:,:))
     EX_sum  = sum   (rhogw_pl (:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
  enddo

  write(ADM_LOG_FID,*) '### Varidation : grid-by-grid diff ###'
  check_rhogw   (:,:,:) = check_rhogw   (:,:,:) - rhogw   (:,:,:)
  check_rhogw_pl(:,:,:) = check_rhogw_pl(:,:,:) - rhogw_pl(:,:,:)

  EX_item =       'check_rhogw'
  EX_max  = maxval(check_rhogw(:,:,:))
  EX_min  = minval(check_rhogw(:,:,:))
  EX_sum  = sum   (check_rhogw(:,:,:))
  write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
  EX_item =       'check_rhogw_pl'
  EX_max  = maxval(check_rhogw_pl(:,:,:))
  EX_min  = minval(check_rhogw_pl(:,:,:))
  EX_sum  = sum   (check_rhogw_pl(:,:,:))
  write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum

  call DEBUG_rapend('DC_vi_rhow_solver_kernel')
  write(*,*) "*** Finish kernel"

  call DEBUG_rapreport

  stop
end program dckernel_vi_rhow_solver
!-------------------------------------------------------------------------------
