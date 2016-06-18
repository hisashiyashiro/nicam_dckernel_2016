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
     OPRT_diffusion_putcoef, &
     OPRT_diffusion
  !-----------------------------------------------------------------------------
  implicit none

  real(DP), allocatable :: ORG_dscl          (:,:,:)
  real(DP), allocatable :: ORG_dscl_pl       (:,:,:)
  real(DP), allocatable :: ORG_scl           (:,:,:)
  real(DP), allocatable :: ORG_scl_pl        (:,:,:)
  real(DP), allocatable :: ORG_kh            (:,:,:)
  real(DP), allocatable :: ORG_kh_pl         (:,:,:)
  real(DP), allocatable :: ORG_cinterp_TN    (:,:,:,:)
  real(DP), allocatable :: ORG_cinterp_TN_pl (:,:,:,:)
  real(DP), allocatable :: ORG_cinterp_HN    (:,:,:,:)
  real(DP), allocatable :: ORG_cinterp_HN_pl (:,:,:)
  real(DP), allocatable :: ORG_cinterp_TRA   (:,:,:)
  real(DP), allocatable :: ORG_cinterp_TRA_pl(:,:)
  real(DP), allocatable :: ORG_cinterp_PRA   (:,:)
  real(DP), allocatable :: ORG_cinterp_PRA_pl(:,:)

  real(RP), allocatable :: dscl          (:,:,:)
  real(RP), allocatable :: dscl_pl       (:,:,:)
  real(RP), allocatable :: scl           (:,:,:)
  real(RP), allocatable :: scl_pl        (:,:,:)
  real(RP), allocatable :: kh            (:,:,:)
  real(RP), allocatable :: kh_pl         (:,:,:)
  real(RP), allocatable :: cinterp_TN    (:,:,:,:)
  real(RP), allocatable :: cinterp_TN_pl (:,:,:,:)
  real(RP), allocatable :: cinterp_HN    (:,:,:,:)
  real(RP), allocatable :: cinterp_HN_pl (:,:,:)
  real(RP), allocatable :: cinterp_TRA   (:,:,:)
  real(RP), allocatable :: cinterp_TRA_pl(:,:)
  real(RP), allocatable :: cinterp_PRA   (:,:)
  real(RP), allocatable :: cinterp_PRA_pl(:,:)

  real(RP), allocatable :: check_dscl   (:,:,:)
  real(RP), allocatable :: check_dscl_pl(:,:,:)

  integer :: iteration
  !=============================================================================

  write(*,*) "[KERNEL] dckernel_diffusion"

  !###############################################################################

  write(*,*) "*** Start  initialize"

  allocate( dscl          (ADM_gall   ,ADM_kall,ADM_lall   )       )
  allocate( dscl_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)       )
  allocate( scl           (ADM_gall   ,ADM_kall,ADM_lall   )       )
  allocate( scl_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)       )
  allocate( kh            (ADM_gall   ,ADM_kall,ADM_lall   )       )
  allocate( kh_pl         (ADM_gall_pl,ADM_kall,ADM_lall_pl)       )
  allocate( cinterp_TN    (ADM_gall   ,ADM_lall   ,AI:AJ,ADM_nxyz) )
  allocate( cinterp_TN_pl (ADM_gall_pl,ADM_lall_pl,2    ,ADM_nxyz) )
  allocate( cinterp_HN    (ADM_gall   ,ADM_lall   ,AI:AJ,ADM_nxyz) )
  allocate( cinterp_HN_pl (ADM_gall_pl,ADM_lall_pl,      ADM_nxyz) )
  allocate( cinterp_TRA   (ADM_gall   ,ADM_lall   ,TI:TJ         ) )
  allocate( cinterp_TRA_pl(ADM_gall_pl,ADM_lall_pl               ) )
  allocate( cinterp_PRA   (ADM_gall   ,ADM_lall                  ) )
  allocate( cinterp_PRA_pl(ADM_gall_pl,ADM_lall_pl               ) )
  ! Todo : first touch with considering NUMA
  dscl          (:,:,:)   = 0.0_RP
  dscl_pl       (:,:,:)   = 0.0_RP
  scl           (:,:,:)   = 0.0_RP
  scl_pl        (:,:,:)   = 0.0_RP
  kh            (:,:,:)   = 0.0_RP
  kh_pl         (:,:,:)   = 0.0_RP
  cinterp_TN    (:,:,:,:) = 0.0_RP
  cinterp_TN_pl (:,:,:,:) = 0.0_RP
  cinterp_HN    (:,:,:,:) = 0.0_RP
  cinterp_HN_pl (:,:,  :) = 0.0_RP
  cinterp_TRA   (:,:,:  ) = 0.0_RP
  cinterp_TRA_pl(:,:    ) = 0.0_RP
  cinterp_PRA   (:,:    ) = 0.0_RP
  cinterp_PRA_pl(:,:    ) = 0.0_RP

  allocate( check_dscl    (ADM_gall   ,ADM_kall,ADM_lall   ) )
  allocate( check_dscl_pl (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
  check_dscl    (:,:,:)   = 0.0_RP
  check_dscl_pl (:,:,:)   = 0.0_RP

  !###############################################################################

  !---< read input data >---
  allocate( ORG_dscl          (ADM_gall   ,ADM_kall,ADM_lall   )       )
  allocate( ORG_dscl_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl)       )
  allocate( ORG_scl           (ADM_gall   ,ADM_kall,ADM_lall   )       )
  allocate( ORG_scl_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl)       )
  allocate( ORG_kh            (ADM_gall   ,ADM_kall,ADM_lall   )       )
  allocate( ORG_kh_pl         (ADM_gall_pl,ADM_kall,ADM_lall_pl)       )
  allocate( ORG_cinterp_TN    (ADM_gall   ,ADM_lall   ,AI:AJ,ADM_nxyz) )
  allocate( ORG_cinterp_TN_pl (ADM_gall_pl,ADM_lall_pl,2    ,ADM_nxyz) )
  allocate( ORG_cinterp_HN    (ADM_gall   ,ADM_lall   ,AI:AJ,ADM_nxyz) )
  allocate( ORG_cinterp_HN_pl (ADM_gall_pl,ADM_lall_pl,      ADM_nxyz) )
  allocate( ORG_cinterp_TRA   (ADM_gall   ,ADM_lall   ,TI:TJ         ) )
  allocate( ORG_cinterp_TRA_pl(ADM_gall_pl,ADM_lall_pl               ) )
  allocate( ORG_cinterp_PRA   (ADM_gall   ,ADM_lall                  ) )
  allocate( ORG_cinterp_PRA_pl(ADM_gall_pl,ADM_lall_pl               ) )

  call MISC_make_idstr(EX_fname,'snapshot.dc_diffusion','pe',SET_prc_me)
  EX_fid = MISC_get_available_fid()
  open( unit   = EX_fid,         &
        file   = trim(EX_fname), &
        form   = 'unformatted',  &
        access = 'sequential',   &
        status = 'old'           )

     read(EX_fid) ORG_dscl          (:,:,:)
     read(EX_fid) ORG_dscl_pl       (:,:,:)
     read(EX_fid) ORG_scl           (:,:,:)
     read(EX_fid) ORG_scl_pl        (:,:,:)
     read(EX_fid) ORG_kh            (:,:,:)
     read(EX_fid) ORG_kh_pl         (:,:,:)
     read(EX_fid) ORG_cinterp_TN    (:,:,:,:)
     read(EX_fid) ORG_cinterp_TN_pl (:,:,:,:)
     read(EX_fid) ORG_cinterp_HN    (:,:,:,:)
     read(EX_fid) ORG_cinterp_HN_pl (:,:,  :)
     read(EX_fid) ORG_cinterp_TRA   (:,:,:  )
     read(EX_fid) ORG_cinterp_TRA_pl(:,:    )
     read(EX_fid) ORG_cinterp_PRA   (:,:    )
     read(EX_fid) ORG_cinterp_PRA_pl(:,:    )
  close(EX_fid)

  dscl          (:,:,:)   = real( ORG_dscl          (:,:,:)  , kind=RP )
  dscl_pl       (:,:,:)   = real( ORG_dscl_pl       (:,:,:)  , kind=RP )
  scl           (:,:,:)   = real( ORG_scl           (:,:,:)  , kind=RP )
  scl_pl        (:,:,:)   = real( ORG_scl_pl        (:,:,:)  , kind=RP )
  kh            (:,:,:)   = real( ORG_kh            (:,:,:)  , kind=RP )
  kh_pl         (:,:,:)   = real( ORG_kh_pl         (:,:,:)  , kind=RP )
  cinterp_TN    (:,:,:,:) = real( ORG_cinterp_TN    (:,:,:,:), kind=RP )
  cinterp_TN_pl (:,:,:,:) = real( ORG_cinterp_TN_pl (:,:,:,:), kind=RP )
  cinterp_HN    (:,:,:,:) = real( ORG_cinterp_HN    (:,:,:,:), kind=RP )
  cinterp_HN_pl (:,:,  :) = real( ORG_cinterp_HN_pl (:,:,  :), kind=RP )
  cinterp_TRA   (:,:,:  ) = real( ORG_cinterp_TRA   (:,:,:  ), kind=RP )
  cinterp_TRA_pl(:,:    ) = real( ORG_cinterp_TRA_pl(:,:    ), kind=RP )
  cinterp_PRA   (:,:    ) = real( ORG_cinterp_PRA   (:,:    ), kind=RP )
  cinterp_PRA_pl(:,:    ) = real( ORG_cinterp_PRA_pl(:,:    ), kind=RP )

  deallocate( ORG_dscl           )
  deallocate( ORG_dscl_pl        )
  deallocate( ORG_scl            )
  deallocate( ORG_scl_pl         )
  deallocate( ORG_kh             )
  deallocate( ORG_kh_pl          )
  deallocate( ORG_cinterp_TN     )
  deallocate( ORG_cinterp_TN_pl  )
  deallocate( ORG_cinterp_HN     )
  deallocate( ORG_cinterp_HN_pl  )
  deallocate( ORG_cinterp_TRA    )
  deallocate( ORG_cinterp_TRA_pl )
  deallocate( ORG_cinterp_PRA    )
  deallocate( ORG_cinterp_PRA_pl )

  !###############################################################################

  !---< operator module setup >---
  call OPRT_diffusion_putcoef( cinterp_TN , cinterp_TN_pl , & ![IN]
                               cinterp_HN , cinterp_HN_pl , & ![IN]
                               cinterp_TRA, cinterp_TRA_pl, & ![IN]
                               cinterp_PRA, cinterp_PRA_pl  ) ![IN]

  check_dscl   (:,:,:) = dscl   (:,:,:)
  check_dscl_pl(:,:,:) = dscl_pl(:,:,:)

  write(*,*) "*** Finish initialize"

  !###############################################################################

  write(*,*) "*** Start  kernel"
  call DEBUG_rapstart('DC_diffusion_kernel')

  do iteration = 1, SET_iteration



     call OPRT_diffusion( dscl(:,:,:), dscl_pl(:,:,:), & !--- OUT
                          scl (:,:,:), scl_pl (:,:,:), & !--- IN
                          kh  (:,:,:), kh_pl  (:,:,:)  ) !--- IN



     write(ADM_LOG_FID,*) '### Input ###'
     call DEBUG_valuecheck( 'check_dscl   ', check_dscl   (:,:,:) )
     call DEBUG_valuecheck( 'check_dscl_pl', check_dscl_pl(:,:,:) )
     call DEBUG_valuecheck( 'scl          ', scl          (:,:,:) )
     call DEBUG_valuecheck( 'scl_pl       ', scl_pl       (:,:,:) )
     call DEBUG_valuecheck( 'kh           ', kh           (:,:,:) )
     call DEBUG_valuecheck( 'kh_pl        ', kh_pl        (:,:,:) )
     write(ADM_LOG_FID,*) '### Output ###'
     call DEBUG_valuecheck( 'dscl         ', dscl         (:,:,:) )
     call DEBUG_valuecheck( 'dscl_pl      ', dscl_pl      (:,:,:) )
  enddo

  write(ADM_LOG_FID,*) '### Varidation by diff ###'
  check_dscl   (:,:,:) = check_dscl   (:,:,:) - dscl   (:,:,:)
  check_dscl_pl(:,:,:) = check_dscl_pl(:,:,:) - dscl_pl(:,:,:)
  call DEBUG_valuecheck( 'check_dscl   ', check_dscl   (:,:,:) )
  call DEBUG_valuecheck( 'check_dscl_pl', check_dscl_pl(:,:,:) )

  call DEBUG_rapend('DC_diffusion_kernel')
  write(*,*) "*** Finish kernel"

  !###############################################################################

  call DEBUG_rapreport

  stop
end program dckernel_diffusion
!-------------------------------------------------------------------------------
