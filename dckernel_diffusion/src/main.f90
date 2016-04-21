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

  real(RP), allocatable :: dscl   (:,:,:)
  real(RP), allocatable :: dscl_pl(:,:,:)
  real(RP), allocatable :: scl    (:,:,:)
  real(RP), allocatable :: scl_pl (:,:,:)
  real(RP), allocatable :: kh     (:,:,:)
  real(RP), allocatable :: kh_pl  (:,:,:)

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

  allocate( dscl          (ADM_gall   ,ADM_kall,ADM_lall   ) )
  allocate( dscl_pl       (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
  allocate( scl           (ADM_gall   ,ADM_kall,ADM_lall   ) )
  allocate( scl_pl        (ADM_gall_pl,ADM_kall,ADM_lall_pl) )
  allocate( kh            (ADM_gall   ,ADM_kall,ADM_lall   ) )
  allocate( kh_pl         (ADM_gall_pl,ADM_kall,ADM_lall_pl) )

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

  call MISC_make_idstr(EX_fname,'snapshot.dc_diffusion','pe',SET_prc_me)
  EX_fid = MISC_get_available_fid()
  open( unit   = EX_fid,         &
        file   = trim(EX_fname), &
        form   = 'unformatted',  &
        access = 'sequential',   &
        status = 'old'           )

     read(EX_fid) dscl          (:,:,:)
     read(EX_fid) dscl_pl       (:,:,:)
     read(EX_fid) scl           (:,:,:)
     read(EX_fid) scl_pl        (:,:,:)
     read(EX_fid) kh            (:,:,:)
     read(EX_fid) kh_pl         (:,:,:)
     read(EX_fid) cinterp_TN    (:,:,:,:)
     read(EX_fid) cinterp_TN_pl (:,:,:,:)
     read(EX_fid) cinterp_HN    (:,:,:,:)
     read(EX_fid) cinterp_HN_pl (:,:,  :)
     read(EX_fid) cinterp_TRA   (:,:,:  )
     read(EX_fid) cinterp_TRA_pl(:,:    )
     read(EX_fid) cinterp_PRA   (:,:    )
     read(EX_fid) cinterp_PRA_pl(:,:    )
  close(EX_fid)

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
     EX_item =       'check_dscl   '
     EX_max  = maxval(check_dscl   (:,:,:))
     EX_min  = minval(check_dscl   (:,:,:))
     EX_sum  = sum   (check_dscl   (:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
     EX_item =       'check_dscl_pl'
     EX_max  = maxval(check_dscl_pl(:,:,:))
     EX_min  = minval(check_dscl_pl(:,:,:))
     EX_sum  = sum   (check_dscl_pl(:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
     EX_item =       'scl    '
     EX_max  = maxval(scl    (:,:,:))
     EX_min  = minval(scl    (:,:,:))
     EX_sum  = sum   (scl    (:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
     EX_item =       'scl_pl '
     EX_max  = maxval(scl_pl (:,:,:))
     EX_min  = minval(scl_pl (:,:,:))
     EX_sum  = sum   (scl_pl (:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
     EX_item =       'kh     '
     EX_max  = maxval(kh     (:,:,:))
     EX_min  = minval(kh     (:,:,:))
     EX_sum  = sum   (kh     (:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
     EX_item =       'kh_pl  '
     EX_max  = maxval(kh_pl  (:,:,:))
     EX_min  = minval(kh_pl  (:,:,:))
     EX_sum  = sum   (kh_pl  (:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum

     write(ADM_LOG_FID,*) '### Output ###'
     EX_item =       'dscl   '
     EX_max  = maxval(dscl   (:,:,:))
     EX_min  = minval(dscl   (:,:,:))
     EX_sum  = sum   (dscl   (:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
     EX_item =       'dscl_pl'
     EX_max  = maxval(dscl_pl(:,:,:))
     EX_min  = minval(dscl_pl(:,:,:))
     EX_sum  = sum   (dscl_pl(:,:,:))
     write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
  enddo

  write(ADM_LOG_FID,*) '### Varidation by diff ###'
  check_dscl   (:,:,:) = check_dscl   (:,:,:) - dscl   (:,:,:)
  check_dscl_pl(:,:,:) = check_dscl_pl(:,:,:) - dscl_pl(:,:,:)

  EX_item    =       'check_dscl   '
  EX_max = maxval(check_dscl   (:,:,:))
  EX_min = minval(check_dscl   (:,:,:))
  EX_sum = sum   (check_dscl   (:,:,:))
  write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum
  EX_item    =       'check_dscl_pl'
  EX_max = maxval(check_dscl_pl(:,:,:))
  EX_min = minval(check_dscl_pl(:,:,:))
  EX_sum = sum   (check_dscl_pl(:,:,:))
  write(ADM_LOG_FID,'(1x,A,A16,3(A,ES24.16))') '+check[',EX_item,'] max=',EX_max,',min=',EX_min,',sum=',EX_sum

  call DEBUG_rapend('DC_diffusion_kernel')
  write(*,*) "*** Finish kernel"

  call DEBUG_rapreport

  stop
end program dckernel_diffusion
!-------------------------------------------------------------------------------
