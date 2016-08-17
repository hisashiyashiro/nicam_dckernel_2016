!-------------------------------------------------------------------------------
!
!+  Program dynamics kernel driver (setup operator)
!
!-------------------------------------------------------------------------------
program dckernel_setup
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_misc
  use mod_gmtr, only: &
     GMTR_p_setup, &
     GMTR_t_setup, &
     GMTR_a_setup
  use mod_oprt, only: &
     OPRT_divergence_setup, &
     OPRT_rotation_setup,   &
     OPRT_gradient_setup,   &
     OPRT_laplacian_setup,  &
     OPRT_diffusion_setup
  !-----------------------------------------------------------------------------
  implicit none

  real(DP), allocatable :: ORG_GRD_x            (:,:,:,:)
  real(DP), allocatable :: ORG_GRD_x_pl         (:,:,:,:)
  real(DP), allocatable :: ORG_GRD_xt           (:,:,:,:,:)
  real(DP), allocatable :: ORG_GRD_xt_pl        (:,:,:,:)
  real(DP), allocatable :: ORG_GRD_s            (:,:,:,:)
  real(DP), allocatable :: ORG_GRD_s_pl         (:,:,:,:)
  real(DP), allocatable :: ORG_GMTR_p           (:,:,:,:)
  real(DP), allocatable :: ORG_GMTR_p_pl        (:,:,:,:)
  real(DP), allocatable :: ORG_GMTR_t           (:,:,:,:,:)
  real(DP), allocatable :: ORG_GMTR_t_pl        (:,:,:,:)
  real(DP), allocatable :: ORG_GMTR_a           (:,:,:,:,:)
  real(DP), allocatable :: ORG_GMTR_a_pl        (:,:,:,:)
  real(DP), allocatable :: ORG_OPRT_coef_div    (:,:,:,:,:)
  real(DP), allocatable :: ORG_OPRT_coef_div_pl (:,:,:)
  real(DP), allocatable :: ORG_OPRT_coef_rot    (:,:,:,:,:)
  real(DP), allocatable :: ORG_OPRT_coef_rot_pl (:,:,:)
  real(DP), allocatable :: ORG_OPRT_coef_grad   (:,:,:,:,:)
  real(DP), allocatable :: ORG_OPRT_coef_grad_pl(:,:,:)
  real(DP), allocatable :: ORG_OPRT_coef_lap    (:,:,:,:)
  real(DP), allocatable :: ORG_OPRT_coef_lap_pl (:,:)
  real(DP), allocatable :: ORG_OPRT_coef_intp   (:,:,:,:,:,:)
  real(DP), allocatable :: ORG_OPRT_coef_intp_pl(:,:,:,:)
  real(DP), allocatable :: ORG_OPRT_coef_diff   (:,:,:,:,:)
  real(DP), allocatable :: ORG_OPRT_coef_diff_pl(:,:,:)

  real(RP), allocatable :: GRD_x            (:,:,:,:)
  real(RP), allocatable :: GRD_x_pl         (:,:,:,:)
  real(RP), allocatable :: GRD_xt           (:,:,:,:,:)
  real(RP), allocatable :: GRD_xt_pl        (:,:,:,:)
  real(RP), allocatable :: GRD_s            (:,:,:,:)
  real(RP), allocatable :: GRD_s_pl         (:,:,:,:)
  real(RP), allocatable :: GMTR_p           (:,:,:,:)
  real(RP), allocatable :: GMTR_p_pl        (:,:,:,:)
  real(RP), allocatable :: GMTR_t           (:,:,:,:,:)
  real(RP), allocatable :: GMTR_t_pl        (:,:,:,:)
  real(RP), allocatable :: GMTR_a           (:,:,:,:,:)
  real(RP), allocatable :: GMTR_a_pl        (:,:,:,:)
  real(RP), allocatable :: OPRT_coef_div    (:,:,:,:,:)
  real(RP), allocatable :: OPRT_coef_div_pl (:,:,:)
  real(RP), allocatable :: OPRT_coef_rot    (:,:,:,:,:)
  real(RP), allocatable :: OPRT_coef_rot_pl (:,:,:)
  real(RP), allocatable :: OPRT_coef_grad   (:,:,:,:,:)
  real(RP), allocatable :: OPRT_coef_grad_pl(:,:,:)
  real(RP), allocatable :: OPRT_coef_lap    (:,:,:,:)
  real(RP), allocatable :: OPRT_coef_lap_pl (:,:)
  real(RP), allocatable :: OPRT_coef_intp   (:,:,:,:,:,:)
  real(RP), allocatable :: OPRT_coef_intp_pl(:,:,:,:)
  real(RP), allocatable :: OPRT_coef_diff   (:,:,:,:,:)
  real(RP), allocatable :: OPRT_coef_diff_pl(:,:,:)

  real(RP), allocatable :: check_GMTR_p           (:,:,:,:)
  real(RP), allocatable :: check_GMTR_p_pl        (:,:,:,:)
  real(RP), allocatable :: check_GMTR_t           (:,:,:,:,:)
  real(RP), allocatable :: check_GMTR_t_pl        (:,:,:,:)
  real(RP), allocatable :: check_GMTR_a           (:,:,:,:,:)
  real(RP), allocatable :: check_GMTR_a_pl        (:,:,:,:)
  real(RP), allocatable :: check_OPRT_coef_div    (:,:,:,:,:)
  real(RP), allocatable :: check_OPRT_coef_div_pl (:,:,:)
  real(RP), allocatable :: check_OPRT_coef_rot    (:,:,:,:,:)
  real(RP), allocatable :: check_OPRT_coef_rot_pl (:,:,:)
  real(RP), allocatable :: check_OPRT_coef_grad   (:,:,:,:,:)
  real(RP), allocatable :: check_OPRT_coef_grad_pl(:,:,:)
  real(RP), allocatable :: check_OPRT_coef_lap    (:,:,:,:)
  real(RP), allocatable :: check_OPRT_coef_lap_pl (:,:)
  real(RP), allocatable :: check_OPRT_coef_intp   (:,:,:,:,:,:)
  real(RP), allocatable :: check_OPRT_coef_intp_pl(:,:,:,:)
  real(RP), allocatable :: check_OPRT_coef_diff   (:,:,:,:,:)
  real(RP), allocatable :: check_OPRT_coef_diff_pl(:,:,:)

  real(RP), allocatable :: IxJ_GMTR_p   (:,:,:,:,:)
  real(RP), allocatable :: IxJ_GMTR_p_pl(:,:,:,:)
  real(RP), allocatable :: IxJ_GMTR_t   (:,:,:,:,:,:)
  real(RP), allocatable :: IxJ_GMTR_t_pl(:,:,:,:)
  real(RP), allocatable :: IxJ_GMTR_a   (:,:,:,:,:,:)
  real(RP), allocatable :: IxJ_GMTR_a_pl(:,:,:,:)

  integer :: iteration
  integer :: i, j, g
  !=============================================================================

  write(*,*) "[KERNEL] dckernel_setup"

  !###############################################################################

  write(*,*) "*** Start  initialize"

  allocate( GRD_x            (ADM_gall   ,K0,ADM_lall   ,      ADM_nxyz)       )
  allocate( GRD_x_pl         (ADM_gall_pl,K0,ADM_lall_pl,      ADM_nxyz)       )
  allocate( GRD_xt           (ADM_gall   ,K0,ADM_lall   ,TI:TJ,ADM_nxyz)       )
  allocate( GRD_xt_pl        (ADM_gall_pl,K0,ADM_lall_pl,      ADM_nxyz)       )
  allocate( GRD_s            (ADM_gall   ,K0,ADM_lall   ,      2)              )
  allocate( GRD_s_pl         (ADM_gall_pl,K0,ADM_lall_pl,      2)              )
  allocate( GMTR_p           (ADM_gall   ,K0,ADM_lall   ,      GMTR_p_nmax   ) )
  allocate( GMTR_p_pl        (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_p_nmax   ) )
  allocate( GMTR_t           (ADM_gall   ,K0,ADM_lall   ,TI:TJ,GMTR_t_nmax   ) )
  allocate( GMTR_t_pl        (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_t_nmax   ) )
  allocate( GMTR_a           (ADM_gall   ,K0,ADM_lall   ,AI:AJ,GMTR_a_nmax   ) )
  allocate( GMTR_a_pl        (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_a_nmax_pl) )
  allocate( OPRT_coef_div    (ADM_iall,ADM_jall,0:6        ,ADM_nxyz,      ADM_lall   ) )
  allocate( OPRT_coef_div_pl (                  0:ADM_vlink,ADM_nxyz,      ADM_lall_pl) )
  allocate( OPRT_coef_rot    (ADM_iall,ADM_jall,0:6        ,ADM_nxyz,      ADM_lall   ) )
  allocate( OPRT_coef_rot_pl (                  0:ADM_vlink,ADM_nxyz,      ADM_lall_pl) )
  allocate( OPRT_coef_grad   (ADM_iall,ADM_jall,0:6        ,ADM_nxyz,      ADM_lall   ) )
  allocate( OPRT_coef_grad_pl(                  0:ADM_vlink,ADM_nxyz,      ADM_lall_pl) )
  allocate( OPRT_coef_lap    (ADM_iall,ADM_jall,0:6        ,               ADM_lall   ) )
  allocate( OPRT_coef_lap_pl (                  0:ADM_vlink,               ADM_lall_pl) )
  allocate( OPRT_coef_intp   (ADM_iall,ADM_jall,1:3        ,ADM_nxyz,TI:TJ,ADM_lall   ) )
  allocate( OPRT_coef_intp_pl(ADM_gall_pl      ,1:3        ,ADM_nxyz,      ADM_lall_pl) )
  allocate( OPRT_coef_diff   (ADM_iall,ADM_jall,1:6        ,ADM_nxyz,      ADM_lall   ) )
  allocate( OPRT_coef_diff_pl(                  1:ADM_vlink,ADM_nxyz,      ADM_lall_pl) )
  ! Todo : first touch with considering NUMA
  GRD_x            (:,:,:,:)     = 0.0_RP
  GRD_x_pl         (:,:,:,:)     = 0.0_RP
  GRD_xt           (:,:,:,:,:)   = 0.0_RP
  GRD_xt_pl        (:,:,:,:)     = 0.0_RP
  GRD_s            (:,:,:,:)     = 0.0_RP
  GRD_s_pl         (:,:,:,:)     = 0.0_RP
  GMTR_p           (:,:,:,:)     = 0.0_RP
  GMTR_p_pl        (:,:,:,:)     = 0.0_RP
  GMTR_t           (:,:,:,:,:)   = 0.0_RP
  GMTR_t_pl        (:,:,:,:)     = 0.0_RP
  GMTR_a           (:,:,:,:,:)   = 0.0_RP
  GMTR_a_pl        (:,:,:,:)     = 0.0_RP
  OPRT_coef_div    (:,:,:,:,:)   = 0.0_RP
  OPRT_coef_div_pl (:,:,:)       = 0.0_RP
  OPRT_coef_rot    (:,:,:,:,:)   = 0.0_RP
  OPRT_coef_rot_pl (:,:,:)       = 0.0_RP
  OPRT_coef_grad   (:,:,:,:,:)   = 0.0_RP
  OPRT_coef_grad_pl(:,:,:)       = 0.0_RP
  OPRT_coef_lap    (:,:,:,:)     = 0.0_RP
  OPRT_coef_lap_pl (:,:)         = 0.0_RP
  OPRT_coef_intp   (:,:,:,:,:,:) = 0.0_RP
  OPRT_coef_intp_pl(:,:,:,:)     = 0.0_RP
  OPRT_coef_diff   (:,:,:,:,:)   = 0.0_RP
  OPRT_coef_diff_pl(:,:,:)       = 0.0_RP

  allocate( check_GMTR_p           (ADM_gall   ,K0,ADM_lall   ,      GMTR_p_nmax   ) )
  allocate( check_GMTR_p_pl        (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_p_nmax   ) )
  allocate( check_GMTR_t           (ADM_gall   ,K0,ADM_lall   ,TI:TJ,GMTR_t_nmax   ) )
  allocate( check_GMTR_t_pl        (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_t_nmax   ) )
  allocate( check_GMTR_a           (ADM_gall   ,K0,ADM_lall   ,AI:AJ,GMTR_a_nmax   ) )
  allocate( check_GMTR_a_pl        (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_a_nmax_pl) )
  allocate( check_OPRT_coef_div    (ADM_iall,ADM_jall,0:6        ,ADM_nxyz,      ADM_lall   ) )
  allocate( check_OPRT_coef_div_pl (                  0:ADM_vlink,ADM_nxyz,      ADM_lall_pl) )
  allocate( check_OPRT_coef_rot    (ADM_iall,ADM_jall,0:6        ,ADM_nxyz,      ADM_lall   ) )
  allocate( check_OPRT_coef_rot_pl (                  0:ADM_vlink,ADM_nxyz,      ADM_lall_pl) )
  allocate( check_OPRT_coef_grad   (ADM_iall,ADM_jall,0:6        ,ADM_nxyz,      ADM_lall   ) )
  allocate( check_OPRT_coef_grad_pl(                  0:ADM_vlink,ADM_nxyz,      ADM_lall_pl) )
  allocate( check_OPRT_coef_lap    (ADM_iall,ADM_jall,0:6        ,               ADM_lall   ) )
  allocate( check_OPRT_coef_lap_pl (                  0:ADM_vlink,               ADM_lall_pl) )
  allocate( check_OPRT_coef_intp   (ADM_iall,ADM_jall,1:3        ,ADM_nxyz,TI:TJ,ADM_lall   ) )
  allocate( check_OPRT_coef_intp_pl(ADM_gall_pl      ,1:3        ,ADM_nxyz,      ADM_lall_pl) )
  allocate( check_OPRT_coef_diff   (ADM_iall,ADM_jall,1:6        ,ADM_nxyz,      ADM_lall   ) )
  allocate( check_OPRT_coef_diff_pl(                  1:ADM_vlink,ADM_nxyz,      ADM_lall_pl) )
  check_GMTR_p           (:,:,:,:)     = 0.0_RP
  check_GMTR_p_pl        (:,:,:,:)     = 0.0_RP
  check_GMTR_t           (:,:,:,:,:)   = 0.0_RP
  check_GMTR_t_pl        (:,:,:,:)     = 0.0_RP
  check_GMTR_a           (:,:,:,:,:)   = 0.0_RP
  check_GMTR_a_pl        (:,:,:,:)     = 0.0_RP
  check_OPRT_coef_div    (:,:,:,:,:)   = 0.0_RP
  check_OPRT_coef_div_pl (:,:,:)       = 0.0_RP
  check_OPRT_coef_rot    (:,:,:,:,:)   = 0.0_RP
  check_OPRT_coef_rot_pl (:,:,:)       = 0.0_RP
  check_OPRT_coef_grad   (:,:,:,:,:)   = 0.0_RP
  check_OPRT_coef_grad_pl(:,:,:)       = 0.0_RP
  check_OPRT_coef_lap    (:,:,:,:)     = 0.0_RP
  check_OPRT_coef_lap_pl (:,:)         = 0.0_RP
  check_OPRT_coef_intp   (:,:,:,:,:,:) = 0.0_RP
  check_OPRT_coef_intp_pl(:,:,:,:)     = 0.0_RP
  check_OPRT_coef_diff   (:,:,:,:,:)   = 0.0_RP
  check_OPRT_coef_diff_pl(:,:,:)       = 0.0_RP

  !###############################################################################

  !---< read input data >---)
  allocate( ORG_GRD_x            (ADM_gall   ,K0,ADM_lall   ,      ADM_nxyz)       )
  allocate( ORG_GRD_x_pl         (ADM_gall_pl,K0,ADM_lall_pl,      ADM_nxyz)       )
  allocate( ORG_GRD_xt           (ADM_gall   ,K0,ADM_lall   ,TI:TJ,ADM_nxyz)       )
  allocate( ORG_GRD_xt_pl        (ADM_gall_pl,K0,ADM_lall_pl,      ADM_nxyz)       )
  allocate( ORG_GRD_s            (ADM_gall   ,K0,ADM_lall   ,      2)              )
  allocate( ORG_GRD_s_pl         (ADM_gall_pl,K0,ADM_lall_pl,      2)              )
  allocate( ORG_GMTR_p           (ADM_gall   ,K0,ADM_lall   ,      GMTR_p_nmax   ) )
  allocate( ORG_GMTR_p_pl        (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_p_nmax   ) )
  allocate( ORG_GMTR_t           (ADM_gall   ,K0,ADM_lall   ,TI:TJ,GMTR_t_nmax   ) )
  allocate( ORG_GMTR_t_pl        (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_t_nmax   ) )
  allocate( ORG_GMTR_a           (ADM_gall   ,K0,ADM_lall   ,AI:AJ,GMTR_a_nmax   ) )
  allocate( ORG_GMTR_a_pl        (ADM_gall_pl,K0,ADM_lall_pl,      GMTR_a_nmax_pl) )
  allocate( ORG_OPRT_coef_div    (ADM_iall,ADM_jall,0:6        ,ADM_nxyz,      ADM_lall   ) )
  allocate( ORG_OPRT_coef_div_pl (                  0:ADM_vlink,ADM_nxyz,      ADM_lall_pl) )
  allocate( ORG_OPRT_coef_rot    (ADM_iall,ADM_jall,0:6        ,ADM_nxyz,      ADM_lall   ) )
  allocate( ORG_OPRT_coef_rot_pl (                  0:ADM_vlink,ADM_nxyz,      ADM_lall_pl) )
  allocate( ORG_OPRT_coef_grad   (ADM_iall,ADM_jall,0:6        ,ADM_nxyz,      ADM_lall   ) )
  allocate( ORG_OPRT_coef_grad_pl(                  0:ADM_vlink,ADM_nxyz,      ADM_lall_pl) )
  allocate( ORG_OPRT_coef_lap    (ADM_iall,ADM_jall,0:6        ,               ADM_lall   ) )
  allocate( ORG_OPRT_coef_lap_pl (                  0:ADM_vlink,               ADM_lall_pl) )
  allocate( ORG_OPRT_coef_intp   (ADM_iall,ADM_jall,1:3        ,ADM_nxyz,TI:TJ,ADM_lall   ) )
  allocate( ORG_OPRT_coef_intp_pl(ADM_gall_pl      ,1:3        ,ADM_nxyz,      ADM_lall_pl) )
  allocate( ORG_OPRT_coef_diff   (ADM_iall,ADM_jall,1:6        ,ADM_nxyz,      ADM_lall   ) )
  allocate( ORG_OPRT_coef_diff_pl(                  1:ADM_vlink,ADM_nxyz,      ADM_lall_pl) )

  call dumpio_syscheck
  call dumpio_mk_fname(EX_fname,'snapshot.dc_setup','pe',SET_prc_me-1,6)
  call dumpio_fopen(EX_fid,EX_fname,IO_FREAD)

  call dumpio_read_data( EX_fid, ADM_gall   *K0*ADM_lall   *  ADM_nxyz               , ORG_GRD_x            (:,:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall_pl*K0*ADM_lall_pl*  ADM_nxyz               , ORG_GRD_x_pl         (:,:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall   *K0*ADM_lall   *2*ADM_nxyz               , ORG_GRD_xt           (:,:,:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_gall_pl*K0*ADM_lall_pl*  ADM_nxyz               , ORG_GRD_xt_pl        (:,:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall   *K0*ADM_lall   *  2                      , ORG_GRD_s            (:,:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall_pl*K0*ADM_lall_pl*  2                      , ORG_GRD_s_pl         (:,:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall   *K0*ADM_lall   *  GMTR_p_nmax            , ORG_GMTR_p           (:,:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall_pl*K0*ADM_lall_pl*  GMTR_p_nmax            , ORG_GMTR_p_pl        (:,:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall   *K0*ADM_lall   *2*GMTR_t_nmax            , ORG_GMTR_t           (:,:,:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_gall_pl*K0*ADM_lall_pl*  GMTR_t_nmax            , ORG_GMTR_t_pl        (:,:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_gall   *K0*ADM_lall   *3*GMTR_a_nmax            , ORG_GMTR_a           (:,:,:,:,:)   )
  call dumpio_read_data( EX_fid, ADM_gall_pl*K0*ADM_lall_pl*  GMTR_a_nmax_pl         , ORG_GMTR_a_pl        (:,:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_iall*ADM_jall*7          *ADM_nxyz*  ADM_lall   , ORG_OPRT_coef_div    (:,:,:,:,:)   )
  call dumpio_read_data( EX_fid,                 (ADM_vlink+1)*ADM_nxyz*  ADM_lall_pl, ORG_OPRT_coef_div_pl (:,:,:)       )
  call dumpio_read_data( EX_fid, ADM_iall*ADM_jall*7          *ADM_nxyz*  ADM_lall   , ORG_OPRT_coef_rot    (:,:,:,:,:)   )
  call dumpio_read_data( EX_fid,                 (ADM_vlink+1)*ADM_nxyz*  ADM_lall_pl, ORG_OPRT_coef_rot_pl (:,:,:)       )
  call dumpio_read_data( EX_fid, ADM_iall*ADM_jall*7          *ADM_nxyz*  ADM_lall   , ORG_OPRT_coef_grad   (:,:,:,:,:)   )
  call dumpio_read_data( EX_fid,                 (ADM_vlink+1)*ADM_nxyz*  ADM_lall_pl, ORG_OPRT_coef_grad_pl(:,:,:)       )
  call dumpio_read_data( EX_fid, ADM_iall*ADM_jall*7          *           ADM_lall   , ORG_OPRT_coef_lap    (:,:,:,:)     )
  call dumpio_read_data( EX_fid,                 (ADM_vlink+1)*           ADM_lall_pl, ORG_OPRT_coef_lap_pl (:,:)         )
  call dumpio_read_data( EX_fid, ADM_iall*ADM_jall*3          *ADM_nxyz*2*ADM_lall   , ORG_OPRT_coef_intp   (:,:,:,:,:,:) )
  call dumpio_read_data( EX_fid, ADM_gall_pl      *3          *ADM_nxyz*  ADM_lall_pl, ORG_OPRT_coef_intp_pl(:,:,:,:)     )
  call dumpio_read_data( EX_fid, ADM_iall*ADM_jall*6          *ADM_nxyz*  ADM_lall   , ORG_OPRT_coef_diff   (:,:,:,:,:)   )
  call dumpio_read_data( EX_fid,                   ADM_vlink  *ADM_nxyz*  ADM_lall_pl, ORG_OPRT_coef_diff_pl(:,:,:)       )

  call dumpio_fclose(EX_fid)

  GRD_x            (:,:,:,:)     = real( ORG_GRD_x            (:,:,:,:)    , kind=RP )
  GRD_x_pl         (:,:,:,:)     = real( ORG_GRD_x_pl         (:,:,:,:)    , kind=RP )
  GRD_xt           (:,:,:,:,:)   = real( ORG_GRD_xt           (:,:,:,:,:)  , kind=RP )
  GRD_xt_pl        (:,:,:,:)     = real( ORG_GRD_xt_pl        (:,:,:,:)    , kind=RP )
  GRD_s            (:,:,:,:)     = real( ORG_GRD_s            (:,:,:,:)    , kind=RP )
  GRD_s_pl         (:,:,:,:)     = real( ORG_GRD_s_pl         (:,:,:,:)    , kind=RP )
  GMTR_p           (:,:,:,:)     = real( ORG_GMTR_p           (:,:,:,:)    , kind=RP )
  GMTR_p_pl        (:,:,:,:)     = real( ORG_GMTR_p_pl        (:,:,:,:)    , kind=RP )
  GMTR_t           (:,:,:,:,:)   = real( ORG_GMTR_t           (:,:,:,:,:)  , kind=RP )
  GMTR_t_pl        (:,:,:,:)     = real( ORG_GMTR_t_pl        (:,:,:,:)    , kind=RP )
  GMTR_a           (:,:,:,:,:)   = real( ORG_GMTR_a           (:,:,:,:,:)  , kind=RP )
  GMTR_a_pl        (:,:,:,:)     = real( ORG_GMTR_a_pl        (:,:,:,:)    , kind=RP )
  OPRT_coef_div    (:,:,:,:,:)   = real( ORG_OPRT_coef_div    (:,:,:,:,:)  , kind=RP )
  OPRT_coef_div_pl (:,:,:)       = real( ORG_OPRT_coef_div_pl (:,:,:)      , kind=RP )
  OPRT_coef_rot    (:,:,:,:,:)   = real( ORG_OPRT_coef_rot    (:,:,:,:,:)  , kind=RP )
  OPRT_coef_rot_pl (:,:,:)       = real( ORG_OPRT_coef_rot_pl (:,:,:)      , kind=RP )
  OPRT_coef_grad   (:,:,:,:,:)   = real( ORG_OPRT_coef_grad   (:,:,:,:,:)  , kind=RP )
  OPRT_coef_grad_pl(:,:,:)       = real( ORG_OPRT_coef_grad_pl(:,:,:)      , kind=RP )
  OPRT_coef_lap    (:,:,:,:)     = real( ORG_OPRT_coef_lap    (:,:,:,:)    , kind=RP )
  OPRT_coef_lap_pl (:,:)         = real( ORG_OPRT_coef_lap_pl (:,:)        , kind=RP )
  OPRT_coef_intp   (:,:,:,:,:,:) = real( ORG_OPRT_coef_intp   (:,:,:,:,:,:), kind=RP )
  OPRT_coef_intp_pl(:,:,:,:)     = real( ORG_OPRT_coef_intp_pl(:,:,:,:)    , kind=RP )
  OPRT_coef_diff   (:,:,:,:,:)   = real( ORG_OPRT_coef_diff   (:,:,:,:,:)  , kind=RP )
  OPRT_coef_diff_pl(:,:,:)       = real( ORG_OPRT_coef_diff_pl(:,:,:)      , kind=RP )

  deallocate( ORG_GRD_x             )
  deallocate( ORG_GRD_x_pl          )
  deallocate( ORG_GRD_xt            )
  deallocate( ORG_GRD_xt_pl         )
  deallocate( ORG_GRD_s             )
  deallocate( ORG_GRD_s_pl          )
  deallocate( ORG_GMTR_p            )
  deallocate( ORG_GMTR_p_pl         )
  deallocate( ORG_GMTR_t            )
  deallocate( ORG_GMTR_t_pl         )
  deallocate( ORG_GMTR_a            )
  deallocate( ORG_GMTR_a_pl         )
  deallocate( ORG_OPRT_coef_div     )
  deallocate( ORG_OPRT_coef_div_pl  )
  deallocate( ORG_OPRT_coef_rot     )
  deallocate( ORG_OPRT_coef_rot_pl  )
  deallocate( ORG_OPRT_coef_grad    )
  deallocate( ORG_OPRT_coef_grad_pl )
  deallocate( ORG_OPRT_coef_lap     )
  deallocate( ORG_OPRT_coef_lap_pl  )
  deallocate( ORG_OPRT_coef_intp    )
  deallocate( ORG_OPRT_coef_intp_pl )
  deallocate( ORG_OPRT_coef_diff    )
  deallocate( ORG_OPRT_coef_diff_pl )

  !###############################################################################

  call CONST_setup ! set PI and Epsilon

  check_GMTR_p           (:,:,:,:)     = GMTR_p           (:,:,:,:)
  check_GMTR_p_pl        (:,:,:,:)     = GMTR_p_pl        (:,:,:,:)
  check_GMTR_t           (:,:,:,:,:)   = GMTR_t           (:,:,:,:,:)
  check_GMTR_t_pl        (:,:,:,:)     = GMTR_t_pl        (:,:,:,:)
  check_GMTR_a           (:,:,:,:,:)   = GMTR_a           (:,:,:,:,:)
  check_GMTR_a_pl        (:,:,:,:)     = GMTR_a_pl        (:,:,:,:)
  check_OPRT_coef_div    (:,:,:,:,:)   = OPRT_coef_div    (:,:,:,:,:)
  check_OPRT_coef_div_pl (:,:,:)       = OPRT_coef_div_pl (:,:,:)
  check_OPRT_coef_rot    (:,:,:,:,:)   = OPRT_coef_rot    (:,:,:,:,:)
  check_OPRT_coef_rot_pl (:,:,:)       = OPRT_coef_rot_pl (:,:,:)
  check_OPRT_coef_grad   (:,:,:,:,:)   = OPRT_coef_grad   (:,:,:,:,:)
  check_OPRT_coef_grad_pl(:,:,:)       = OPRT_coef_grad_pl(:,:,:)
  check_OPRT_coef_lap    (:,:,:,:)     = OPRT_coef_lap    (:,:,:,:)
  check_OPRT_coef_lap_pl (:,:)         = OPRT_coef_lap_pl (:,:)
  check_OPRT_coef_intp   (:,:,:,:,:,:) = OPRT_coef_intp   (:,:,:,:,:,:)
  check_OPRT_coef_intp_pl(:,:,:,:)     = OPRT_coef_intp_pl(:,:,:,:)
  check_OPRT_coef_diff   (:,:,:,:,:)   = OPRT_coef_diff   (:,:,:,:,:)
  check_OPRT_coef_diff_pl(:,:,:)       = OPRT_coef_diff_pl(:,:,:)

  allocate( IxJ_GMTR_p   (ADM_iall,ADM_jall,K0,ADM_lall   ,      GMTR_p_nmax   ) )
  allocate( IxJ_GMTR_p_pl(ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_p_nmax   ) )
  allocate( IxJ_GMTR_t   (ADM_iall,ADM_jall,K0,ADM_lall   ,TI:TJ,GMTR_t_nmax   ) )
  allocate( IxJ_GMTR_t_pl(ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_t_nmax   ) )
  allocate( IxJ_GMTR_a   (ADM_iall,ADM_jall,K0,ADM_lall   ,AI:AJ,GMTR_a_nmax   ) )
  allocate( IxJ_GMTR_a_pl(ADM_gall_pl      ,K0,ADM_lall_pl,      GMTR_a_nmax_pl) )

  IxJ_GMTR_p    = reshape(GMTR_p,shape(IxJ_GMTR_p))
  IxJ_GMTR_p_pl = GMTR_p_pl
  IxJ_GMTR_t    = reshape(GMTR_t,shape(IxJ_GMTR_t))
  IxJ_GMTR_t_pl = GMTR_t_pl
  IxJ_GMTR_a    = reshape(GMTR_a,shape(IxJ_GMTR_a))
  IxJ_GMTR_a_pl = GMTR_a_pl

  write(*,*) "*** Finish initialize"

  !###############################################################################

  write(*,*) "*** Start  kernel"
  call DEBUG_rapstart('DC_setup_kernel')

  do iteration = 1, SET_iteration



     !--- calc geometrical information for cell point
     call GMTR_p_setup( GRD_x (:,:,:,:),   GRD_x_pl (:,:,:,:), & ! [IN]
                        GRD_xt(:,:,:,:,:), GRD_xt_pl(:,:,:,:), & ! [IN]
                        GRD_s (:,:,:,:),   GRD_s_pl (:,:,:,:), & ! [IN]
                        GMTR_p(:,:,:,:),   GMTR_p_pl(:,:,:,:), & ! [OUT]
                        GRD_rscale                             ) ! [IN]

     ! dummy for MPI communication
     do j = 1, ADM_jall
     do i = 1, ADM_iall
        g = (j-1)*ADM_gall_1d + i
        if (      i < ADM_imin .OR. i > ADM_imax &
             .OR. j < ADM_jmin .OR. j > ADM_jmax ) then

           GMTR_p(g,:,:,:) = check_GMTR_p(g,:,:,:)

        endif
     enddo
     enddo

     do g = ADM_gmin_pl, ADM_gmax_pl
        GMTR_p_pl(g,:,:,:) = check_GMTR_p_pl(g,:,:,:)
     enddo

     !--- calc geometrical information for cell vertex (triangle)
     call GMTR_t_setup( GRD_x (:,:,:,:),   GRD_x_pl (:,:,:,:), & ! [IN]
                        GRD_xt(:,:,:,:,:), GRD_xt_pl(:,:,:,:), & ! [IN]
                        GMTR_t(:,:,:,:,:), GMTR_t_pl(:,:,:,:), & ! [OUT]
                        GRD_rscale                             ) ! [IN]

     !--- calc geometrical information for cell arc
     call GMTR_a_setup( GRD_x (:,:,:,:),   GRD_x_pl (:,:,:,:), & ! [IN]
                        GRD_xt(:,:,:,:,:), GRD_xt_pl(:,:,:,:), & ! [IN]
                        GMTR_a(:,:,:,:,:), GMTR_a_pl(:,:,:,:), & ! [OUT]
                        GRD_rscale                             ) ! [IN]

     ! NOTE: The arrays for GMTR is reshaped after setup. This is tentative treatment

     call OPRT_divergence_setup( IxJ_GMTR_p    (:,:,:,:,:),   IxJ_GMTR_p_pl    (:,:,:,:), & ! [IN]
                                 IxJ_GMTR_t    (:,:,:,:,:,:), IxJ_GMTR_t_pl    (:,:,:,:), & ! [IN]
                                 IxJ_GMTR_a    (:,:,:,:,:,:), IxJ_GMTR_a_pl    (:,:,:,:), & ! [IN]
                                 OPRT_coef_div (:,:,:,:,:),   OPRT_coef_div_pl (:,:,:)    ) ! [OUT]

     call OPRT_rotation_setup  ( IxJ_GMTR_p    (:,:,:,:,:),   IxJ_GMTR_p_pl    (:,:,:,:), & ! [IN]
                                 IxJ_GMTR_t    (:,:,:,:,:,:), IxJ_GMTR_t_pl    (:,:,:,:), & ! [IN]
                                 IxJ_GMTR_a    (:,:,:,:,:,:), IxJ_GMTR_a_pl    (:,:,:,:), & ! [IN]
                                 OPRT_coef_rot (:,:,:,:,:),   OPRT_coef_rot_pl (:,:,:)    ) ! [OUT]

     call OPRT_gradient_setup  ( IxJ_GMTR_p    (:,:,:,:,:),   IxJ_GMTR_p_pl    (:,:,:,:), & ! [IN]
                                 IxJ_GMTR_t    (:,:,:,:,:,:), IxJ_GMTR_t_pl    (:,:,:,:), & ! [IN]
                                 IxJ_GMTR_a    (:,:,:,:,:,:), IxJ_GMTR_a_pl    (:,:,:,:), & ! [IN]
                                 OPRT_coef_grad(:,:,:,:,:),   OPRT_coef_grad_pl(:,:,:)    ) ! [OUT]

     call OPRT_laplacian_setup ( IxJ_GMTR_p    (:,:,:,:,:),   IxJ_GMTR_p_pl    (:,:,:,:), & ! [IN]
                                 IxJ_GMTR_t    (:,:,:,:,:,:), IxJ_GMTR_t_pl    (:,:,:,:), & ! [IN]
                                 IxJ_GMTR_a    (:,:,:,:,:,:), IxJ_GMTR_a_pl    (:,:,:,:), & ! [IN]
                                 OPRT_coef_lap (:,:,:,:),     OPRT_coef_lap_pl (:,:)      ) ! [OUT]

     call OPRT_diffusion_setup ( IxJ_GMTR_p    (:,:,:,:,:),   IxJ_GMTR_p_pl    (:,:,:,:), & ! [IN]
                                 IxJ_GMTR_t    (:,:,:,:,:,:), IxJ_GMTR_t_pl    (:,:,:,:), & ! [IN]
                                 IxJ_GMTR_a    (:,:,:,:,:,:), IxJ_GMTR_a_pl    (:,:,:,:), & ! [IN]
                                 OPRT_coef_intp(:,:,:,:,:,:), OPRT_coef_intp_pl(:,:,:,:), & ! [OUT]
                                 OPRT_coef_diff(:,:,:,:,:),   OPRT_coef_diff_pl(:,:,:)    ) ! [OUT]



     write(ADM_LOG_FID,*) '### Input ###'
     call DEBUG_valuecheck( 'GRD_x                  ', GRD_x                  (:,:,:,:)     )
     call DEBUG_valuecheck( 'GRD_x_pl               ', GRD_x_pl               (:,:,:,:)     )
     call DEBUG_valuecheck( 'GRD_xt                 ', GRD_xt                 (:,:,:,:,:)   )
     call DEBUG_valuecheck( 'GRD_xt_pl              ', GRD_xt_pl              (:,:,:,:)     )
     call DEBUG_valuecheck( 'GRD_s                  ', GRD_s                  (:,:,:,:)     )
     call DEBUG_valuecheck( 'GRD_s_pl               ', GRD_s_pl               (:,:,:,:)     )
     call DEBUG_valuecheck( 'check_GMTR_p           ', check_GMTR_p           (:,:,:,:)     )
     call DEBUG_valuecheck( 'check_GMTR_p_pl        ', check_GMTR_p_pl        (:,:,:,:)     )
     call DEBUG_valuecheck( 'check_GMTR_t           ', check_GMTR_t           (:,:,:,:,:)   )
     call DEBUG_valuecheck( 'check_GMTR_t_pl        ', check_GMTR_t_pl        (:,:,:,:)     )
     call DEBUG_valuecheck( 'check_GMTR_a           ', check_GMTR_a           (:,:,:,:,:)   )
     call DEBUG_valuecheck( 'check_GMTR_a_pl        ', check_GMTR_a_pl        (:,:,:,:)     )
     call DEBUG_valuecheck( 'check_OPRT_coef_div    ', check_OPRT_coef_div    (:,:,:,:,:)   )
     call DEBUG_valuecheck( 'check_OPRT_coef_div_pl ', check_OPRT_coef_div_pl (:,:,:)       )
     call DEBUG_valuecheck( 'check_OPRT_coef_rot    ', check_OPRT_coef_rot    (:,:,:,:,:)   )
     call DEBUG_valuecheck( 'check_OPRT_coef_rot_pl ', check_OPRT_coef_rot_pl (:,:,:)       )
     call DEBUG_valuecheck( 'check_OPRT_coef_grad   ', check_OPRT_coef_grad   (:,:,:,:,:)   )
     call DEBUG_valuecheck( 'check_OPRT_coef_grad_pl', check_OPRT_coef_grad_pl(:,:,:)       )
     call DEBUG_valuecheck( 'check_OPRT_coef_lap    ', check_OPRT_coef_lap    (:,:,:,:)     )
     call DEBUG_valuecheck( 'check_OPRT_coef_lap_pl ', check_OPRT_coef_lap_pl (:,:)         )
     call DEBUG_valuecheck( 'check_OPRT_coef_intp   ', check_OPRT_coef_intp   (:,:,:,:,:,:) )
     call DEBUG_valuecheck( 'check_OPRT_coef_intp_pl', check_OPRT_coef_intp_pl(:,:,:,:)     )
     call DEBUG_valuecheck( 'check_OPRT_coef_diff   ', check_OPRT_coef_diff   (:,:,:,:,:)   )
     call DEBUG_valuecheck( 'check_OPRT_coef_diff_pl', check_OPRT_coef_diff_pl(:,:,:)       )
     write(ADM_LOG_FID,*) '### Output ###'
     call DEBUG_valuecheck( 'GMTR_p           ', GMTR_p           (:,:,:,:)     )
     call DEBUG_valuecheck( 'GMTR_p_pl        ', GMTR_p_pl        (:,:,:,:)     )
     call DEBUG_valuecheck( 'GMTR_t           ', GMTR_t           (:,:,:,:,:)   )
     call DEBUG_valuecheck( 'GMTR_t_pl        ', GMTR_t_pl        (:,:,:,:)     )
     call DEBUG_valuecheck( 'GMTR_a           ', GMTR_a           (:,:,:,:,:)   )
     call DEBUG_valuecheck( 'GMTR_a_pl        ', GMTR_a_pl        (:,:,:,:)     )
     call DEBUG_valuecheck( 'OPRT_coef_div    ', OPRT_coef_div    (:,:,:,:,:)   )
     call DEBUG_valuecheck( 'OPRT_coef_div_pl ', OPRT_coef_div_pl (:,:,:)       )
     call DEBUG_valuecheck( 'OPRT_coef_rot    ', OPRT_coef_rot    (:,:,:,:,:)   )
     call DEBUG_valuecheck( 'OPRT_coef_rot_pl ', OPRT_coef_rot_pl (:,:,:)       )
     call DEBUG_valuecheck( 'OPRT_coef_grad   ', OPRT_coef_grad   (:,:,:,:,:)   )
     call DEBUG_valuecheck( 'OPRT_coef_grad_pl', OPRT_coef_grad_pl(:,:,:)       )
     call DEBUG_valuecheck( 'OPRT_coef_lap    ', OPRT_coef_lap    (:,:,:,:)     )
     call DEBUG_valuecheck( 'OPRT_coef_lap_pl ', OPRT_coef_lap_pl (:,:)         )
     call DEBUG_valuecheck( 'OPRT_coef_intp   ', OPRT_coef_intp   (:,:,:,:,:,:) )
     call DEBUG_valuecheck( 'OPRT_coef_intp_pl', OPRT_coef_intp_pl(:,:,:,:)     )
     call DEBUG_valuecheck( 'OPRT_coef_diff   ', OPRT_coef_diff   (:,:,:,:,:)   )
     call DEBUG_valuecheck( 'OPRT_coef_diff_pl', OPRT_coef_diff_pl(:,:,:)       )
  enddo

  write(ADM_LOG_FID,*) '### Varidation by diff ###'
  check_GMTR_p           (:,:,:,:)     = check_GMTR_p           (:,:,:,:)     - GMTR_p           (:,:,:,:)
  check_GMTR_p_pl        (:,:,:,:)     = check_GMTR_p_pl        (:,:,:,:)     - GMTR_p_pl        (:,:,:,:)
  check_GMTR_t           (:,:,:,:,:)   = check_GMTR_t           (:,:,:,:,:)   - GMTR_t           (:,:,:,:,:)
  check_GMTR_t_pl        (:,:,:,:)     = check_GMTR_t_pl        (:,:,:,:)     - GMTR_t_pl        (:,:,:,:)
  check_GMTR_a           (:,:,:,:,:)   = check_GMTR_a           (:,:,:,:,:)   - GMTR_a           (:,:,:,:,:)
  check_GMTR_a_pl        (:,:,:,:)     = check_GMTR_a_pl        (:,:,:,:)     - GMTR_a_pl        (:,:,:,:)
  check_OPRT_coef_div    (:,:,:,:,:)   = check_OPRT_coef_div    (:,:,:,:,:)   - OPRT_coef_div    (:,:,:,:,:)
  check_OPRT_coef_div_pl (:,:,:)       = check_OPRT_coef_div_pl (:,:,:)       - OPRT_coef_div_pl (:,:,:)
  check_OPRT_coef_rot    (:,:,:,:,:)   = check_OPRT_coef_rot    (:,:,:,:,:)   - OPRT_coef_rot    (:,:,:,:,:)
  check_OPRT_coef_rot_pl (:,:,:)       = check_OPRT_coef_rot_pl (:,:,:)       - OPRT_coef_rot_pl (:,:,:)
  check_OPRT_coef_grad   (:,:,:,:,:)   = check_OPRT_coef_grad   (:,:,:,:,:)   - OPRT_coef_grad   (:,:,:,:,:)
  check_OPRT_coef_grad_pl(:,:,:)       = check_OPRT_coef_grad_pl(:,:,:)       - OPRT_coef_grad_pl(:,:,:)
  check_OPRT_coef_lap    (:,:,:,:)     = check_OPRT_coef_lap    (:,:,:,:)     - OPRT_coef_lap    (:,:,:,:)
  check_OPRT_coef_lap_pl (:,:)         = check_OPRT_coef_lap_pl (:,:)         - OPRT_coef_lap_pl (:,:)
  check_OPRT_coef_intp   (:,:,:,:,:,:) = check_OPRT_coef_intp   (:,:,:,:,:,:) - OPRT_coef_intp   (:,:,:,:,:,:)
  check_OPRT_coef_intp_pl(:,:,:,:)     = check_OPRT_coef_intp_pl(:,:,:,:)     - OPRT_coef_intp_pl(:,:,:,:)
  check_OPRT_coef_diff   (:,:,:,:,:)   = check_OPRT_coef_diff   (:,:,:,:,:)   - OPRT_coef_diff   (:,:,:,:,:)
  check_OPRT_coef_diff_pl(:,:,:)       = check_OPRT_coef_diff_pl(:,:,:)       - OPRT_coef_diff_pl(:,:,:)
  call DEBUG_valuecheck( 'check_GMTR_p           ', check_GMTR_p           (:,:,:,:)     )
  call DEBUG_valuecheck( 'check_GMTR_p_pl        ', check_GMTR_p_pl        (:,:,:,:)     )
  call DEBUG_valuecheck( 'check_GMTR_t           ', check_GMTR_t           (:,:,:,:,:)   )
  call DEBUG_valuecheck( 'check_GMTR_t_pl        ', check_GMTR_t_pl        (:,:,:,:)     )
  call DEBUG_valuecheck( 'check_GMTR_a           ', check_GMTR_a           (:,:,:,:,:)   )
  call DEBUG_valuecheck( 'check_GMTR_a_pl        ', check_GMTR_a_pl        (:,:,:,:)     )
  call DEBUG_valuecheck( 'check_OPRT_coef_div    ', check_OPRT_coef_div    (:,:,:,:,:)   )
  call DEBUG_valuecheck( 'check_OPRT_coef_div_pl ', check_OPRT_coef_div_pl (:,:,:)       )
  call DEBUG_valuecheck( 'check_OPRT_coef_rot    ', check_OPRT_coef_rot    (:,:,:,:,:)   )
  call DEBUG_valuecheck( 'check_OPRT_coef_rot_pl ', check_OPRT_coef_rot_pl (:,:,:)       )
  call DEBUG_valuecheck( 'check_OPRT_coef_grad   ', check_OPRT_coef_grad   (:,:,:,:,:)   )
  call DEBUG_valuecheck( 'check_OPRT_coef_grad_pl', check_OPRT_coef_grad_pl(:,:,:)       )
  call DEBUG_valuecheck( 'check_OPRT_coef_lap    ', check_OPRT_coef_lap    (:,:,:,:)     )
  call DEBUG_valuecheck( 'check_OPRT_coef_lap_pl ', check_OPRT_coef_lap_pl (:,:)         )
  call DEBUG_valuecheck( 'check_OPRT_coef_intp   ', check_OPRT_coef_intp   (:,:,:,:,:,:) )
  call DEBUG_valuecheck( 'check_OPRT_coef_intp_pl', check_OPRT_coef_intp_pl(:,:,:,:)     )
  call DEBUG_valuecheck( 'check_OPRT_coef_diff   ', check_OPRT_coef_diff   (:,:,:,:,:)   )
  call DEBUG_valuecheck( 'check_OPRT_coef_diff_pl', check_OPRT_coef_diff_pl(:,:,:)       )

  call DEBUG_rapend('DC_setup_kernel')
  write(*,*) "*** Finish kernel"

  !###############################################################################

  call DEBUG_rapreport

  stop
end program dckernel_setup
!-------------------------------------------------------------------------------
