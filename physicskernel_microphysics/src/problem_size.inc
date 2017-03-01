!-------------------------------------------------------------------------------
!
! Problem size and global parameters
!
!-------------------------------------------------------------------------------

  integer,  public, parameter :: H_SHORT      = 32
  integer,  public, parameter :: H_LONG       = 1024
  integer,  public            :: IO_FID_LOG   = 6
  integer,  public, parameter :: ADM_prc_all  = 1

  integer,  public, parameter :: IO_FREAD = 0

  integer,  public, parameter :: ADM_KNONE    = 1

!  integer,  public            :: ADM_gall         = 16900 ! horizontal grid number
  integer,  public            :: ADM_gall_in_orig = 16641 ! horizontal grid number (for physics)
  integer,  public            :: ADM_gall_in      = 16641 ! horizontal grid number (for physics)
  integer,  public            :: ADM_vlayer   = 94
  integer,  public            :: ADM_kall     = 96   ! vertical   grid number
  integer,  public            :: ADM_kmin     = 2
  integer,  public            :: ADM_kmax     = 95
  integer,  public            :: ADM_lall     = 1
  integer,  public            :: TRC_VMAX     = 6

  integer,  public            :: kdim         = 96
  integer,  public            :: kmin         = 2
  integer,  public            :: kmax         = 95
  integer,  public            :: knone        = 1
  integer,  public            :: nqmax        = 6

  real(RP), public            :: TIME_DTL     = 6.0_RP ! [sec]
  real(RP), public            :: TIME_CTIME   = 0.0_RP
  integer,  public            :: TIME_CSTEP   = 49
  integer,  public            :: SET_l        = 1
  integer,  public            :: SET_rgnid    = 4

  ! const
  real(RP), public            :: CONST_PI     = 3.141592653589793_RP !< pi
  real(RP), public            :: CONST_EPS    = 2.220446E-16_RP      !< small number

  real(RP), public, parameter :: CONST_UNDEF  = -9.9999E+30_RP
  real(RP), public, parameter :: CONST_GRAV   =   9.79764_RP
  real(RP), public, parameter :: CONST_STB    =   5.670373E-8_RP

  real(RP), public, parameter :: CONST_Rdry   =  287.04_RP
  real(RP), public, parameter :: CONST_CPdry  = 1004.64_RP
  real(RP), public, parameter :: CONST_CVdry  =  717.56_RP
  real(RP), public, parameter :: CONST_Rvap   =  461.46_RP
  real(RP), public, parameter :: CONST_CPvap  = 1846.0_RP
  real(RP), public, parameter :: CONST_CVvap  = 1384.54_RP
  real(RP), public, parameter :: CONST_CL     = 4218.0_RP
  real(RP), public, parameter :: CONST_CI     = 2006.0_RP

  real(RP), public, parameter :: CONST_KAPPA  = 0.2857256619550070_RP
  real(RP), public, parameter :: CONST_EPSvap = 0.6219718309859156_RP
  real(RP), public, parameter :: CONST_PSAT0  =     610.7_RP
  real(RP), public, parameter :: CONST_LHV0   = 2501000.0_RP
  real(RP), public, parameter :: CONST_LHS0   = 2834000.0_RP
  real(RP), public, parameter :: CONST_LHF0   =  333000.0_RP
  real(RP), public, parameter :: CONST_LHV00  = 3148911.8_RP
  real(RP), public, parameter :: CONST_LHS00  = 2877704.0_RP
  real(RP), public, parameter :: CONST_LHF00  = -271207.8_RP
  real(RP), public, parameter :: CONST_Pstd   =  101325.0_RP
  real(RP), public, parameter :: CONST_PRE00  =  100000.0_RP
  real(RP), public, parameter :: CONST_Tstd   =    288.15_RP
  real(RP), public, parameter :: CONST_TEM00  =    273.15_RP
  real(RP), public, parameter :: CONST_TMELT  =    273.15_RP
  real(RP), public, parameter :: CONST_DWATR  =    1000.0_RP
  real(RP), public, parameter :: CONST_DICE   =     916.8_RP
  real(RP), public            :: CONST_LHV    = CONST_LHV0
  real(RP), public            :: CONST_LHS    = CONST_LHS0
  real(RP), public            :: CONST_LHF    = CONST_LHF0

  real(RP), public            :: PI    = 3.141592653589793_RP !< pi
  real(RP), public            :: EPS   = 2.220446E-16_RP      !< small number
  real(RP), public, parameter :: UNDEF = CONST_UNDEF
  real(RP), public, parameter :: GRAV  = CONST_GRAV
  real(RP), public, parameter :: Rdry  = CONST_Rdry
  real(RP), public, parameter :: Rvap  = CONST_Rvap
  real(RP), public, parameter :: CL    = CONST_CL
  real(RP), public            :: LHV   = CONST_LHV0
  real(RP), public            :: LHS   = CONST_LHS0
  real(RP), public            :: LHF   = CONST_LHF0
  real(RP), public, parameter :: rho_w = CONST_DWATR
  real(RP), public, parameter :: rhow  = CONST_DWATR
  real(RP), public, parameter :: DWATR = CONST_DWATR
  real(RP), public, parameter :: DICE  = CONST_DICE
  real(RP), public, parameter :: Pstd  = CONST_Pstd
  real(RP), public, parameter :: Tstd  = CONST_Tstd
  real(RP), public, parameter :: TEM00 = CONST_TEM00

  ! grd, gmtr, vmtr
  character(len=H_LONG), public :: vgrid_fname ='./vgrid94.dat'

  real(RP), public, allocatable :: GRD_gz   (:)
  real(RP), public, allocatable :: GRD_gzh  (:)
  real(RP), public, allocatable :: GRD_dgz  (:)
  real(RP), public, allocatable :: GRD_dgzh (:)
  real(RP), public, allocatable :: GRD_rdgz (:)
  real(RP), public, allocatable :: GRD_rdgzh(:)
  real(RP), public, allocatable :: GRD_afact(:)
  real(RP), public, allocatable :: GRD_bfact(:)
  real(RP), public, allocatable :: GRD_cfact(:)
  real(RP), public, allocatable :: GRD_dfact(:)

  ! runconf
  character(len=H_SHORT), public, parameter :: EIN_TYPE            = 'SIMPLE'
  character(len=H_SHORT), public, parameter :: RAIN_TYPE           = 'COLD'
  character(len=H_SHORT), public, parameter :: MP_TYPE             = 'NSW6'
  integer,                public, parameter :: MP_DIV_NUM          = 1
  logical,                public, parameter :: opt_2moment_water   = .false.
  logical,                public, parameter :: ISOTOPE             = .false.
  logical,                public, parameter :: opt_offline_aerosol = .false.
  logical,                public, parameter :: opt_aerosol_forcing = .false.
  character(len=H_SHORT), public, parameter :: PRCIP_TRN_ECORRECT  = 'KIN2EIN'

  character(len=H_SHORT), public, parameter :: RD_TYPE             = 'MSTRNX'
  character(len=H_SHORT), public, parameter :: AE_TYPE             = 'NONE'
  integer,                public, parameter :: KAPCL               = 7
  integer,                public, parameter :: NCRF                = 2
  integer,                public, parameter :: NRBND               = 3
  integer,                public, parameter :: NRBND_VIS           = 1
  integer,                public, parameter :: NRBND_NIR           = 2
  integer,                public, parameter :: NRBND_IR            = 3
  integer,                public, parameter :: NRDIR               = 2
  integer,                public, parameter :: NRDIR_DIRECT        = 1
  integer,                public, parameter :: NRDIR_DIFFUSE       = 2
  integer,                public, parameter :: NTAU                = 7
  integer,                public, parameter :: NPRES               = 7
  integer,                public, parameter :: HYDRO_MAX           = 7

  integer, public :: NQW_STR = 2
  integer, public :: NQW_END = 6
  integer, public :: I_QV =  1
  integer, public :: I_QC =  2
  integer, public :: I_QR =  3
  integer, public :: I_QI =  4
  integer, public :: I_QS =  5
  integer, public :: I_QG =  6
  integer, public :: I_QH =  7
  integer, public :: I_NC =  7
  integer, public :: I_NR =  8
  integer, public :: I_NI =  9
  integer, public :: I_NS = 10
  integer, public :: I_NG = 11

  character(len=16), public :: TRC_name(6) = (/ 'QV', &
                                                'QC', &
                                                'QR', &
                                                'QI', &
                                                'QS', &
                                                'QG'  /)

  real(RP), public :: CVW(6) = (/ CONST_CVdry, &
                                  CONST_CVdry, &
                                  CONST_CVdry, &
                                  CONST_CVdry, &
                                  CONST_CVdry, &
                                  CONST_CVdry  /)
  real(RP), public :: CPW(6) = (/ CONST_CPdry, &
                                  CONST_CVdry, &
                                  CONST_CVdry, &
                                  CONST_CVdry, &
                                  CONST_CVdry, &
                                  CONST_CVdry  /)

  integer, public :: SET_iteration = 1
  logical, public :: SET_check     = .true.
