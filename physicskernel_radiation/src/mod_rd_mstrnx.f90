!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Radiative Transfer
!!
!! @par Description
!!          2-stream, k-distribution broadband radiative transfer scheme mstrnX
!!          Reference : Nakajima and Tanaka(1986)   : J.Quant.Spectrosc.Radiat.Transfer,35,pp.13–21
!!                      Nakajima et al.(2000)       : Appl.Opt.,39,pp.4869–4878
!!                      Sekiguchi and Nakajima(2008): J.Quant.Spectrosc.Radiat.Transfer,109,pp.2779–2793
!!
!! @author NICAM developers
!!
!<
!-------------------------------------------------------------------------------
module mod_rd_mstrnx
  !-----------------------------------------------------------------------------
  !
  !++ Note : radiation grid
  !
  !      KMAX_RAD = kmax - kmin + 1 = kdim - 2
  !      - - - : full level
  !      ===== : half level
  !
  !      - - - - - - - - - - - kmax+1=kdim              : NA
  !      (top)================             kmax+1       :          1
  !      - - - - - - - - - - - kmax                     : 1
  !      ---------------------                          :
  !      - - - - - - - - - - -                          :
  !      ---------------------                          :
  !      - - - - - - - - - - - k                        : kmax+1-k
  !      ---------------------             k            :          kmax+2-k
  !      - - - - - - - - - - -                          :
  !      ---------------------                          :
  !      - - - - - - - - - - - kmin                     : KMAX_RAD
  !      (surface)============             kmin         :          KMAX_RAD+1
  !      - - - - - - - - - - - kmin-1=1                 : NA
  !      (not used)-----------             1            :          NA
  !
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
  public :: rd_mstrnx_init
  public :: rd_mstrnx

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: OPPARM2
  private :: get_cot_550nm
  private :: DTRN31
  private :: PTFIT2
  private :: CNTCFC2
  private :: RMDIDX
  private :: SCATAE
  private :: SCATRY
  private :: SCATCL
  private :: PLKEXP
  private :: PLANKS
  private :: PLANKF
  private :: TWST
  private :: ADDING
  private :: RTS_MR
  private :: ADDING_MR
  private :: BCVR
  private :: MULTMR
  private :: M3X2
  private :: M3X3

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,  private, parameter :: NRBFLX =  4  ! 3:SR,ALB=0 4:LR,GTMP+1
  integer,  private, parameter :: KMOL   =  6  ! 1:H2O, 2:CO2, 3:O3, 4:N2O, 5:CH4, 6:O2

  integer,  private, parameter :: KPLK   =  2
  integer,  private, parameter :: KDA    =  1
  integer,  private, parameter :: KCLD   =  3

  integer,  private, parameter :: KSFC   =  7
  integer,  private, parameter :: KWNB   = 29
  integer,  private, parameter :: KFLG   =  8
  integer,  private, parameter :: KCH    = 13
  integer,  private, parameter :: KAO3   =  3
  integer,  private, parameter :: KPG    = 26
  integer,  private, parameter :: KTG    =  3
  integer,  private, parameter :: KDMAX  =  6
  ! [comment] 08/05/30 T.Mitsui
  !  KDMAX is Number of Moment of Phase Function and two-stream aprrox. never used high order.
  !  1: Ce/V, 2: Ca/V, 3: Cs*G(2)/V,..., 6: Cs*G(5)/V with G(1)=1
  !  See. Nakajima et al.(2000), Applied Optics, eq.(15)

!ESC!  integer,  private, parameter :: KAPCL = 7
  ! 1: Soil Dust
  ! 2: Carbonaceous (BC/OC=0.3)
  ! 3: Carbonaceous (BC/OC=0.15)
  ! 4: Carbonaceous (BC/OC=0)
  ! 5: Black Carbon (external mixture)
  ! 6: Sulfate
  ! 7: Sea Salt

  integer,  private, parameter :: KPLNK  =  5
  integer,  private, parameter :: KH2O   =  3
  integer,  private, parameter :: KKDT   =  3
  integer,  private, parameter :: KKDP   =  8
  integer,  private, parameter :: KCFC   = 28

  integer,  private, parameter :: KCTYP  =  2 ! cloud type (St,Cu)

!ESC!  integer,  private, parameter :: NCRF =  2
!ESC!  integer,  private, parameter :: NRDIR = 2
!ESC!  integer,  private, parameter :: NRBND = 3
  integer,  private, parameter :: ICMAX = 3

  integer,  private :: KMAX_RAD ! radiation grid number (downward)

  ! Lambert Parameter
  real(RP), private            :: AMUI(2)
  real(RP), private            :: AMUA(2)
  real(RP), private            :: AMX (2)
  real(RP), private            :: WMP (2)
  real(RP), private            :: WMM (2)
  real(RP), private            :: WMX (2)
  real(RP), private            :: PA0 (2)
  real(RP), private            :: PA1 (2)
  real(RP), private            :: P00 (2)
  real(RP), private            :: P01F(2)
  real(RP), private            :: EPS (2)
  real(RP), private            :: EPST

  integer,  private            :: NWNB
  integer,  private            :: NMOL ( KWNB )
  integer,  private            :: MID ( KMOL, KWNB ) !! Molecular ID
  real(RP), private            :: WNBND ( KWNB+1 )
  integer,  private            :: iflgb ( KFLG, KWNB )
  real(RP), private            :: APLNK ( KPLNK, KWNB )
  real(RP), private            :: FSOL ( KWNB ) !! solar spectra reference
  real(RP), private            :: FSTOT !! solar constant reference
  real(RP), private            :: RY ( KWNB )
  real(RP), private            :: QMOL ( KDMAX )
  integer,  private            :: NCH ( KWNB )
  real(RP), private            :: WGT0 ( KCH, KWNB )
  integer,  private            :: NTG, NPG
  real(RP), private            :: TG ( KTG )
  real(RP), private            :: PLG ( KPG )
  real(RP), private            :: AKD ( KPG, KTG, KCH, KMOL, KWNB )
  real(RP), private            :: SKD ( KPG, KTG, KCH, KWNB )
  real(RP), private            :: ACFC ( KCFC, KWNB )

  !  Q is Moment of Phase Function and two-stream aprrox.
  !  1=Extinction Coefficient
  !  2=Absorption Coefiicient
  !  3=Asymmetry factor
  !  4=Truncation factor( to use Delta-Function Adjustment )
  real(RP), private, allocatable :: Q    (:,:,:,:) ! (KMODE,KTPCL,KDA*2+2,KWNB)
  real(RP), private, allocatable :: RMODE(:,:,:)   ! (0:KMODE+1,KTPCL,2)
  integer,  private, allocatable :: NMODE(:)       ! (KTPCL)

  ! NAMELIST : NM_RD_MSTRN_TR
  character(len=H_LONG),  private :: PARA         = 'PARA.bnd29ch111sp'  ! parameter file
  real(RP),               private :: PPMCO2       = 345.0_RP             ! CO2 concentration
  real(RP),               private :: PPMN2O       =   0.3_RP             ! N2O concentration
  real(RP),               private :: PPMCH4       =   1.7_RP             ! CH4 concentration
  real(RP),               private :: PPMCFC(KCFC)                        ! CFC concentration
  real(RP),               private :: PPMO2        =  2.1E+5_RP           ! O2  concentration
  real(RP),               private :: QSTR         =  2.5E-6_RP           ! q of stratosphere
  real(RP),               private :: PSTRMX       = 300.E+2_RP           ! max p of stratopause
  real(RP),               private :: PSTRMN       =  50.E+2_RP           ! min p of stratopause
  real(RP),               private :: GCRSTR       =   1.E-4_RP           ! crit. dT/dz tropopause

  ! NAMELIST : NM_RD_MSTRN_CLOUD
  logical,                private :: MAX_RND_CLD_OVERLAP = .false.       ! use maximal/random overlapping method?
  logical,                private :: OPT_RE_CLOUDMP      = .false.       ! use cloud effective radius calculated in microphysics scheme?
  logical,                private :: use_ice_dens        = .false.       ! use density of ice for ice cloud, instead of density of liquid water?
  real(RP),               private :: radius_cumulus      = 40.E-6_RP     ! cloud drop radius (cumulus)
  real(RP),               private :: radius_stratus      =  8.E-6_RP     ! cloud drop radius (stratus)

  real(RP),               private :: Q2PPM                               ! conversion factor for QV (kg/kg->PPM)
  real(RP),               private :: C2PPM                               ! conversion factor for QC (kg/kg->PPM)
  real(RP),               private :: I2PPM                               ! conversion factor for QI (kg/kg->PPM)

  ! NAMELIST : NM_RD_MSTRN_AEROSOL
  logical,                private :: OVARRA = .true.                     ! use aerosol forcing
  logical,                private :: OAEROF = .false.                    ! compute aerosol radiative forcing (aero-noaero)
  real(RP),               private :: ARHMAX = 0.99_RP                    ! maximum RH
  real(RP),               private :: RCMIN  = 1.E-6_RP                   ! minimum of radius
  ! [FIX] 2016/05/13 T.Seiki
!!$  real(RP),            private :: RCMAX  = 40.E-6_RP                  ! maximum of radius
  real(RP),               private :: RCMAX  = 1000.E-6_RP                ! maximum of radius [m]

  ! NAMELIST : NM_RD_MSTRN_ISCCP
  character(len=H_SHORT), private :: isccp_version  = 'NONE'
  character(len=H_LONG),  private :: TAUTBL         = 'tautab.formatted' ! tau tbl file
  character(len=H_LONG),  private :: INVTBL         = 'invtau.formatted' ! inversion file
  integer,                private :: ncol           = 1                  ! number of subcolumns in ISCCP
  integer,                private :: overlap        = 1                  ! cloud overlap type (1=max,2=rand,3=max/rand)
  character(len=H_SHORT), private :: DFQ_TYPE       = 'UNDEF'            ! nighttime values of dfq_isccp
                                                                         ! 'UNDEF' :  undef
                                                                         ! 'FQ'    :  fq_isccp
                                                                         ! '-FQ'   : -fq_isccp
  logical,                private :: calc_nighttime = .false.            ! calc nighttime (no sunlit) region?
  integer,                private :: debug          = 0
  integer,                private :: debugcol       = 0
  integer,                private :: ncolprint      = 0                  ! print out column decomposition

  real(4),                private :: tautab(0:255)                       ! tau table
  integer,                private :: invtau(-20:45000)                   ! inversion table

  real(RP), private, allocatable :: boxtau  (:,:)                        ! optical thickness       in each column
  real(RP), private, allocatable :: boxptop (:,:)                        ! cloud top pressure (mb) in each column
  real(RP), private, allocatable :: frac_out(:,:,:)

  ! NAMELIST : NM_RD_MSTRN_INIT
  integer,                private :: KMODE                      ! mode radius grid
  integer,                private :: KTPCL
  integer,                private :: KCPCL             = 2
  logical,                private :: opt_cloud_onoff   = .true.
  logical,                private :: opt_rain_onoff    = .true.
  logical,                private :: opt_ice_onoff     = .true.
  logical,                private :: opt_snow_onoff    = .true.
  logical,                private :: opt_graupel_onoff = .true.

  integer,                private :: I_QC_RD           = -999
  integer,                private :: I_QR_RD           = -999
  integer,                private :: I_QI_RD           = -999
  integer,                private :: I_QS_RD           = -999
  integer,                private :: I_QG_RD           = -999

  integer,  private, allocatable :: nq2IP(:)
  integer,  private, allocatable :: IP2nq(:)
  logical,  private, allocatable :: opt_cloud_type_onoff(:)

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> setup
  subroutine rd_mstrnx_init( &
       ijdim )
!ESC!    use mod_adm, only: &
!ESC!       ADM_proc_stop,    &
!ESC!       kmin => ADM_kmin, &
!ESC!       kmax => ADM_kmax
!ESC!    use mod_const, only: &
!ESC!       Rdry  => CONST_Rdry,  &
!ESC!       RVAP  => CONST_RVAP,  &
!ESC!       DWATR => CONST_DWATR, &
!ESC!       DICE  => CONST_DICE,  &
!ESC!       Tstd  => CONST_TEM00, &
!ESC!       Pstd  => CONST_Pstd
!ESC!    use mod_runconf, only: &
!ESC!       AE_TYPE,             &
!ESC!       opt_aerosol_forcing, &
!ESC!       opt_offline_aerosol, &
!ESC!       I_QC, I_QR, I_QI,    &
!ESC!       I_QS, I_QG, I_QH,    &
!ESC!       HYDRO_MAX
    implicit none

    integer,  intent(in) :: ijdim

    NAMELIST / NM_RD_MSTRN_TR / &
       PARA,   &
       PPMCO2, &
       PPMN2O, &
       PPMCH4, &
       PPMCFC, &
       PPMO2,  &
       QSTR,   &
       PSTRMX, &
       PSTRMN, &
       GCRSTR

    NAMELIST / NM_RD_MSTRN_CLOUD / &
       MAX_RND_CLD_OVERLAP, &
       OPT_RE_CLOUDMP,      &
       use_ice_dens,        &
       radius_cumulus,      &
       radius_stratus

    NAMELIST / NM_RD_MSTRN_AEROSOL / &
       OVARRA, &
       OAEROF, &
       ARHMAX, &
       RCMIN,  &
       RCMAX

    NAMELIST / NM_RD_MSTRN_ISCCP / &
       isccp_version,  &
       TAUTBL,         &
       INVTBL,         &
       ncol,           &
       OVERLAP,        &
       DFQ_TYPE,       &
       calc_nighttime, &
       debug,          &
       debugcol,       &
       ncolprint

    NAMELIST / NM_RD_MSTRN_INIT / &
       KCPCL,            &
       opt_cloud_onoff,  &
       opt_rain_onoff,   &
       opt_ice_onoff,    &
       opt_snow_onoff,   &
       opt_graupel_onoff

    integer  :: IW, IP, nq
    integer  :: ifid, ierr
    !---------------------------------------------------------------------------

    KMAX_RAD = kmax - kmin + 1

    !--- read parameters
    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '+++ Module[RD mstrnX AR5]/Category[nhm physics]'

    PPMCFC(:) = 1.E-4_RP

!ESC!    rewind(IO_FID_CONF)
!ESC!    read(IO_FID_CONF,nml=NM_RD_MSTRN_TR,iostat=ierr)
!ESC!    if ( ierr < 0 ) then
!ESC!       write(IO_FID_LOG,*) '*** NM_RD_MSTRN_TR is not specified. use default.'
!ESC!    elseif( ierr > 0 ) then
!ESC!       write(*,         *) 'xxx Not appropriate names in namelist NM_RD_MSTRN_TR. STOP.'
!ESC!       write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist NM_RD_MSTRN_TR. STOP.'
!ESC!       call ADM_proc_stop
!ESC!    endif
!ESC!    write(IO_FID_LOG,nml=NM_RD_MSTRN_TR)
!ESC!
!ESC!    rewind(IO_FID_CONF)
!ESC!    read(IO_FID_CONF,nml=NM_RD_MSTRN_CLOUD,iostat=ierr)
!ESC!    if ( ierr < 0 ) then
!ESC!       write(IO_FID_LOG,*) '*** NM_RD_MSTRN_CLOUD is not specified. use default.'
!ESC!    elseif( ierr > 0 ) then
!ESC!       write(*,         *) 'xxx Not appropriate names in namelist NM_RD_MSTRN_CLOUD. STOP.'
!ESC!       write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist NM_RD_MSTRN_CLOUD. STOP.'
!ESC!       call ADM_proc_stop
!ESC!    endif
!ESC!    write(IO_FID_LOG,nml=NM_RD_MSTRN_CLOUD)

PPMCO2 =  393.8_RP
PPMCH4 =  1.886_RP
PPMN2O = 0.3247_RP

MAX_RND_CLD_OVERLAP = .false.
OPT_RE_CLOUDMP      = .true.
use_ice_dens        = .true.

    Q2PPM = 1.E+6_RP * RVAP / Rdry                     ! [kg/kg] => 1.E+6*[m/m] (parts per million)
    C2PPM = 1.E+6_RP * Pstd  / ( Rdry*Tstd ) / DWATR   ! [kg/kg] => 1.E+6*[m/m] (parts per million)
    if ( use_ice_dens ) then
       I2PPM = 1.E+6_RP * Pstd / ( Rdry*Tstd ) / DICE  ! [kg/kg] => 1.E+6*[m/m] (parts per million)
    else
       I2PPM = 1.E+6_RP * Pstd / ( Rdry*Tstd ) / DWATR ! [kg/kg] => 1.E+6*[m/m] (parts per million)
    endif

    if ( opt_aerosol_forcing ) then
       OAEROF = .true.
    else
       OAEROF = .false.
    endif

!ESC!    rewind(IO_FID_CONF)
!ESC!    read(IO_FID_CONF,nml=NM_RD_MSTRN_AEROSOL,iostat=ierr)
!ESC!    if ( ierr < 0 ) then
!ESC!       write(IO_FID_LOG,*) '*** NM_RD_MSTRN_AEROSOL is not specified. use default.'
!ESC!    elseif( ierr > 0 ) then
!ESC!       write(*,         *) 'xxx Not appropriate names in namelist NM_RD_MSTRN_AEROSOL. STOP.'
!ESC!       write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist NM_RD_MSTRN_AEROSOL. STOP.'
!ESC!       call ADM_proc_stop
!ESC!    endif
!ESC!    write(IO_FID_LOG,nml=NM_RD_MSTRN_AEROSOL)

    if (      AE_TYPE == 'SPRINTARS'     &
         .OR. AE_TYPE == 'SPRINTARS_CRM' &
         .OR. opt_offline_aerosol        ) then
       ! do_nothing
    else
       OVARRA = .false. ! force false
       OAEROF = .false. ! force false
    endif

    if ( OAEROF ) then
       OVARRA = .true. ! force true
    endif

!ESC!    rewind(IO_FID_CONF)
!ESC!    read(IO_FID_CONF,nml=NM_RD_MSTRN_ISCCP,iostat=ierr)
!ESC!    if ( ierr < 0 ) then
!ESC!       write(IO_FID_LOG,*) '*** NM_RD_MSTRN_ISCCP is not specified. use default.'
!ESC!    elseif( ierr > 0 ) then
!ESC!       write(*,         *) 'xxx Not appropriate names in namelist NM_RD_MSTRN_ISCCP. STOP.'
!ESC!       write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist NM_RD_MSTRN_ISCCP. STOP.'
!ESC!       call ADM_proc_stop
!ESC!    endif
!ESC!    write(IO_FID_LOG,nml=NM_RD_MSTRN_ISCCP)

    allocate( boxtau  (ijdim,ncol)          )
    allocate( boxptop (ijdim,ncol)          )
    allocate( frac_out(ijdim,ncol,KMAX_RAD) )

!ESC!    if ( isccp_version == 'old' ) then
!ESC!       write(IO_FID_LOG,*) '*** Reading TAUTBL for ISCCP simulator : ', trim(TAUTBL)
!ESC!       ifid = IO_get_available_fid()
!ESC!       open(ifid,file=trim(TAUTBL),form='formatted',status='old')
!ESC!          read(ifid,*) tautab
!ESC!       close(ifid)
!ESC!
!ESC!       write(IO_FID_LOG,*) '*** Reading INVTBL for ISCCP simulator : ', trim(INVTBL)
!ESC!       ifid = IO_get_available_fid()
!ESC!       open(ifid,file=trim(INVTBL),form='formatted',status='old')
!ESC!          read(ifid,*) invtau
!ESC!       close(ifid)
!ESC!    endif

!ESC!    rewind(IO_FID_CONF)
!ESC!    read(IO_FID_CONF,nml=NM_RD_MSTRN_INIT,iostat=ierr)
!ESC!    if ( ierr < 0 ) then
!ESC!       write(IO_FID_LOG,*) '*** NM_RD_MSTRN_INIT is not specified. use default.'
!ESC!    elseif( ierr > 0 ) then
!ESC!       write(*,         *) 'xxx Not appropriate names in namelist NM_RD_MSTRN_INIT. STOP.'
!ESC!       write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist NM_RD_MSTRN_INIT. STOP.'
!ESC!       call ADM_proc_stop
!ESC!    endif
!ESC!    write(IO_FID_LOG,nml=NM_RD_MSTRN_INIT)

    KTPCL = KCPCL + KAPCL

    allocate( nq2IP(HYDRO_MAX) )
    allocate( IP2nq(KCPCL)     )
    nq2IP(:) = -999
    IP2nq(:) = -999

    allocate( opt_cloud_type_onoff(KCPCL) )
    opt_cloud_type_onoff(:) = .false.

    if ( KCPCL == 2 ) then ! liquid and ice

       KMODE   = 8

       I_QC_RD = 1
       I_QI_RD = 2

       nq2IP(I_QC)    = I_QC_RD
       nq2IP(I_QI)    = I_QI_RD

       IP2nq(I_QC_RD) = I_QC
       IP2nq(I_QI_RD) = I_QI

       opt_cloud_type_onoff(I_QC_RD) = opt_cloud_onoff
       opt_cloud_type_onoff(I_QI_RD) = opt_ice_onoff

    elseif( KCPCL == 5 ) then ! cloud, rain, ice, snow ,graupel

       KMODE   = 20

       I_QC_RD = 1
       I_QR_RD = 2
       I_QI_RD = 3
       I_QS_RD = 4
       I_QG_RD = 5

       nq2IP(I_QC)    = I_QC_RD
       nq2IP(I_QR)    = I_QR_RD
       nq2IP(I_QI)    = I_QI_RD
       nq2IP(I_QS)    = I_QS_RD
       nq2IP(I_QG)    = I_QG_RD

       IP2nq(I_QC_RD) = I_QC
       IP2nq(I_QR_RD) = I_QR
       IP2nq(I_QI_RD) = I_QI
       IP2nq(I_QS_RD) = I_QS
       IP2nq(I_QG_RD) = I_QG

       opt_cloud_type_onoff(I_QC_RD) = opt_cloud_onoff
       opt_cloud_type_onoff(I_QR_RD) = opt_rain_onoff
       opt_cloud_type_onoff(I_QI_RD) = opt_ice_onoff
       opt_cloud_type_onoff(I_QS_RD) = opt_snow_onoff
       opt_cloud_type_onoff(I_QG_RD) = opt_graupel_onoff

    endif

    write(IO_FID_LOG,*) '*** KCPCL=',KCPCL

    do IP = 1, KCPCL
       write(IO_FID_LOG,'(1x,2(A,I5))') '*** IP2nq, IP:', IP, ' => nq:', IP2nq(IP)
    enddo

    do nq = I_QC, I_QG
       write(IO_FID_LOG,'(1x,2(A,I5))') '*** nq2IP, nq:', nq, ' => IP:', nq2IP(nq)
    enddo

    do IP = 1, KCPCL
       write(IO_FID_LOG,'(1x,A,I5,A,L2)') &
       '*** opt_cloud_type_onoff, IP:', IP, ' = ', opt_cloud_type_onoff(IP)
    enddo



    allocate( Q    (KMODE,KTPCL,KDA*2+2,KWNB) )
    allocate( RMODE(0:KMODE+1,KTPCL,2) )
    allocate( NMODE(KTPCL) )

    call OPPARM2    !--- reference solar constant

    FSTOT = 0.0_RP
    do IW = 1, NWNB
       FSTOT = FSTOT + FSOL(IW)
    enddo
    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '*** Reference solar constant : FSTOT = ', FSTOT

    AMUI(1) = sqrt( 3.0_RP )
    AMUA(1) = 1.0_RP / AMUI(1)
    AMX (1) = 1.0_RP
    WMP (1) = sqrt( AMUA(1) )
    WMM (1) = 1.0_RP / WMP(1)
    WMX (1) = WMP(1) / AMX(1)
    PA0 (1) = WMM(1)**2 / 2.0_RP
    PA1 (1) = 3.0_RP / 2.0_RP * ( WMM(1)*AMUA(1) )**2
    P00 (1) = WMM(1) / 2.0_RP
    P01F(1) = 3.0_RP / 2.0_RP * ( WMM(1)*AMUA(1) )

    AMUI(2) = 1.66_RP
    AMUA(2) = 1.0_RP / AMUI(2)
    AMX (2) = AMUA(2) * 2.0_RP
    WMP (2) = sqrt( AMUA(2) )
    WMM (2) = 1.0_RP / WMP(2)
    WMX (2) = WMP(2) / AMX(2)
    PA0 (2) = WMM(2)**2 / 2.0_RP
    PA1 (2) = 3.0_RP / 2.0_RP * ( WMM(2)*AMUA(2) )**2
    P00 (2) = WMM(2) / 2.0_RP
    P01F(2) = 3.0_RP / 2.0_RP * ( WMM(2)*AMUA(2) )

    EPS(1) = 1.E-5_RP
    EPS(2) = 1.E-10_RP
    EPST   = 1.E-4_RP

    return
  end subroutine rd_mstrnx_init

  !-----------------------------------------------------------------------------
  ! read radiation parameter
  subroutine OPPARM2
!ESC!    use mod_adm, only: &
!ESC!       ADM_proc_stop
    implicit none

    real(RP) :: SR(KWNB,KSFC) !! Parameter for surface
    real(RP) :: QZ(KMODE)
    integer  :: NSFC, NPCL, NDA, NPLK, NFLG, NAPLNK, NACFC
    integer  :: NWNBP, NPCLP
    integer  :: IW, IPG, ITG, IFLG, IPL, ISFC, IDA, IPCL, IMOL, ICH, ICFC
    integer  :: IMODE, MMDZ, NMDZ

    character(len=1) :: CDUM

    integer  :: ifid, ierr
    !---------------------------------------------------------------------------

    write(IO_FID_LOG,*) '*** READING PARA ', trim(PARA)

    ifid = IO_get_available_fid()
    open(ifid,file=trim(PARA),form='formatted',status='old', iostat=ierr)
    if ( ierr /= 0 ) then
       write(IO_FID_LOG,*) 'xxx FILE NOT FOUND : ', trim(PARA)
       call ADM_proc_stop
    endif

    read(ifid,*) NWNB, NDA, NPG, NTG, NFLG, NACFC

    if ( NWNB > KWNB .OR. NFLG > KFLG .OR.&
         NPG /= KPG .OR. NTG /= KTG .OR. &
         NACFC /= KCFC ) then
       write(IO_FID_LOG,*) 'xxx INVALID RADIATION WORK SIZE (1) ', &
            NWNB, NFLG, NPG, NTG, NACFC
       write(IO_FID_LOG,*) 'xxx FOR                         ', &
            KWNB, KFLG, KPG, KTG, KCFC
       call ADM_proc_stop
    endif

    ! BAND BOUNDARIES
    read(ifid,'(A)') CDUM
    read(ifid,*) ( WNBND(IW), IW=1,NWNB+1 )

    ! LOG(PRESSURE) GRIDS
    read(ifid,'(A)') CDUM
    read(ifid,*) ( PLG( IPG ), IPG=1,NPG )

    ! TEMPERATURE GRIDS
    read(ifid,'(A)') CDUM
    read(ifid,*) ( TG( ITG ), ITG=1,NTG )

    ! QUANTITIES FOR EACH BAND
    do IW = 1, NWNB
       !* OPTICAL PROPERTY FLAG
       read(ifid,'(A)') CDUM
       read(ifid,*) ( iflgb( ifLG,IW ), ifLG=1,NFLG )

       if ( KFLG > NFLG ) then
          if ( WNBND(IW) >= 10000.0_RP/0.8_RP ) then
             iflgb(KFLG,IW) = 1 !! PAR
          else
             iflgb(KFLG,IW) = 0
          endif
          ! ISCCP simulator
          if    ( WNBND(IW)   <= 1.493E+4_RP .AND. &
                  WNBND(IW+1) >  1.493E+4_RP ) then
             iflgb(KFLG,IW) = 2 !! 0.67 micron
          elseif( WNBND(IW)   <= 2.0E+4_RP .AND. &
                  WNBND(IW+1) >  2.0E+4_RP ) then
             iflgb(KFLG,IW) = 3 !! 0.5 micron
          elseif( ( WNBND(IW)   <=  901.0_RP .AND. &
                    WNBND(IW+1) >   901.0_RP       ) .AND.&
                  ( WNBND(IW)   >  1042.0_RP .OR.  & !! exclude 9.6 micron
                    WNBND(IW+1) <  1042.0_RP       ) ) then
             iflgb(KFLG,IW) = -1 !! 11.1 micron
          elseif( ( WNBND(IW)   <= 952.0_RP .AND. &
                    WNBND(IW+1) >  952.0_RP       ) .AND.&
                  ( WNBND(IW)   >  1042.0_RP .OR. &!! exclude 9.6 micron
                    WNBND(IW+1) <  1042.0_RP      ) ) then
             iflgb(KFLG,IW) = -1 !! 10.5 micron
          endif

          write( IO_FID_LOG,'(A,I3,A,I3,A,2E12.5)' )&
               ' *** IW=', IW, &
               ' iflgb(KFLG,IW)=', iflgb(KFLG,IW),&
               ' WNBND(IW,IW+1)=', WNBND(IW), WNBND(IW+1)
       endif

       !* GAS ABSORPTION

       !* NUMBER OF CHANNELS
       read(ifid,'(A)') CDUM
       read(ifid,*) NCH(IW)
       if ( NCH(IW) > KCH ) then
          write(IO_FID_LOG,*) 'xxx INVALID RADIATION WORK KCH ', NCH(IW)
          call ADM_proc_stop
       endif

       !* WEIGHTS FOR CHANNELS
       WGT0(:,IW) = 0.0_RP
       if ( NCH(IW) > 0 ) then
          read(ifid,'(A)') CDUM
          read(ifid,*) ( WGT0( ICH,IW ), ICH=1,NCH(IW) )
       endif

       !* MOLECULES
       AKD(:,:,:,:,IW) = 0.0_RP
       if ( NCH(IW) > 0 ) then
          read(ifid,'(A)') CDUM
          read(ifid,*) NMOL(IW)
          do IMOL = 1, NMOL(IW)
             read(ifid,*) MID( IMOL,IW )
             if ( MID( IMOL,IW ) == 6 ) MID( IMOL,IW ) = 5
             if ( MID( IMOL,IW ) == 7 ) MID( IMOL,IW ) = 6
             !*** K-WIDTH
             do ITG = 1, NTG
                do IPG = 1, NPG
                   read(ifid,*) ( AKD( IPG,ITG,ICH,IMOL,IW ),ICH=1,NCH(IW) )
                enddo
             enddo
          enddo
       endif

       !*** SELF CONTINUUM
       SKD(:,:,:,IW) = 0.0_RP
       if ( iflgb( 5, IW ) > 0 ) then
          read(ifid,'(A)') CDUM
          do ITG = 1, NTG
             do IPG = 1, NPG
                read(ifid,*) ( SKD( IPG,ITG,ICH,IW ),ICH=1,NCH(IW) )
             enddo
          enddo
       endif

       !*** CFCs
       ACFC(:,IW) = 0.0_RP
       if ( iflgb( 7, IW ) > 0 ) then
          read(ifid,'(A)') CDUM
          read(ifid,*) ( ACFC( ICFC,IW ), ICFC=1,NACFC )
       endif
    enddo

    ! < VARDATA >

    RMODE=1.E+5_RP
    write(IO_FID_LOG,*) '*** READING VARDATA : ', trim(PARA)
    read(ifid,*) NPCL

    if ( NPCL /= KTPCL ) then
       write(IO_FID_LOG,*) 'xxx [OPPARM] Invalid number of particle type! STOP.'
       write(IO_FID_LOG,*) 'xxx KTPCL = KCPCL + KAPCL (input,requested)=', NPCL, KTPCL
       call ADM_proc_stop
    endif

    do IPCL = 1, NPCL
       read(ifid,'(A)') CDUM
       read(ifid,*) MMDZ, NMDZ
       ! [comment] 08/05/30, RMODE has unit[cm]
       read(ifid,*) ( RMODE( IMODE,IPCL,1 ), IMODE=1,NMDZ )
       NMODE( IPCL ) = NMDZ
       RMODE( 0, IPCL,1 ) = -1.E+5_RP
       RMODE( NMDZ+1,IPCL,1 ) = 1.E+5_RP

       do IMODE = 0, NMDZ
          if (  RMODE( IMODE+1,IPCL,1 ) > RMODE( IMODE, IPCL,1 )  ) then
             RMODE( IMODE,IPCL,2 ) = 1.0_RP/( RMODE( IMODE+1,IPCL,1 ) -RMODE( IMODE, IPCL,1 ) )
          endif
       enddo
    enddo

    ! < PARAPC >

    write(IO_FID_LOG,*) '*** READING PARAPC : ', trim(PARA)
    read(ifid,*) NWNBP,NSFC,NPCLP,NDA,NPLK,NAPLNK

    if ( NWNBP /= NWNB .OR. NSFC > KSFC .OR.&
         NPCL > KTPCL .OR. NPCLP /= NPCL .OR.&
         NDA > KDA .OR. NPLK > KPLK+1 .OR. &
         NAPLNK /= KPLNK ) then
       write(IO_FID_LOG,*) 'xxx INVALID RADIATION WORK SIZE (2) ',&
            NWNBP,NSFC,NPCL,NPCLP,NDA,NPLK,NAPLNK
       write(IO_FID_LOG,*) 'xxx FOR                         ',&
            NWNB,KSFC,KTPCL,NPCLP,KDA,KPLK+1,KPLNK
       call ADM_proc_stop
    endif

    ! BAND BOUNDARIES
    read(ifid,'(A)') CDUM
    read(ifid,*) ( WNBND(IW), IW=1,NWNB+1 )

    ! QUANTITIES FOR EACH BAND
    do IW = 1, NWNB
       !* PLANK FUNCTIONS
       read(ifid,'(A)') CDUM
       read(ifid,*) ( APLNK( IPL,IW ), IPL=1,NAPLNK )
       !* SOLAR INSOLATION
       read(ifid,'(A)') CDUM
       read(ifid,*) FSOL(IW)
       !* SURFACE PROPERTIES
       read(ifid,'(A)') CDUM
       read(ifid,*) ( SR( IW,ISFC ), ISFC=1,KSFC )
       !* RAYLEIGH SCATTERING
       read(ifid,'(A)') CDUM
       read(ifid,*) RY(IW)
       read(ifid,'(A)') CDUM

       do IDA = 1, KDMAX
          !** MOMENTS FOR RAYLEIGH SCATTERING PHASE FUNCTION
          read(ifid,*) QMOL( IDA )
          !** MOMENTS FOR PARTICLE SCATTERING PHASE FUNCTION
          do IPCL = 1, NPCL
             read(ifid,*) ( QZ( IMODE ), IMODE=1,NMODE( IPCL ) )
             if ( IDA <= KDA*2+2 ) then ! two-stream: IDA<=4, four-stream: IDA<=6
                do IMODE = 1, NMODE( IPCL )
                   Q( IMODE,IPCL,IDA,IW ) = QZ( IMODE )
                enddo
             endif
          enddo
       enddo
    enddo

    close(ifid)

    write(IO_FID_LOG,*) '*** READING PARA DONE.'

    return
  end subroutine OPPARM2

  !-----------------------------------------------------------------------------
  subroutine get_cot_550nm( &
       ijdim,  &
       kdim,   &
       kmin,   &
       kmax,   &
       kcpcl,  &
       volume, &
       re_cm,  &
       zh,     &
       tauc,   &
       taucl,  &
       tauclk, &
       tauci,  &
       taucik  )
!ESC!    use mod_adm, only: &
!ESC!       ADM_proc_stop
!ESC!    use mod_runconf, only: &
!ESC!       I_QC, I_QR, I_QI, &
!ESC!       I_QS, I_QG, I_QH, &
!ESC!       HYDRO_MAX
    implicit none

    integer,  intent(in)  :: ijdim                              ! horizontal grid index
    integer,  intent(in)  :: kdim                               ! vertical grid index
    integer,  intent(in)  :: kmin                               ! atmosphere bottom layer
    integer,  intent(in)  :: kmax                               ! ...           top layer
    integer,  intent(in)  :: kcpcl                              ! number of cloud species
    real(RP), intent(in)  :: volume(ijdim,kcpcl,KCTYP,KMAX_RAD) ! [cm3/m3]
    real(RP), intent(in)  :: re_cm (ijdim,kcpcl,KCTYP,KMAX_RAD) ! [cm]
    real(RP), intent(in)  :: zh    (ijdim,kdim)                 ! half level altitude [m]
    real(RP), intent(out) :: tauc  (ijdim,HYDRO_MAX)            ! each components
    real(RP), intent(out) :: taucl (ijdim)
    real(RP), intent(out) :: tauclk(ijdim,kdim)
    real(RP), intent(out) :: tauci (ijdim)
    real(RP), intent(out) :: taucik(ijdim,kdim)

    real(RP) :: re   (ijdim,kcpcl,KMAX_RAD)
    real(RP) :: wtauc(ijdim,kdim)
    real(RP) :: wtaur(ijdim,kdim)
    real(RP) :: wtaui(ijdim,kdim)
    real(RP) :: wtaus(ijdim,kdim)
    real(RP) :: wtaug(ijdim,kdim)
    real(RP) :: dz   (ijdim)
    real(RP) :: cext_c, cext_r, cext_i, cext_s, cext_g
    real(RP) :: wt

    ! from Mie Calculation
    integer, parameter :: irmax = 20
    ! effective radius[m]
    real(RP), parameter :: re_lut(0:irmax+1)=(/&
         -1.E+5_RP, &
         1.E-6_RP, 1.39E-6_RP, 1.95E-6_RP, 2.71E-6_RP, 3.79E-6_RP, 5.28E-6_RP, 7.37E-6_RP, 10.28E-6_RP, &
         14.34E-6_RP, 20.00E-6_RP, 29.58E-6_RP, 43.73E-6_RP, 64.67E-6_RP, 95.64E-6_RP, 141.42E-6_RP, &
         209.13E-6_RP, 309.25E-6_RP, 457.31E-6_RP, 676.24E-6_RP, 1000.00E-6_RP, &
         1.E+5_RP /) ! maximum
    real(RP)            :: rre_lut(0:irmax)

    ! extinction coefficient per volume [cm2/cm3] @ 550nm
    ! cloud
    real(RP), save :: cext_c_lut(irmax)=(/ &
         17203.7782172935_RP, 12560.8897330232_RP, 8664.97131126605_RP, 6096.58605245161_RP, &
         4286.19571545505_RP, 3026.73088387572_RP, 2147.86267679752_RP, 1520.10259526331_RP, &
         1084.32846277628_RP, 774.052096156682_RP, 518.198215125097_RP, 348.116506755950_RP, &
         234.796516493421_RP, 158.338581955142_RP, 106.794126006500_RP, 72.0971796691872_RP, &
         48.7097412229701_RP, 32.9088928255357_RP, 22.2381034027089_RP, 15.0293555468837_RP  /)
    ! rain
    real(RP), save :: cext_r_lut(irmax)=(/ &
         19624.0191781032_RP, 13115.7959643712_RP, 8938.17095854434_RP, 6194.48391104529_RP, &
         4336.66459501264_RP, 3054.32348938201_RP, 2159.28115737525_RP, 1531.21959222514_RP, &
         1088.57858703676_RP, 774.800922626182_RP, 519.903786179962_RP, 349.360568281822_RP, &
         235.182085819604_RP, 158.510877659205_RP, 106.916379863338_RP, 72.1638480073628_RP, &
         48.7336857486411_RP, 32.9218057324168_RP, 22.2449551268957_RP, 15.0333744259211_RP  /)
    ! Parameters in Fu(1996) are not dependent on particle types.
    ! Ice particles type is specified as hexagonal plate.
    real(RP), save :: cext_i_lut(1:irmax)=(/ &
         14989.4683103716_RP, 10745.5176019459_RP, 7703.15171578400_RP, 5522.16734032808_RP, &
         3958.68253147633_RP, 2837.86535597540_RP, 2034.38384225318_RP, 1458.39111390763_RP, &
         1045.43480762711_RP, 749.363397962046_RP, 506.662046731739_RP, 342.536874885672_RP, &
         231.548320858464_RP, 156.525207264162_RP, 105.849128611131_RP, 71.5797680358731_RP, &
         48.4053412559754_RP, 32.7337895385909_RP, 22.1360070966238_RP, 14.9693273247235_RP  /)
    real(RP), save :: cext_s_lut(1:irmax)=(/ &
         14989.4683103716_RP, 10745.5176019459_RP, 7703.15171578400_RP, 5522.16734032808_RP, &
         3958.68253147633_RP, 2837.86535597540_RP, 2034.38384225318_RP, 1458.39111390763_RP, &
         1045.43480762711_RP, 749.363397962046_RP, 506.662046731739_RP, 342.536874885672_RP, &
         231.548320858464_RP, 156.525207264162_RP, 105.849128611131_RP, 71.5797680358731_RP, &
         48.4053412559754_RP, 32.7337895385909_RP, 22.1360070966238_RP, 14.9693273247235_RP  /)
    real(RP), save :: cext_g_lut(1:irmax)=(/ &
         14989.4683103716_RP, 10745.5176019459_RP, 7703.15171578400_RP, 5522.16734032808_RP, &
         3958.68253147633_RP, 2837.86535597540_RP, 2034.38384225318_RP, 1458.39111390763_RP, &
         1045.43480762711_RP, 749.363397962046_RP, 506.662046731739_RP, 342.536874885672_RP, &
         231.548320858464_RP, 156.525207264162_RP, 105.849128611131_RP, 71.5797680358731_RP, &
         48.4053412559754_RP, 32.7337895385909_RP, 22.1360070966238_RP, 14.9693273247235_RP  /)

    logical, save :: opt_ice_fu98 = .true.  ! consider non-spherical particles in detail

    namelist / NM_RD_MSTRNX_COT / &
       opt_ice_fu98

    logical, save :: flag_first = .true.

    integer  :: ij, k, kk, ir, ir1, ir2
    integer  :: ierr
    !---------------------------------------------------------------------------

    call PROF_rapstart('____RD_COT550NM')

    if ( flag_first ) then
       flag_first = .false.

!ESC!       rewind(IO_FID_CONF)
!ESC!       read(IO_FID_CONF,nml=NM_RD_MSTRNX_COT,iostat=ierr)
!ESC!       if ( ierr < 0 ) then
!ESC!          write(IO_FID_LOG,*) '*** NM_RD_MSTRNX_COT is not specified. use default.'
!ESC!       elseif( ierr > 0 ) then
!ESC!          write(*,         *) 'xxx Not appropriate names in namelist NM_RD_MSTRNX_COT. STOP.'
!ESC!          write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist NM_RD_MSTRNX_COT. STOP.'
!ESC!          call ADM_proc_stop
!ESC!       endif
!ESC!       write(IO_FID_LOG,nml=NM_RD_MSTRNX_COT)

       if( .not.opt_ice_fu98 )then
          ! These are based on Mie theory with equivalent effective size
          cext_i_lut(1:irmax)=(/ &
               21032.9691919490_RP, 13570.9971797305_RP, 9586.40361387653_RP, 6669.64408400313_RP, &
               4671.71534035534_RP, 3285.71300607513_RP, 2323.42181584591_RP, 1646.99953562399_RP, &
               1167.61026685840_RP, 828.358877418047_RP, 557.507438705250_RP, 375.893134098444_RP, &
               253.678920901014_RP, 171.137091320702_RP, 115.447093291965_RP, 77.8786440752527_RP, &
               52.5827790795412_RP, 35.5199687311447_RP, 24.0010208367202_RP, 16.2180136603262_RP  /)
          cext_s_lut(1:irmax)=(/ &
               21103.4637044187_RP, 13605.1500696995_RP, 9660.87839471139_RP, 6712.97470761242_RP, &
               4702.93498700681_RP, 3307.17088296459_RP, 2338.76437114754_RP, 1657.94884686484_RP, &
               1175.47429176329_RP, 833.659101767424_RP, 561.182535807979_RP, 378.351433999661_RP, &
               255.357301190004_RP, 172.263748680847_RP, 116.209101543876_RP, 78.3901270643259_RP, &
               52.9281179495022_RP, 35.7534264548752_RP, 24.1589287041989_RP, 16.3245657264746_RP  /)
          cext_g_lut(1:irmax)=(/ &
               20546.5974158939_RP, 13326.4116808347_RP, 9286.10101969333_RP, 6475.22253341383_RP, &
               4535.07828585620_RP, 3191.07215197754_RP, 2256.19565301852_RP, 1599.04999034174_RP, &
               1133.48300038223_RP, 804.934543746440_RP, 541.523034526722_RP, 365.138962393542_RP, &
               246.378155768906_RP, 166.220951701169_RP, 112.127420861654_RP, 75.6475100400376_RP, &
               51.0766625430279_RP, 34.5024489239109_RP, 23.3129977808318_RP, 15.7534988139508_RP  /)
       endif
    endif

    tauc  (:,:) = 0.0_RP
    taucl (:)   = 0.0_RP
    tauci (:)   = 0.0_RP
    tauclk(:,:) = 0.0_RP
    taucik(:,:) = 0.0_RP

    re(:,:,:) = re_cm(:,:,1,:) * 1.E-2_RP

    do ir=0, irmax
       rre_lut(ir) = 1.0_RP/(re_lut(ir+1) - re_lut(ir))
    enddo

    if ( kcpcl == 2 ) then ! Cloud and Cloud Ice
       do ir = 0, irmax
          ir1 = max(1    ,ir  )
          ir2 = min(irmax,ir+1)

          do k = kmin, kmax
             kk = kmax+1-k
             dz(:) = (zh(:,k+1)-zh(:,k))
             do ij = 1, ijdim
                ! linear interpolation
                if(  (re_lut(ir+1) > re(ij,I_QC_RD,kk)) .and. (re_lut(ir) <= re(ij,I_QC_RD,kk)) )then
                   wt             = (re(ij,I_QC_RD,kk)-re_lut(ir))*rre_lut(ir)     ! [no unit]
                   cext_c         = cext_c_lut(ir1)*(1.0_RP-wt) + cext_c_lut(ir2)*wt ! [cm2/cm3]
                   tauclk( ij,k ) = cext_c*1.E+2_RP &                ! [cm2/cm3] => [m2/m3]
                        *           volume(ij,I_QC_RD,1,kk)*1.E-6_RP& ! [cm3]=>[m3]
                        *           dz(ij)                       ! [m]
                   taucl( ij )  = taucl( ij ) + tauclk( ij,k )
                endif
                ! linear interpolation
                if( (re_lut(ir+1) > re(ij,I_QI_RD,kk)) .and. (re_lut(ir) <= re(ij,I_QI_RD,kk)) )then
                   wt             = (re(ij,I_QI_RD,kk)-re_lut(ir))*rre_lut(ir)
                   cext_i         = cext_g_lut(ir1)*(1.0_RP-wt) + cext_g_lut(ir2)*wt
                   taucik( ij,k ) = cext_i*1.E+2_RP*volume(ij,I_QI_RD,1,kk)*1.E-6_RP*dz(ij)
                   tauci( ij )  = tauci( ij ) + taucik( ij,k )
                endif
             enddo
          enddo
       enddo
       tauc(:,I_QC) = taucl(:)
       tauc(:,I_QI) = tauci(:)

    elseif( kcpcl == 5 ) then ! Cloud, Rain, Cloud Ice, Snow, Graupel

       wtauc(:,:) = 0.0_RP
       wtaur(:,:) = 0.0_RP
       wtaui(:,:) = 0.0_RP
       wtaus(:,:) = 0.0_RP
       wtaug(:,:) = 0.0_RP
       do ir=0,irmax
          ir1=max(1,ir)
          ir2=min(irmax,ir+1)
          do k = kmin, kmax
             kk = kmax+1-k
             dz(:) = (zh(:,k+1)-zh(:,k))
             do ij = 1, ijdim
                ! Cloud
                if( (re_lut(ir+1) > re(ij,I_QC_RD,kk)) .and. (re_lut(ir) <= re(ij,I_QC_RD,kk)) )then
                   wt          = (re(ij,I_QC_RD,kk)-re_lut(ir))*rre_lut(ir)
                   cext_c      = cext_c_lut(ir1)*(1.0_RP-wt) + cext_c_lut(ir2)*wt
                   wtauc(ij,k) = cext_c*1.E+2_RP*volume(ij,I_QC_RD,1,kk)*1.E-6_RP*dz(ij)
                endif
                ! Rain
                if( (re_lut(ir+1) > re(ij,I_QR_RD,kk)) .and. (re_lut(ir) <= re(ij,I_QR_RD,kk)) )then
                   wt          = (re(ij,I_QR_RD,kk)-re_lut(ir))*rre_lut(ir)
                   cext_r      = cext_r_lut(ir1)*(1.0_RP-wt) + cext_r_lut(ir2)*wt
                   wtaur(ij,k) = cext_r*1.E+2_RP*volume(ij,I_QR_RD,1,kk)*1.E-6_RP*dz(ij)
                endif
                ! Ice
                if( (re_lut(ir+1) > re(ij,I_QI_RD,kk)) .and. (re_lut(ir) <= re(ij,I_QI_RD,kk)) )then
                   wt          = (re(ij,I_QI_RD,kk)-re_lut(ir))*rre_lut(ir)
                   cext_i      = cext_i_lut(ir1)*(1.0_RP-wt) + cext_i_lut(ir2)*wt
                   wtaui(ij,k) = cext_i*1.E+2_RP*volume(ij,I_QI_RD,1,kk)*1.E-6_RP*dz(ij)
                endif
                ! Snow
                if( (re_lut(ir+1) > re(ij,I_QS_RD,kk)) .and. (re_lut(ir) <= re(ij,I_QS_RD,kk)) )then
                   wt          = (re(ij,I_QS_RD,kk)-re_lut(ir))*rre_lut(ir)
                   cext_s      = cext_s_lut(ir1)*(1.0_RP-wt) + cext_s_lut(ir2)*wt
                   wtaus(ij,k) = cext_s*1.E+2_RP*volume(ij,I_QS_RD,1,kk)*1.E-6_RP*dz(ij)
                endif
                ! Graupel
                if( (re_lut(ir+1) > re(ij,I_QG_RD,kk)) .and. (re_lut(ir) <= re(ij,I_QG_RD,kk)) )then
                   wt          = (re(ij,I_QG_RD,kk)-re_lut(ir))*rre_lut(ir)
                   cext_g      = cext_g_lut(ir1)*(1.0_RP-wt) + cext_g_lut(ir2)*wt
                   wtaug(ij,k) = cext_g*1.E+2_RP*volume(ij,I_QG_RD,1,kk)*1.E-6_RP*dz(ij)
                endif

             enddo
          enddo
       enddo

       tauclk(:,:) = wtauc(:,:) + wtaur(:,:)
       taucik(:,:) = wtaui(:,:) + wtaus(:,:) + wtaug(:,:)

       ! Optical Thickness[m2/m2]
       do k = kmin, kmax
          taucl(:)      = taucl(:)      + tauclk(:,k)
          tauci(:)      = tauci(:)      + taucik(:,k)
          tauc (:,I_QC) = tauc (:,I_QC) + wtauc(:,k)
          tauc (:,I_QR) = tauc (:,I_QR) + wtaur(:,k)
          tauc (:,I_QI) = tauc (:,I_QI) + wtaui(:,k)
          tauc (:,I_QS) = tauc (:,I_QS) + wtaus(:,k)
          tauc (:,I_QG) = tauc (:,I_QG) + wtaug(:,k)
       enddo

    endif

    call PROF_rapend  ('____RD_COT550NM')

    return
  end subroutine get_cot_550nm

  !-----------------------------------------------------------------------------
  subroutine rd_mstrnx( &
       ijdim,          &
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
       GDTM,           &
       GDPM,           &
       GDZM,           &
       GDT,            &
       GDP,            &
       GDQ,            &
       GDO3,           &
       cloud_volume,   &
       cumulus_volume, &
       RE_all,         &
       GDCFRC,         &
       CUMFRC,         &
       SINS,           &
       COSZ,           &
       GRALB,          &
       GDTG,           &
       OUTQLD,         &
       AERDFS,         &
       AERDFL,         &
       AERDFS_TRP,     &
       AERDFL_TRP,     &
       tbb_11um,       &
       aot_ext_vis,    & ! [Add] 2016/05/18 T.Seiki
       aot_abs_vis     ) ! [Add] 2016/05/18 T.Seiki
!ESC!    use mod_adm, only: &
!ESC!       kdim  => ADM_kall, &
!ESC!       ADM_kmin,          &
!ESC!       ADM_kmax
!ESC!    use mod_const, only: &
!ESC!       CONST_UNDEF, &
!ESC!       CONST_Rdry,  &
!ESC!       CONST_TEM00, &
!ESC!       CONST_Pstd
!ESC!    use mod_runconf, only: &
!ESC!       NTAU  => NTAU_ISCCP,  &
!ESC!       NPRES => NPRES_ISCCP, &
!ESC!       RAIN_TYPE,            &
!ESC!       AE_TYPE,              &
!ESC!       opt_offline_aerosol,  &
!ESC!       opt_2moment_water,    &
!ESC!       HYDRO_MAX
    use mod_satadjust, only: &
       SATURATION_qsat_liq
!ESC!    use mod_history, only: &
!ESC!       history_in
    implicit none

    integer,  intent(in)  :: ijdim
    real(RP), intent(out) :: RFLXSD        (ijdim,kdim,NCRF)        ! downward SW radiation flux
    real(RP), intent(out) :: RFLXSU        (ijdim,kdim,NCRF)        ! upward   SW radiation flux
    real(RP), intent(out) :: RFLXLD        (ijdim,kdim,NCRF)        ! downward LW radiation flux
    real(RP), intent(out) :: RFLXLU        (ijdim,kdim,NCRF)        ! upward   LW radiation flux
    real(RP), intent(out) :: DRFLXS        (ijdim,kdim,NCRF)        ! SW derivative
    real(RP), intent(out) :: DRFLXL        (ijdim,kdim,NCRF)        ! LW derivative
    real(RP), intent(out) :: RFSFCD        (ijdim,NRDIR,NRBND,NCRF) ! surface downward SW radiation flux
    real(RP), intent(out) :: dfq_isccp     (ijdim,NTAU,NPRES)       ! Occurrence of cloud category in ISCCP
    real(RP), intent(out) :: TAUC_ALL      (ijdim,HYDRO_MAX)        ! cloud optical thickness (550nm,all particles)
    real(RP), intent(out) :: TAUCL         (ijdim)                  ! cloud optical thickness (550nm,liquid)
    real(RP), intent(out) :: TAUCI         (ijdim)                  ! cloud optical thickness (550nm,ice)
    real(RP), intent(out) :: TAUCLK        (ijdim,kdim)             ! COT in each layer (liquid)
    real(RP), intent(out) :: TAUCIK        (ijdim,kdim)             ! COT in each layer (ice)
    real(RP), intent(in)  :: GDTM          (ijdim,kdim)             ! temperature (cell face)   [K]
    real(RP), intent(in)  :: GDPM          (ijdim,kdim)             ! pressure    (cell face)   [Pa]
    real(RP), intent(in)  :: GDZM          (ijdim,kdim)             ! altitude    (cell face)   [m]
    real(RP), intent(in)  :: GDT           (ijdim,kdim)             ! temperature (cell center) [K]
    real(RP), intent(in)  :: GDP           (ijdim,kdim)             ! pressure    (cell center) [Pa]
    real(RP), intent(in)  :: GDQ           (ijdim,kdim)             ! mixing ratio of water vapor [kg/kg]
    real(RP), intent(in)  :: GDO3          (ijdim,kdim)             ! mixing ratio of O3          [ppmv]
    real(RP), intent(in)  :: cloud_volume  (ijdim,kdim,HYDRO_MAX)   ! grid-resolved cloud volume concentration [m3/m3]
    real(RP), intent(in)  :: cumulus_volume(ijdim,kdim,1:2)         ! cumulus       cloud volume concentration [m3/m3]
    real(RP), intent(in)  :: RE_all        (ijdim,kdim,HYDRO_MAX)   ! Effective Radius [m]
    real(RP), intent(in)  :: GDCFRC        (ijdim,kdim)             ! area fraction of grid-resolved cloud
    real(RP), intent(in)  :: CUMFRC        (ijdim)                  ! area fraction of cumulus       cloud
    real(RP), intent(in)  :: SINS          (ijdim)                  ! solar insolation
    real(RP), intent(in)  :: COSZ          (ijdim)                  ! cosine(solar zenith angle)
    real(RP), intent(in)  :: GRALB         (ijdim,NRDIR,NRBND)      ! surface albedo
    real(RP), intent(in)  :: GDTG          (ijdim)                  ! surface temperature
    real(RP), intent(in)  :: OUTQLD        (ijdim,kdim,KAPCL)       ! mixing ratio of aerosol mass [kg/kg]
    real(RP), intent(out) :: AERDFS        (ijdim,kdim,NCRF)        ! aerosol radiative forcing (shortwave)
    real(RP), intent(out) :: AERDFL        (ijdim,kdim,NCRF)        ! aerosol radiative forcing (longwave)
    real(RP), intent(out) :: AERDFS_TRP    (ijdim,NCRF)             ! aerosol radiative forcing (shortwave, tropopause)
    real(RP), intent(out) :: AERDFL_TRP    (ijdim,NCRF)             ! aerosol radiative forcing (longwave,  tropopause)
    real(RP), intent(out) :: tbb_11um      (ijdim)                  ! black body temperature @ 11um
    ! [Add] 2016/05/18 T.Seiki
    real(RP), intent(out) :: aot_ext_vis   (ijdim,KAPCL)            ! aerosol optical thickness (absorption+scatter)
    real(RP), intent(out) :: aot_abs_vis   (ijdim,KAPCL)            ! aerosol optical thickness (only absorption)

    real(RP) :: FU     (ijdim,KMAX_RAD+1,NRBFLX,NCRF) ! upward   SW/LW radiation flux
    real(RP) :: FD     (ijdim,KMAX_RAD+1,NRBFLX,NCRF) ! downward SW/LW radiation flux
    real(RP) :: TAUV   (ijdim,KMAX_RAD,4)             ! mean 0.55 & 0.67 micron optical depth

    ! reversed vertical order
    real(RP) :: PB     (ijdim,KMAX_RAD+1)             ! pressure    (cell face)   [hPa]
    real(RP) :: PL     (ijdim,KMAX_RAD  )             ! pressure    (cell center) [hPa]
    real(RP) :: TB     (ijdim,KMAX_RAD+1)             ! temperature (cell face)   [K]
    real(RP) :: TL     (ijdim,KMAX_RAD  )             ! temperature (cell center) [K]
    real(RP) :: DZ     (ijdim,KMAX_RAD  )             ! layer thickness [m]
    real(RP) :: RH     (ijdim,KMAX_RAD  )             ! relative humidity

    ! KCPCL=1:liquid, 2:ice
    ! KCTYP=1:stratus(small), 2:cumulus(big)
    real(RP) :: CGAS   (ijdim,KMOL,       KMAX_RAD)
    real(RP) :: CCFC   (      KCFC,       KMAX_RAD)
    real(RP) :: CCPCL  (ijdim,KCPCL,KCTYP,KMAX_RAD)
    real(RP) :: CAPCL  (ijdim,KAPCL,      KMAX_RAD)
    real(RP) :: RCPCL  (ijdim,KCPCL,KCTYP,KMAX_RAD)   ! cloud effective radius [cm]
    real(RP) :: RAPCL  (ijdim,KAPCL,      KMAX_RAD)
    real(RP) :: FCLD   (ijdim,            KMAX_RAD)

    integer  :: KSTRT  (ijdim)                        ! layer number of tropopause
    integer  :: KCTOP

    ! for aerosol radiative forcing
    real(RP) :: RFLXSU0(ijdim,kdim,NCRF)              ! upward SW rad. (w/o aerosol)
    real(RP) :: RFLXSD0(ijdim,kdim,NCRF)              ! down.  SW rad. (w/o aerosol)
    real(RP) :: RFLXLU0(ijdim,kdim,NCRF)              ! upward LW rad. (w/o aerosol)
    real(RP) :: RFLXLD0(ijdim,kdim,NCRF)              ! down.  LW rad. (w/o aerosol)
    logical, parameter :: OVARRAoff = .false.        ! force not to calc aerosol optical thickness

    ! for ISCCP
    integer,  parameter :: top_height           = 1    ! definition cloud top, 1:adjusted by IR+VIS, 2:not adjusted, 3:adjusted by IR only
    integer,  parameter :: top_height_direction = 2    ! direction for finding atmosphere pressure level
    real(RP), parameter :: emsfc_lw             = 1.0_RP ! 10.5 micron emissivity of surface (fraction)

    real(RP) :: dem_s        (ijdim,KMAX_RAD  ) ! 10.5 micron longwave emissivity of stratiform clouds
    real(RP) :: dem_c        (ijdim,KMAX_RAD  ) ! 10.5 micron longwave emissivity of convective clouds
    real(RP) :: PB100        (ijdim,KMAX_RAD+1) ! pressure of full model levels [Pa]
    real(RP) :: PL100        (ijdim,KMAX_RAD  ) ! pressure of half model levels [Pa]
    real(RP) :: QL           (ijdim,KMAX_RAD  ) ! water vapor specific humidity (kg vapor/ kg air)
    real(RP) :: conv         (ijdim,KMAX_RAD  ) ! convective cloud cover
    real(RP) :: rsunlit      (ijdim)            ! solar insolation
    integer  :: sunlit       (ijdim)            ! 1 for day points, 0 for night time
    integer  :: seed         (ijdim)            ! seed value for random number generator

    real(RP) :: fq_isccp     (ijdim,NTAU,NPRES) ! the fraction of the model grid box covered by each of the 49 ISCCP D level cloud types
    real(RP) :: totalcldarea (ijdim)            ! the fraction of model grid box columns
    real(RP) :: meanptop     (ijdim)            ! mean cloud top pressure (mb), linear averaging in cloud top pressure
    real(RP) :: meantaucld   (ijdim)            ! mean optical thickness,       linear averaging in albedo performed
    real(RP) :: meanalbedocld(ijdim)            ! mean cloud albedo,            linear averaging in albedo performed
    real(RP) :: meantb       (ijdim)            ! mean all-sky   10.5 micron brightness temperature
    real(RP) :: meantbclr    (ijdim)            ! mean clear-sky 10.5 micron brightness temperature
    integer  :: m1, m2

    real(RP) :: AERDEN(KAPCL) ! aerosol density
    real(RP) :: AERRD (KAPCL) ! aerosol radius

    data AERDEN / 2.5E+3_RP, 1.43E+3_RP, 1.46E+3_RP,  1.5E+3_RP, 1.25E+3_RP, 1.769E+3_RP,  2.2E+3_RP /
    data AERRD  /  4.E-6_RP,    -1.0_RP,    -1.0_RP,    -1.0_RP,   4.E-8_RP,     -1.0_RP,   2.E-6_RP / ! in [m]

    real(RP) :: C2PPMA
    real(RP) :: GAM
    real(RP) :: QSAT(ijdim,kdim)

    integer  :: kmin, kmax, kdim_
    real(RP) :: UNDEF, Rdry, Tstd, Pstd

    integer  :: ij, k, kk, ic, ip, nq
    !---------------------------------------------------------------------------

if ( TIME_CSTEP == EX_STEP ) then
   call PROF_rapstart('_kernel_partR1')
endif

    kdim_ = kdim
    kmin  = ADM_kmin
    kmax  = ADM_kmax

    UNDEF = CONST_UNDEF
    Rdry  = CONST_Rdry
    Tstd  = CONST_TEM00
    Pstd  = CONST_Pstd

    !---< change index of vertical coordinate into downward direction half level >---

    call SATURATION_qsat_liq( ijdim, kdim, GDT(:,:), GDP(:,:), QSAT(:,:) )

    !$omp parallel default(none),private(ij,k,kk), &
    !$omp shared(ijdim,kmin,kmax,PB,TB,GDPM,GDTM)
    do k = kmin, kmax+1
       kk = kmax+1 - k + 1

       ! cell wall value
       !$omp do
       do ij = 1, ijdim
          PB(ij,kk) = GDPM(ij,k) * 1.E-2_RP ! [Pa->hPa]
          TB(ij,kk) = GDTM(ij,k)
       enddo
       !$omp end do
    enddo
    !$omp end parallel

    !$omp parallel default(none),private(ij,k,kk,nq,ip), &
    !$omp shared(ijdim,kmin,kmax,KCPCL,I_QC_RD,I_QI_RD,                                                           &
    !$omp        PL,TL,QL,DZ,RH,CGAS,CCFC,CCPCL,FCLD,CONV,GDP,GDT,GDQ,GDZM,QSAT,GDO3,cloud_volume,cumulus_volume, &
    !$omp        CUMFRC,GDCFRC,ARHMAX,Q2PPM,PPMCO2,PPMN2O,PPMCH4,PPMO2,PPMCFC,nq2ip,opt_cloud_type_onoff          )
    do k = kmin, kmax
       kk = kmax - k + 1

       ! cell center value
       !$omp do
       do ij = 1, ijdim
          PL(ij,kk) = GDP(ij,k) * 1.E-2_RP ! [Pa->hPa]
          TL(ij,kk) = GDT(ij,k)
          QL(ij,kk) = GDQ(ij,k)

          DZ(ij,kk) = GDZM(ij,k+1) - GDZM(ij,k)

          RH(ij,kk) = max( min( GDQ(ij,k)/QSAT(ij,k), ARHMAX ), 0.0_RP )

          ! mixing ratio of gases
          CGAS(ij,1,kk) = max(GDQ(ij,k),1.E-7_RP) * Q2PPM
          CGAS(ij,2,kk) = PPMCO2
          CGAS(ij,3,kk) = max(GDO3(ij,k),1.E-10_RP)
          CGAS(ij,4,kk) = PPMN2O
          CGAS(ij,5,kk) = PPMCH4
          CGAS(ij,6,kk) = PPMO2
       enddo
       !$omp end do nowait

       ! mixing ratio of CFCs
       !$omp do
       do nq = 1, KCFC
          CCFC(nq,kk) = PPMCFC(nq)
       enddo
       !$omp end do nowait

       !$omp do
       do ip = 1, KCPCL
       do ij = 1, ijdim
          CCPCL(ij,ip,1,kk) = 0.0_RP
          CCPCL(ij,ip,2,kk) = 0.0_RP
       enddo
       enddo
       !$omp end do

       ! mixing ratio of grid resolved clouds
       do nq = 1, HYDRO_MAX
          ip = nq2ip(nq)

          if ( ip > 0 ) then
             if ( opt_cloud_type_onoff(ip) ) then ! true=on, false=off
                !$omp do
                do ij = 1, ijdim
                   CCPCL(ij,ip,1,kk) = CCPCL(ij,ip,1,kk) + cloud_volume(ij,k,nq) * 1.E+6_RP ! [m3/m3=>cm3/m3]
                enddo
                !$omp end do
             endif
          endif
       enddo

       !$omp do
       do ij = 1, ijdim
          ! mixing ratio of grid unresolved (cumulus) clouds
          CCPCL(ij,I_QC_RD,2,kk) = cumulus_volume(ij,k,1) * 1.E+6_RP ! [m3/m3=>cm3/m3]
          CCPCL(ij,I_QI_RD,2,kk) = cumulus_volume(ij,k,2) * 1.E+6_RP ! [m3/m3=>cm3/m3]

          ! cloud fraction
          ! [note] consider cumulus cloud fraction even if the bluk microphysics schemes are selected.
          FCLD(ij,kk) = ( 1.0_RP-CUMFRC(ij) ) * GDCFRC(ij,k) + CUMFRC(ij) ! total
          CONV(ij,kk) = (        CUMFRC(ij) )                             ! cumulus only
       enddo
       !$omp end do
    enddo
    !$omp end parallel

    ! mixing ratio of aerosols
    if (      AE_TYPE == 'SPRINTARS'     &
         .OR. AE_TYPE == 'SPRINTARS_CRM' &
         .OR. opt_offline_aerosol        ) then

       !$omp parallel default(none),private(ij,k,kk,ip,C2PPMA), &
       !$omp shared(ijdim,kmin,kmax,CAPCL,AERDEN,OUTQLD,Pstd,Tstd,Rdry)
       do ip = 1, KAPCL
          C2PPMA = 1.E+6_RP * Pstd / ( Rdry * Tstd * AERDEN(ip) )

          do k = kmin, kmax
             kk = kmax - k + 1

             !$omp do
             do ij = 1, ijdim
                CAPCL(ij,ip,kk) = max( OUTQLD(ij,k,ip), 0.0_RP ) * C2PPMA
             enddo
             !$omp end do
          enddo
       enddo
       !$omp end parallel
    else
       !$omp parallel workshare
       CAPCL(:,:,:) = 0.0_RP
       !$omp end parallel workshare
    endif

    ! stratosphere
    !$omp parallel workshare
    KSTRT(:) = KMAX
    !$omp end parallel workshare
    KCTOP    = 1
    !$omp parallel default(none),private(ij,k,kk,GAM), &
    !$omp shared(ijdim,kmin,kmax,KSTRT,KCTOP,CGAS,CCPCL,FCLD, &
    !$omp        GDTM,GDZM,GDP,PSTRMN,PSTRMX,GCRSTR,QSTR,Q2PPM)
    do k = kmin, kmax
       kk = kmax - k + 1

       !$omp do
       do ij = 1, ijdim
          GAM = ( GDTM(ij,k+1)-GDTM(ij,k) ) / ( GDZM(ij,k+1)-GDZM(ij,k) )

          if (    ( GDP(ij,k) < PSTRMX .AND. GAM > GCRSTR ) &
               .OR. GDP(ij,k) < PSTRMN                      ) then
             KSTRT(ij) = min(k,KSTRT(ij))
          else
             KCTOP = min(kk,KCTOP)
          endif

          if ( k >= KSTRT(ij) ) then
             if( QSTR > 0.0_RP ) CGAS(ij,1,kk) = QSTR * Q2PPM
             CCPCL(ij,:,:,kk) = 0.0_RP
             FCLD (ij,    kk) = 0.0_RP
          endif
       enddo
       !$omp end do
    enddo
    !$omp end parallel

    ! cloud mode radius
    if ( OPT_RE_CLOUDMP ) then ! from cloud microphysics module, including aerosol indirect effect
       !$omp parallel default(none),private(ij,k,kk,ip,nq), &
       !$omp shared(ijdim,kmin,kmax,KCPCL,CCPCL,RCPCL,ip2nq,re_all,RCMIN,RCMAX)
       do k = kmin, kmax
          kk = kmax - k + 1

          do ip = 1, KCPCL
             nq = ip2nq(ip)

             !$omp do
             do ij = 1, ijdim
                if    ( re_all(ij,k,nq) < RCMIN+1.E-30_RP ) then ! too small to effect, neglect radiative forcing
                   CCPCL(ij,ip,1,kk) = 0.0_RP
                   RCPCL(ij,ip,1,kk) = RCMIN * 100.0_RP ! [m=>cm]
                elseif( re_all(ij,k,nq) > RCMAX-1.E-30_RP ) then ! too large to effect, neglect radiative forcing
                   CCPCL(ij,ip,1,kk) = 0.0_RP
                   RCPCL(ij,ip,1,kk) = RCMAX * 100.0_RP ! [m=>cm]
                else                                          ! suitable range
                   RCPCL(ij,ip,1,kk) = re_all(ij,k,nq) * 100.0_RP ! [m=>cm]
                endif
             enddo
             !$omp end do
          enddo
       enddo
       !$omp end parallel
    else
       !$omp parallel workshare
       RCPCL(:,:,1,:) = radius_stratus * 100.0_RP ! [m=>cm]
       !$omp end parallel workshare
    endif
    !$omp parallel workshare
    RCPCL(:,:,2,:) = radius_cumulus * 100.0_RP ! [m=>cm]
    !$omp end parallel workshare

    ! aerosol mode radius
    do ip = 1, KAPCL
       if ( AERRD(ip) < 0.0_RP ) then
          !$omp parallel workshare
          RAPCL(:,ip,:) = RH(:,:)
          !$omp end parallel workshare
       else
          !$omp parallel workshare
          RAPCL(:,ip,:) = AERRD(ip) * 100.0_RP ! [m=>cm]
          !$omp end parallel workshare
       endif
    enddo

    ! cloud optical thickness (only for stratiform or resolved cloud)
    call get_cot_550nm( ijdim,             & ! [IN]
                        kdim,              & ! [IN]
                        kmin,              & ! [IN]
                        kmax,              & ! [IN]
                        KCPCL,             & ! [IN]
                        CCPCL   (:,:,:,:), & ! [IN]
                        RCPCL   (:,:,:,:), & ! [IN]
                        GDZM    (:,:),     & ! [IN]
                        TAUC_ALL(:,:),     & ! [OUT]
                        TAUCL   (:),       & ! [OUT]
                        TAUCLK  (:,:),     & ! [OUT]
                        TAUCI   (:),       & ! [OUT]
                        TAUCIK  (:,:)      ) ! [OUT]

    ! convert rho(z) => standard atmosphere

    !$omp parallel do default(none),private(ij,k,ip), &
    !$omp shared(ijdim,KMAX_RAD,KCPCL,CCPCL,PL,TL,Pstd,Tstd)
    do k  = 1, KMAX_RAD
    do ip = 1, KCPCL
    do ij = 1, ijdim
       CCPCL(ij,ip,1,k) = CCPCL(ij,ip,1,k) * ( Pstd*1.E-2_RP / PL(ij,k) ) / ( Tstd / TL(ij,k) )
       CCPCL(ij,ip,2,k) = CCPCL(ij,ip,2,k) * ( Pstd*1.E-2_RP / PL(ij,k) ) / ( Tstd / TL(ij,k) )
    enddo
    enddo
    enddo
    !$omp end parallel do

    !---< main routine of radiation (Nakajima) >---

    if ( OAEROF ) then ! additional calculation for estimation of aerosol forcing

       call DTRN31( ijdim,             & ! [IN]
                    KMAX_RAD,          & ! [IN]
                    KCTOP,             & ! [IN]
                    SINS    (:),       & ! [IN]
                    COSZ    (:),       & ! [IN]
                    PB      (:,:),     & ! [IN]
                    PL      (:,:),     & ! [IN]
                    TB      (:,:),     & ! [IN]
                    TL      (:,:),     & ! [IN]
                    DZ      (:,:),     & ! [IN]
                    GDTG    (:),       & ! [IN]
                    CGAS    (:,:,:),   & ! [IN]
                    CCFC    (:,:),     & ! [IN]
                    GRALB   (:,:,:),   & ! [IN]
                    CCPCL   (:,:,:,:), & ! [IN]
                    CAPCL   (:,:,:),   & ! [IN]
                    FCLD    (:,:),     & ! [IN]
                    CUMFRC  (:),       & ! [IN]
                    RCPCL   (:,:,:,:), & ! [IN]
                    RAPCL   (:,:,:),   & ! [IN]
                    OVARRAoff,         & ! [IN] force false
                    FD      (:,:,:,:), & ! [OUT]
                    FU      (:,:,:,:), & ! [OUT]
                    RFSFCD  (:,:,:,:), & ! [OUT]
                    TAUV    (:,:,:),   & ! [OUT]
                    dem_s   (:,:),     & ! [OUT]
                    dem_c   (:,:),     & ! [OUT]
                    tbb_11um(:),       & ! [OUT]
                    aot_ext_vis(:,:),  & ! [OUT] [add] 2016/05/18 T.Seiki
                    aot_abs_vis(:,:)   ) ! [OUT] [add] 2016/05/18 T.Seiki

       !--- reverse vertical order

       !$omp parallel default(none),private(ij,k,kk,ic), &
       !$omp shared(ijdim,kmin,kmax,RFLXSU0,RFLXSD0,RFLXLU0,RFLXLD0,FU,FD)
       do ic = 1, NCRF
          do k  = kmin, kmax+1 ! half level
             kk = kmax+1 - k + 1

             !$omp do
             do ij = 1, ijdim
                RFLXSU0(ij,k,ic) = FU(ij,kk,1,ic)
                RFLXSD0(ij,k,ic) = FD(ij,kk,1,ic)
                RFLXLU0(ij,k,ic) = FU(ij,kk,2,ic)
                RFLXLD0(ij,k,ic) = FD(ij,kk,2,ic)
             enddo
             !$omp end do
          enddo

          !$omp do
          do ij = 1, ijdim
             RFLXSU0(ij,kmin-1,ic) = 0.0_RP
             RFLXSD0(ij,kmin-1,ic) = 0.0_RP
             RFLXLU0(ij,kmin-1,ic) = 0.0_RP
             RFLXLD0(ij,kmin-1,ic) = 0.0_RP
          enddo
          !$omp end do
       enddo
       !$omp end parallel

    endif

if ( TIME_CSTEP == EX_STEP ) then
   call PROF_rapend  ('_kernel_partR1')
   call PROF_rapstart('_kernel_partR2')
endif

    call DTRN31( ijdim,             & ! [IN]
                 KMAX_RAD,          & ! [IN]
                 KCTOP,             & ! [IN]
                 SINS    (:),       & ! [IN]
                 COSZ    (:),       & ! [IN]
                 PB      (:,:),     & ! [IN]
                 PL      (:,:),     & ! [IN]
                 TB      (:,:),     & ! [IN]
                 TL      (:,:),     & ! [IN]
                 DZ      (:,:),     & ! [IN]
                 GDTG    (:),       & ! [IN]
                 CGAS    (:,:,:),   & ! [IN]
                 CCFC    (:,:),     & ! [IN]
                 GRALB   (:,:,:),   & ! [IN]
                 CCPCL   (:,:,:,:), & ! [IN]
                 CAPCL   (:,:,:),   & ! [IN]
                 FCLD    (:,:),     & ! [IN]
                 CUMFRC  (:),       & ! [IN]
                 RCPCL   (:,:,:,:), & ! [IN]
                 RAPCL   (:,:,:),   & ! [IN]
                 OVARRA,            & ! [IN]
                 FD      (:,:,:,:), & ! [OUT]
                 FU      (:,:,:,:), & ! [OUT]
                 RFSFCD  (:,:,:,:), & ! [OUT]
                 TAUV    (:,:,:),   & ! [OUT]
                 dem_s   (:,:),     & ! [OUT]
                 dem_c   (:,:),     & ! [OUT]
                 tbb_11um(:),       & ! [OUT]
                 aot_ext_vis(:,:),  & ! [OUT] [add] 2016/05/18 T.Seiki
                 aot_abs_vis(:,:)   ) ! [OUT] [add] 2016/05/18 T.Seiki

if ( TIME_CSTEP == EX_STEP ) then
   call PROF_rapend  ('_kernel_partR2')
   call PROF_rapstart('_kernel_partR3')
endif

    !$omp parallel default(none),private(ij,k,kk,ic), &
    !$omp shared(ijdim,kmin,kmax,RFLXSU,RFLXSD,RFLXLU,RFLXLD,DRFLXS,DRFLXL,FU,FD)
    do ic = 1, NCRF
       ! reverse vertical order
       do k  = kmin, kmax+1 ! half level
          kk = kmax+1 - k + 1

          !$omp do
          do ij = 1, ijdim
             RFLXSU(ij,k,ic) = FU(ij,kk,1,ic)
             RFLXSD(ij,k,ic) = FD(ij,kk,1,ic)
             RFLXLU(ij,k,ic) = FU(ij,kk,2,ic)
             RFLXLD(ij,k,ic) = FD(ij,kk,2,ic)
             DRFLXS(ij,k,ic) = FU(ij,kk,3,ic) - FU(ij,kk,1,ic)
             DRFLXL(ij,k,ic) = FU(ij,kk,4,ic) - FU(ij,kk,2,ic)
          enddo
          !$omp end do
       enddo

       !$omp do
       do ij = 1, ijdim
          RFLXSU(ij,kmin-1,ic) = 0.0_RP
          RFLXSD(ij,kmin-1,ic) = 0.0_RP
          RFLXLU(ij,kmin-1,ic) = 0.0_RP
          RFLXLD(ij,kmin-1,ic) = 0.0_RP
          DRFLXS(ij,kmin-1,ic) = 0.0_RP
          DRFLXL(ij,kmin-1,ic) = 0.0_RP
       enddo
       !$omp end do

    enddo
    !$omp end parallel

    if ( OAEROF ) then ! calc aerosol forcing

       !$omp parallel default(none),private(ij,k,ic), &
       !$omp shared(ijdim,kdim_,AERDFS,AERDFL,AERDFS_TRP,AERDFL_TRP, &
       !$omp        RFLXSU0,RFLXSD0,RFLXSU,RFLXSD,RFLXLU0,RFLXLD0,RFLXLU,RFLXLD,KSTRT)
       do ic = 1, NCRF
          !$omp do
          do k  = 1, kdim_
          do ij = 1, ijdim
             AERDFS(:,k,ic) = ( RFLXSU0(:,k,ic) - RFLXSD0(:,k,ic) ) &
                            - ( RFLXSU (:,k,ic) - RFLXSD (:,k,ic) )
             AERDFL(:,k,ic) = ( RFLXLU0(:,k,ic) - RFLXLD0(:,k,ic) ) &
                            - ( RFLXLU (:,k,ic) - RFLXLD (:,k,ic) )
          enddo
          enddo
          !$omp end do

          !$omp do
          do ij = 1, ijdim
             AERDFS_TRP(ij,ic) = AERDFS(ij,KSTRT(ij),ic)
             AERDFL_TRP(ij,ic) = AERDFL(ij,KSTRT(ij),ic)
          enddo
          !$omp end do
       enddo
       !$omp end parallel
    else
       !$omp parallel workshare
       AERDFS    (:,:,:) = UNDEF
       AERDFL    (:,:,:) = UNDEF
       AERDFS_TRP(:,:)   = UNDEF
       AERDFL_TRP(:,:)   = UNDEF
       !$omp end parallel workshare
    endif

!ESC!    call PROF_rapstart('____RD_ISCCP')
!ESC!
!ESC!    !---< ISCCP simulator >---
!ESC!
!ESC!    if ( isccp_version /= 'NONE' ) then
!ESC!
!ESC!       ! sunlit = 1 if sins > 0
!ESC!       ! sunlit = 0 if sins <= 0
!ESC!       rsunlit(:) = 0.5_RP - sign(0.5_RP,-sins(:)+1.E-8_RP)
!ESC!
!ESC!       if ( calc_nighttime ) then ! prevent ignoring nighttime region
!ESC!          sunlit(:) = 1
!ESC!       else
!ESC!          sunlit(:) = nint( rsunlit(:) )
!ESC!       endif
!ESC!
!ESC!       seed(:) = int( ( GDP(:,kmin+1) - real(int(GDP(:,kmin+1)),kind=RP) ) * 100.0_RP ) + 1
!ESC!
!ESC!       if ( isccp_version == 'old' ) then
!ESC!
!ESC!          call ISCCP_CLOUD_TYPES_old( debug,        & ! in
!ESC!                                      debugcol,     & ! in
!ESC!                                      ijdim,        & ! in
!ESC!                                      sunlit,       & ! in
!ESC!                                      KMAX_RAD,     & ! in
!ESC!                                      ncol,         & ! in
!ESC!                                      seed,         & ! in
!ESC!                                      PL,           & ! inout PL <= PL*100
!ESC!                                      PB,           & ! inout PB <= PB*100
!ESC!                                      QL,           & ! in
!ESC!                                      FCLD,         & ! in
!ESC!                                      CONV,         & ! in
!ESC!                                      TAUV(1,1,3),  & ! in
!ESC!                                      TAUV(1,1,4),  & ! in
!ESC!                                      top_height,   & ! in
!ESC!                                      overlap,      & ! in
!ESC!                                      tautab,       & ! in
!ESC!                                      invtau,       & ! in
!ESC!                                      GDTG,         & ! in
!ESC!                                      emsfc_lw,     & ! in
!ESC!                                      TL,           & ! in
!ESC!                                      dem_s,        & ! in
!ESC!                                      dem_c,        & ! in
!ESC!                                      fq_isccp,     & ! out
!ESC!                                      totalcldarea, & ! out
!ESC!                                      meanptop,     & ! out
!ESC!                                      meantaucld,   & ! out
!ESC!                                      boxtau,       & ! out
!ESC!                                      boxptop,      & ! out
!ESC!                                      ncolprint     ) ! in
!ESC!
!ESC!       elseif(      isccp_version == '4_1' &
!ESC!               .OR. isccp_version == '4.1' ) then
!ESC!
!ESC!          PL100 = PL * 100.0_RP
!ESC!          PB100 = PB * 100.0_RP
!ESC!
!ESC!          call ISCCP_CLOUD_TYPES_icarus4_1( debug,                & ! in
!ESC!                                            debugcol,             & ! in
!ESC!                                            ijdim,                & ! in
!ESC!                                            sunlit,               & ! in
!ESC!                                            KMAX_RAD,             & ! in
!ESC!                                            ncol,                 & ! in
!ESC!                                            seed,                 & ! in
!ESC!                                            PL100,                & ! inout  (Pa)
!ESC!                                            PB100,                & ! inout  (Pa)
!ESC!                                            QL,                   & ! in
!ESC!                                            FCLD,                 & ! in
!ESC!                                            CONV,                 & ! in
!ESC!                                            TAUV(:,:,3),          & ! in
!ESC!                                            TAUV(:,:,4),          & ! in
!ESC!                                            top_height,           & ! in
!ESC!                                            top_height_direction, & ! top height direction
!ESC!                                            overlap,              & ! in
!ESC!                                            frac_out(:,:,:),      &
!ESC!                                            GDTG,                 & ! in
!ESC!                                            emsfc_lw,             & ! in
!ESC!                                            TL,                   & ! in
!ESC!                                            dem_s,                & ! in
!ESC!                                            dem_c,                & ! in
!ESC!                                            fq_isccp,             & ! out
!ESC!                                            totalcldarea,         & ! out
!ESC!                                            meanptop,             & ! out
!ESC!                                            meantaucld,           & ! out
!ESC!                                            meanalbedocld,        & ! out
!ESC!                                            meantb,               & !?
!ESC!                                            meantbclr,            &!?
!ESC!                                            boxtau,               & ! out
!ESC!                                            boxptop               ) ! out
!ESC!       else
!ESC!
!ESC!          !$omp parallel workshare
!ESC!          fq_isccp(:,:,:) = 0.0_RP
!ESC!          !$omp end parallel workshare
!ESC!
!ESC!       endif
!ESC!
!ESC!       if ( DFQ_TYPE == 'FQ' ) then
!ESC!
!ESC!          !$omp parallel do default(none),private(ij,m1,m2), &
!ESC!          !$omp shared(ijdim,dfq_isccp,fq_isccp)
!ESC!          do m2 = 1, NPRES
!ESC!          do m1 = 1, NTAU
!ESC!          do ij = 1, ijdim
!ESC!             dfq_isccp(ij,m1,m2) = fq_isccp(ij,m1,m2)
!ESC!          enddo
!ESC!          enddo
!ESC!          enddo
!ESC!          !$omp end parallel do
!ESC!
!ESC!       elseif( DFQ_TYPE == 'UNDEF' ) then
!ESC!
!ESC!          !$omp parallel do default(none),private(ij,m1,m2), &
!ESC!          !$omp shared(ijdim,dfq_isccp,fq_isccp,rsunlit,UNDEF)
!ESC!          do m2 = 1, NPRES
!ESC!          do m1 = 1, NTAU
!ESC!          do ij = 1, ijdim
!ESC!             dfq_isccp(ij,m1,m2) = (        rsunlit(ij) ) * fq_isccp(ij,m1,m2) & ! if sins >  0
!ESC!                                 + ( 1.0_RP-rsunlit(ij) ) * UNDEF                ! if sins <= 0
!ESC!          enddo
!ESC!          enddo
!ESC!          enddo
!ESC!          !$omp end parallel do
!ESC!
!ESC!       elseif( DFQ_TYPE == '-FQ' ) then
!ESC!
!ESC!          !$omp parallel do default(none),private(ij,m1,m2), &
!ESC!          !$omp shared(ijdim,dfq_isccp,fq_isccp,rsunlit)
!ESC!          do m2 = 1, NPRES
!ESC!          do m1 = 1, NTAU
!ESC!          do ij = 1, ijdim
!ESC!             dfq_isccp(ij,m1,m2) = (        rsunlit(ij) ) * fq_isccp(ij,m1,m2) & ! if sins >  0
!ESC!                                 - ( 1.0_RP-rsunlit(ij) ) * fq_isccp(ij,m1,m2)   ! if sins <= 0
!ESC!          enddo
!ESC!          enddo
!ESC!          enddo
!ESC!          !$omp end parallel do
!ESC!
!ESC!       else
!ESC!
!ESC!          !$omp parallel workshare
!ESC!          dfq_isccp(:,:,:) = 0.0_RP
!ESC!          !$omp end parallel workshare
!ESC!
!ESC!       endif
!ESC!
!ESC!    endif ! isccp?
!ESC!
!ESC!    call PROF_rapend  ('____RD_ISCCP')

if ( TIME_CSTEP == EX_STEP ) then
   call PROF_rapend  ('_kernel_partR3')
endif

    return
  end subroutine rd_mstrnx

  !-----------------------------------------------------------------------------
  subroutine DTRN31( &
       ijdim,       &
       klev,        &
       KCTOP,       &
       SINS,        &
       AMS0,        &
       PB,          &
       PL,          &
       TB,          &
       TL,          &
       DZ,          &
       GTMP,        &
       CGAS,        &
       CCFC,        &
       GRALB,       &
       CCPCL,       &
       CAPCL,       &
       FCLD,        &
       FCUM,        &
       RCPCL,       &
       RAPCL,       &
       OVARRA,      &
       FD,          &
       FU,          &
       FDS,         &
       TAUVIS,      &
       dem_s,       &
       dem_c,       &
       tbb_11um,    &
       aot_ext_vis, & ! [Add] 2016/05/18 T.Seiki
       aot_abs_vis  ) ! [Add] 2016/05/18 T.Seiki
!ESC!    use mod_adm, only: &
!ESC!       ADM_proc_stop
!ESC!    use mod_const, only: &
!ESC!       CONST_PI,   &
!ESC!       CONST_Pstd, &
!ESC!       CONST_TEM00
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: klev                                      ! No. of levels (usually =KMAX_RAD)
    integer,  intent(in)  :: KCTOP                                     ! cloud top level
    real(RP), intent(in)  :: SINS       (ijdim)                        ! solar insolation
    real(RP), intent(in)  :: AMS0       (ijdim)                        ! COS(Solar zenith angle)
    real(RP), intent(in)  :: PB         (ijdim,KMAX_RAD+1)             ! Pressure interface (hPa)
    real(RP), intent(in)  :: PL         (ijdim,KMAX_RAD  )             ! Pressure layer (hPa)
    real(RP), intent(in)  :: TB         (ijdim,KMAX_RAD+1)             ! Temp. interface (K)
    real(RP), intent(in)  :: TL         (ijdim,KMAX_RAD  )             ! Temp. layer (K)
    real(RP), intent(in)  :: DZ         (ijdim,KMAX_RAD  )             ! dz [m]
    real(RP), intent(in)  :: GTMP       (ijdim)                        ! Ground temp. (K)
    real(RP), intent(in)  :: CGAS       (ijdim,KMOL,KMAX_RAD)          ! Gas concentration [PPMV] 1: H2O, 2: CO2, 3: O3, 4: N2O, 5: CH4, 6: O2
    real(RP), intent(in)  :: CCFC       (KCFC, KMAX_RAD )              ! CFC concentration
    real(RP), intent(in)  :: GRALB      (ijdim,NRDIR,NRBND)            ! surface albedo(1:VIS,2:NIR,3:IR)
    real(RP), intent(in)  :: CCPCL      (ijdim,KCPCL,KCTYP,KMAX_RAD)   ! cloud   concentration [ppmv]
    real(RP), intent(in)  :: CAPCL      (ijdim,KAPCL      ,KMAX_RAD)   ! aerosol concentration [ppmv]
    real(RP), intent(in)  :: FCLD       (ijdim,KMAX_RAD)               ! fractional cloudiness
    real(RP), intent(in)  :: FCUM       (ijdim)                        ! cumulus fractional area
    real(RP), intent(in)  :: RCPCL      (ijdim,KCPCL,KCTYP,KMAX_RAD)   ! cloud   mode Radius
    real(RP), intent(in)  :: RAPCL      (ijdim,KAPCL,      KMAX_RAD)   ! aerosol mode Radius
    logical,  intent(in)  :: OVARRA                                    ! use aerosol radius
    real(RP), intent(out) :: FD         (ijdim,KMAX_RAD+1,NRBFLX,NCRF) ! downward
    real(RP), intent(out) :: FU         (ijdim,KMAX_RAD+1,NRBFLX,NCRF) ! Upward
    real(RP), intent(out) :: FDS        (ijdim,NRDIR,NRBND,NCRF)       ! downward rad sfc.
    real(RP), intent(out) :: TAUVIS     (ijdim,KMAX_RAD,4)             ! 0.5 & 0.67 micron tau
    real(RP), intent(out) :: dem_s      (ijdim,KMAX_RAD)               ! 10.5 micron longwave emissivity
    real(RP), intent(out) :: dem_c      (ijdim,KMAX_RAD)               ! 10.5 micron longwave emissivity
    real(RP), intent(out) :: tbb_11um   (ijdim,1)                      ! black body temperature @ 11um ir
    real(RP), intent(out) :: aot_ext_vis(ijdim,KAPCL)                  ! aerosol optical thickness (absorption+scatter) [Add] 2016/05/18 T.Seiki
    real(RP), intent(out) :: aot_abs_vis(ijdim,KAPCL)                  ! aerosol optical thickness (only absorption)

    real(RP) :: aot_ext_wrk(ijdim,KAPCL)
    real(RP) :: aot_abs_wrk(ijdim,KAPCL)

    real(RP) :: DPRE  (ijdim,KMAX_RAD) ! d(pressure)_z
    real(RP) :: RDZ   (ijdim,KMAX_RAD) ! rho * dz [kg/m3*km]
    real(RP) :: AMS   (ijdim)
    real(RP) :: AMSINV(ijdim)
    real(RP) :: CUMWGT(ijdim,KCLD)

    ! variables for each broadband range
    integer  :: IRGN           ! 1:SW 2:LW
    integer  :: NPLK           ! 0 for SW, 2 for LW
    integer  :: IBND           ! 1:VIS 2:NIR 3:IR
    integer  :: IDIR           ! 1:DIRECT 2:DIFFUSE
    real(RP) :: FSOLW(ijdim)   ! solar flux with modification
    real(RP) :: WNAVR
    real(RP) :: WGT(KCH)

    real(RP) :: BP  (ijdim,KMAX_RAD+1)
    real(RP) :: BL  (ijdim,KMAX_RAD  )
    real(RP) :: BGND(ijdim,2)

    real(RP) :: TKD  (ijdim,KCH,    KMAX_RAD     ) ! tau (GAS)
    real(RP) :: TCON (ijdim,        KMAX_RAD     ) ! tau (CFC)
    real(RP) :: TAUP (ijdim,        KMAX_RAD,KCLD) ! tau (Cloud,Aerosol,Rayleigh scattering)
    real(RP) :: SCAP (ijdim,        KMAX_RAD,KCLD) ! scattering cross-section
    real(RP) :: G    (ijdim,KDA*2+1,KMAX_RAD,KCLD) ! asymmetry factor
    real(RP) :: QQ   (ijdim,KDA*2+2              ) ! 1:Extinction 2:Absorption 3:Asymmetry factor 4:Truncation factor

    integer  :: IDXCL(ijdim,KCPCL,KCTYP,KMAX_RAD)  ! index    for LUT (cloud)
    integer  :: IDXAE(ijdim,KAPCL,      KMAX_RAD)  ! index    for LUT (aerosol)
    real(RP) :: FXCL (ijdim,KCPCL,KCTYP,KMAX_RAD)  ! fraction for LUT (cloud)
    real(RP) :: FXAE (ijdim,KAPCL,      KMAX_RAD)  ! fraction for LUT (aerosol)
    real(RP) :: TAUCV(ijdim)                       ! for output of tau in specific band
    integer  :: IFCLD, M

    ! variables for each sub-channel
    real(RP) :: TAU (ijdim)
    real(RP) :: CPLK(ijdim,KPLK+1)
    real(RP) :: EXPD(ijdim,KMAX_RAD+1,KCLD)
    real(RP) :: RE  (ijdim,KMAX_RAD+1,KCLD) ! reflection
    real(RP) :: TE  (ijdim,KMAX_RAD+1,KCLD) ! transmission
    real(RP) :: SER (ijdim,KMAX_RAD+1,KCLD) ! upward   source
    real(RP) :: SET (ijdim,KMAX_RAD+1,KCLD) ! downward source
    real(RP) :: SERM(ijdim,KMAX_RAD+1,KCLD) ! upward   source, for NO-MAXR
    real(RP) :: SETM(ijdim,KMAX_RAD+1,KCLD) ! downward source, for NO-MAXR
    real(RP) :: EXPT(ijdim,KMAX_RAD  ,KCLD) ! cumulative direct solar transmission

    ! for maximal/random cloud overlapping method
    real(RP) :: B     (ijdim,    KMAX_RAD  ,4) ! cloud cover interaction matrix
    real(RP) :: EXPDMR(ijdim,    KMAX_RAD+1,2) ! cumulative direct solar transmission
    real(RP) :: TRNS0 (ijdim,    KMAX_RAD+1,2) ! direct solar transmission for solar source terms
    real(RP) :: RET   (ijdim,2,2,KMAX_RAD+1  ) ! reflection   matrices of sublayers(0,1)
    real(RP) :: RER   (ijdim,2,2,KMAX_RAD+1  ) ! reflection   matrices of sublayers(1,0)
    real(RP) :: TET   (ijdim,2,2,KMAX_RAD+1  ) ! transmission matrices of sublayers(0,1)
    real(RP) :: TER   (ijdim,2,2,KMAX_RAD+1  ) ! transmission matrices of sublayers(1,0)
    real(RP) :: SETT  (ijdim,2,  KMAX_RAD+1  ) ! downward source matrices
    real(RP) :: SERR  (ijdim,2,  KMAX_RAD+1  ) ! upward   source matrices

    ! Adding method
    real(RP) :: ALB  (ijdim,NRDIR)
    real(RP) :: FLXU (ijdim,KMAX_RAD+1)
    real(RP) :: FLXD (ijdim,KMAX_RAD+1)
    real(RP) :: FLXDD(ijdim,KMAX_RAD+1)
    real(RP) :: FDIR (ijdim)
    real(RP) :: FX, FY
    integer  :: IRBF

    real(RP) :: ir_11um(ijdim)    ! spectrum irradiance @ 11um
    integer  :: iw_11um           ! wave-length index of 11um
    real(RP) :: wl_11um, dwl_11um ! wave-length, wave-length range
    real(DP) :: c2_tbb, c1_tbb
    logical  :: flag_first_ir_11um

    real(DP), parameter :: kb =   1.38062E-23_DP ! Boltzmann constant [J/K]
    real(DP), parameter :: hp =     6.626E-34_DP ! Planck constant    [J*sec]
    real(DP), parameter :: Cs = 2.99792458E+8_DP  ! velocity of light  [m/s]

    real(RP) :: PI, Tstd, Pstd

    integer  :: ICL, ITG1, ICMX, ICX
    integer  :: ij, K, IW, ICH, IP, IDA
    !---------------------------------------------------------------------------

    call PROF_rapstart('____RD_DTRN31')

if ( TIME_CSTEP == EX_STEP ) then
   call PROF_rapstart('_kernel_partR2s1')
endif

    PI   = CONST_PI
    Tstd = CONST_TEM00
    Pstd = CONST_Pstd

    !$omp parallel default(none),private(ij,k,IRBF,IDIR,IBND,ICX), &
    !$omp shared(ijdim,KMAX_RAD,FD,FU,FDS)

    !$omp do
    do ICX  = 1, NCRF
    do IRBF = 1, NRBFLX
    do k    = 1, KMAX_RAD+1
    do ij   = 1, ijdim
       FD(ij,k,IRBF,ICX) = 0.0_RP
       FU(ij,k,IRBF,ICX) = 0.0_RP
    enddo
    enddo
    enddo
    enddo
    !$omp end do

    !$omp do
    do ICX  = 1, NCRF
    do IBND = 1, NRBND
    do IDIR = 1, NRDIR
    do ij   = 1, ijdim
       FDS(ij,IDIR,IBND,ICX) = 0.0_RP
    enddo
    enddo
    enddo
    enddo
    !$omp end do

    !$omp end parallel

    ! atmospheric layer thickness
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,klev,DPRE,RDZ,PB,PL,TL,DZ,Pstd,Tstd)
    do k  = 1, klev
    do ij = 1, ijdim
       DPRE(ij,k) = abs( PB(ij,k) - PB(ij,k+1) )
       ! normalized dz[km] from hydrostatic approx.
       !RDZ(ij,k) = DPRE(ij,k) * (Tstd*Rdry) / (Pstd*GRAV) * 1.E-3_RP                   ! diagnosed dz
       RDZ(ij,k) = (PL(ij,k)/(Pstd*1.E-2_RP)) / (TL(ij,k)/Tstd) * DZ(ij,k) * 1.E-3_RP ! model dz
    enddo
    enddo
    !$omp end parallel do

    !$omp parallel do default(none),private(ij), &
    !$omp shared(ijdim,AMS,AMSINV,CUMWGT,ir_11um,AMS0,FCUM)
    do ij = 1, ijdim
       ! cos(solar zenith angle)
       AMS   (ij) = max( AMS0(ij), 1.E-3_RP )
       AMSINV(ij) = 1.0_RP / AMS(ij)

       ! weight for cumulus
       CUMWGT(ij,1) = 1.0_RP - FCUM(ij)
       CUMWGT(ij,2) = 1.0_RP
       CUMWGT(ij,3) = FCUM(ij)
    enddo
    !$omp end parallel do

    ! B matrix for maximal/random overlap
    if ( MAX_RND_CLD_OVERLAP ) then
       call BCVR( ijdim,      & ! [IN]
                  klev,       & ! [IN]
                  FCLD(:,:),  & ! [IN]
                  B   (:,:,:) ) ! [OUT]
    endif

    ! prepare lookup table index and fraction
    call RMDIDX( ijdim,          & ! [IN]
                 klev,           & ! [IN]
                 KCPCL,          & ! [IN]
                 KCTYP,          & ! [IN]
                 KAPCL,          & ! [IN]
                 RMODE(:,:,:),   & ! [IN]
                 RCPCL(:,:,:,:), & ! [IN]
                 RAPCL(:,:,:),   & ! [IN]
                 IDXCL(:,:,:,:), & ! [OUT]
                 IDXAE(:,:,:),   & ! [OUT]
                 FXCL (:,:,:,:), & ! [OUT]
                 FXAE (:,:,:)    ) ! [OUT]

    flag_first_ir_11um = .true.

    !$omp parallel default(none),private(ij,IP), &
    !$omp shared(ijdim,ir_11um,aot_ext_vis,aot_abs_vis)

    !$omp do
    do ij = 1, ijdim
       ir_11um(ij) = 0.0_RP
    enddo
    !$omp end do nowait

    !$omp do
    do IP = 1, KAPCL
    do ij = 1, ijdim
       aot_ext_vis(ij,IP) = 0.0_RP ! [Add] 2016/05/18 T.Seiki
       aot_abs_vis(ij,IP) = 0.0_RP ! [Add] 2016/05/18 T.Seiki
    enddo
    enddo
    !$omp end do

    !$omp end parallel

if ( TIME_CSTEP == EX_STEP ) then
   call PROF_rapend  ('_kernel_partR2s1')
endif

    do IW = 1, NWNB

if ( TIME_CSTEP == EX_STEP ) then
   call PROF_rapstart('_kernel_partR2s2')
endif

       if (       iflgb(3,IW) >= 1 &
            .AND. iflgb(4,IW) >= 1 ) then
          write(IO_FID_LOG,*) 'xxx iflgb ERROR IN DTRN31'
          call ADM_proc_stop
       endif

       if ( iflgb(4,IW) > 0 ) then ! Solar band
          IRGN = 1 ! SW
          NPLK = 0 ! SW
       else
          IRGN = 2 ! LW
          NPLK = 2 ! LW
       endif

       if ( iflgb(KFLG,IW) > 0 ) then
          IBND = 1 ! visible
       elseif( iflgb(4,IW) > 0 ) then
          IBND = 2 ! near-IR
       else
          IBND = 3 ! IR
       endif

       !$omp parallel do default(none),private(ij), &
       !$omp shared(IW,ijdim,FSOLW,SINS,FSOL,FSTOT)
       do ij = 1, ijdim
          FSOLW(ij) = SINS(ij) * FSOL(IW) / FSTOT
       enddo
       !$omp end parallel do

       WNAVR = sqrt( WNBND(IW)*WNBND(IW+1) )

       ! channel weight
       do ICH = 1, NCH(IW)
          WGT(ICH) = WGT0(ICH,IW)
       enddo

       if ( iflgb(1,IW) <= 0 ) then ! no gas absorption (debug)
          WGT(:) = 0.0_RP
          WGT(1) = 1.0_RP
       endif

       ! PLANK FUNCTIONS
       call PLANKS( ijdim,       & ! [IN]
                    TB   (:,:),  & ! [IN]
                    TL   (:,:),  & ! [IN]
                    GTMP (:),    & ! [IN]
                    WNAVR,       & ! [IN]
                    APLNK(:,IW), & ! [IN]
                    NPLK,        & ! [IN]
                    klev,        & ! [IN]
                    BP   (:,:),  & ! [OUT]
                    BL   (:,:),  & ! [OUT]
                    BGND (:,1),  & ! [OUT]
                    BGND (:,2)   ) ! [OUT]

       call PROF_rapstart('____RD_DTRN_GAS')

       ! P-T-FITTING
       call PTFIT2( ijdim,             & ! [IN]
                    klev,              & ! [IN]
                    PL   (:,:),        & ! [IN]
                    TL   (:,:),        & ! [IN]
                    RDZ  (:,:),        & ! [IN]
                    NMOL (IW),         & ! [IN]
                    MID  (:,IW),       & ! [IN]
                    AKD  (:,:,:,:,IW), & ! [IN]
                    SKD  (:,:,:,IW),   & ! [IN]
                    CGAS (:,:,:),      & ! [IN]
                    iflgb(:,IW),       & ! [IN]
                    NCH  (IW),         & ! [IN]
                    TG   (:),          & ! [IN]
                    PLG  (:),          & ! [IN]
                    TKD  (:,:,:)       ) ! [OUT]

       ! CONTINUM & CFC
       call CNTCFC2( ijdim,       & ! [IN]
                     klev,        & ! [IN]
                     RDZ  (:,:),  & ! [IN]
                     ACFC (:,IW), & ! [IN]
                     CCFC (:,:),  & ! [IN]
                     iflgb(7,IW), & ! [IN]
                     TCON (:,:)   ) ! [OUT]

       call PROF_rapend  ('____RD_DTRN_GAS')
       call PROF_rapstart('____RD_DTRN_PTCL')

       do k = 1, klev

          ! Phase function, sigle scattering albedo, etc.
          ! ICL= 1:CLOUDY(MIXED), 2:CLEAR, 3:CUMULUS

          ! moments of scattering (aerosols)
          if ( OVARRA ) then
             call SCATAE( ijdim,                 & ! [IN]
                          CAPCL      (:,:,k),    & ! [IN]
                          RDZ        (:,k),      & ! [IN]
                          Q          (:,:,:,IW), & ! [IN]
                          iflgb      (2,IW),     & ! [IN]
                          IDXAE      (:,:,k),    & ! [IN]
                          FXAE       (:,:,k),    & ! [IN]
                          NMODE      (:),        & ! [IN]
                          aot_ext_wrk(:,:),      & ! [out] [Add] 2016/05/18 T.Seiki
                          aot_abs_wrk(:,:),      & ! [out] [Add] 2016/05/18 T.Seiki
                          QQ         (:,:)       ) ! [OUT]
          else
             !$omp parallel default(none),private(ij,IDA,IP), &
             !$omp shared(ijdim,QQ,aot_ext_wrk,aot_abs_wrk)

             !$omp do
             do IDA = 1, KDA*2+2
             do ij  = 1, ijdim
                QQ(ij,IDA) = 0.0_RP
             enddo
             enddo
             !$omp end do nowait

             !$omp do
             do IP = 1, KAPCL
             do ij = 1, ijdim
                aot_ext_wrk(ij,IP) = 0.0_RP ! [Add] 2016/05/18 T.Seiki
                aot_abs_wrk(ij,IP) = 0.0_RP ! [Add] 2016/05/18 T.Seiki
             enddo
             enddo
             !$omp end do

             !$omp end parallel
          endif

          ![Add] 2016/05/18 T.Seiki
          ! get aot in the wave-length rage [500nm < wave length < 600nm]
          if (       10000.0_RP / WNAVR > 0.5_RP &
               .AND. 10000.0_RP / WNAVR < 0.6_RP ) then
             !$omp parallel do default(none),private(ij,IP), &
             !$omp shared(ijdim,aot_ext_vis,aot_abs_vis,aot_ext_wrk,aot_abs_wrk)
             do IP = 1, KAPCL
             do ij = 1, ijdim
                aot_ext_vis(ij,IP) = aot_ext_vis(ij,IP) + aot_ext_wrk(ij,IP)
                aot_abs_vis(ij,IP) = aot_abs_vis(ij,IP) + aot_abs_wrk(ij,IP)
             enddo
             enddo
             !$omp end parallel do
          endif

          ! moments of scattering (Rayleigh)
          call SCATRY( ijdim,       & ! [IN]
                       RY   (IW),   & ! [IN]
                       DPRE (:,k),  & ! [IN]
                       Pstd,        & ! [IN]
                       QMOL (:),    & ! [IN]
                       iflgb(2,IW), & ! [IN]
                       QQ   (:,:)   ) ! [INOUT]

          ICMX = ICMAX
          if ( K < KCTOP ) ICMX = 1 ! no clouds in stratosphere

          do ICL = 1, ICMX

             IFCLD = 1 ! cloudy sky
             M     = 1 ! stratiform cloud
             if( ICL == 2 .OR. ICMX == 1 ) IFCLD = 0 ! clear sky
             if( ICL == 3 )                M     = 2 ! convective cloud

             ! moments of scattering (clouds)
             call SCATCL( ijdim,            & ! [IN]
                          CCPCL(:,:,M,k),   & ! [IN]
                          RDZ  (:,k),       & ! [IN]
                          Q    (:,:,:,IW),  & ! [IN]
                          iflgb(2,IW),      & ! [IN]
                          IFCLD,            & ! [IN]
                          IDXCL(:,:,M,k),   & ! [IN]
                          FXCL (:,:,M,k),   & ! [IN]
                          NMODE(:),         & ! [IN]
                          QQ   (:,:),       & ! [IN]
                          TAUP (:,k,ICL),   & ! [OUT]
                          SCAP (:,k,ICL),   & ! [OUT]
                          G    (:,:,k,ICL), & ! [OUT]
                          TAUCV(:)          ) ! [OUT]

             if ( iflgb(KFLG,IW) == 3 ) then ! 0.5 micron tau
                if ( ICL == 1 ) then
                   !$omp parallel do default(none),private(ij), &
                   !$omp shared(k,ijdim,TAUVIS,TAUCV)
                   do ij = 1, ijdim
                      TAUVIS(ij,k,1) = TAUCV(ij)
                   enddo
                   !$omp end parallel do
                endif
                if ( ICL == 3 ) then
                   !$omp parallel do default(none),private(ij), &
                   !$omp shared(k,ijdim,TAUVIS,TAUCV)
                   do ij = 1, ijdim
                      TAUVIS(ij,k,2) = TAUCV(ij)
                   enddo
                   !$omp end parallel do
                endif
             endif

             if ( iflgb(KFLG,IW) == 2 ) then !! 0.67 micron tau
                if ( ICL == 1 ) then
                   !$omp parallel do default(none),private(ij), &
                   !$omp shared(k,ijdim,TAUVIS,TAUCV)
                   do ij = 1, ijdim
                      TAUVIS(ij,k,3) = TAUCV(ij)
                   enddo
                   !$omp end parallel do
                endif
                if ( ICL == 3 ) then
                   !$omp parallel do default(none),private(ij), &
                   !$omp shared(k,ijdim,TAUVIS,TAUCV)
                   do ij = 1, ijdim
                      TAUVIS(ij,k,4) = TAUCV(ij)
                   enddo
                   !$omp end parallel do
                endif
             endif

          enddo ! ICL loop

       enddo ! K loop

       call PROF_rapend  ('____RD_DTRN_PTCL')

if ( TIME_CSTEP == EX_STEP ) then
   call PROF_rapend  ('_kernel_partR2s2')
endif

       do ICH = 1, max(NCH(IW),1)

if ( TIME_CSTEP == EX_STEP ) then
   call PROF_rapstart('_kernel_partR2s34')
endif
          call PROF_rapstart('____RD_DTRN_TWST')

          if ( MAX_RND_CLD_OVERLAP ) then
             !$omp parallel default(none),private(ij,k,icl), &
             !$omp shared(ijdim,klev,EXPD,EXPDMR)

             !$omp do
             do ICL = 2, ICMAX
             do k   = 1, klev+1
             do ij  = 1, ijdim
                EXPD(ij,k,ICL) = 1.0_RP
             enddo
             enddo
             enddo
             !$omp end do nowait

             !$omp do
             do ij  = 1, ijdim
                EXPDMR(ij,1,1) = 0.0_RP
                EXPDMR(ij,1,2) = 1.0_RP
             enddo
             !$omp end do

             !$omp end parallel
          else
             !$omp parallel do default(none),private(ij,icl), &
             !$omp shared(ijdim,EXPD)
             do ICL = 1, ICMAX
             do ij  = 1, ijdim
                EXPD(ij,1,ICL) = 1.0_RP
             enddo
             enddo
             !$omp end parallel do
          endif

          do k = 1, klev

             ! ICL= 1: CLOUDY(MIXED), 2: CLEAR, (3: CUMULUS)
             ICMX = ICMAX
             if ( K < KCTOP ) ICMX = 1 ! reduce calculation

             do ICL = 1, ICMX
                ! TOTAL TAU AND OMEGA

                !$omp parallel do default(none),private(ij), &
                !$omp shared(k,ICH,ICL,ijdim,TAU,TKD,TAUP,TCON)
                do ij  = 1, ijdim
                   TAU(ij) = TKD(ij,ICH,k) + TAUP(ij,k,ICL) + TCON(ij,k)
                enddo
                !$omp end parallel do

                ! PLANK EXPANSION COEFFICIENTS
                call PLKEXP( ijdim, BP(:,k:k+1), BL(:,k), TAU(:), NPLK, CPLK(:,:) )

                call TWST( ijdim,             & ! [IN]
                           IRGN,              & ! [IN]
                           AMS   (:),         & ! [IN]
                           AMSINV(:),         & ! [IN]
                           TAU   (:),         & ! [IN]
                           SCAP  (:,k,ICL),   & ! [IN]
                           G     (:,:,k,ICL), & ! [IN]
                           FSOLW (:),         & ! [IN]
                           CPLK  (:,:),       & ! [IN]
                           EXPD  (:,k,ICL),   & ! [IN]
                           EXPD  (:,k,1  ),   & ! [IN]
                           RE    (:,k,ICL),   & ! [OUT]
                           TE    (:,k,ICL),   & ! [OUT]
                           SER   (:,k,ICL),   & ! [OUT]
                           SET   (:,k,ICL),   & ! [OUT]
                           SERM  (:,k,ICL),   & ! [OUT]
                           SETM  (:,k,ICL),   & ! [OUT]
                           EXPT  (:,k,ICL)    ) ! [OUT]
             enddo ! ICL loop

             if ( MAX_RND_CLD_OVERLAP ) then
                !$omp parallel do default(none),private(ij,ICL), &
                !$omp shared(k,ijdim,ICMX,RE,TE,SER,SET,EXPT)
                do ICL = ICMX+1, ICMAX
                do ij  = 1, ijdim
                   RE  (ij,k,ICL) = RE  (ij,k,1)
                   TE  (ij,k,ICL) = TE  (ij,k,1)
                   SER (ij,k,ICL) = SER (ij,k,1)
                   SET (ij,k,ICL) = SET (ij,k,1)
                   EXPT(ij,k,ICL) = EXPT(ij,k,1)
                enddo
                enddo
                !$omp end parallel do
             else
                if ( ICMX >= 2 ) then
                   !$omp parallel do default(none),private(ij), &
                   !$omp shared(k,ijdim,RE,TE,SER,SET,EXPT,SERM,SETM,FCLD)
                   do ij  = 1, ijdim
                      RE  (ij,k,1) = FCLD(ij,k) * RE  (ij,k,1) + ( 1.0_RP-FCLD(ij,k) ) * RE  (ij,k,2)
                      TE  (ij,k,1) = FCLD(ij,k) * TE  (ij,k,1) + ( 1.0_RP-FCLD(ij,k) ) * TE  (ij,k,2)
                      SER (ij,k,1) = FCLD(ij,k) * SERM(ij,k,1) + ( 1.0_RP-FCLD(ij,k) ) * SERM(ij,k,2)
                      SET (ij,k,1) = FCLD(ij,k) * SETM(ij,k,1) + ( 1.0_RP-FCLD(ij,k) ) * SETM(ij,k,2)
                      EXPT(ij,k,1) = FCLD(ij,k) * EXPT(ij,k,1) + ( 1.0_RP-FCLD(ij,k) ) * EXPT(ij,k,2)
                   enddo
                   !$omp end parallel do
                endif

                !$omp parallel do default(none),private(ij,ICL), &
                !$omp shared(k,ijdim,ICMX,RE,TE,SER,SET,EXPT)
                do ICL = ICMX+1, ICMAX
                do ij  = 1, ijdim
                   RE  (ij,k,ICL) = RE  (ij,k,1)
                   TE  (ij,k,ICL) = TE  (ij,k,1)
                   SER (ij,k,ICL) = SER (ij,k,1)
                   SET (ij,k,ICL) = SET (ij,k,1)
                   EXPT(ij,k,ICL) = EXPT(ij,k,1)
                enddo
                enddo
                !$omp end parallel do

                ! semi-random overlap
                !$omp parallel default(none),private(ij,icl), &
                !$omp shared(k,ijdim,EXPD,EXPT)
                do ICL = 1, ICMAX
                   !$omp do
                   do ij = 1, ijdim
                      EXPD(ij,k+1,ICL) = EXPD(ij,k,ICL) * EXPT(ij,k,ICL)
                   enddo
                   !$omp end do
                enddo
                !$omp end parallel
             endif

          enddo ! k loop

          if ( MAX_RND_CLD_OVERLAP ) then

             !$omp parallel default(none),private(ij,k,icl), &
             !$omp shared(IW,ICH,ijdim,klev,EXPD,EXPT)
             do ICL = 2, ICMAX
             do k   = 1, klev
                !$omp do
                do ij = 1, ijdim
                   EXPD(ij,k+1,ICL) = EXPD(ij,k,ICL) * EXPT(ij,k,ICL)
                enddo
                !$omp end do
             enddo
             enddo
             !$omp end parallel

             !$omp parallel default(none),private(ij,k), &
             !$omp shared(IW,ICH,ijdim,klev,TRNS0,EXPDMR,B,EXPT)
             do k  = 1, klev
                !$omp do
                do ij = 1, ijdim
                   TRNS0 (ij,k  ,1) = EXPDMR(ij,k,1) * (        B(ij,k,3) ) &
                                    + EXPDMR(ij,k,2) * ( 1.0_RP-B(ij,k,1) )
                   TRNS0 (ij,k  ,2) = EXPDMR(ij,k,1) * ( 1.0_RP-B(ij,k,3) ) &
                                    + EXPDMR(ij,k,2) * (        B(ij,k,1) )

                   EXPDMR(ij,k+1,1) = TRNS0 (ij,k,1) * EXPT(ij,k,1)
                   EXPDMR(ij,k+1,2) = TRNS0 (ij,k,2) * EXPT(ij,k,2)
                enddo
                !$omp end do
             enddo
             !$omp end parallel

             ! make T,R,S matrixes for maximal/random overlap
             call RTS_MR( ijdim,          & ! [IN]
                          klev,           & ! [IN]
                          irgn,           & ! [IN]
                          RE   (:,:,:),   & ! [IN]
                          TE   (:,:,:),   & ! [IN]
                          SET  (:,:,:),   & ! [IN]
                          SER  (:,:,:),   & ! [IN]
                          TRNS0(:,:,:),   & ! [IN]
                          FCLD (:,:),     & ! [IN]
                          B    (:,:,:),   & ! [IN]
                          RET  (:,:,:,:), & ! [OUT]
                          RER  (:,:,:,:), & ! [OUT]
                          TET  (:,:,:,:), & ! [OUT]
                          TER  (:,:,:,:), & ! [OUT]
                          SETT (:,:,:),   & ! [OUT]
                          SERR (:,:,:)    ) ! [OUT]
          endif

          ! ISCCP simulator
          if ( iflgb(KFLG,IW) == -1 ) then ! 10.5 micron emissivity
             !$omp parallel do default(none),private(ij,k), &
             !$omp shared(IW,ICH,ijdim,klev,dem_s,dem_c,RE,TE)
             do k  = 1, klev
             do ij = 1, ijdim
                dem_s(ij,k) = 1.0_RP - RE(ij,k,1) - TE(ij,k,1)
                dem_c(ij,k) = 1.0_RP - RE(ij,k,3) - TE(ij,k,3)
             enddo
             enddo
             !$omp end parallel do
          endif

          call PROF_rapend  ('____RD_DTRN_TWST')

if ( TIME_CSTEP == EX_STEP ) then
   call PROF_rapend  ('_kernel_partR2s34')
   call PROF_rapstart('_kernel_partR2s5')
endif

          call PROF_rapstart('____RD_DTRN_ADDING')

          ! ICX =1 : all sky,    ICX =2 : clear sky
          ! ITG1=1 : normal,     ITG1=2 :  (no albedo for sw, T(sfc) + 1 [K] for lw)
          ! FY(FDIR): direct ray, FX(FLXD): direct ray + diffusion ray

          do ICL = 1, ICMAX

             ICX = ICL
             if( ICL == 3 ) ICX = 1

             do ITG1 = 1, 2

                if ( IRGN == 1 .AND. ITG1 == 2 ) then
                   !$omp parallel do default(none),private(ij,IDIR), &
                   !$omp shared(ijdim,ALB)
                   do IDIR = 1, NRDIR
                   do ij   = 1, ijdim
                      ALB(ij,IDIR) = 0.0_RP
                   enddo
                   enddo
                   !$omp end parallel do
                else
                   !$omp parallel do default(none),private(ij,IDIR), &
                   !$omp shared(IBND,ijdim,ALB,GRALB)
                   do IDIR = 1, NRDIR
                   do ij   = 1, ijdim
                      ALB(ij,IDIR) = GRALB(ij,IDIR,IBND)
                   enddo
                   enddo
                   !$omp end parallel do
                endif

                if (       MAX_RND_CLD_OVERLAP &
                     .AND. ICL == 1            ) then
                   ! calc maximal/random flux by the adding method (for all-sky only)
                   call ADDING_MR( ijdim,           & ! [IN]
                                   klev,            & ! [IN]
                                   irgn,            & ! [IN]
                                   RER   (:,:,:,:), & ! [INOUT]
                                   RET   (:,:,:,:), & ! [INOUT]
                                   TER   (:,:,:,:), & ! [INOUT]
                                   TET   (:,:,:,:), & ! [INOUT]
                                   SERR  (:,:,:),   & ! [INOUT]
                                   SETT  (:,:,:),   & ! [INOUT]
                                   FCLD  (:,:),     & ! [IN]
                                   EXPDMR(:,:,:),   & ! [IN]
                                   AMS   (:),       & ! [IN]
                                   ALB   (:,:),     & ! [IN]
                                   BGND  (:,ITG1),  & ! [IN]
                                   FSOLW (:),       & ! [IN]
                                   FLXU  (:,:),     & ! [OUT]
                                   FLXD  (:,:),     & ! [OUT]
                                   FDIR  (:)        ) ! [OUT]
                else
                   ! boundary condition
                   !$omp parallel do default(none),private(ij,k), &
                   !$omp shared(IW,ICH,ijdim,klev,FLXDD,AMS,FSOLW,EXPD,ICL)
                   do k  = 1, klev+1
                   do ij = 1, ijdim
                      FLXDD(ij,k) = AMS(ij) * FSOLW(ij) * EXPD(ij,k,ICL)
                   enddo
                   enddo
                   !$omp end parallel do

                   !$omp parallel do default(none),private(ij), &
                   !$omp shared(IW,ICH,ijdim,klev,FDIR,FLXDD)
                   do ij = 1, ijdim
                      FDIR(ij) = FLXDD(ij,klev+1)
                   enddo
                   !$omp end parallel do

                   if (       MAX_RND_CLD_OVERLAP &
                        .AND. IRGN == 1           ) then
                      !$omp parallel do default(none),private(ij,k), &
                      !$omp shared(IW,ICH,ijdim,klev,SER,SET,EXPD,ICL)
                      do k  = 1, klev
                      do ij = 1, ijdim
                         SER(ij,k,ICL) = SER(ij,k,ICL) * EXPD(ij,k,ICL)
                         SET(ij,k,ICL) = SET(ij,k,ICL) * EXPD(ij,k,ICL)
                      enddo
                      enddo
                      !$omp end parallel do
                   endif

                   !$omp parallel do default(none),private(ij), &
                   !$omp shared(IW,ICH,ijdim,klev,RE,TE,SER,SET,ALB,WMP,AMUI,FDIR,BGND,PI,irgn,ICL,ITG1)
                   do ij = 1, ijdim
                      RE (ij,klev+1,ICL) = ALB(ij,2)
                      TE (ij,klev+1,ICL) = 0.0_RP
                      SER(ij,klev+1,ICL) = WMP(irgn) * ( (        ALB(ij,1) ) * AMUI(irgn) * FDIR(ij)       & ! direct downward
                                                       + ( 1.0_RP-ALB(ij,2) ) * 2.0_RP * PI * BGND(ij,ITG1) ) ! surface emission
                      SET(ij,klev+1,ICL) = 0.0_RP
                   enddo
                   !$omp end parallel do

                   call ADDING( ijdim,          & ! [IN]
                                klev,           & ! [IN]
                                irgn,           & ! [IN]
                                RE   (:,:,ICL), & ! [IN]
                                TE   (:,:,ICL), & ! [IN]
                                SER  (:,:,ICL), & ! [IN]
                                SET  (:,:,ICL), & ! [IN]
                                FLXDD(:,:),     & ! [IN]
                                FLXD (:,:),     & ! [OUT]
                                FLXU (:,:)      ) ! [OUT]
                endif

                IRBF = IRGN + (ITG1-1)*2
                if ( ICX == 1 .OR. ICX == 2 ) then
                   !$omp parallel do default(none),private(ij,k), &
                   !$omp shared(IW,ICH,ijdim,klev,FD,FU,CUMWGT,WGT,FLXD,FLXU,IRBF,ICX,ICL)
                   do k  = 1, klev+1
                   do ij = 1, ijdim
                      FD(ij,k,IRBF,ICX) = FD(ij,k,IRBF,ICX) + CUMWGT(ij,ICL) * WGT(ICH) * FLXD(ij,k)
                      FU(ij,k,IRBF,ICX) = FU(ij,k,IRBF,ICX) + CUMWGT(ij,ICL) * WGT(ICH) * FLXU(ij,k)
                   enddo
                   enddo
                   !$omp end parallel do

                   if ( ITG1 == 1 ) then
                      !$omp parallel do default(none),private(ij,FX,FY), &
                      !$omp shared(IW,ICH,ijdim,FDS,CUMWGT,WGT,FLXD,FDIR,KMAX_RAD,IBND,ICX,ICL)
                      do ij = 1, ijdim
                         FX = CUMWGT(ij,ICL) * WGT(ICH) * FLXD(ij,KMAX_RAD+1)
                         FY = CUMWGT(ij,ICL) * WGT(ICH) * FDIR(ij)
                         FDS(ij,1,IBND,ICX) = FDS(ij,1,IBND,ICX) + FY
                         FDS(ij,2,IBND,ICX) = FDS(ij,2,IBND,ICX) + FX-FY
                      enddo
                      !$omp end parallel do
                   endif
                endif

                ! preparation for output to compare with IR satellite
                if (       flag_first_ir_11um   & ! if 11 micron
                     .AND. IFLGB(KFLG,IW) == -1 &
                     .AND. ICL            ==  1 &
                     .AND. ITG1           ==  1 ) then ! only 1 band

                   iw_11um = IW

                   !$omp parallel do default(none),private(ij), &
                   !$omp shared(IW,ICH,ijdim,ir_11um,CUMWGT,WGT,FLXU,ICL)
                   do ij = 1, ijdim
                      ir_11um(ij) = ir_11um(ij) + CUMWGT(ij,ICL) * WGT(ICH) * FLXU(ij,1) ! Here ir_11um is defined as B(lambda,T)*dlambda
                   enddo
                   !$omp end parallel do

                   if ( ICH == max(NCH(IW),1) ) then ! true until summation of all sub-channel
                      flag_first_ir_11um = .false.
                   endif
                endif

             enddo ! ITG1 loop

          enddo ! ICL loop

          call PROF_rapend  ('____RD_DTRN_ADDING')

if ( TIME_CSTEP == EX_STEP ) then
   call PROF_rapend  ('_kernel_partR2s5')
endif

       enddo ! ICH loop

    enddo ! IW loop

    ! check wavelength
    wl_11um  = 0.01_RP / sqrt( WNBND(iw_11um)*WNBND(iw_11um+1) ) ! [m]
    dwl_11um = 0.01_RP / WNBND(iw_11um) - 0.01_RP / WNBND(iw_11um+1)

    c1_tbb   = ( 2.0_DP * hp * Cs * Cs ) / ( wl_11um**5 )
    c2_tbb   = hp * Cs / ( kb * wl_11um )

    ! spectrum irradiance[W/m2/m] = pi*B(T,lambda)
    !$omp parallel do default(none),private(ij), &
    !$omp shared(ijdim,c1_tbb,c2_tbb,ir_11um,tbb_11um,dwl_11um,PI)
    do ij = 1, ijdim
       ir_11um (ij)   = max(ir_11um(ij),1.E-4_RP) / (dwl_11um*PI)
       tbb_11um(ij,1) = c2_tbb / ( log( 1.0_DP + c1_tbb/ir_11um(ij) ) )
    enddo
    !$omp end parallel do

    call PROF_rapend  ('____RD_DTRN31')

    return
  end subroutine DTRN31

  !-----------------------------------------------------------------------------
  ! k-distribution fitted in T,P
  subroutine PTFIT2( &
       ijdim, &
       klev,  &
       P,     &
       T,     &
       RDZ,   &
       NMOL,  &
       MID,   &
       AKD,   &
       SKD,   &
       CGAS,  &
       iflgb, &
       NCH,   &
       TG,    &
       PLG,   &
       TKD    )
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: klev
    real(RP), intent(in)  :: P    (ijdim,klev)       ! P-interpolation point
    real(RP), intent(in)  :: T    (ijdim,klev)       ! T-interpolation point
    real(RP), intent(in)  :: RDZ  (ijdim,klev)       ! layer thickness * density
    integer,  intent(in)  :: NMOL                    ! Number of molecular spiecies
    integer,  intent(in)  :: MID  (KMOL)             ! Molecular ID
    real(RP), intent(in)  :: AKD  (KPG,KTG,KCH,KMOL) ! Fitting for Abs.X-section
    real(RP), intent(in)  :: SKD  (KPG,KTG,KCH)      ! Fitting for Abs.X-section
    real(RP), intent(in)  :: CGAS (ijdim,KMOL,klev)  ! Gas concentration in PPMV 1:H2O, 2:CO2, 3:O3, 4:N2O, 5:CH4, 6:O2
    integer,  intent(in)  :: iflgb(KFLG)             ! Optical flag
    integer,  intent(in)  :: NCH                     ! Number of channels
    real(RP), intent(in)  :: TG   (KTG)              ! Grid of Temperature
    real(RP), intent(in)  :: PLG  (KPG)              ! Grid of Log10(Prs)
    real(RP), intent(out) :: TKD  (ijdim,KCH,klev)   ! Optical thickness

    integer  :: NNP (ijdim,klev)
    real(RP) :: DPRE(ijdim,klev)
    real(RP) :: PL  (ijdim,klev)
    real(RP) :: TL  (ijdim,klev)
    real(RP) :: GASM(ijdim,klev)

    real(RP) :: AKT1, AKT2, AKT3
    real(RP) :: ACFT

    real(RP), save :: TLG(KTG)
    real(RP), save :: RTLG32, RTLG21, RTG31
    logical,  save :: OFIRST = .true.

    integer  :: ij, k, im, ic, ipg
    !---------------------------------------------------------------------------

    if (ofirst) then
       ofirst = .false.

       TLG(:) = log10( TG(:) )
       RTLG32 = 1.0_RP / ( TLG(3)-TLG(2) )
       RTLG21 = 1.0_RP / ( TLG(2)-TLG(1) )
       RTG31  = 1.0_RP / ( TG(3) -TG(1)  )
    endif

    !$omp parallel do default(none),private(ij,k,ic), &
    !$omp shared(ijdim,klev,NCH,TKD)
    do k  = 1, klev
    do ic = 1, NCH
    do ij = 1, ijdim
       TKD(ij,ic,k) = 0.0_RP
    enddo
    enddo
    enddo
    !$omp end parallel do

    if ( iflgb(1) <= 0 ) then ! no gas absorption
       return
    endif

    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,klev,PL,TL,NNP,P,T)
    do k   = 1, klev
    do ij  = 1, ijdim
       PL (ij,k) = log10( P(ij,k) )
       TL (ij,k) = log10( T(ij,k) )
       NNP(ij,k) = KPG
    enddo
    enddo
    !$omp end parallel do

    !$omp parallel do default(none),private(ij,k,ipg), &
    !$omp shared(ijdim,klev,PL,PLG,NNP)
    do k   = 1, klev
    do ipg = KPG, 2, -1
    do ij  = 1, ijdim
       if( PL(ij,k) >= PLG(ipg) ) NNP(ij,k) = ipg
    enddo
    enddo
    enddo
    !$omp end parallel do

    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,klev,DPRE,PL,PLG,NNP)
    do k  = 1, klev
    do ij = 1, ijdim
       DPRE(ij,k) = ( PL(ij,k)-PLG(NNP(ij,k)-1) ) / ( PLG(NNP(ij,k))-PLG(NNP(ij,k)-1) )
    enddo
    enddo
    !$omp end parallel do

    !$omp parallel default(none),private(ij,k,ic,im,AKT1,AKT2,AKT3,ACFT), &
    !$omp shared(ijdim,klev,NCH,NMOL,GASM,CGAS,RDZ,MID,TKD,AKD,NNP,DPRE,T,TG,TL,TLG,RTLG32,RTLG21,RTG31)
    do im = 1, NMOL

       !$omp do
       do k  = 1, klev
       do ij = 1, ijdim
          GASM(ij,k) = CGAS(ij,MID(im),k) * RDZ(ij,k) * 1.E-1_RP ! [mol/cm2/km]
       enddo
       enddo
       !$omp end do

       !$omp do
       do k  = 1, klev
       do ic = 1, NCH
       do ij = 1, ijdim

          AKT1 = AKD(NNP(ij,k)-1,1,ic,im) &
               + ( AKD(NNP(ij,k),1,ic,im) - AKD(NNP(ij,k)-1,1,ic,im) ) * DPRE(ij,k)
          AKT2 = AKD(NNP(ij,k)-1,2,ic,im) &
               + ( AKD(NNP(ij,k),2,ic,im) - AKD(NNP(ij,k)-1,2,ic,im) ) * DPRE(ij,k)
          AKT3 = AKD(NNP(ij,k)-1,3,ic,im) &
               + ( AKD(NNP(ij,k),3,ic,im) - AKD(NNP(ij,k)-1,3,ic,im) ) * DPRE(ij,k)

          ACFT = AKT2                                           &
               + ( ( AKT3-AKT2 ) * ( T(ij,k)-TG(1) ) * RTLG32   &
                 + ( AKT2-AKT1 ) * ( TG(3)-T(ij,k) ) * RTLG21 ) &
               * ( TL(ij,k)-TLG(2) ) * RTG31

          TKD(ij,ic,k) = TKD(ij,ic,k) + 10**ACFT * GASM(ij,k)

       enddo
       enddo
       enddo
       !$omp end do

    enddo
    !$omp end parallel

    ! self broadening
    if ( iflgb(5) > 0 ) then
       !$omp parallel default(none),private(ij,k,ic,AKT1,AKT2,AKT3,ACFT), &
       !$omp shared(ijdim,klev,NCH,GASM,CGAS,RDZ,MID,TKD,SKD,NNP,DPRE,T,TG,TL,TLG,RTLG32,RTLG21,RTG31)

       !$omp do
       do k  = 1, klev
       do ij = 1, ijdim
          GASM(ij,k) = CGAS(ij,MID(1),k) * RDZ(ij,k) * 1.E-1_RP ! [mol/cm2/km]
       enddo
       enddo
       !$omp end do

       !$omp do
       do k  = 1, klev
       do ic = 1, NCH
       do ij = 1, ijdim

          AKT1 = SKD(NNP(ij,k)-1,1,ic) &
               + ( SKD(NNP(ij,k),1,ic) - SKD(NNP(ij,k)-1,1,ic) ) * DPRE(ij,k)
          AKT2 = SKD(NNP(ij,k)-1,2,ic) &
               + ( SKD(NNP(ij,k),2,ic) - SKD(NNP(ij,k)-1,2,ic) ) * DPRE(ij,k)
          AKT3 = SKD(NNP(ij,k)-1,3,ic) &
               + ( SKD(NNP(ij,k),3,ic) - SKD(NNP(ij,k)-1,3,ic) ) * DPRE(ij,k)

          ACFT = AKT2                                           &
               + ( ( AKT3-AKT2 ) * ( T(ij,k)-TG(1) ) * RTLG32   &
                 + ( AKT2-AKT1 ) * ( TG(3)-T(ij,k) ) * RTLG21 ) &
               * ( TL(ij,k)-TLG(2) ) * RTG31

          TKD(ij,ic,k) = TKD(ij,ic,k) &
                       + 10**ACFT * GASM(ij,k)**2 / ( GASM(ij,k) + RDZ(ij,k)*1.E+5_RP )

       enddo
       enddo
       enddo
       !$omp end do

       !$omp end parallel
    endif

    return
  end subroutine PTFIT2

  !-----------------------------------------------------------------------------
  ! optical thickness for cfc
  subroutine CNTCFC2( &
       ijdim, &
       klev,  &
       RDZ,   &
       ACFC,  &
       CCFC,  &
       IFLG,  &
       TCON   )
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: klev
    real(RP), intent(in)  :: RDZ (ijdim,klev) ! layer thickness * density
    real(RP), intent(in)  :: ACFC(KCFC)       ! Fitting for Abs.X-section
    real(RP), intent(in)  :: CCFC(KCFC,klev)  ! Gas concentration in PPMV
    integer,  intent(in)  :: IFLG             ! Optical flag
    real(RP), intent(out) :: TCON(ijdim,klev) ! Optical thickness

    real(RP) :: GASM

    integer  :: ij, k, icfc
    !---------------------------------------------------------------------------

    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,klev,TCON)
    do k  = 1, klev
    do ij = 1, ijdim
       TCON(ij,k) = 0.0_RP
    enddo
    enddo
    !$omp end parallel do

    if ( IFLG > 0 ) then
       !$omp parallel default(none),private(ij,k,icfc,GASM), &
       !$omp shared(ijdim,klev,TCON,CCFC,RDZ,ACFC)
       do icfc = 1, KCFC
          !$omp do
          do k    = 1, klev
          do ij   = 1, ijdim
             GASM = CCFC(icfc,k) * RDZ(ij,k) * 1.E-1_RP ! [mol/cm2/km]

             TCON(ij,k) = TCON(ij,k) + 10.0_RP**ACFC(icfc) * GASM
          enddo
          enddo
          !$omp end do
       enddo
       !$omp end parallel
    endif

    return
  end subroutine CNTCFC2

  !-----------------------------------------------------------------------------
  ! create idx table.
  subroutine RMDIDX( &
       ijdim, klev,         &
       KCPCL, KCTYP, KAPCL, &
       RMODE,               &
       RCPCL, RAPCL,        &
       IDXCL, IDXAE,        &
       FXCL,  FXAE          )
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: klev
    integer,  intent(in)  :: KCPCL
    integer,  intent(in)  :: KCTYP
    integer,  intent(in)  :: KAPCL
    real(RP), intent(in)  :: RMODE(0:KMODE+1,KTPCL,2)
    real(RP), intent(in)  :: RCPCL(ijdim,KCPCL,KCTYP,klev)
    real(RP), intent(in)  :: RAPCL(ijdim,KAPCL,      klev)
    integer,  intent(out) :: IDXCL(ijdim,KCPCL,KCTYP,klev)
    integer,  intent(out) :: IDXAE(ijdim,KAPCL,      klev)
    real(RP), intent(out) :: FXCL (ijdim,KCPCL,KCTYP,klev)
    real(RP), intent(out) :: FXAE (ijdim,KAPCL,      klev)

    real(RP) :: rpls, rmns
    integer  :: ij,k,IP,kt,im,idx
    !---------------------------------------------------------------------------

    !$omp parallel default(none),private(ij,k,IP,im,rpls,rmns), &
    !$omp shared(ijdim,klev,KCPCL,KAPCL,KMODE,IDXAE,RAPCL,RMODE)
    do k  = 1, klev
    do IP = 1, KAPCL ! aerosol species
       !$omp do
       do ij = 1, ijdim
          IDXAE(ij,IP,k) = 0
       enddo
       !$omp end do

       do im = 0, KMODE ! mode radius/ mode relative humidity
          rpls = RMODE(im+1,IP+KCPCL,1)
          rmns = RMODE(im  ,IP+KCPCL,1)

          !$omp do
          do ij = 1, ijdim
             if (       RAPCL(ij,IP,k) <  rpls &
                  .AND. RAPCL(ij,IP,k) >= rmns ) then
                IDXAE(ij,IP,k) = im
             endif
          enddo
          !$omp end do
       enddo
    enddo
    enddo
    !$omp end parallel

    !$omp parallel default(none),private(ij,k,kt,IP,im,rpls,rmns), &
    !$omp shared(ijdim,klev,KCPCL,KCTYP,KMODE,IDXCL,RCPCL,RMODE)
    do k  = 1, klev
    do kt = 1, KCTYP ! stratiform/cumulus
    do IP = 1, KCPCL ! cloud type(cloud, rain, ice, snow, graupel....)
       !$omp do
       do ij = 1, ijdim
          IDXCL(ij,IP,kt,k) = 0
       enddo
       !$omp end do

       do im = 0, KMODE ! mode radius
          rpls = RMODE(im+1,IP,1)
          rmns = RMODE(im  ,IP,1)

          !$omp do
          do ij = 1, ijdim
             if (       RCPCL(ij,IP,kt,k) <  rpls &
                  .AND. RCPCL(ij,IP,kt,k) >= rmns ) then
                IDXCL(ij,IP,kt,k) = im
             endif
          enddo
          !$omp end do
       enddo
    enddo
    enddo
    enddo
    !$omp end parallel

    !$omp parallel do default(none),private(ij,k,IP,idx), &
    !$omp shared(ijdim,klev,KCPCL,KAPCL,FXAE,IDXAE,RAPCL,RMODE)
    do k  = 1, klev
    do IP = 1, KAPCL
    do ij = 1, ijdim
       idx = IDXAE(ij,IP,k)
       FXAE(ij,IP,k) = ( RAPCL(ij,IP,k) - RMODE(idx,KCPCL+IP,1) ) &
                     * RMODE(idx,KCPCL+IP,2)
    enddo
    enddo
    enddo
    !$omp end parallel do

    !$omp parallel do default(none),private(ij,k,kt,IP,idx), &
    !$omp shared(ijdim,klev,KCPCL,KCTYP,FXCL,IDXCL,RCPCL,RMODE)
    do k  = 1, klev
    do kt = 1, KCTYP
    do IP = 1, KCPCL
    do ij = 1, ijdim
       idx = IDXCL(ij,IP,kt,k)
       FXCL(ij,IP,kt,k) = ( RCPCL(ij,IP,kt,k) - RMODE(idx,IP,1) ) &
                        * RMODE(idx,IP,2)
    enddo
    enddo
    enddo
    enddo
    !$omp end parallel do

    return
  end subroutine RMDIDX

  !-----------------------------------------------------------------------------
  ! moments of scattering (aerosols)
  subroutine SCATAE( &
       ijdim,   &
       CAPCL,   &
       RDZ,     &
       Q,       &
       IFLG,    &
       IDXAE,   &
       FXAE,    &
       NMODE,   &
       aot_ext, & ! intent out [add] 2016/05/18 T.Seiki
       aot_abs, & ! intent out [add] 2016/05/18 T.Seiki
       QQ       )
    implicit none

    integer,  intent(in)  :: ijdim
    real(RP), intent(in)  :: CAPCL  (ijdim,KAPCL)
    real(RP), intent(in)  :: RDZ    (ijdim)
    real(RP), intent(in)  :: Q      (KMODE,KTPCL,KDA*2+2)
    integer,  intent(in)  :: IFLG
    integer,  intent(in)  :: IDXAE  (ijdim,KAPCL)
    real(RP), intent(in)  :: FXAE   (ijdim,KAPCL)
    integer,  intent(in)  :: NMODE  (KTPCL)
    real(RP), intent(out) :: aot_ext(ijdim,KAPCL) ! [add] 2016/05/18 T.Seiki
    real(RP), intent(out) :: aot_abs(ijdim,KAPCL) ! [add] 2016/05/18 T.Seiki
    real(RP), intent(out) :: QQ     (ijdim,KDA*2+2)

    real(RP) :: QA(ijdim,KAPCL,KDA*2+2)
    real(RP) :: QB(ijdim,KAPCL,KDA*2+2)
    integer  :: IMA, IMB

    integer  :: ij, IDA, IP
    !---------------------------------------------------------------------------

    !$omp parallel default(none),private(ij,IDA,IP), &
    !$omp shared(ijdim,QQ,aot_ext,aot_abs)

    !$omp do
    do IDA = 1, KDA*2+2
    do ij  = 1, ijdim
       QQ(ij,IDA) = 0.0_RP
    enddo
    enddo
    !$omp end do nowait

    !$omp do
    do IP = 1, KAPCL
    do ij = 1, ijdim
       aot_ext(ij,IP) = 0.0_RP ! [Add] 2016/05/18 T.Seiki
       aot_abs(ij,IP) = 0.0_RP ! [Add] 2016/05/18 T.Seiki
    enddo
    enddo
    !$omp end do

    !$omp end parallel

    if( IFLG <= 0 ) return

    ! [comment] this means we never use extra-interpolation beyond NMODE.
    !$omp parallel do default(none),private(ij,IDA,IP,IMA,IMB), &
    !$omp shared(ijdim,KCPCL,QA,QB,IDXAE,NMODE,Q)
    do ij  = 1, ijdim
       do IDA = 1, KDA*2+2
       do IP  = 1, KAPCL
          IMA = max( IDXAE(ij,IP)  , 1               )
          IMB = min( IDXAE(ij,IP)+1, NMODE(IP+KCPCL) )

          QA(ij,IP,IDA) = Q(IMA,IP+KCPCL,IDA)
          QB(ij,IP,IDA) = Q(IMB,IP+KCPCL,IDA)
       enddo
       enddo
    enddo
    !$omp end parallel do

    !$omp parallel default(none),private(ij,IDA,IP), &
    !$omp shared(ijdim,QQ,aot_ext,aot_abs,FXAE,QA,QB,CAPCL,RDZ)
    do IP = 1, KAPCL

       !$omp do
       do IDA = 1, KDA*2+2
       do ij  = 1, ijdim
          QQ(ij,IDA) = QQ(ij,IDA) &
                     + ( ( 1.0_RP-FXAE(ij,IP) ) * QA(ij,IP,IDA) &
                       + (        FXAE(ij,IP) ) * QB(ij,IP,IDA) ) * CAPCL(ij,IP) * RDZ(ij) * 0.1_RP
       enddo
       enddo
       !$omp end do nowait

       ! [Add] 2016/05/18 T.Seiki
       !$omp do
       do ij = 1, ijdim
          aot_ext(ij,IP) = ( ( 1.0_RP-FXAE(ij,IP) ) * QA(ij,IP,1) &
                           + (        FXAE(ij,IP) ) * QB(ij,IP,1) ) * CAPCL(ij,IP) * RDZ(ij) * 0.1_RP
          aot_abs(ij,IP) = ( ( 1.0_RP-FXAE(ij,IP) ) * QA(ij,IP,2) &
                           + (        FXAE(ij,IP) ) * QB(ij,IP,2) ) * CAPCL(ij,IP) * RDZ(ij) * 0.1_RP
       enddo
       !$omp end do
    enddo
    !$omp end parallel

    return
  end subroutine SCATAE

  !-----------------------------------------------------------------------------
  ! moments of scattering (Rayleigh)
  subroutine SCATRY( &
       ijdim, &
       RY,    &
       DPRE,  &
       Pstd,  &
       QMOL,  &
       IFLG,  &
       QQ     )
    implicit none

    integer,  intent(in)    :: ijdim
    real(RP), intent(in)    :: RY
    real(RP), intent(in)    :: DPRE(ijdim)
    real(RP), intent(in)    :: Pstd         ! [Pa]
    real(RP), intent(in)    :: QMOL(KDMAX)
    integer,  intent(in)    :: IFLG
    real(RP), intent(inout) :: QQ  (ijdim,KDA*2+2)

    integer  :: ij, IDA
    !---------------------------------------------------------------------------

    if( IFLG <= 0 ) return

    !$omp parallel do default(none),private(ij,IDA), &
    !$omp shared(ijdim,QQ,RY,DPRE,Pstd,QMOL)
    do IDA = 1, KDA*2
    do ij  = 1, ijdim
       QQ(ij,IDA) = QQ(ij,IDA) + RY * DPRE(ij) / (Pstd*1.E-2_RP) * QMOL(IDA)
    enddo
    enddo
    !$omp end parallel do

    return
  end subroutine SCATRY

  !-----------------------------------------------------------------------------
  ! moments of scattering (clouds)
  subroutine SCATCL( &
       ijdim, &
       CCPCL, &
       RDZ,   &
       Q,     &
       IFLG,  &
       IFCLD, &
       IDXCL, &
       FXCL,  &
       NMODE, &
       QQ,    &
       TAUP,  &
       SCAP,  &
       G,     &
       TAUCV  )
    implicit none

    integer,  intent(in)  :: ijdim
    real(RP), intent(in)  :: CCPCL(ijdim,KCPCL)         ! volume ratio of particles [m/m]*1.E+6_RP
    real(RP), intent(in)  :: RDZ  (ijdim)               ! rho(ij,k)*dz(ij,k)
    real(RP), intent(in)  :: Q    (KMODE,KTPCL,KDA*2+2) ! moment derived by parameter table
    integer,  intent(in)  :: IFLG
    integer,  intent(in)  :: IFCLD
    integer,  intent(in)  :: IDXCL(ijdim,KCPCL)         ! index of table
    real(RP), intent(in)  :: FXCL (ijdim,KCPCL)         ! weight coefficient = (r-r1)/(r2-r1)
    integer,  intent(in)  :: NMODE(KTPCL)               ! num. of index of parameter table
    real(RP), intent(in)  :: QQ   (ijdim,KDA*2+2)       ! accumlated ext.coef.(rayleigh+aerosol mie)
    real(RP), intent(out) :: TAUP (ijdim)               ! extinction cross-section of particles
    real(RP), intent(out) :: SCAP (ijdim)               ! scattering cross-section of particles
    real(RP), intent(out) :: G    (ijdim,KDA*2+1)       ! asymmetry factor and high order moment.
    real(RP), intent(out) :: TAUCV(ijdim)               ! 0.5 micron tau

    real(RP) :: QX(ijdim,KDA*2+2)
    real(RP) :: QA(ijdim,KCPCL,KDA*2+2)
    real(RP) :: QB(ijdim,KCPCL,KDA*2+2)
    integer  :: IMA, IMB

    integer  :: ij, IDA, IP
    !---------------------------------------------------------------------------

    !$omp parallel do default(none),private(ij), &
    !$omp shared(ijdim,G)
    do ij = 1, ijdim
       G(ij,1) = 1.0_RP
    enddo
    !$omp end parallel do

    if ( IFLG <= 0 ) then

       !$omp parallel default(none),private(ij,IDA), &
       !$omp shared(ijdim,TAUP,SCAP,TAUCV,G)

       !$omp do
       do ij = 1, ijdim
          TAUP (ij) = 0.0_RP
          SCAP (ij) = 0.0_RP
          TAUCV(ij) = 0.0_RP
       enddo
       !$omp end do nowait

       !$omp do
       do IDA = 2, KDA*2+1
       do ij  = 1, ijdim
          G(ij,IDA) = 0.0_RP
       enddo
       enddo
       !$omp end do

       !$omp end parallel

       return
    endif

    !$omp parallel do default(none),private(ij,IDA), &
    !$omp shared(ijdim,QX)
    do IDA = 1, KDA*2+2
    do ij  = 1, ijdim
       QX(ij,IDA) = 0.0_RP
    enddo
    enddo
    !$omp end parallel do

    if ( IFCLD >= 1 ) then
       ! [comment] this means we never use extra-interpolation beyond NMODE.
       !$omp parallel do default(none),private(ij,IDA,IP,IMA,IMB), &
       !$omp shared(ijdim,KCPCL,QA,QB,IDXCL,NMODE,Q)
       do ij  = 1, ijdim
       do IDA = 1, KDA*2+2
       do IP  = 1, KCPCL ! 1=liquid, 2=ice
          IMA = max( IDXCL(ij,IP)  , 1         )
          IMB = min( IDXCL(ij,IP)+1, NMODE(IP) )

          QA(ij,IP,IDA) = Q(IMA,IP,IDA)
          QB(ij,IP,IDA) = Q(IMB,IP,IDA)
       enddo
       enddo
       enddo
       !$omp end parallel do

       ! [comment] 08/05/30 T.Mitsui
       ! Q     is extinction coefficient(or other coeffcients) / volume.
       ! CCPCL is volume of particles ([m/m]*1.E+6_RP).
       ! FXCL  is weighting coefficient of interpolation between two radius.

       !$omp parallel default(none),private(ij,IDA,IP), &
       !$omp shared(ijdim,KCPCL,QX,FXCL,QA,QB,CCPCL)
       do IP  = 1, KCPCL ! 1=liquid, 2=ice
          !$omp do
          do IDA = 1, KDA*2+2
          do ij  = 1, ijdim
             QX(ij,IDA) = QX(ij,IDA) + ( ( 1.0_RP-FXCL(ij,IP) ) * QA(ij,IP,IDA) &
                                       + (        FXCL(ij,IP) ) * QB(ij,IP,IDA) ) * CCPCL(ij,IP)
          enddo
          enddo
          !$omp end do
       enddo
       !$omp end parallel
    endif


    !$omp parallel default(none),private(ij,IDA), &
    !$omp shared(ijdim,G,TAUCV,TAUP,SCAP,QX,QQ,RDZ)

    !$omp do
    do ij = 1, ijdim
       ! [comment] 08/05/30 T.Mitsui, What's 0.1 ???
       ! TAUP is sum of optical thickness(extinction) of all particles.
       ! QQ(2) and QX(2) is absorption part of extinction so
       ! (scattering) = (extinction) - (absorption)
       TAUCV(ij) = QX(ij,1) * RDZ(ij) * 0.1_RP ! 0.5 micron tau
       TAUP (ij) = TAUCV(ij) + QQ(ij,1)
       SCAP (ij) = ( QX(ij,1) * RDZ(ij) * 0.1_RP + QQ(ij,1) ) &
                 - ( QX(ij,2) * RDZ(ij) * 0.1_RP + QQ(ij,2) )
    enddo
    !$omp end do

    !$omp do
    do IDA = 3, KDA*2+2
    do ij  = 1, ijdim
       G(ij,IDA-1) = ( QX(ij,IDA) * RDZ(ij) * 0.1_RP + QQ(ij,IDA) ) / SCAP(ij)
    enddo
    enddo
    !$omp end do

    !$omp end parallel

    return
  end subroutine SCATCL

  !-----------------------------------------------------------------------------
  ! plank function expansion (KPLK=2)
  subroutine PLKEXP( &
       ijdim, &
       BP,    &
       BL,    &
       TAU,   &
       IFLG,  &
       CPLK   )
    implicit none

    integer,  intent(in)  :: ijdim
    real(RP), intent(in)  :: BP  (ijdim,2)
    real(RP), intent(in)  :: BL  (ijdim)
    real(RP), intent(in)  :: TAU (ijdim)
    integer,  intent(in)  :: IFLG
    real(RP), intent(out) :: CPLK(ijdim,KPLK+1)

    integer  :: ij, k
    !---------------------------------------------------------------------------

    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,CPLK)
    do k  = 1, KPLK+1
    do ij = 1, ijdim
       CPLK(ij,k) = 0.0_RP
    enddo
    enddo
    !$omp end parallel do

    if ( IFLG > 0 ) then
      !$omp parallel do default(none),private(ij), &
      !$omp shared(ijdim,CPLK,TAU,BP,BL)
       do ij = 1, ijdim
          if ( TAU(ij) > 0.0_RP )then
             CPLK(ij,1) = BP(ij,1)
             CPLK(ij,2) = ( 4.0_RP*BL(ij)-BP(ij,2)-3.0_RP*BP(ij,1) ) / TAU(ij)
             CPLK(ij,3) = (      BP(ij,2)+BP(ij,1)-2.0_RP*BL(ij)   ) / TAU(ij)**2 * 2.0_RP
          endif
       enddo
       !$omp end parallel do
    endif

    return
  end subroutine PLKEXP

  !-----------------------------------------------------------------------------
  ! PLANK FUNCTIONS
  subroutine PLANKS( &
       ijdim, &
       TB,    &
       TL,    &
       GTMP,  &
       WNAVR, &
       APLNK, &
       IFLG,  &
       klev,  &
       BP,    &
       BL,    &
       BGND,  &
       BGND1  )
    implicit none

    integer,  intent(in)  :: ijdim
    real(RP), intent(in)  :: TB   (ijdim,KMAX_RAD+1)
    real(RP), intent(in)  :: TL   (ijdim,KMAX_RAD  )
    real(RP), intent(in)  :: GTMP (ijdim)
    real(RP), intent(in)  :: WNAVR
    real(RP), intent(in)  :: APLNK(KPLNK)
    integer,  intent(in)  :: IFLG
    integer,  intent(in)  :: klev !! No. of levels (usually =KMAX_RAD)
    real(RP), intent(out) :: BP   (ijdim,KMAX_RAD+1)
    real(RP), intent(out) :: BL   (ijdim,KMAX_RAD  )
    real(RP), intent(out) :: BGND (ijdim)
    real(RP), intent(out) :: BGND1(ijdim)

    real(RP) :: GTMP1(ijdim)

    integer  :: ij, k
    !---------------------------------------------------------------------------

    !--- PLANK FUNCTIONS
    if ( IFLG > 0 ) then

       !$omp parallel do default(none),private(ij), &
       !$omp shared(ijdim,GTMP1,GTMP)
       do ij = 1, ijdim
          GTMP1(ij) = GTMP(ij) + 1.0_RP
       enddo
       !$omp end parallel do

       call PLANKF(ijdim,klev+1,WNAVR,APLNK,TB   ,BP   )
       call PLANKF(ijdim,klev  ,WNAVR,APLNK,TL   ,BL   )
       call PLANKF(ijdim,1     ,WNAVR,APLNK,GTMP ,BGND )
       call PLANKF(ijdim,1     ,WNAVR,APLNK,GTMP1,BGND1)
    else
       !$omp parallel default(none),private(ij,k), &
       !$omp shared(ijdim,KMAX_RAD,BP,BL,BGND,BGND1)

       !$omp do
       do k  = 1, KMAX_RAD+1
       do ij = 1, ijdim
          BP(ij,k) = 0.0_RP
       enddo
       enddo
       !$omp end do nowait

       !$omp do
       do k  = 1, KMAX_RAD
       do ij = 1, ijdim
          BL(ij,k) = 0.0_RP
       enddo
       enddo
       !$omp end do nowait

       !$omp do
       do ij = 1, ijdim
          BGND (ij) = 0.0_RP
          BGND1(ij) = 0.0_RP
       enddo
       !$omp end do

       !$omp end parallel
    endif

    return
  end subroutine PLANKS

  !-----------------------------------------------------------------------------
  ! PLANK FUNCTION
  subroutine PLANKF( &
       ijdim, &
       kdim,  &
       WNAVR, &
       APLNK, &
       T,     &
       BP     )
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(RP), intent(in)  :: WNAVR
    real(RP), intent(in)  :: APLNK(KPLNK)
    real(RP), intent(in)  :: T (ijdim,kdim)
    real(RP), intent(out) :: BP(ijdim,kdim)

    real(RP) :: WL, BP0
    real(RP) :: X, BPX
    integer  :: ij, K
    !---------------------------------------------------------------------------

    WL  = 10000.0_RP / WNAVR
    BP0 = 1.0_RP / WL**3

    !$omp parallel do default(none),private(ij,k,X,BPX), &
    !$omp shared(ijdim,kdim,BP,T,APLNK,WL,BP0)
    do k  = 1, kdim
    do ij = 1, ijdim
       X   = 1.0_RP / ( WL*T(ij,K) )
       BPX = ( APLNK(5)*X + APLNK(4) )*X + APLNK(3)
       BPX = ( BPX     *X + APLNK(2) )*X + APLNK(1)

       BP(ij,K) = BP0 * exp(-BPX) / X
    enddo
    enddo
    !$omp end parallel do

    return
  end subroutine PLANKF

  !-----------------------------------------------------------------------------
  ! two-stream
  subroutine TWST( &
       ijdim, &
       IRGN,  &
       AM0,   &
       rAM0,  &
       THK,   &
       SCAP,  &
       G,     &
       FSOL,  &
       CPLK,  &
       EXPD,  &
       EXPDM, &
       R,     &
       T,     &
       ER,    &
       ET,    &
       ERM,   &
       ETM,   &
       EXPT   )
!ESC!    use mod_const, only: &
!ESC!       CONST_PI
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: IRGN                ! 1 for SHORT, 2 for LONG
    real(RP), intent(in)  :: AM0  (ijdim)        ! COS(zenith angle)
    real(RP), intent(in)  :: rAM0 (ijdim)        ! COS(zenith angle), inverse
    real(RP), intent(in)  :: THK  (ijdim)        ! Optical thickness
    real(RP), intent(in)  :: SCAP (ijdim)        ! scattering tau
    real(RP), intent(in)  :: G    (ijdim,3)      ! phase function. g(i)
    real(RP), intent(in)  :: FSOL (ijdim)        ! Solar irradiance
    real(RP), intent(in)  :: CPLK (ijdim,KPLK+1) ! Thermal emission
    real(RP), intent(in)  :: EXPD (ijdim)        ! Direct Solar Trans.
    real(RP), intent(in)  :: EXPDM(ijdim)        ! Direct Solar Trans., Averaged
    real(RP), intent(out) :: R    (ijdim)        ! Transmission Matrix
    real(RP), intent(out) :: T    (ijdim)
    real(RP), intent(out) :: ER   (ijdim)
    real(RP), intent(out) :: ET   (ijdim)
    real(RP), intent(out) :: ERM  (ijdim)
    real(RP), intent(out) :: ETM  (ijdim)
    real(RP), intent(out) :: EXPT (ijdim)

    real(RP) :: OMG(ijdim)
    real(RP) :: TAUsw
    real(RP) :: LWsw
    real(RP) :: SWsw

    real(RP) :: rEPS, PI
    real(RP) :: FACT
    real(RP) :: TAU
    real(RP) :: OMGT, GT, BB
    real(RP) :: X, Y
    real(RP) :: ZEIG
    real(RP) :: E

    real(RP) :: CPLF
    real(RP) :: CPL0, CPL1, CPL2
    real(RP) :: QC1, QC2, QC3, QC4

    real(RP) :: SPls, SMns
    real(RP) :: DPls, DMns, DP0, DM0, DPX
    real(RP) :: SUM2
    real(RP) :: ERL, ETL

    real(RP) :: QGA, QGB, VP0, VM0
    real(RP) :: SQEXP
    real(RP) :: ERS, ETS

    integer  :: ij
    !---------------------------------------------------------------------------

    PI   = CONST_PI
    rEPS = 1.0_RP-EPS(IRGN)

    ! solar source switch
    if ( IRGN == 1 ) then
       SWsw = 1.0_RP
    else
       SWsw = 0.0_RP
    endif
    ! thermal source switch
    if ( IRGN == 2 ) then
       LWsw = 1.0_RP
    else
       LWsw = 0.0_RP
    endif

#if defined _OPENMP && _OPENMP >= 201511
    !$omp parallel do simd schedule(simd:static), default(none), &
#else
    !$omp parallel do default(none), &
#endif
    !$omp private(FACT,TAU,TAUsw,OMGT,GT,X,Y,ZEIG,E,SMns,SPls,BB,CPLF,CPL0,CPL1,CPL2,                &
    !$omp         QC3,QC4,QC1,QC2,DP0,DM0,DPX,DPls,DMns,SUM2,ERL,ETL,QGA,QGB,VP0,VM0,SQEXP,ERS,ETS), &
    !$omp shared(IRGN,ijdim,R,T,ER,ET,ERM,ETM,EXPT,OMG,AM0,rAM0,THK,SCAP,G,FSOL,CPLK,EXPD,EXPDM,     &
    !$omp        rEPS,SWsw,LWsw,EPST,AMUI,PA1,PA0,P01F,P00,WMM,PI,AMUA)
    do ij = 1, ijdim

       if ( THK(ij) > 0.0_RP ) then
          OMG(ij) = SCAP(ij) / THK(ij)
       else
          OMG(ij) = 1.0_RP
       endif

       FACT  = G(ij,1) / ( G(ij,1) - G(ij,3)*OMG(ij) )
       TAU   = THK(ij) / FACT
       ! TAU switch
       TAUsw = 0.5_RP - sign( 0.5_RP, EPST-TAU )

       OMGT  = ( G(ij,1)-G(ij,3) )*OMG(ij) / ( G(ij,1)-G(ij,3)*OMG(ij) )
       OMGT  = min( OMGT, rEPS )
       GT    = ( G(ij,2)-G(ij,3) ) / ( G(ij,1)-G(ij,3) )

       X     = AMUI(IRGN) - 2.0_RP * OMGT * GT * PA1(IRGN)
       Y     = AMUI(IRGN) - 2.0_RP * OMGT      * PA0(IRGN)

       ZEIG  = sqrt( max( X*Y, 0.0_RP ) )

       E     = exp( -TAU * ZEIG )

       SMns  = 2.0_RP * OMGT * GT * P01F(IRGN) * AM0(ij)
       SPls  = 2.0_RP * OMGT *      P00 (IRGN)

       BB    = 1.0_RP / ( (X+ZEIG)*(X+ZEIG) - E*(X-ZEIG)*E*(X-ZEIG) )

       R(ij) = ( X*X - ZEIG*ZEIG ) * ( 1.0_RP - E*E ) * BB
       T(ij) = 4.0_RP * E * X * ZEIG * BB

       !--- thermal source
       CPLF  = WMM(IRGN) * 2.0_RP * PI * ( 1.0_RP-OMGT )
       CPL0  = CPLF * CPLK(ij,1)
       CPL1  = CPLF * CPLK(ij,2) * FACT
       CPL2  = CPLF * CPLK(ij,3) * FACT * FACT

       ! if ( TAU(ij) >  EPST )
       QC3   = CPL1 / Y
       QC4   = CPL2 / Y * 2.0_RP / X
       QC1   = ( CPL0+QC4 ) / Y
       QC2   = QC3 / X
       DP0   = QC1 - QC2
       DM0   = QC1 + QC2
       DPX   = CPL2 / Y * TAU * TAU
       DPls  = DP0 + (QC3-QC4)*TAU + DPX
       DMns  = DM0 + (QC3+QC4)*TAU + DPX

       ! if ( TAU(ij) <= EPST )
       SUM2  = ( CPL0*AMUA(IRGN)*AMUA(IRGN) - CPL1*2.0_RP*AMUA(IRGN) + 2.0_RP*CPL2 ) * TAU / 6.0_RP &
             - ( CPL0*AMUA(IRGN)            - CPL1                                 )       / 2.0_RP

       ERL   = LWsw * ( (          TAUsw ) * ( DM0  - R(ij)*DP0 - T(ij)*DMns ) &
                      + ( 1.0_RP - TAUsw ) * ( CPL0*TAU + SUM2*TAU*TAU )     )

       ETL   = LWsw * ( (          TAUsw ) * ( DPls - T(ij)*DP0 - R(ij)*DMns ) &
                      + ( 1.0_RP - TAUsw ) * ( CPL0*TAU + SUM2*TAU*TAU )     )

       !--- solar source
       EXPT(ij) = (          SWsw ) * exp( -TAU * rAM0(ij) ) &
                + ( 1.0_RP - SWsw ) * 1.0_RP

       ! if ( TAU(ij) > EPST )
       QGA   = ( X*SPls*AM0(ij) + SMns     ) &
             / ( X*Y   *AM0(ij) - rAM0(ij) )
       QGB   = ( QGA*rAM0(ij) + SMns ) / X
       VP0   = 0.5_RP * ( QGA + QGB )
       VM0   = 0.5_RP * ( QGA - QGB )

       ! if ( TAU(ij) <= EPST )
       SQEXP = sqrt( EXPT(ij) )

       ERS   = SWsw * ( (          TAUsw ) * ( VM0          - R(ij)*VP0 - T(ij)*VM0*EXPT(ij) ) &
                      + ( 1.0_RP - TAUsw ) * ( 0.5_RP * (SPls-SMns) * TAU * SQEXP            ) )

       ETS   = SWsw * ( (          TAUsw ) * ( VP0*EXPT(ij) - T(ij)*VP0 - R(ij)*VM0*EXPT(ij) ) &
                      + ( 1.0_RP - TAUsw ) * ( 0.5_RP * (SPls+SMns) * TAU * SQEXP            ) )

       ERM(ij) = ERL + ERS * EXPDM(ij) * FSOL(ij)
       ETM(ij) = ETL + ETS * EXPDM(ij) * FSOL(ij)
       ER (ij) = ERL + ERS * EXPD (ij) * FSOL(ij)
       ET (ij) = ETL + ETS * EXPD (ij) * FSOL(ij)
    enddo
#if defined _OPENMP && _OPENMP >= 201511
    !$omp end parallel do simd
#else
    !$omp end parallel do
#endif

    return
  end subroutine TWST

  !-----------------------------------------------------------------------------
  ! Adding-doubling method
  subroutine ADDING( &
       ijdim, &
       klev,  &
       irgn,  &
       RE,    &
       TE,    &
       SER,   &
       SET,   &
       RDDN,  &
       RDN,   &
       RUP    )
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: klev                   ! # of levels (usually=kmax)
    integer,  intent(in)  :: irgn                   ! 1 for SHORT, 2 for LONG
    real(RP), intent(in)  :: RE  (ijdim,KMAX_RAD+1)
    real(RP), intent(in)  :: TE  (ijdim,KMAX_RAD+1)
    real(RP), intent(in)  :: SER (ijdim,KMAX_RAD+1)
    real(RP), intent(in)  :: SET (ijdim,KMAX_RAD+1)
    real(RP), intent(in)  :: RDDN(ijdim,KMAX_RAD+1) ! direct downward flux (W/m2)
    real(RP), intent(out) :: RDN (ijdim,KMAX_RAD+1) ! downward flux (W/m2)
    real(RP), intent(out) :: RUP (ijdim,KMAX_RAD+1) ! upward   flux (W/m2)

    real(RP) :: RL(ijdim,KMAX_RAD+1)
    real(RP) :: SL(ijdim,KMAX_RAD+1)
    real(RP) :: RT(ijdim,KMAX_RAD+1)
    real(RP) :: ST(ijdim,KMAX_RAD+1)
    real(RP) :: temp, recip

    integer, parameter :: simdlen = 16
    integer  :: blk, vec, veclen

    integer  :: ij, k
    !---------------------------------------------------------------------------

    !$omp parallel default(none) &
    !$omp private(ij,k,blk,vec,veclen,recip,temp) &
    !$omp shared(ijdim,klev,RL,SL,RT,ST,RE,SER,SET,irgn,RUP,RDN,WMX,RDDN,TE)

    !$omp do
#if defined _OPENMP && _OPENMP >= 201511
    do blk = 1, ijdim, simdlen
       veclen = min( ijdim-blk+1, simdlen )
       !$omp simd simdlen(16)
       do vec = 1, veclen
          ij = blk + vec - 1
#else
    do ij = 1, ijdim
#endif
          ! boundary condition
          RL(ij,klev+1) = RE (ij,klev+1)
          SL(ij,klev+1) = SER(ij,klev+1)
          RT(ij,1)      = RE (ij,1)
          ST(ij,1)      = SET(ij,1)

          ! upward adding
          do k  = klev, 1, -1
             recip = 1.0_RP / ( 1.0_RP - RL(ij,k+1)*RE(ij,k) )

             RL(ij,k) = RE (ij,k) + TE(ij,k) * ( RL(ij,k+1)*TE (ij,k)              ) * recip
             SL(ij,k) = SER(ij,k) + TE(ij,k) * ( RL(ij,k+1)*SET(ij,k) + SL(ij,k+1) ) * recip
          enddo

          ! boundary condition
          RUP(ij,1) = WMX(irgn) * SL(ij,1)
          RDN(ij,1) = RDDN(ij,1)

          do k  = 2, klev+1
             ! downward adding
             recip = 1.0_RP / ( 1.0_RP - RT(ij,k-1)*RE(ij,k) )

             RT(ij,k) = RE (ij,k) + TE(ij,k) * ( RT(ij,k-1)*TE (ij,k)              ) * recip
             ST(ij,k) = SET(ij,k) + TE(ij,k) * ( RT(ij,k-1)*SER(ij,k) + ST(ij,k-1) ) * recip

             ! flux
             temp = ( ST(ij,k-1) + RT(ij,k-1)*SL(ij,k) ) / ( 1.0_RP - RT(ij,k-1)*RL(ij,k) )

             RUP(ij,k) = WMX(irgn) * ( SL(ij,k) + RL(ij,k) * temp )
             RDN(ij,k) = WMX(irgn) * (                       temp ) + RDDN(ij,k)
          enddo
#if defined _OPENMP && _OPENMP >= 201511
       enddo
#endif
    enddo
    !$omp end do

    !$omp end parallel

    return
  end subroutine ADDING

  !-----------------------------------------------------------------------------
  ! make T,R,S matrixes for maximal/random approx.
  subroutine RTS_MR( &
       ijdim, &
       klev,  &
       IRGN,  &
       R,     &
       T,     &
       ST,    &
       SR,    &
       TRNS0, &
       FCLD,  &
       B,     &
       RET,   &
       RER,   &
       TET,   &
       TER,   &
       SET,   &
       SER    )
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: klev                            ! No. of levels (=KMAX_RAD)
    integer,  intent(in)  :: IRGN                            ! Optical flag
    real(RP), intent(in)  :: R    (ijdim,KMAX_RAD+1,KCLD)    ! reflection matrices of sublayers.
    real(RP), intent(in)  :: T    (ijdim,KMAX_RAD+1,KCLD)    ! transmission matrices of sublayers.
    real(RP), intent(in)  :: ST   (ijdim,KMAX_RAD+1,KCLD)    ! downward source matrices.
    real(RP), intent(in)  :: SR   (ijdim,KMAX_RAD+1,KCLD)    ! upward source matrices.
    real(RP), intent(in)  :: TRNS0(ijdim,KMAX_RAD+1,2)
    real(RP), intent(in)  :: FCLD (ijdim,KMAX_RAD)           ! fractional cloudiness
    real(RP), intent(in)  :: B    (ijdim,KMAX_RAD,4)         ! cloud cover interaction matrix
    real(RP), intent(out) :: RET  (ijdim,2,2,KMAX_RAD+1) ! reflection matrices of sublayers(1,0).
    real(RP), intent(out) :: RER  (ijdim,2,2,KMAX_RAD+1) ! reflection matrices of sublayers(0,1).
    real(RP), intent(out) :: TET  (ijdim,2,2,KMAX_RAD+1) ! transmission matrices of sublayers(1,0).
    real(RP), intent(out) :: TER  (ijdim,2,2,KMAX_RAD+1) ! transmission matrices of sublayers(0,1).
    real(RP), intent(out) :: SET  (ijdim,2,  KMAX_RAD+1) ! upward source matrices.
    real(RP), intent(out) :: SER  (ijdim,2,  KMAX_RAD+1) ! downward source matrices.

    integer  :: ij, k
    !---------------------------------------------------------------------------

    if ( IRGN == 1 ) then

       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(irgn,ijdim,klev,SER,SET,SR,ST,TRNS0)
       do k  = 1, klev
       do ij = 1, ijdim
          SER(ij,1,k) = SR(ij,k,1) * TRNS0(ij,k,1)
          SER(ij,2,k) = SR(ij,k,2) * TRNS0(ij,k,2)
          SET(ij,1,k) = ST(ij,k,1) * TRNS0(ij,k,1)
          SET(ij,2,k) = ST(ij,k,2) * TRNS0(ij,k,2)
       enddo
       enddo
       !$omp end parallel do

    elseif( IRGN == 2 ) then

       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(irgn,ijdim,klev,SER,SET,SR,ST,FCLD)
       do k  = 1, klev
       do ij = 1, ijdim
          SER(ij,1,k) = SR(ij,k,1) * (        FCLD(ij,k) )
          SER(ij,2,k) = SR(ij,k,2) * ( 1.0_RP-FCLD(ij,k) )
          SET(ij,1,k) = ST(ij,k,1) * (        FCLD(ij,k) )
          SET(ij,2,k) = ST(ij,k,2) * ( 1.0_RP-FCLD(ij,k) )
       enddo
       enddo
       !$omp end parallel do

    endif

    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(irgn,ijdim,klev,RET,RER,TET,TER,R,T,B)
    do k  = 1, klev
    do ij = 1, ijdim
       RET(ij,1,1,k) = R(ij,k,1) * (        B(ij,k,3) )
       RET(ij,1,2,k) = R(ij,k,1) * ( 1.0_RP-B(ij,k,1) )
       RET(ij,2,1,k) = R(ij,k,2) * ( 1.0_RP-B(ij,k,3) )
       RET(ij,2,2,k) = R(ij,k,2) * (        B(ij,k,1) )

       RER(ij,1,1,k) = R(ij,k,1) * (        B(ij,k,4) )
       RER(ij,1,2,k) = R(ij,k,1) * ( 1.0_RP-B(ij,k,2) )
       RER(ij,2,1,k) = R(ij,k,2) * ( 1.0_RP-B(ij,k,4) )
       RER(ij,2,2,k) = R(ij,k,2) * (        B(ij,k,2) )

       TET(ij,1,1,k) = T(ij,k,1) * (        B(ij,k,3) )
       TET(ij,1,2,k) = T(ij,k,1) * ( 1.0_RP-B(ij,k,1) )
       TET(ij,2,1,k) = T(ij,k,2) * ( 1.0_RP-B(ij,k,3) )
       TET(ij,2,2,k) = T(ij,k,2) * (        B(ij,k,1) )

       TER(ij,1,1,k) = T(ij,k,1) * (        B(ij,k,4) )
       TER(ij,1,2,k) = T(ij,k,1) * ( 1.0_RP-B(ij,k,2) )
       TER(ij,2,1,k) = T(ij,k,2) * ( 1.0_RP-B(ij,k,4) )
       TER(ij,2,2,k) = T(ij,k,2) * (        B(ij,k,2) )
    enddo
    enddo
    !$omp end parallel do

    return
  end subroutine RTS_MR

  !-----------------------------------------------------------------------------
  !--- cal. maximal/random flux by the adding method.
  subroutine ADDING_MR( &
       ijdim, &
       klev,  &
       IRGN,  &
       RER,   &
       RET,   &
       TER,   &
       TET,   &
       SER,   &
       SET,   &
       FCLD,  &
       EXPD,  &
       AM0,   &
       GALB,  &
       BGND,  &
       FSOL,  &
       FU,    &
       FD,    &
       FDI    )
!ESC!    use mod_const, only: &
!ESC!       CONST_PI
    implicit none

    integer,  intent(in)    :: ijdim
    integer,  intent(in)    :: klev                       ! No. of levels (=KMAX_RAD)
    integer,  intent(in)    :: IRGN                       ! Optical flag
    real(RP), intent(inout) :: RER (ijdim,2,2,KMAX_RAD+1) ! REFLECTION MATRICES OF SUBLAYERS(1,0).
    real(RP), intent(inout) :: RET (ijdim,2,2,KMAX_RAD+1) ! REFLECTION MATRICES OF SUBLAYERS(0,1).
    real(RP), intent(inout) :: TER (ijdim,2,2,KMAX_RAD+1) ! TRANSMISSION MATRICES OF SUBLAYERS(1,0).
    real(RP), intent(inout) :: TET (ijdim,2,2,KMAX_RAD+1) ! TRANSMISSION MATRICES OF SUBLAYERS(0,1).
    real(RP), intent(inout) :: SER (ijdim,2,  KMAX_RAD+1) ! UPWARD SOURCE MATRICES.
    real(RP), intent(inout) :: SET (ijdim,2,  KMAX_RAD+1) ! downward SOURCE MATRICES.
    real(RP), intent(in)    :: FCLD(ijdim,KMAX_RAD)       ! fractional cloudiness
    real(RP), intent(in)    :: EXPD(ijdim,KMAX_RAD+1,2)
    real(RP), intent(in)    :: AM0 (ijdim)                ! COS(zenith angle)
    real(RP), intent(in)    :: GALB(ijdim,NRDIR)          ! Ground albedo (dir,diff)
    real(RP), intent(in)    :: BGND(ijdim)                ! emission from ground
    real(RP), intent(in)    :: FSOL(ijdim)                ! Solar irradiance
    real(RP), intent(out)   :: FU  (ijdim,KMAX_RAD+1)     ! UPWARD FLUXES.
    real(RP), intent(out)   :: FD  (ijdim,KMAX_RAD+1)     ! downward FLUXES.
    real(RP), intent(out)   :: FDI (ijdim)                ! SFC Direct downward flux

    real(RP) :: RUL0(ijdim,2,2,KMAX_RAD+1)
    real(RP) :: SD0L(ijdim,2,  KMAX_RAD+1)
    real(RP) :: RUP (ijdim,2,  KMAX_RAD+1) !! UPWARD INTENSITIES.
    real(RP) :: RDN (ijdim,2,  KMAX_RAD+1) !! downward INTENSITIES.

    real(RP) :: WK1 (ijdim,2,2)
    real(RP) :: WK2 (ijdim,2,2)
    real(RP) :: WK3 (ijdim,2,2)
    real(RP) :: WK5 (ijdim,2)
    real(RP) :: WK6 (ijdim,2)
    real(RP) :: WK7 (ijdim,2)
    real(RP) :: SE
    real(RP) :: PI

    integer  :: ij, K, N, k1, k2, k3
    !-----------------------------------------------------------------------------

    PI = CONST_PI

    ! LAMBERT SURFACE
    !$omp parallel workshare
    RER(:,1,1,klev+1) = GALB(:,2)
    RER(:,1,2,klev+1) = 0.0_RP
    RER(:,2,1,klev+1) = 0.0_RP
    RER(:,2,2,klev+1) = GALB(:,2)
    RET(:,1,1,klev+1) = GALB(:,2)
    RET(:,1,2,klev+1) = 0.0_RP
    RET(:,2,1,klev+1) = 0.0_RP
    RET(:,2,2,klev+1) = GALB(:,2)
    TER(:,:,:,klev+1) = 0.0_RP
    TET(:,:,:,klev+1) = 0.0_RP
    SET(:,:,  klev+1) = 0.0_RP
    !$omp end parallel workshare

    if ( IRGN == 1 ) then

       !$omp parallel do default(none),private(ij,SE), &
       !$omp shared(ijdim,klev,SER,GALB,AM0,FSOL,EXPD,WMP,AMUA)
       do ij = 1, ijdim
          SE = WMP(1) * GALB(ij,1) * AM0(ij) / AMUA(1) * FSOL(ij)
          SER(ij,1,klev+1) = SE * EXPD(ij,klev+1,1)
          SER(ij,2,klev+1) = SE * EXPD(ij,klev+1,2)
       enddo
       !$omp end parallel do

    elseif( IRGN == 2 ) then

       !$omp parallel do default(none),private(ij,SE), &
       !$omp shared(ijdim,klev,SER,GALB,BGND,FCLD,WMP,AMUA,PI)
       do ij = 1, ijdim
          SE = WMP(2) * 2.0_RP * PI * ( 1.0_RP-GALB(ij,2) ) * BGND(ij)
          SER(ij,1,klev+1) = SE * (          FCLD(ij,klev) )
          SER(ij,2,klev+1) = SE * ( 1.0_RP - FCLD(ij,klev) )
       enddo
       !$omp end parallel do

    endif

    ! ADDING RE<L,0>, SET<0,L>
    !$omp parallel workshare
    SD0L(:,:,  1) = SET(:,:,  1)
    RUL0(:,:,:,1) = RER(:,:,:,1)
    !$omp end parallel workshare

    do k = 2, klev
       k1 = K-1
       k2 = K

       call M3X3  ( WK1(:,:,:), RUL0(:,:,:,k1), RET (:,:,:,k2), ijdim )
       call MULTMR( WK2(:,:,:), WK1 (:,:,:),                    ijdim )
       call M3X3  ( WK3(:,:,:), TET (:,:,:,k2), WK2 (:,:,:),    ijdim )
       ! ADDING 'SD0L'
       call M3X2  ( WK5(:,:),   RUL0(:,:,:,k1), SER (:,:,k2),   ijdim )
       !$omp parallel workshare
       WK6 (:,:)    = WK5(:,:) + SD0L(:,:,k1)
       !$omp end parallel workshare
       call M3X2  ( WK7(:,:),   WK3 (:,:,:),    WK6 (:,:),      ijdim )
       !$omp parallel workshare
       SD0L(:,:,k2) = WK7(:,:) + SET (:,:,k2)
       !$omp end parallel workshare
       ! ADDING 'RUL0'
       call M3X3  ( WK1(:,:,:), WK3 (:,:,:),    RUL0(:,:,:,k1), ijdim )
       call M3X3  ( WK2(:,:,:), WK1 (:,:,:),    TER (:,:,:,k2), ijdim )
       !$omp parallel workshare
       RUL0(:,:,:,k2) = WK2(:,:,:) + RER(:,:,:,k2)
       !$omp end parallel workshare
    enddo

    ! ADDING RADIANCE
    ! RUP,RDN : k1 (LAYER BOUNDARY NUMBER)
    k1 = klev
    k2 = klev+1

    call M3X3  ( WK1(:,:,:),  RUL0(:,:,:,k1), RET (:,:,:,k2), ijdim )
    call MULTMR( WK2(:,:,:),  WK1 (:,:,:),                    ijdim )
    call M3X3  ( WK1(:,:,:),  RET (:,:,:,k2), RUL0(:,:,:,k1), ijdim )
    call MULTMR( WK3(:,:,:),  WK1 (:,:,:),                    ijdim )
    ! RDN
    call M3X2  ( WK5(:,:),    RUL0(:,:,:,k1), SER (:,:,k2),   ijdim )
    !$omp parallel workshare
    WK6(:,:) = WK5(:,:) + SD0L(:,:,k1)
    !$omp end parallel workshare
    call M3X2  ( RDN(:,:,k2), WK2 (:,:,:),    WK6 (:,:),      ijdim )
    ! RUP
    call M3X2  ( WK5(:,:),    RET (:,:,:,k2), SD0L(:,:,k1),   ijdim )
    !$omp parallel workshare
    WK6(:,:) = WK5(:,:) + SER (:,:,k2)
    !$omp end parallel workshare
    call M3X2  ( RUP(:,:,k2), WK3 (:,:,:),    WK6 (:,:),      ijdim )

    do k = klev-1, 1, -1
       k1 = K
       k2 = K+1
       k3 = K+2

       call M3X3  ( WK1(:,:,:),  RUL0(:,:,:,k1), RET (:,:,:,k2), ijdim )
       call MULTMR( WK2(:,:,:),  WK1 (:,:,:),                    ijdim )
       call M3X3  ( WK1(:,:,:),  RET (:,:,:,k2), RUL0(:,:,:,k1), ijdim )
       call MULTMR( WK3(:,:,:),  WK1 (:,:,:),                    ijdim )
       ! RDN
       call M3X3  ( WK1(:,:,:),  RUL0(:,:,:,k1), TER (:,:,:,k2), ijdim )
       call M3X2  ( WK5(:,:),    WK1 (:,:,:),    RUP (:,:,k3),   ijdim )
       call M3X2  ( WK6(:,:),    RUL0(:,:,:,k1), SER (:,:,k2),   ijdim )
       !$omp parallel workshare
       WK7(:,:) = WK5(:,:) + WK6(:,:) + SD0L(:,:,k1)
       !$omp end parallel workshare
       call M3X2  ( RDN(:,:,k2), WK2 (:,:,:),    WK7 (:,:),      ijdim )
       ! RUP
       call M3X2  ( WK5(:,:),    TER (:,:,:,k2), RUP (:,:,k3),   ijdim )
       call M3X2  ( WK6(:,:),    RET (:,:,:,k2), SD0L(:,:,k1),   ijdim )
       !$omp parallel workshare
       WK7(:,:) = WK5(:,:) + WK6(:,:) + SER(:,:,k2)
       !$omp end parallel workshare
       call M3X2  ( RUP(:,:,k2), WK3 (:,:,:),    WK7 (:,:),      ijdim )
    enddo

    k1 = 1
    k2 = 2

    ! RDN
    !$omp parallel workshare
    RDN(:,:,k1) = 0.0_RP
    !$omp end parallel workshare
    ! RUP
    call M3X2( WK5(:,:), TER(:,:,:,k1), RUP(:,:,k2), ijdim )
    !$omp parallel workshare
    RUP(:,:,k1) = WK5(:,:) + SER(:,:,k1)
    !$omp end parallel workshare

    ! FLUX
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(IRGN,ijdim,klev,FU,FD,RUP,RDN,AM0,FSOL,EXPD,WMX)
    do k  = 1, klev+1
    do ij = 1, ijdim
       FU(ij,k) = WMX(IRGN) * ( RUP(ij,1,k) + RUP(ij,2,k) )
       FD(ij,k) = WMX(IRGN) * ( RDN(ij,1,k) + RDN(ij,2,k) ) &
                + AM0(ij) * FSOL(ij) * ( EXPD(ij,k,1) + EXPD(ij,k,2) )
    enddo
    enddo
    !$omp end parallel do

    !$omp parallel workshare
    FDI(:) = AM0(:) * FSOL(:) * ( EXPD(:,klev+1,1) + EXPD(:,klev+1,2) )
    !$omp end parallel workshare

    return
  end subroutine ADDING_MR

  !-----------------------------------------------------------------------------
  ! cal. B matrix for maximal/random approx.
  subroutine BCVR( &
       ijdim, &
       klev,  &
       FCLD,  &
       B      )
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: klev                   ! No. of levels (usually =KMAX_RAD )
    real(RP), intent(in)  :: FCLD(ijdim,KMAX_RAD)   ! fractional cloudiness
    real(RP), intent(out) :: B   (ijdim,KMAX_RAD,4) ! cloud cover interaction matrix

    real(RP) :: FCLR(ijdim,KMAX_RAD) ! 1 - FCLD
    real(RP) :: zerosw
    real(RP), parameter :: EPS = 1.E-14_RP

    integer  :: ij, k
    !---------------------------------------------------------------------------

    !$omp parallel workshare
    FCLR(:,:) = 1.0_RP - FCLD(:,:)
    !$omp end parallel workshare

    k = 1
    !$omp parallel do default(none),private(ij,zerosw), &
    !$omp shared(k,ijdim,B,FCLR,FCLD)
    do ij = 1, ijdim
       B(ij,k,1) = FCLR(ij,k)

       zerosw = 0.5_RP - sign(0.5_RP,FCLR(ij,k+1)-EPS)
       B(ij,k,2) = (1.0_RP-zerosw) * ( 1.0_RP-max(FCLD(ij,k),FCLD(ij,k+1)) ) / ( FCLR(ij,k+1)+zerosw ) &
                 + (       zerosw)

       B(ij,k,3) = 1.0_RP

       zerosw = 0.5_RP - sign(0.5_RP,FCLD(ij,k+1)-EPS)
       B(ij,k,4) = (1.0_RP-zerosw) * (        min(FCLD(ij,k+1),FCLD(ij,k)) ) / ( FCLD(ij,k+1)+zerosw ) &
                 + (       zerosw)
    enddo
    !$omp end parallel do

    !$omp parallel do default(none),private(ij,k,zerosw), &
    !$omp shared(ijdim,klev,B,FCLR,FCLD)
    do k  = 2, klev-1
    do ij = 1, ijdim
       zerosw = 0.5_RP - sign(0.5_RP,FCLR(ij,k-1)-EPS)
       B(ij,k,1) = (1.0_RP-zerosw) * ( 1.0_RP-max(FCLD(ij,k),FCLD(ij,k-1)) ) / ( FCLR(ij,k-1)+zerosw ) &
                 + (       zerosw)

       zerosw = 0.5_RP - sign(0.5_RP,FCLR(ij,k+1)-EPS)
       B(ij,k,2) = (1.0_RP-zerosw) * ( 1.0_RP-max(FCLD(ij,k),FCLD(ij,k+1)) ) / ( FCLR(ij,k+1)+zerosw ) &
                 + (       zerosw)

       zerosw = 0.5_RP - sign(0.5_RP,FCLD(ij,k-1)-EPS)
       B(ij,k,3) = (1.0_RP-zerosw) * (        min(FCLD(ij,k-1),FCLD(ij,k)) ) / ( FCLD(ij,k-1)+zerosw ) &
                 + (       zerosw)

       zerosw = 0.5_RP - sign(0.5_RP,FCLD(ij,k+1)-EPS)
       B(ij,k,4) = (1.0_RP-zerosw) * (        min(FCLD(ij,k+1),FCLD(ij,k)) ) / ( FCLD(ij,k+1)+zerosw ) &
                 + (       zerosw)
    enddo
    enddo
    !$omp end parallel do

    k = klev
    !$omp parallel do default(none),private(ij,zerosw), &
    !$omp shared(k,ijdim,B,FCLR,FCLD)
    do ij = 1, ijdim
       zerosw = 0.5_RP - sign(0.5_RP,FCLR(ij,k-1)-EPS)
       B(ij,k,1) = (1.0_RP-zerosw) * ( 1.0_RP-max(FCLD(ij,k),FCLD(ij,k-1)) ) / ( FCLR(ij,k-1)+zerosw ) &
                 + (       zerosw)

       B(ij,k,2) = 1.0_RP

       zerosw = 0.5_RP - sign(0.5_RP,FCLD(ij,k-1)-EPS)
       B(ij,k,3) = (1.0_RP-zerosw) * (        min(FCLD(ij,k-1),FCLD(ij,k)) ) / ( FCLD(ij,k-1)+zerosw ) &
                 + (       zerosw)

       B(ij,k,4) = 1.0_RP
    enddo
    !$omp end parallel do

    return
  end subroutine BCVR

  !-----------------------------------------------------------------------------
  ! calculation of multiple reflection between two layers
  subroutine MULTMR( &
       CC, &
       CD, &
       NN  )
    implicit none

    integer,  intent(in)  :: NN
    real(RP), intent(out) :: CC(NN,2,2)
    real(RP), intent(in)  :: CD(NN,2,2) ! R = R1 * R2.
    !---------------------------------------------------------------------------

    ! CC = ( 1 - CD )**-1 = 1 + CD + CD**2 + CD**3 + ...

    ! 1 / DET(1-CD)
    !$omp parallel default(none),shared(CC,CD)

    !$omp workshare
    CC(:,2,2) = 1.0_RP / ( ( 1.0_RP-CD(:,1,1) ) * ( 1.0_RP-CD(:,2,2) ) - CD(:,2,1) * CD(:,1,2) )
    !$omp end workshare

    !$omp workshare
    CC(:,1,1) = ( 1.0_RP - CD(:,2,2) ) * CC(:,2,2)
    !$omp end workshare
    !$omp workshare
    CC(:,2,1) = (          CD(:,2,1) ) * CC(:,2,2)
    !$omp end workshare
    !$omp workshare
    CC(:,1,2) = (          CD(:,1,2) ) * CC(:,2,2)
    !$omp end workshare
    !$omp workshare
    CC(:,2,2) = ( 1.0_RP - CD(:,1,1) ) * CC(:,2,2)
    !$omp end workshare

    !$omp end parallel

    return
  end subroutine MULTMR

  !-----------------------------------------------------------------------------
  ! C=A*B 3(2)-DIM = 3(2)-DIM * 2(1)-DIM
  subroutine M3X2( &
       C,  &
       A,  &
       B,  &
       NN  )
    implicit none

    integer,  intent(in)  :: NN
    real(RP), intent(out) :: C(NN,2)   ! A*B.
    real(RP), intent(in)  :: A(NN,2,2) ! 3-DIM ARRAY A.
    real(RP), intent(in)  :: B(NN,2)   ! 2-DIM ARRAY B.
    !---------------------------------------------------------------------------

    !$omp parallel default(none),shared(C,A,B)

    !$omp workshare
    C(:,1) = A(:,1,1)*B(:,1) + A(:,1,2)*B(:,2)
    !$omp end workshare
    !$omp workshare
    C(:,2) = A(:,2,1)*B(:,1) + A(:,2,2)*B(:,2)
    !$omp end workshare

    !$omp end parallel

    return
  end subroutine M3X2

  !-----------------------------------------------------------------------------
  ! C=A*B 3(2)-DIM = 3(2)-DIM * 3(2)-DIM
  subroutine M3X3( &
       C,  &
       A,  &
       B,  &
       NN  )
    implicit none

    integer,  intent(in)  :: NN
    real(RP), intent(out) :: C(NN,2,2) ! A*B.
    real(RP), intent(in)  :: A(NN,2,2) ! 3-DIM ARRAY A.
    real(RP), intent(in)  :: B(NN,2,2) ! 3-DIM ARRAY B.
    !---------------------------------------------------------------------------

    !$omp parallel default(none),shared(C,A,B)

    !$omp workshare
    C(:,1,1) = A(:,1,1)*B(:,1,1) + A(:,1,2)*B(:,2,1)
    !$omp end workshare
    !$omp workshare
    C(:,2,1) = A(:,2,1)*B(:,1,1) + A(:,2,2)*B(:,2,1)
    !$omp end workshare
    !$omp workshare
    C(:,1,2) = A(:,1,1)*B(:,1,2) + A(:,1,2)*B(:,2,2)
    !$omp end workshare
    !$omp workshare
    C(:,2,2) = A(:,2,1)*B(:,1,2) + A(:,2,2)*B(:,2,2)
    !$omp end workshare

    !$omp end parallel

    return
  end subroutine M3X3

end module mod_rd_mstrnx
