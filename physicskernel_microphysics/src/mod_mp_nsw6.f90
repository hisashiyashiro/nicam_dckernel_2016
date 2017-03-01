!-------------------------------------------------------------------------------
!> module ATMOSPHERE / Physics Cloud Microphysics
!!
!! @par Description
!!          6-category one-moment bulk cloud microphysics NSW6
!!          Reference: Tomita(2008)
!!
!!          Reference : Berry(1968)                : Proc.Phys.Soc.,B66,pp.688-694
!!                      Hsie et al.(1980)          : J.Appl.Meteor.,19,pp.950-977
!!                      Grabowski(1998)            : J.Atmos.Sci.,55,pp.3283-3298
!!                      Hon and Lim(2006)          : J.Korean Meteor.Soc.,42,pp.129-151
!!                      Lin et al.(1983)           : J.Appl.Meteor.,22,pp.1065-1092
!!                      Prupacher and Klett(1978)  : Micophysics of Clouds and
!!                                                   Precipitation.
!!                                                   Kluwer Academic Publishers
!!                      Ruttledge and Hobbs(1983)  : J.Atmos.Sci.,40,pp.1185-1206
!!                      Ruttledge and Hobbs(1984)  : J.Atmos.Sci.,40,pp.2949-2977
!!                      Heymsfield and Doner(1990) : J.Atmos.Sci.,??,pp.1865-1877
!!                      Khairoutdinov and Kogan (2000)
!!                      Roh and Satoh(2014)        : J.Atmos.Sci.,71,pp.2654-2673
!!                      Aerosol Indirect Effect:
!!                      Suzuki et al. (2004)       : J.Atmos.Sci.,61,pp.179-194
!!
!! @author NICAM developers
!!
!<
!-------------------------------------------------------------------------------
module mod_mp_nsw6
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
  !++ Public procedure
  !
  public :: mp_nsw6_init
  public :: mp_nsw6
!ESC!  public :: mp_nsw6_terminal_velocity
!ESC!  public :: mp_nsw6_effective_radius

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: Bergeron_param

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, parameter :: dens00 = 1.28_RP   !< standard density [kg/m3]

  ! Parameter for Marshall-Palmer distribution
  real(RP), private            :: N0r    = 8.E+6_RP  !< intercept parameter for rain    [1/m4]
  real(RP), private            :: N0s    = 3.E+6_RP  !< intercept parameter for snow    [1/m4]
  real(RP), private            :: N0g    = 4.E+6_RP  !< intercept parameter for graupel [1/m4]

  real(RP), private            :: rho_s  = 100.0_RP  !< density of snow    [kg/m3]
  real(RP), private            :: rho_g  = 400.0_RP  !< density of graupel [kg/m3]
                                                     !   graupel : 400
                                                     !   hail    : 917

  real(RP), private            :: C_d    =   0.6_RP  !< drag coefficient for graupel
  real(RP), private            :: Cr     = 130.0_RP
  real(RP), private            :: Cs     =  4.84_RP

  ! Empirical parameter
  real(RP), private            :: Ar, As, Ag
  real(RP), private            :: Br, Bs, Bg
  real(RP), private            :: Cg
  real(RP), private            :: Dr, Ds, Dg

  ! GAMMA function
  real(RP), private            :: GAM, GAM_2, GAM_3

  real(RP), private            :: GAM_1br, GAM_2br, GAM_3br
  real(RP), private            :: GAM_3dr
  real(RP), private            :: GAM_6dr
  real(RP), private            :: GAM_1brdr
  real(RP), private            :: GAM_5dr_h

  real(RP), private            :: GAM_1bs, GAM_2bs, GAM_3bs
  real(RP), private            :: GAM_3ds
  real(RP), private            :: GAM_1bsds
  real(RP), private            :: GAM_5ds_h

  real(RP), private            :: GAM_1bg, GAM_3dg
  real(RP), private            :: GAM_1bgdg
  real(RP), private            :: GAM_5dg_h

  !---< Khairoutdinov and Kogan (2000) >---
  real(RP), private            :: sw_kk2000  = 0.0_RP      !< switch for k-k scheme

  !---< Roh and Satoh (2014) >---
  real(RP), private            :: sw_roh2014 = 0.0_RP      !< switch for Roh scheme
  real(RP), private            :: ln10                 !< log(10)
  real(RP), private            :: coef_a(10) = (/ 5.065339_RP, -0.062659_RP, -3.032362_RP, 0.029469_RP, -0.000285_RP, &
                                                  0.31255_RP,   0.000204_RP,  0.003199_RP, 0.0_RP,      -0.015952_RP  /)
  real(RP), private            :: coef_b(10) = (/ 0.476221_RP, -0.015896_RP,  0.165977_RP, 0.007468_RP, -0.000141_RP, &
                                                  0.060366_RP,  0.000079_RP,  0.000594_RP, 0.0_RP,      -0.003577_RP  /)

  ! Accretion parameter
  real(RP), private            :: Eiw        = 1.0_RP      !< collection efficiency of cloud ice for cloud water
  real(RP), private            :: Erw        = 1.0_RP      !< collection efficiency of rain    for cloud water
  real(RP), private            :: Esw        = 1.0_RP      !< collection efficiency of snow    for cloud water
  real(RP), private            :: Egw        = 1.0_RP      !< collection efficiency of graupel for cloud water
  real(RP), private            :: Eri        = 1.0_RP      !< collection efficiency of rain    for cloud ice
  real(RP), private            :: Esi        = 1.0_RP      !< collection efficiency of snow    for cloud ice
  real(RP), private            :: Egi        = 0.1_RP      !< collection efficiency of graupel for cloud ice
  real(RP), private            :: Esr        = 1.0_RP      !< collection efficiency of snow    for rain
  real(RP), private            :: Egr        = 1.0_RP      !< collection efficiency of graupel for rain
  real(RP), private            :: Egs        = 1.0_RP      !< collection efficiency of graupel for snow
  real(RP), private            :: gamma_sacr = 25.E-3_RP   !< effect of low temperature for Esi
  real(RP), private            :: gamma_gacs = 90.E-3_RP   !< effect of low temperature for Egs
  real(RP), private            :: mi         = 4.19E-13_RP !< mass of one cloud ice crystal [kg]

  ! Auto-conversion parameter
  real(RP), private, parameter :: Nc_lnd     = 2000.0_RP   !< number concentration of cloud water (land)  [1/cc]
  real(RP), private, parameter :: Nc_ocn     =   50.0_RP   !< number concentration of cloud water (ocean) [1/cc]
  real(RP), private            :: Nc_def                   !< number concentration of cloud water         [1/cc]

  real(RP), private            :: beta_saut  =  1.E-3_RP   !< auto-conversion factor beta  for ice
  real(RP), private            :: gamma_saut = 25.E-3_RP   !< auto-conversion factor gamma for ice
  real(RP), private            :: beta_gaut  =  1.E-3_RP   !< auto-conversion factor beta  for snow
  real(RP), private            :: gamma_gaut = 90.E-3_RP   !< auto-conversion factor gamma for snow
  real(RP), private            :: qicrt_saut =  0.0_RP     !< mixing ratio threshold for Psaut [kg/kg]
  real(RP), private            :: qscrt_gaut =  6.E-4_RP   !< mixing ratio threshold for Pgaut [kg/kg]

  ! Evaporation, Sublimation parameter
  real(RP), private, parameter :: Ka0        = 2.428E-2_RP !< thermal diffusion coefficient of air at 0C,1atm [J/m/s/K]
  real(RP), private, parameter :: dKa_dT     =  7.47E-5_RP !< Coefficient of Ka depending on temperature      [J/m/s/K/K]
  real(RP), private, parameter :: Kd0        = 2.222E-5_RP !< diffusion coefficient of water vapor in the air at 0C,1atm [m2/s]
  real(RP), private, parameter :: dKd_dT     =  1.37E-7_RP !< Coefficient of Dw depending on temperature                 [m2/s/K]
  real(RP), private, parameter :: nu0        = 1.718E-5_RP !< kinematic viscosity of air at 0C,1atm      [m2/s*kg/m3]
  real(RP), private, parameter :: dnu_dT     =  5.28E-8_RP !< Coefficient of mu depending on temperature [m2/s/K*kg/m3]

  real(RP), private            :: f1r        = 0.78_RP     !< ventilation factor 1 for rain
  real(RP), private            :: f2r        = 0.27_RP     !< ventilation factor 2 for rain
  real(RP), private            :: f1s        = 0.65_RP     !< ventilation factor 1 for snow
  real(RP), private            :: f2s        = 0.39_RP     !< ventilation factor 2 for snow
  real(RP), private            :: f1g        = 0.78_RP     !< ventilation factor 1 for graupel
  real(RP), private            :: f2g        = 0.27_RP     !< ventilation factor 2 for graupel

  ! Freezing parameter
  real(RP), private            :: A_frz      = 0.66_RP     !< freezing factor [/K]
  real(RP), private            :: B_frz      = 100.0_RP    !< freezing factor [/m3/s]

  ! Bergeron process parameter
  real(RP), private            :: mi40       = 2.46E-10_RP !< mass              of a 40 micron ice crystal [kg]
  real(RP), private            :: mi50       = 4.80E-10_RP !< mass              of a 50 micron ice crystal [kg]
  real(RP), private            :: vti50      = 1.0_RP      !< terminal velocity of a 50 micron ice crystal [m/s]
  real(RP), private            :: Ri50       = 5.E-5_RP    !< radius            of a 50 micron ice crystal [m]

  ! Explicit ice generation
  real(RP), private            :: sw_expice  = 0.0_RP      !< switch for explicit ice generation
  real(RP), private, parameter :: Nc_ihtr    = 300.0_RP    !< cloud number concentration for heterogeneous ice nucleation [1/cc]
  real(RP), private, parameter :: Di_max     = 500.E-6_RP
  real(RP), private, parameter :: Di_a       = 11.9_RP

  integer,  private, parameter :: wk_nmax = 49
  integer,  private, parameter :: I_dqv_dt  =  1 !
  integer,  private, parameter :: I_dqc_dt  =  2 !
  integer,  private, parameter :: I_dqr_dt  =  3 !
  integer,  private, parameter :: I_dqi_dt  =  4 !
  integer,  private, parameter :: I_dqs_dt  =  5 !
  integer,  private, parameter :: I_dqg_dt  =  6 !
  integer,  private, parameter :: I_delta1  =  7 ! separation switch for r->s,g
  integer,  private, parameter :: I_delta2  =  8 ! separation switch for s->g
  integer,  private, parameter :: I_spsati  =  9 ! separation switch for ice sublimation
  integer,  private, parameter :: I_iceflg  = 10 ! separation switch for T > 0
  integer,  private, parameter :: I_RLMDr   = 11
  integer,  private, parameter :: I_RLMDs   = 12
  integer,  private, parameter :: I_RLMDg   = 13
  integer,  private, parameter :: I_Piacr   = 14 ! r->s,g
  integer,  private, parameter :: I_Psacr   = 15 ! r->s,g
  integer,  private, parameter :: I_Praci   = 16 ! i->s,g
  integer,  private, parameter :: I_Pigen   = 17 ! v->i
  integer,  private, parameter :: I_Pidep   = 18 ! v->i
  integer,  private, parameter :: I_Psdep   = 19 ! v->s
  integer,  private, parameter :: I_Pgdep   = 20 ! v->g
  integer,  private, parameter :: I_Praut   = 21 ! c->r
  integer,  private, parameter :: I_Pracw   = 22 ! c->r
  integer,  private, parameter :: I_Pihom   = 23 ! c->i
  integer,  private, parameter :: I_Pihtr   = 24 ! c->i
  integer,  private, parameter :: I_Psacw   = 25 ! c->s
  integer,  private, parameter :: I_Psfw    = 26 ! c->s
  integer,  private, parameter :: I_Pgacw   = 27 ! c->g
  integer,  private, parameter :: I_Prevp   = 28 ! r->v
  integer,  private, parameter :: I_Piacr_s = 29 ! r->s
  integer,  private, parameter :: I_Psacr_s = 30 ! r->s
  integer,  private, parameter :: I_Piacr_g = 31 ! r->g
  integer,  private, parameter :: I_Psacr_g = 32 ! r->g
  integer,  private, parameter :: I_Pgacr   = 33 ! r->g
  integer,  private, parameter :: I_Pgfrz   = 34 ! r->g
  integer,  private, parameter :: I_Pisub   = 35 ! i->v
  integer,  private, parameter :: I_Pimlt   = 36 ! i->c
  integer,  private, parameter :: I_Psaut   = 37 ! i->s
  integer,  private, parameter :: I_Praci_s = 38 ! i->s
  integer,  private, parameter :: I_Psaci   = 39 ! i->s
  integer,  private, parameter :: I_Psfi    = 40 ! i->s
  integer,  private, parameter :: I_Praci_g = 41 ! i->g
  integer,  private, parameter :: I_Pgaci   = 42 ! i->g
  integer,  private, parameter :: I_Pssub   = 43 ! s->v
  integer,  private, parameter :: I_Psmlt   = 44 ! s->r
  integer,  private, parameter :: I_Pgaut   = 45 ! s->g
  integer,  private, parameter :: I_Pracs   = 46 ! s->g
  integer,  private, parameter :: I_Pgacs   = 47 ! s->g
  integer,  private, parameter :: I_Pgsub   = 48 ! g->v
  integer,  private, parameter :: I_Pgmlt   = 49 ! g->r

  character(len=H_SHORT), private :: w_name(wk_nmax)

  data w_name / 'dqv_dt ', &
                'dqc_dt ', &
                'dqr_dt ', &
                'dqi_dt ', &
                'dqs_dt ', &
                'dqg_dt ', &
                'delta1 ', &
                'delta2 ', &
                'spsati ', &
                'iceflg ', &
                'RLMDr  ', &
                'RLMDs  ', &
                'RLMDg  ', &
                'Piacr  ', &
                'Psacr  ', &
                'Praci  ', &
                'Pigen  ', &
                'Pidep  ', &
                'Psdep  ', &
                'Pgdep  ', &
                'Praut  ', &
                'Pracw  ', &
                'Pihom  ', &
                'Pihtr  ', &
                'Psacw  ', &
                'Psfw   ', &
                'Pgacw  ', &
                'Prevp  ', &
                'Piacr_s', &
                'Psacr_s', &
                'Piacr_g', &
                'Psacr_g', &
                'Pgacr  ', &
                'Pgfrz  ', &
                'Pisub  ', &
                'Pimlt  ', &
                'Psaut  ', &
                'Praci_s', &
                'Psaci  ', &
                'Psfi   ', &
                'Praci_g', &
                'Pgaci  ', &
                'Pssub  ', &
                'Psmlt  ', &
                'Pgaut  ', &
                'Pracs  ', &
                'Pgacs  ', &
                'Pgsub  ', &
                'Pgmlt  '  /

  character(len=H_SHORT), private :: precip_transport_type = '3WATER'
  character(len=H_SHORT), private :: precip_scheme_type    = 'Default' !
                                                           ! 'Default' = SL_UPSTREAM
                                                           ! 'Flux-Semilag_new', same as default but fast
                                                           ! 'Upwind-Euler'

  logical,                private :: OPT_EXPLICIT_ICEGEN   = .false.   ! enable explicit ice generation?
  logical,                private :: OPT_INDIR             = .false.   ! enable aerosol indirect effect?
  logical,                private :: Roh_flag              = .false.   ! enable setting by Roh and Satoh (2014)?

  real(RP),               private :: sw_constVti           = 0.0_RP
  real(RP),               private :: CONST_Vti                         ! force constant terminal velocity for ice

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> setup
  subroutine mp_nsw6_init
!ESC!    use mod_misc, only: &
!ESC!       MISC_gammafunc
!ESC!    use mod_adm, only: &
!ESC!       ADM_proc_stop
!ESC!    use mod_const, only: &
!ESC!       UNDEF => CONST_UNDEF, &
!ESC!       PI    => CONST_PI,    &
!ESC!       GRAV  => CONST_GRAV,  &
!ESC!       rho_w => CONST_DWATR
    implicit none

    logical  :: INC_PGAUT   = .false.
    real(RP) :: PSAUT_BETA0

    character(len=H_SHORT) :: qr_aut_acc_type = 'Default' !
                                              ! 'Default' ! Berry-type scheme
                                              ! 'KK2000'  ! Khairoutdinov and Kogan (2000)

    namelist / nm_mp_nsw6 / &
         N0r,                   &
         N0s,                   &
         N0g,                   &
         rho_g,                 &
         rho_s,                 &
         C_d,                   &
         Cr,                    &
         Cs,                    &
         Eiw,                   &
         Erw,                   &
         Esw,                   &
         Egw,                   &
         Eri,                   &
         Esi,                   &
         Egi,                   &
         Esr,                   &
         Egr,                   &
         Egs,                   &
         gamma_sacr,            &
         gamma_gacs,            &
         mi,                    &
         Nc_def,                &
         PSAUT_BETA0,           & ! 10/10/03 noda
         gamma_saut,            &
         qicrt_saut,            &
         beta_gaut,             &
         gamma_gaut,            &
         qscrt_gaut,            &
         f1r,                   &
         f2r,                   &
         f1s,                   &
         f2s,                   &
         f1g,                   &
         f2g,                   &
         A_frz,                 &
         B_frz,                 &
         precip_scheme_type,    &
         precip_transport_type, &
         INC_PGAUT,             &
         OPT_EXPLICIT_ICEGEN,   &
         OPT_INDIR,             &
         CONST_Vti,             & ! 10/09/03 noda
         Roh_flag,              & ! [add] WS. Roh 20141011
         qr_aut_acc_type          ! [Add] T.Seiki 2015/02/26

    integer :: ierr
    !---------------------------------------------------------------------------

    Nc_def      = Nc_ocn
    CONST_Vti   = UNDEF
    PSAUT_BETA0 = beta_saut

!ESC!    !--- read parameters
!ESC!    write(IO_FID_LOG,*)
!ESC!    write(IO_FID_LOG,*) '+++ Module[MP NSW6]/Category[nhm physics]'
!ESC!    rewind(IO_FID_CONF)
!ESC!    read(IO_FID_CONF,nml=nm_mp_nsw6,iostat=ierr)
!ESC!    if(ierr<0) then
!ESC!       write(IO_FID_LOG,*) '*** nm_mp_nsw6 is not specified. use default.'
!ESC!    else if(ierr>0) then
!ESC!       write(*,         *) 'xxx Not appropriate names in namelist nm_mp_nsw6. STOP.'
!ESC!       write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist nm_mp_nsw6. STOP.'
!ESC!       call ADM_proc_stop
!ESC!    endif
!ESC!    write(IO_FID_LOG,nml=nm_mp_nsw6)
gamma_sacr=6.E-2_RP
psaut_beta0=5.E-3_RP
gamma_saut=6.E-2_RP
precip_scheme_type="Flux-Semilag_new"
roh_flag=.true.
qr_aut_acc_type="KK2000"

    beta_saut = PSAUT_BETA0

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '*** Calculation flag of sedimentation:'
    write(IO_FID_LOG,*) '*** QV => NO'
    write(IO_FID_LOG,*) '*** QC => NO'
    write(IO_FID_LOG,*) '*** QR => YES'
    if    ( precip_transport_type == '3WATER' ) then
       write(IO_FID_LOG,*) '*** QI => NO'
    elseif( precip_transport_type == '4WATER' ) then
       write(IO_FID_LOG,*) '*** QI => YES'
    endif
    write(IO_FID_LOG,*) '*** QS => YES'
    write(IO_FID_LOG,*) '*** QG => YES'

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '*** Precipitation(sedimentation) scheme:'
    if    ( precip_scheme_type == 'Upwind-Euler' ) then
       write(IO_FID_LOG,*) '*** => Upwind-Euler'
    elseif( precip_scheme_type == 'Flux-Semilag_new' ) then
       write(IO_FID_LOG,*) '*** => Flux-Semilag_new'
    else
       write(IO_FID_LOG,*) '*** => Default(Flux-Semilag)'
    endif

    !--- empirical coefficients A, B, C, D
    Ar = PI * rho_w / 6.0_RP
    As = PI * rho_s / 6.0_RP
    Ag = PI * rho_g / 6.0_RP

    Br = 3.0_RP
    Bs = 3.0_RP
    Bg = 3.0_RP

    Cg = sqrt( ( 4.0_RP * rho_g * GRAV ) / ( 3.0_RP * dens00 * C_d ) )

    Dr = 0.50_RP
    Ds = 0.25_RP
    Dg = 0.50_RP

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '*** Use setting of Roh and Satoh(2014)?:'
    if ( Roh_flag ) then ! overwrite parameters
       write(IO_FID_LOG,*) '*** => Yes'
       OPT_EXPLICIT_ICEGEN = .true.

       sw_roh2014 = 1.0_RP
       N0g        = 4.E+8_RP
       As         = 0.069_RP
       Bs         = 2.0_RP
       Esi        = 0.25_RP
       Egi        = 0.0_RP
       Egs        = 0.0_RP
    else
       write(IO_FID_LOG,*) '*** => No : default'
    endif

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '*** Use explicit ice generation scheme?:'
    if ( OPT_EXPLICIT_ICEGEN ) then
       write(IO_FID_LOG,*) '*** => Yes'
       sw_expice = 1.0_RP
    else
       write(IO_FID_LOG,*) '*** => No : default'
       sw_expice = 0.0_RP
    endif

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '*** Autoconversion & Accretion scheme for QC->Qr:'
    if    ( qr_aut_acc_type == "Default" ) then
       write(IO_FID_LOG,*) '*** => Berry(1968) : default'
       sw_kk2000 = 0.0_RP
    elseif( qr_aut_acc_type == "KK2000" ) then
       write(IO_FID_LOG,*) '*** => Khairoutdinov and Kogan(2000)'
       sw_kk2000 = 1.0_RP
    else
       write(*,         *) 'xxx Not appropriate qr_aut_acc_type. STOP', trim(qr_aut_acc_type)
       write(IO_FID_LOG,*) 'xxx Not appropriate qr_aut_acc_type. STOP', trim(qr_aut_acc_type)
       call ADM_proc_stop
    endif

    if ( .NOT. INC_PGAUT ) then
       beta_gaut = 0.0_RP
    endif

    if ( CONST_Vti /= UNDEF ) then
       CONST_Vti   = abs(CONST_Vti)
       sw_constVti = 1.0_RP
    else
       sw_constVti = 0.0_RP
    endif

    GAM       = 1.0_RP ! =0!
    GAM_2     = 1.0_RP ! =1!
    GAM_3     = 2.0_RP ! =2!

    GAM_1br   = MISC_gammafunc( 1.0_RP + Br ) ! = 4!
    GAM_2br   = MISC_gammafunc( 2.0_RP + Br ) ! = 5!
    GAM_3br   = MISC_gammafunc( 3.0_RP + Br ) ! = 6!
    GAM_3dr   = MISC_gammafunc( 3.0_RP + Dr )
    GAM_6dr   = MISC_gammafunc( 6.0_RP + Dr )
    GAM_1brdr = MISC_gammafunc( 1.0_RP + Br + Dr )
    GAM_5dr_h = MISC_gammafunc( 0.5_RP * (5.0_RP+Dr) )

    GAM_1bs   = MISC_gammafunc( 1.0_RP + Bs ) ! = 4!
    GAM_2bs   = MISC_gammafunc( 2.0_RP + Bs ) ! = 5!
    GAM_3bs   = MISC_gammafunc( 3.0_RP + Bs ) ! = 6!
    GAM_3ds   = MISC_gammafunc( 3.0_RP + Ds )
    GAM_1bsds = MISC_gammafunc( 1.0_RP + Bs + Ds )
    GAM_5ds_h = MISC_gammafunc( 0.5_RP * (5.0_RP+Ds) )

    GAM_1bg   = MISC_gammafunc( 1.0_RP + Bg ) ! = 4!
    GAM_3dg   = MISC_gammafunc( 3.0_RP + Dg )
    GAM_1bgdg = MISC_gammafunc( 1.0_RP + Bg + Dg)
    GAM_5dg_h = MISC_gammafunc( 0.5_RP * (5.0_RP+Dg) )

    ln10 = log(10.0_RP)

    return
  end subroutine mp_nsw6_init

  !-----------------------------------------------------------------------------
  subroutine mp_nsw6( &
       ijdim,          &
       l_region,       &
       rhog,           &
       rhogvx,         &
       rhogvy,         &
       rhogvz,         &
       rhogw,          &
       rhoge,          &
       rhogq,          &
       vx,             &
       vy,             &
       vz,             &
       w,              &
       UNCCN,          &
       rho,            &
       tem,            &
       pre,            &
       q,              &
       qd,             &
       precip,         &
       precip_rhoe,    &
       precip_lh_heat, &
       precip_rhophi,  &
       precip_rhokin,  &
       gprec,          &
       rceff,          &
       rctop,          &
       rwtop,          &
       tctop,          &
       rceff_cld,      &
       rctop_cld,      &
       rwtop_cld,      &
       tctop_cld,      &
       gsgam2,         &
       gsgam2h,        &
       gam2,           &
       gam2h,          &
       ix,             &
       iy,             &
       iz,             &
       jx,             &
       jy,             &
       jz,             &
       z,              &
       dt              )
!ESC!    use mod_misc, only: &
!ESC!       MISC_gammafunc
!ESC!    use mod_adm, only: &
!ESC!       kdim => ADM_kall, &
!ESC!       kmin => ADM_kmin, &
!ESC!       kmax => ADM_kmax
!ESC!    use mod_const, only: &
!ESC!       CONST_UNDEF,          &
!ESC!       CONST_EPS,            &
!ESC!       CONST_PI,             &
!ESC!       CONST_Rvap,           &
!ESC!       CL => CONST_CL,       &
!ESC!       CONST_LHV0,           &
!ESC!       CONST_LHS0,           &
!ESC!       CONST_LHF0,           &
!ESC!       rho_w => CONST_DWATR, &
!ESC!       TEM00 => CONST_TEM00, &
!ESC!       CONST_Pstd
!ESC!    use mod_runconf, only: &
!ESC!       LHV,               &
!ESC!       LHF,               &
!ESC!       nqmax => TRC_VMAX, &
!ESC!       NQW_STR,           &
!ESC!       NQW_END,           &
!ESC!       I_QV,              &
!ESC!       I_QC,              &
!ESC!       I_QR,              &
!ESC!       I_QI,              &
!ESC!       I_QS,              &
!ESC!       I_QG
    use mod_thrmdyn, only: &
       thrmdyn_cv, &
       thrmdyn_qd
    use mod_satadjust, only: &
       SATURATION_psat_liq,  &
       SATURATION_psat_ice,  &
       SATURATION_setrange,  &
       SATURATION_adjustment
    use mod_precip_transport, only: &
       precip_transport_new
!ESC!       precip_transport_nwater, &
!ESC!       precip_transport_new,    &
!ESC!       precip_transport_euler
!ESC!    use mod_history, only: &
!ESC!       history_in
    implicit none

    integer,  intent(in)    :: ijdim
    integer,  intent(in)    :: l_region
    real(RP), intent(inout) :: rhog          (ijdim,kdim)
    real(RP), intent(inout) :: rhogvx        (ijdim,kdim)
    real(RP), intent(inout) :: rhogvy        (ijdim,kdim)
    real(RP), intent(inout) :: rhogvz        (ijdim,kdim)
    real(RP), intent(inout) :: rhogw         (ijdim,kdim)
    real(RP), intent(inout) :: rhoge         (ijdim,kdim)
    real(RP), intent(inout) :: rhogq         (ijdim,kdim,nqmax)
    real(RP), intent(in)    :: vx            (ijdim,kdim)
    real(RP), intent(in)    :: vy            (ijdim,kdim)
    real(RP), intent(in)    :: vz            (ijdim,kdim)
    real(RP), intent(in)    :: w             (ijdim,kdim)
    real(RP), intent(in)    :: UNCCN         (ijdim,kdim)
    real(RP), intent(inout) :: rho           (ijdim,kdim)
    real(RP), intent(inout) :: tem           (ijdim,kdim)
    real(RP), intent(inout) :: pre           (ijdim,kdim)
    real(RP), intent(inout) :: q             (ijdim,kdim,nqmax)
    real(RP), intent(out)   :: qd            (ijdim,kdim)
    real(RP), intent(out)   :: precip        (ijdim,2)
    real(RP), intent(out)   :: precip_rhoe   (ijdim)
    real(RP), intent(out)   :: precip_lh_heat(ijdim)
    real(RP), intent(out)   :: precip_rhophi (ijdim)
    real(RP), intent(out)   :: precip_rhokin (ijdim)
    real(RP), intent(out)   :: gprec         (ijdim,kdim)
    real(RP), intent(out)   :: rceff         (ijdim,kdim)
    real(RP), intent(out)   :: rctop         (ijdim,1)
    real(RP), intent(out)   :: rwtop         (ijdim,1)
    real(RP), intent(out)   :: tctop         (ijdim,1)
    real(RP), intent(out)   :: rceff_cld     (ijdim,kdim)
    real(RP), intent(out)   :: rctop_cld     (ijdim,1)
    real(RP), intent(out)   :: rwtop_cld     (ijdim,1)
    real(RP), intent(out)   :: tctop_cld     (ijdim,1)
    real(RP), intent(in)    :: gsgam2        (ijdim,kdim)
    real(RP), intent(in)    :: gsgam2h       (ijdim,kdim)
    real(RP), intent(in)    :: gam2          (ijdim,kdim)
    real(RP), intent(in)    :: gam2h         (ijdim,kdim)
    real(RP), intent(in)    :: ix            (ijdim)
    real(RP), intent(in)    :: iy            (ijdim)
    real(RP), intent(in)    :: iz            (ijdim)
    real(RP), intent(in)    :: jx            (ijdim)
    real(RP), intent(in)    :: jy            (ijdim)
    real(RP), intent(in)    :: jz            (ijdim)
    real(RP), intent(in)    :: z             (ijdim,kdim)
    real(RP), intent(in)    :: dt

    ! working
    real(RP) :: drhogqv(ijdim,kdim)
    real(RP) :: drhogqc(ijdim,kdim)
    real(RP) :: drhogqi(ijdim,kdim)
    real(RP) :: drhogqr(ijdim,kdim)
    real(RP) :: drhogqs(ijdim,kdim)
    real(RP) :: drhogqg(ijdim,kdim)

    real(RP) :: psatl(ijdim,kdim)
    real(RP) :: qsatl(ijdim,kdim)                !< saturated water vapor for liquid water [kg/kg]
    real(RP) :: psati(ijdim,kdim)
    real(RP) :: qsati(ijdim,kdim)                !< saturated water vapor for ice water    [kg/kg]
    real(RP) :: Nc   (ijdim,kdim)                !< Number concentration of cloud water [1/cc]

    real(RP) :: dens                             !< density
    real(RP) :: temp                             !< T [K]
    real(RP) :: qv                               !< mixing ratio of water vapor  [kg/kg]
    real(RP) :: qc                               !< mixing ratio of liquid water [kg/kg]
    real(RP) :: qr                               !< mixing ratio of rain         [kg/kg]
    real(RP) :: qi                               !< mixing ratio of ice water    [kg/kg]
    real(RP) :: qs                               !< mixing ratio of snow         [kg/kg]
    real(RP) :: qg                               !< mixing ratio of graupel      [kg/kg]
    real(RP) :: qv_t                             !< tendency     of water vapor  [kg/kg/s]
    real(RP) :: qc_t                             !< tendency     of liquid water [kg/kg/s]
    real(RP) :: qr_t                             !< tendency     of rain         [kg/kg/s]
    real(RP) :: qi_t                             !< tendency     of ice water    [kg/kg/s]
    real(RP) :: qs_t                             !< tendency     of snow         [kg/kg/s]
    real(RP) :: qg_t                             !< tendency     of graupel      [kg/kg/s]
    real(RP) :: Sliq                             !< saturated ratio S for liquid water [0-1]
    real(RP) :: Sice                             !< saturated ratio S for ice water    [0-1]
    real(RP) :: Rdens                            !< 1 / density
    real(RP) :: rho_fact                         !< density factor
    real(RP) :: temc                             !< T - T0 [K]

    real(RP) :: RLMDr, RLMDr_2, RLMDr_3
    real(RP) :: RLMDs, RLMDs_2, RLMDs_3
    real(RP) :: RLMDg, RLMDg_2, RLMDg_3
    real(RP) :: RLMDr_1br, RLMDr_2br, RLMDr_3br
    real(RP) :: RLMDs_1bs, RLMDs_2bs, RLMDs_3bs
    real(RP) :: RLMDr_dr, RLMDr_3dr, RLMDr_5dr
    real(RP) :: RLMDs_ds, RLMDs_3ds, RLMDs_5ds
    real(RP) :: RLMDg_dg, RLMDg_3dg, RLMDg_5dg
    real(RP) :: RLMDr_7
    real(RP) :: RLMDr_6dr

    !---< Roh and Satoh (2014) >---
    real(RP) :: tems, Xs2
    real(RP) :: MOMs_0, MOMs_1, MOMs_2
    real(RP) :: MOMs_0bs, MOMs_1bs, MOMs_2bs
    real(RP) :: MOMs_2ds, MOMs_5ds_h, RMOMs_Vt
    real(RP) :: coef_at(4), coef_bt(4)
    real(RP) :: loga_, b_, nm

    real(RP) :: Vti, Vtr, Vts, Vtg               !< terminal velocity
    real(RP) :: Esi_mod, Egs_mod                 !< modified accretion efficiency
    real(RP) :: rhoqc                            !< rho * qc
    real(RP) :: Pracw_orig,  Pracw_kk            !< accretion       term by orig  & k-k scheme
    real(RP) :: Praut_berry, Praut_kk            !< auto-conversion term by berry & k-k scheme
    real(RP) :: Dc                               !< relative variance
    real(RP) :: betai, betas                     !< sticky parameter for auto-conversion
    real(RP) :: Ka                               !< thermal diffusion coefficient of air
    real(RP) :: Kd                               !< diffusion coefficient of water vapor in air
    real(RP) :: Nu                               !< kinematic viscosity of air
    real(RP) :: Glv, Giv, Gil                    !< thermodynamic function
    real(RP) :: ventr, vents, ventg              !< ventilation factor
    real(RP) :: net, fac, fac_sw
    real(RP) :: zerosw, tmp

    !---< Bergeron process >---
    real(RP) :: sw_bergeron                      !< if 0C<T<30C, sw=1
    real(RP) :: a1 (ijdim,kdim)                  !<
    real(RP) :: a2 (ijdim,kdim)                  !<
    real(RP) :: ma2(ijdim,kdim)                  !< 1-a2
    real(RP) :: dt1                              !< time during which the an ice particle of 40um grows to 50um
    real(RP) :: Ni50                             !< number concentration of ice particle of 50um

    !---< Explicit ice generation >---
    real(RP) :: sw, rhoqi, XNi, XMi, Di, Ni0, Qi0

    !---< Effective radius >---
    real(RP), parameter :: r2_min = 1.E-10_RP
    real(RP), parameter :: q_min  = 1.E-5_RP
    real(RP)            :: xf_qc                 !< mean mass   of qc [kg]
    real(RP)            :: rf_qc                 !< mean radius of qc [m]
    real(RP)            :: r2_qc                 !< r^2 moment  of qc
    real(RP)            :: r2_qr                 !< r^2 moment  of qr
    real(RP)            :: r3_qc                 !< r^3 moment  of qc
    real(RP)            :: r3_qr                 !< r^3 moment  of qr
    real(RP)            :: dgamma_a, GAM_dgam23, GAM_dgam
    real(RP)            :: coef_dgam, coef_xf

    !---< Precipitation >---
    logical  :: preciptation_flag(nqmax)

    real(RP) :: Vt    (ijdim,kdim,nqmax)
    real(RP) :: cva   (ijdim,kdim)
    real(RP) :: rgs   (ijdim,kdim)
    real(RP) :: rgsh  (ijdim,kdim)

    real(RP) :: wk      (wk_nmax)
!    real(RP) :: ml_wk   (ijdim,kdim,wk_nmax) tentatively remove due to the stack overflow
    real(RP) :: ml_Pconv(ijdim,kdim)
    real(RP) :: ml_Pconw(ijdim,kdim)
    real(RP) :: ml_Pconi(ijdim,kdim)

    real(RP) :: UNDEF, EPS, PI, Rvap, LHV0, LHS0, LHF0, PRE00

    integer, parameter :: simdlen = 8
    integer  :: blk, vec, veclen

    integer  :: ij, k, nq, ip
    !---------------------------------------------------------------------------

    call PROF_rapstart('____MP_NSW6')

    UNDEF = CONST_UNDEF
    EPS   = CONST_EPS
    PI    = CONST_PI
    Rvap  = CONST_Rvap
    LHV0  = CONST_LHV0
    LHS0  = CONST_LHS0
    LHF0  = CONST_LHF0
    PRE00 = CONST_Pstd

    call negative_filter( ijdim,         & ! [IN]
                          rhog  (:,:),   & ! [INOUT]
                          rhoge (:,:),   & ! [INOUT]
                          rhogq (:,:,:), & ! [INOUT]
                          rho   (:,:),   & ! [INOUT]
                          tem   (:,:),   & ! [IN]
                          pre   (:,:),   & ! [INOUT]
                          q     (:,:,:), & ! [INOUT]
                          gsgam2(:,:)    ) ! [IN]

    if ( OPT_INDIR ) then  !  aerosol indirect effect
       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kdim,Nc,UNCCN,Nc_def), &
       !$omp collapse(2)
       do k  = 1, kdim
       do ij = 1, ijdim
          Nc(ij,k) = max( UNCCN(ij,k)*1.E-6_RP, Nc_def ) ! [#/m3]->[#/cc]
       enddo
       enddo
       !$omp end parallel do
    else
       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kdim,Nc,Nc_def), &
       !$omp collapse(2)
       do k  = 1, kdim
       do ij = 1, ijdim
          Nc(ij,k) = Nc_def
       enddo
       enddo
       !$omp end parallel do
    endif

    ! saturation water contensts
    call SATURATION_psat_liq( ijdim, kdim, tem(:,:), psatl(:,:) )
    call SATURATION_psat_ice( ijdim, kdim, tem(:,:), psati(:,:) )

    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,qsatl,qsati,psatl,psati,rho,tem,Rvap), &
    !$omp collapse(2)
    do k  = 1, kdim
    do ij = 1, ijdim
       qsatl(ij,k) = psatl(ij,k) / ( rho(ij,k) * Rvap * tem(ij,k) )
       qsati(ij,k) = psati(ij,k) / ( rho(ij,k) * Rvap * tem(ij,k) )
    enddo
    enddo
    !$omp end parallel do

    ! Bergeron process parameters
    call Bergeron_param( ijdim,    & ! [IN]
                         tem(:,:), & ! [IN]
                         a1 (:,:), & ! [OUT]
                         a2 (:,:), & ! [OUT]
                         ma2(:,:)  ) ! [OUT]

    ! work for Effective Radius of Liquid Water
    ! dgamma_a is parameter of Gamma-Dist.
    ! "Berry and Reinhardt(1974-a) eq.(16-a,b)"
    ! 0.38 * {var1/2x} =~ {var1/2 r}
    !         dgamma_a = 1.0_RP /{var x}
    Dc         = 0.146_RP - 5.964E-2_RP * log( Nc_def / 2000.0_RP )
    dgamma_a   = 0.1444_RP / Dc**2
    GAM_dgam   = MISC_gammafunc( dgamma_a )
    GAM_dgam23 = MISC_gammafunc( dgamma_a + 2.0_RP/3.0_RP )
    coef_dgam  = GAM_dgam23 / GAM_dgam * dgamma_a**(-2.0_RP/3.0_RP)
    coef_xf    = 3.0_RP / 4.0_RP / PI / rho_w

if ( TIME_CSTEP == EX_STEP ) then
   call PROF_rapstart('_kernel_partM12')
endif

    !$omp parallel do default(none), &
    !$omp private(ij,k,blk,vec,veclen)                                                                        &
    !$omp private(ip,dens,temp,qv,qc,qr,qi,qs,qg,Sliq,Sice,Rdens,rho_fact,temc,sw_bergeron,zerosw,tmp,        &
    !$omp         RLMDr,RLMDr_dr,RLMDr_2,RLMDr_3,RLMDr_7,RLMDr_1br,RLMDr_2br,RLMDr_3br,RLMDr_3dr,             &
    !$omp         RLMDr_5dr,RLMDr_6dr,RLMDs,RLMDs_ds,RLMDs_2,RLMDs_3,RLMDs_1bs,RLMDs_2bs,RLMDs_3bs,           &
    !$omp         RLMDs_3ds,RLMDs_5ds,RLMDg,RLMDg_dg,RLMDg_2,RLMDg_3,RLMDg_3dg,RLMDg_5dg,                     &
    !$omp         MOMs_0,MOMs_1,MOMs_2,MOMs_0bs,MOMs_1bs,MOMs_2bs,MOMs_2ds,MOMs_5ds_h,RMOMs_Vt,               &
    !$omp         coef_bt,coef_at,Xs2,tems,loga_,b_,nm)                                                       &
    !$omp private(Vti,Vtr,Vts,Vtg,Esi_mod,Egs_mod,Pracw_orig,Pracw_kk,rhoqc,Dc,Praut_berry,Praut_kk,          &
    !$omp         betai,betas,Ka,Kd,NU,Glv,Giv,Gil,ventr,vents,ventg,dt1,Ni50,sw,rhoqi,XNi,XMi,Di,Ni0,Qi0,    &
    !$omp         xf_qc,rf_qc,r2_qc,r2_qr,r3_qc,r3_qr,                                                        &
    !$omp         net,fac,fac_sw,qv_t,qc_t,qr_t,qi_t,qs_t,qg_t,wk)                                            &
!    !$omp private(ml_wk)                                                                                     &
    !$omp shared(ijdim,kmin,kmax,rho,tem,pre,q,qsatl,qsati,Nc,dt,sw_expice,sw_roh2014,sw_kk2000,sw_constVti,  &
    !$omp        drhogqv,drhogqc,drhogqr,drhogqi,drhogqs,drhogqg,rhog,Vt,CONST_Vti,                           &
    !$omp        coef_dgam,coef_xf,rceff,rctop,tctop,rceff_cld,rctop_cld,tctop_cld,                           &
    !$omp        I_QV,I_QC,I_QR,I_QI,I_QS,I_QG,UNDEF,EPS,PI,Rvap,LHV0,LHS0,LHF0,PRE00)                        &
    !$omp shared(GAM,GAM_2,GAM_3,GAM_1br,GAM_2br,GAM_3br,GAM_3dr,GAM_6dr,GAM_1brdr,GAM_5dr_h,GAM_1bs,GAM_2bs, &
    !$omp        GAM_3bs,GAM_3ds,GAM_1bsds,GAM_5ds_h,GAM_1bg,GAM_3dg,GAM_1bgdg,GAM_5dg_h,coef_a,coef_b,ln10,  &
    !$omp        N0r,N0s,N0g,Ar,As,Ag,Br,Bs,Bg,Cr,Cs,Cg,Dr,Ds,Dg,Eiw,Erw,Esw,Egw,Eri,Esi,Egi,Esr,Egr,Egs,     &
    !$omp        gamma_sacr,gamma_gacs,mi,beta_saut,gamma_saut,beta_gaut,gamma_gaut,qicrt_saut,qscrt_gaut,    &
    !$omp        f1r,f2r,f1s,f2s,f1g,f2g,A_frz,B_frz,mi40,mi50,vti50,Ri50,a1,a2,ma2)                          &
    !$omp collapse(2)
    do k  = kmin, kmax
#if defined _OPENMP && _OPENMP >= 201511
    do blk = 1, ijdim, simdlen
       veclen = min( ijdim-blk+1, simdlen )
       !$omp simd simdlen(8)
       do vec = 1, veclen
          ij = blk + vec - 1
#else
    do ij = 1, ijdim
#endif

          dens     = rho(ij,k)
          temp     = tem(ij,k)
          qv       = max( q(ij,k,I_QV), 0.0_RP )
          qc       = max( q(ij,k,I_QC), 0.0_RP )
          qr       = max( q(ij,k,I_QR), 0.0_RP )
          qi       = max( q(ij,k,I_QI), 0.0_RP )
          qs       = max( q(ij,k,I_QS), 0.0_RP )
          qg       = max( q(ij,k,I_QG), 0.0_RP )

          ! saturation ratio S
          Sliq     = qv / max( qsatl(ij,k), EPS )
          Sice     = qv / max( qsati(ij,k), EPS )

          Rdens    = 1.0_RP / dens
          rho_fact = sqrt( dens00 * Rdens )
          temc     = temp - TEM00

          wk(I_delta1) = ( 0.5_RP + sign(0.5_RP, qr - 1.E-4_RP ) )

          wk(I_delta2) = ( 0.5_RP + sign(0.5_RP, 1.E-4_RP - qr ) ) &
                       * ( 0.5_RP + sign(0.5_RP, 1.E-4_RP - qs ) )

          wk(I_spsati) = 0.5_RP + sign(0.5_RP, Sice - 1.0_RP )

          wk(I_iceflg) = 0.5_RP - sign( 0.5_RP, temc ) ! 0: warm, 1: ice

          wk(I_dqv_dt) = qv / dt
          wk(I_dqc_dt) = qc / dt
          wk(I_dqr_dt) = qr / dt
          wk(I_dqi_dt) = qi / dt
          wk(I_dqs_dt) = qs / dt
          wk(I_dqg_dt) = qg / dt

          sw_bergeron = ( 0.5_RP + sign(0.5_RP, temc + 30.0_RP ) ) &
                      * ( 0.5_RP + sign(0.5_RP, 0.0_RP - temc  ) ) &
                      * ( 1.0_RP - sw_expice                     )

          ! slope parameter lambda (Rain)
          zerosw = 0.5_RP - sign(0.5_RP, qr - 1.E-12_RP )
          RLMDr  = sqrt(sqrt( dens * qr / ( Ar * N0r * GAM_1br ) + zerosw )) * ( 1.0_RP - zerosw )

          RLMDr_dr  = sqrt( RLMDr )       ! **Dr
          RLMDr_2   = RLMDr**2
          RLMDr_3   = RLMDr**3
          RLMDr_7   = RLMDr**7
          RLMDr_1br = RLMDr**4 ! (1+Br)
          RLMDr_2br = RLMDr**5 ! (2+Br)
          RLMDr_3br = RLMDr**6 ! (3+Br)
          RLMDr_3dr = RLMDr**3 * RLMDr_dr
          RLMDr_5dr = RLMDr**5 * RLMDr_dr
          RLMDr_6dr = RLMDr**6 * RLMDr_dr

          ! slope parameter lambda (Snow)
          zerosw = 0.5_RP - sign(0.5_RP, qs - 1.E-12_RP )
          RLMDs  = sqrt(sqrt( dens * qs / ( As * N0s * GAM_1bs ) + zerosw )) * ( 1.0_RP - zerosw )

          RLMDs_ds  = sqrt( sqrt(RLMDs) ) ! **Ds
          RLMDs_2   = RLMDs**2
          RLMDs_3   = RLMDs**3
          RLMDs_1bs = RLMDs**4 ! (1+Bs)
          RLMDs_2bs = RLMDs**5 ! (2+Bs)
          RLMDs_3bs = RLMDs**6 ! (3+Bs)
          RLMDs_3ds = RLMDs**3 * RLMDs_ds
          RLMDs_5ds = RLMDs**5 * RLMDs_ds

          MOMs_0     = N0s * GAM       * RLMDs           ! Ns * 0th moment
          MOMs_1     = N0s * GAM_2     * RLMDs_2         ! Ns * 1st moment
          MOMs_2     = N0s * GAM_3     * RLMDs_3         ! Ns * 2nd moment
          MOMs_0bs   = N0s * GAM_1bs   * RLMDs_1bs       ! Ns * 0+bs
          MOMs_1bs   = N0s * GAM_2bs   * RLMDs_2bs       ! Ns * 1+bs
          MOMs_2bs   = N0s * GAM_3bs   * RLMDs_3bs       ! Ns * 2+bs
          MOMs_2ds   = N0s * GAM_3ds   * RLMDs_3ds       ! Ns * 2+ds
          MOMs_5ds_h = N0s * GAM_5ds_h * sqrt(RLMDs_5ds) ! Ns * (5+ds)/2
          RMOMs_Vt   = GAM_1bsds / GAM_1bs * RLMDs_ds

          !---< modification by Roh and Satoh (2014) >---
          ! bimodal size distribution of snow
          Xs2    = dens * qs / As
          zerosw = 0.5_RP - sign(0.5_RP, Xs2 - 1.E-12_RP )

          tems       = min( -0.1_RP, temc )
          coef_at(1) = coef_a( 1) + tems * ( coef_a( 2) + tems * ( coef_a( 5) + tems * coef_a( 9) ) )
          coef_at(2) = coef_a( 3) + tems * ( coef_a( 4) + tems *   coef_a( 7) )
          coef_at(3) = coef_a( 6) + tems *   coef_a( 8)
          coef_at(4) = coef_a(10)
          coef_bt(1) = coef_b( 1) + tems * ( coef_b( 2) + tems * ( coef_b( 5) + tems * coef_b( 9) ) )
          coef_bt(2) = coef_b( 3) + tems * ( coef_b( 4) + tems *   coef_b( 7) )
          coef_bt(3) = coef_b( 6) + tems *   coef_b( 8)
          coef_bt(4) = coef_b(10)
          ! 0th moment
          loga_ = coef_at(1)
          b_    = coef_bt(1)
          MOMs_0 = (        sw_roh2014 ) * exp(ln10*loga_) * exp(log(Xs2+zerosw)*b_) * ( 1.0_RP-zerosw ) &
                 + ( 1.0_RP-sw_roh2014 ) * MOMs_0
          ! 1st moment
          nm = 1.0_RP
          loga_ = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
          b_    = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
          MOMs_1 = (        sw_roh2014 ) * exp(ln10*loga_) * exp(log(Xs2+zerosw)*b_) * ( 1.0_RP-zerosw ) &
                 + ( 1.0_RP-sw_roh2014 ) * MOMs_1
          ! 2nd moment
          MOMs_2 = (        sw_roh2014 ) * Xs2 &
                 + ( 1.0_RP-sw_roh2014 ) * MOMs_2
          ! 0 + Bs(=2) moment
          nm = 2.0_RP
          loga_ = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
          b_    = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
          MOMs_0bs = (        sw_roh2014 ) * exp(ln10*loga_) * exp(log(Xs2+zerosw)*b_) * ( 1.0_RP-zerosw ) &
                   + ( 1.0_RP-sw_roh2014 ) * MOMs_0bs
          ! 1 + Bs(=2) moment
          nm = 3.0_RP
          loga_ = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
          b_    = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
          MOMs_1bs = (        sw_roh2014 ) * exp(ln10*loga_) * exp(log(Xs2+zerosw)*b_) * ( 1.0_RP-zerosw ) &
                   + ( 1.0_RP-sw_roh2014 ) * MOMs_1bs
          ! 2 + Bs(=2) moment
          nm = 4.0_RP
          loga_ = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
          b_    = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
          MOMs_2bs = (        sw_roh2014 ) * exp(ln10*loga_) * exp(log(Xs2+zerosw)*b_) * ( 1.0_RP-zerosw ) &
               + ( 1.0_RP-sw_roh2014 ) * MOMs_2bs
          ! 2 + Ds(=0.25) moment
          nm = 2.25_RP
          loga_ = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
          b_    = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
          MOMs_2ds = (        sw_roh2014 ) * exp(ln10*loga_) * exp(log(Xs2+zerosw)*b_) * ( 1.0_RP-zerosw ) &
                   + ( 1.0_RP-sw_roh2014 ) * MOMs_2ds
          ! ( 3 + Ds(=0.25) ) / 2  moment
          nm = 1.625_RP
          loga_ = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
          b_    = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
          MOMs_5ds_h = (        sw_roh2014 ) * exp(ln10*loga_) * exp(log(Xs2+zerosw)*b_) * ( 1.0_RP-zerosw ) &
                     + ( 1.0_RP-sw_roh2014 ) * MOMs_5ds_h
          ! Bs(=2) + Ds(=0.25) moment
          nm = 2.25_RP
          loga_ = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
          b_    = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )

          RMOMs_Vt = (        sw_roh2014 ) * exp(ln10*loga_) * exp(log(Xs2+zerosw)*b_) * ( 1.0_RP-zerosw ) / ( MOMs_0bs + zerosw ) &
                   + ( 1.0_RP-sw_roh2014 ) * RMOMs_Vt

          ! slope parameter lambda (Graupel)
          zerosw = 0.5_RP - sign(0.5_RP, qg - 1.E-12_RP )
          RLMDg  = sqrt(sqrt( dens * qg / ( Ag * N0g * GAM_1bg ) + zerosw )) * ( 1.0_RP - zerosw )

          RLMDg_dg  = sqrt( RLMDg )       ! **Dg
          RLMDg_2   = RLMDg**2
          RLMDg_3   = RLMDg**3
          RLMDg_3dg = RLMDg**3 * RLMDg_dg
          RLMDg_5dg = RLMDg**5 * RLMDg_dg

          wk(I_RLMDr) = RLMDr
          wk(I_RLMDs) = RLMDs
          wk(I_RLMDg) = RLMDg

          !---< terminal velocity >---
          zerosw = 0.5_RP - sign(0.5_RP, qi - 1.E-8_RP )
          Vti = ( 1.0_RP-sw_constVti ) * (-3.29_RP) * exp( log( dens*qi+zerosw )*0.16_RP ) * ( 1.0_RP-zerosw ) &
              + (        sw_constVti ) * (-CONST_Vti)
          Vtr = -Cr * rho_fact * GAM_1brdr / GAM_1br * RLMDr_dr
          Vts = -Cs * rho_fact * RMOMs_Vt
          Vtg = -Cg * rho_fact * GAM_1bgdg / GAM_1bg * RLMDg_dg

          !---< Nucleation >---
          ! [Pigen] ice nucleation
          Ni0 = max( exp(-0.1_RP*temc), 1.0_RP ) * 1000.0_RP
          Qi0 = 4.92E-11_RP * exp(log(Ni0)*1.33_RP) * Rdens

          wk(I_Pigen) = max( min( Qi0-qi, qv-qsati(ij,k) ), 0.0_RP ) / dt

          !---< Accretion >---
          Esi_mod = min( Esi, Esi * exp( gamma_sacr * temc ) )
          Egs_mod = min( Egs, Egs * exp( gamma_gacs * temc ) )

          ! [Pracw] accretion rate of cloud water by rain
          Pracw_orig = qc * 0.25_RP * PI * Erw * N0r * Cr * GAM_3dr * RLMDr_3dr * rho_fact

          zerosw     = 0.5_RP - sign(0.5_RP, qc*qr - 1.E-12_RP )
          Pracw_kk   = 67.0_RP * exp( log( qc*qr+zerosw )*1.15_RP ) * ( 1.0_RP-zerosw ) ! eq.(33) in KK(2000)

          ! switch orig / k-k scheme
          wk(I_Pracw) = ( 1.0_RP-sw_kk2000 ) * Pracw_orig &
                      + (        sw_kk2000 ) * Pracw_kk

          ! [Psacw] accretion rate of cloud water by snow
          wk(I_Psacw) = qc * 0.25_RP * PI * Esw       * Cs * MOMs_2ds            * rho_fact

          ! [Pgacw] accretion rate of cloud water by graupel
          wk(I_Pgacw) = qc * 0.25_RP * PI * Egw * N0g * Cg * GAM_3dg * RLMDg_3dg * rho_fact

          ! [Praci] accretion rate of cloud ice by rain
          wk(I_Praci) = qi * 0.25_RP * PI * Eri * N0r * Cr * GAM_3dr * RLMDr_3dr * rho_fact

          ! [Psaci] accretion rate of cloud ice by snow
          wk(I_Psaci) = qi * 0.25_RP * PI * Esi_mod   * Cs * MOMs_2ds            * rho_fact

          ! [Pgaci] accretion rate of cloud ice by grupel
          wk(I_Pgaci) = qi * 0.25_RP * PI * Egi * N0g * Cg * GAM_3dg * RLMDg_3dg * rho_fact

          ! [Piacr] accretion rate of rain by cloud ice
          wk(I_Piacr) = qi * Ar / Mi * 0.25_RP * PI * Eri * N0r * Cr * GAM_6dr * RLMDr_6dr * rho_fact

          ! [Psacr] accretion rate of rain by snow
          wk(I_Psacr) = Ar * 0.25_RP * PI * Rdens * Esr * N0r       * abs(Vtr-Vts) &
                      * (          GAM_1br * RLMDr_1br * MOMs_2          &
                        + 2.0_RP * GAM_2br * RLMDr_2br * MOMs_1          &
                        +          GAM_3br * RLMDr_3br * MOMs_0          )

          ! [Pgacr] accretion rate of rain by graupel
          wk(I_Pgacr) = Ar * 0.25_RP * PI * Rdens * Egr * N0g * N0r * abs(Vtg-Vtr) &
                      * (          GAM_1br * RLMDr_1br * GAM_3 * RLMDg_3 &
                        + 2.0_RP * GAM_2br * RLMDr_2br * GAM_2 * RLMDg_2 &
                        +          GAM_3br * RLMDr_3br * GAM   * RLMDg   )

          ! [Pracs] accretion rate of snow by rain
          wk(I_Pracs) = As * 0.25_RP * PI * Rdens * Esr       *  N0r * abs(Vtr-Vts) &
                      * (          MOMs_0bs            * GAM_3 * RLMDr_3 &
                        + 2.0_RP * MOMs_1bs            * GAM_2 * RLMDr_2 &
                        +          MOMs_2bs            * GAM   * RLMDr   )

          ! [Pgacs] accretion rate of snow by graupel
          wk(I_Pgacs) = As * 0.25_RP * PI * Rdens * Egs_mod   * N0g * abs(Vtg-Vts) &
                      * (          MOMs_0bs            * GAM_3 * RLMDg_3 &
                        + 2.0_RP * MOMs_1bs            * GAM_2 * RLMDg_2 &
                        +          MOMs_2bs            * GAM   * RLMDg   )

          !---< Autoconversion >---

          ! [Praut] auto-conversion rate from cloud water to rain
          rhoqc = dens * qc * 1000.0_RP ! [g/m3]
          Dc    = 0.146_RP - 5.964E-2_RP * log( Nc(ij,k) / 2000.0_RP )
          Praut_berry = Rdens * 1.67E-5_RP * rhoqc * rhoqc / ( 5.0_RP + 3.66E-2_RP * Nc(ij,k) / ( Dc * rhoqc + EPS ) )

          zerosw      = 0.5_RP - sign(0.5_RP, qc - 1.E-12_RP )
          Praut_kk    = 1350.0_RP                                           &
                      * exp( log( qc+zerosw )*2.47_RP ) * ( 1.0_RP-zerosw ) &
                      * exp( log( Nc(ij,k) )*(-1.79_RP) )                     ! eq.(29) in KK(2000)

          Praut_kk    = 1350.0_RP * qc**2.47_RP * Nc(ij,k)**(-1.79_RP)

          ! switch berry / k-k scheme
          wk(I_Praut) = ( 1.0_RP-sw_kk2000 ) * Praut_berry &
                      + (        sw_kk2000 ) * Praut_kk

          ! [Psaut] auto-conversion rate from cloud ice to snow
          betai = min( beta_saut, beta_saut * exp( gamma_saut * temc ) )
          wk(I_Psaut) = max( betai*(qi-qicrt_saut), 0.0_RP )

          ! [Pgaut] auto-conversion rate from snow to graupel
          betas = min( beta_gaut, beta_gaut * exp( gamma_gaut * temc ) )
          wk(I_Pgaut) = max( betas*(qs-qscrt_gaut), 0.0_RP )

          !---< Evaporation, Sublimation, Melting, and Freezing >---

          Ka  = ( Ka0 + dKa_dT * temc )
          Kd  = ( Kd0 + dKd_dT * temc ) * PRE00 / pre(ij,k)
          NU  = ( nu0 + dnu_dT * temc ) * Rdens

          Glv = 1.0_RP / ( LHV0/(Ka*temp) * ( LHV0/(Rvap*temp) - 1.0_RP ) + 1.0_RP/(Kd*dens*QSATL(ij,k)) )
          Giv = 1.0_RP / ( LHS0/(Ka*temp) * ( LHS0/(Rvap*temp) - 1.0_RP ) + 1.0_RP/(Kd*dens*QSATI(ij,k)) )
          Gil = 1.0_RP / ( LHF0/(Ka*temc) )

          ! [Prevp] evaporation rate of rain
          ventr = f1r * GAM_2 * RLMDr_2 + f2r * sqrt( Cr * rho_fact / NU * RLMDr_5dr ) * GAM_5dr_h

          wk(I_Prevp) = 2.0_RP * PI * Rdens * N0r * ( 1.0_RP-min(Sliq,1.0_RP) ) * Glv * ventr

          ! [Pidep,Pisub] deposition/sublimation rate for ice
          rhoqi = max(dens*qi,EPS)
          XNi   = min( max( 5.38E+7_RP * exp( log(rhoqi)*0.75_RP ), 1.E+3_RP ), 1.E+6_RP )
          XMi   = rhoqi / XNi
          Di    = min( Di_a * sqrt(XMi), Di_max )

          tmp = 4.0_RP * Di * XNi * Rdens * ( Sice-1.0_RP ) * Giv

          wk(I_Pidep) = (        wk(I_spsati) ) * ( tmp) ! Sice > 1
          wk(I_Pisub) = ( 1.0_RP-wk(I_spsati) ) * (-tmp) ! Sice < 1

          ! [Pihom] homogenious freezing at T < -40C
          sw = ( 0.5_RP - sign(0.5_RP, temc + 40.0_RP ) ) ! if T < -40C, sw=1

          wk(I_Pihom) = sw * qc / dt

          ! [Pihtr] heteroginous freezing at -40C < T < 0C
          sw = ( 0.5_RP + sign(0.5_RP, temc + 40.0_RP ) ) &
             * ( 0.5_RP - sign(0.5_RP, temc           ) ) ! if -40C < T < 0C, sw=1

          wk(I_Pihtr) = sw * ( dens / rho_w * qc**2 / ( Nc_ihtr * 1.E+6_RP ) ) &
                      * B_frz * ( exp(-A_frz*temc) - 1.0_RP )

          ! [Pimlt] ice melting at T > 0C
          sw = ( 0.5_RP + sign(0.5_RP, temc           ) ) ! if T > 0C, sw=1

          wk(I_Pimlt) = sw * qi / dt

          ! [Psdep,Pssub] deposition/sublimation rate for snow
          vents = f1s * MOMs_1          + f2s * sqrt( Cs * rho_fact / NU             ) * MOMs_5ds_h

          tmp = 2.0_RP * PI * Rdens *       ( Sice-1.0_RP ) * Giv * vents

          wk(I_Psdep) = (        wk(I_spsati) ) * ( tmp) ! Sice > 1
          wk(I_Pssub) = ( 1.0_RP-wk(I_spsati) ) * (-tmp) ! Sice < 1

          ! [Psmlt] melting rate of snow
          wk(I_Psmlt) = 2.0_RP * PI * Rdens *       Gil * vents &
                      + CL * temc / LHF0 * ( wk(I_Psacw) + wk(I_Psacr) )
          wk(I_Psmlt) = max( wk(I_Psmlt), 0.0_RP )

          ! [Pgdep/pgsub] deposition/sublimation rate for graupel
          ventg = f1g * GAM_2 * RLMDg_2 + f2g * sqrt( Cg * rho_fact / NU * RLMDg_5dg ) * GAM_5dg_h

          tmp = 2.0_RP * PI * Rdens * N0g * ( Sice-1.0_RP ) * Giv * ventg

          wk(I_Pgdep) = (        wk(I_spsati) ) * ( tmp) ! Sice > 1
          wk(I_Pgsub) = ( 1.0_RP-wk(I_spsati) ) * (-tmp) ! Sice < 1

          ! [Pgmlt] melting rate of graupel
          wk(I_Pgmlt) = 2.0_RP * PI * Rdens * N0g * Gil * ventg &
                      + CL * temc / LHF0 * ( wk(I_Pgacw) + wk(I_Pgacr) )
          wk(I_Pgmlt) = max( wk(I_Pgmlt), 0.0_RP )

          ! [Pgfrz] freezing rate of graupel
          wk(I_Pgfrz) = 2.0_RP * PI * Rdens * N0r * 60.0_RP * B_frz * Ar * ( exp(-A_frz*temc) - 1.0_RP ) * RLMDr_7

          ! [Psfw,Psfi] ( Bergeron process ) growth rate of snow by Bergeron process from cloud water/ice
          dt1  = ( exp( log(mi50)*ma2(ij,k) ) &
                 - exp( log(mi40)*ma2(ij,k) ) ) / ( a1(ij,k) * ma2(ij,k) )
          Ni50 = qi * dt / ( mi50 * dt1 )

          wk(I_Psfw) = Ni50 * ( a1(ij,k) * mi50**a2(ij,k) &
                     + PI * Eiw * dens * qc * Ri50*Ri50 * vti50 )
          wk(I_Psfi) = qi / dt1

          !---< limiter >---
          wk(I_Pigen) = min( wk(I_Pigen), wk(I_dqv_dt) ) * (        wk(I_iceflg) ) * sw_expice
          wk(I_Pidep) = min( wk(I_Pidep), wk(I_dqv_dt) ) * (        wk(I_iceflg) ) * sw_expice
          wk(I_Psdep) = min( wk(I_Psdep), wk(I_dqv_dt) ) * (        wk(I_iceflg) )
          wk(I_Pgdep) = min( wk(I_Pgdep), wk(I_dqv_dt) ) * (        wk(I_iceflg) )

          wk(I_Pracw) = wk(I_Pracw)                           &
                      + wk(I_Psacw) * ( 1.0_RP-wk(I_iceflg) ) & ! c->r by s
                      + wk(I_Pgacw) * ( 1.0_RP-wk(I_iceflg) )   ! c->r by g

          wk(I_Praut) = min( wk(I_Praut), wk(I_dqc_dt) )
          wk(I_Pracw) = min( wk(I_Pracw), wk(I_dqc_dt) )
          wk(I_Pihom) = min( wk(I_Pihom), wk(I_dqc_dt) ) * (        wk(I_iceflg) ) * sw_expice
          wk(I_Pihtr) = min( wk(I_Pihtr), wk(I_dqc_dt) ) * (        wk(I_iceflg) ) * sw_expice
          wk(I_Psacw) = min( wk(I_Psacw), wk(I_dqc_dt) ) * (        wk(I_iceflg) )
          wk(I_Psfw ) = min( wk(I_Psfw ), wk(I_dqc_dt) ) * (        wk(I_iceflg) ) * sw_bergeron
          wk(I_Pgacw) = min( wk(I_Pgacw), wk(I_dqc_dt) ) * (        wk(I_iceflg) )

          wk(I_Prevp) = min( wk(I_Prevp), wk(I_dqr_dt) )
          wk(I_Piacr) = min( wk(I_Piacr), wk(I_dqr_dt) ) * (        wk(I_iceflg) )
          wk(I_Psacr) = min( wk(I_Psacr), wk(I_dqr_dt) ) * (        wk(I_iceflg) )
          wk(I_Pgacr) = min( wk(I_Pgacr), wk(I_dqr_dt) ) * (        wk(I_iceflg) )
          wk(I_Pgfrz) = min( wk(I_Pgfrz), wk(I_dqr_dt) ) * (        wk(I_iceflg) )

          wk(I_Pisub) = min( wk(I_Pisub), wk(I_dqi_dt) ) * (        wk(I_iceflg) ) * sw_expice
          wk(I_Pimlt) = min( wk(I_Pimlt), wk(I_dqi_dt) ) * ( 1.0_RP-wk(I_iceflg) ) * sw_expice
          wk(I_Psaut) = min( wk(I_Psaut), wk(I_dqi_dt) ) * (        wk(I_iceflg) )
          wk(I_Praci) = min( wk(I_Praci), wk(I_dqi_dt) ) * (        wk(I_iceflg) )
          wk(I_Psaci) = min( wk(I_Psaci), wk(I_dqi_dt) ) * (        wk(I_iceflg) )
          wk(I_Psfi ) = min( wk(I_Psfi ), wk(I_dqi_dt) ) * (        wk(I_iceflg) ) * sw_bergeron
          wk(I_Pgaci) = min( wk(I_Pgaci), wk(I_dqi_dt) ) * (        wk(I_iceflg) )

          wk(I_Pssub) = min( wk(I_Pssub), wk(I_dqs_dt) ) * (        wk(I_iceflg) )
          wk(I_Psmlt) = min( wk(I_Psmlt), wk(I_dqs_dt) ) * ( 1.0_RP-wk(I_iceflg) )
          wk(I_Pgaut) = min( wk(I_Pgaut), wk(I_dqs_dt) ) * (        wk(I_iceflg) )
          wk(I_Pracs) = min( wk(I_Pracs), wk(I_dqs_dt) ) * (        wk(I_iceflg) )
          wk(I_Pgacs) = min( wk(I_Pgacs), wk(I_dqs_dt) )

          wk(I_Pgsub) = min( wk(I_Pgsub), wk(I_dqg_dt) ) * (        wk(I_iceflg) )
          wk(I_Pgmlt) = min( wk(I_Pgmlt), wk(I_dqg_dt) ) * ( 1.0_RP-wk(I_iceflg) )

          wk(I_Piacr_s) = ( 1.0_RP - wk(I_delta1) ) * wk(I_Piacr)
          wk(I_Piacr_g) = (          wk(I_delta1) ) * wk(I_Piacr)
          wk(I_Praci_s) = ( 1.0_RP - wk(I_delta1) ) * wk(I_Praci)
          wk(I_Praci_g) = (          wk(I_delta1) ) * wk(I_Praci)
          wk(I_Psacr_s) = (          wk(I_delta2) ) * wk(I_Psacr)
          wk(I_Psacr_g) = ( 1.0_RP - wk(I_delta2) ) * wk(I_Psacr)
          wk(I_Pracs  ) = ( 1.0_RP - wk(I_delta2) ) * wk(I_Pracs)

          ! [QC]
          net = &
              + wk(I_Pimlt  ) & ! [prod] i->c
              - wk(I_Praut  ) & ! [loss] c->r
              - wk(I_Pracw  ) & ! [loss] c->r
              - wk(I_Pihom  ) & ! [loss] c->i
              - wk(I_Pihtr  ) & ! [loss] c->i
              - wk(I_Psacw  ) & ! [loss] c->s
              - wk(I_Psfw   ) & ! [loss] c->s
              - wk(I_Pgacw  )   ! [loss] c->g

          fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
          fac    = (          fac_sw ) &
                 + ( 1.0_RP - fac_sw ) * min( -wk(I_dqc_dt)/(net-fac_sw), 1.0_RP ) ! loss limiter

          wk(I_Pimlt  ) = wk(I_Pimlt  ) * fac
          wk(I_Praut  ) = wk(I_Praut  ) * fac
          wk(I_Pracw  ) = wk(I_Pracw  ) * fac
          wk(I_Pihom  ) = wk(I_Pihom  ) * fac
          wk(I_Pihtr  ) = wk(I_Pihtr  ) * fac
          wk(I_Psacw  ) = wk(I_Psacw  ) * fac
          wk(I_Psfw   ) = wk(I_Psfw   ) * fac
          wk(I_Pgacw  ) = wk(I_Pgacw  ) * fac

          ! [QI]
          net = &
              + wk(I_Pigen  ) & ! [prod] v->i
              + wk(I_Pidep  ) & ! [prod] v->i
              + wk(I_Pihom  ) & ! [prod] c->i
              + wk(I_Pihtr  ) & ! [prod] c->i
              - wk(I_Pisub  ) & ! [loss] i->v
              - wk(I_Pimlt  ) & ! [loss] i->c
              - wk(I_Psaut  ) & ! [loss] i->s
              - wk(I_Praci_s) & ! [loss] i->s
              - wk(I_Psaci  ) & ! [loss] i->s
              - wk(I_Psfi   ) & ! [loss] i->s
              - wk(I_Praci_g) & ! [loss] i->g
              - wk(I_Pgaci  )   ! [loss] i->g

          fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
          fac    = (          fac_sw ) &
                 + ( 1.0_RP - fac_sw ) * min( -wk(I_dqi_dt)/(net-fac_sw), 1.0_RP ) ! loss limiter

          wk(I_Pigen  ) = wk(I_Pigen  ) * fac
          wk(I_Pidep  ) = wk(I_Pidep  ) * fac
          wk(I_Pihom  ) = wk(I_Pihom  ) * fac
          wk(I_Pihtr  ) = wk(I_Pihtr  ) * fac
          wk(I_Pisub  ) = wk(I_Pisub  ) * fac
          wk(I_Pimlt  ) = wk(I_Pimlt  ) * fac
          wk(I_Psaut  ) = wk(I_Psaut  ) * fac
          wk(I_Praci_s) = wk(I_Praci_s) * fac
          wk(I_Psaci  ) = wk(I_Psaci  ) * fac
          wk(I_Psfi   ) = wk(I_Psfi   ) * fac
          wk(I_Praci_g) = wk(I_Praci_g) * fac
          wk(I_Pgaci  ) = wk(I_Pgaci  ) * fac

          ! [QR]
          net = &
              + wk(I_Praut  ) & ! [prod] c->r
              + wk(I_Pracw  ) & ! [prod] c->r
              + wk(I_Psmlt  ) & ! [prod] s->r
              + wk(I_Pgmlt  ) & ! [prod] g->r
              - wk(I_Prevp  ) & ! [loss] r->v
              - wk(I_Piacr_s) & ! [loss] r->s
              - wk(I_Psacr_s) & ! [loss] r->s
              - wk(I_Piacr_g) & ! [loss] r->g
              - wk(I_Psacr_g) & ! [loss] r->g
              - wk(I_Pgacr  ) & ! [loss] r->g
              - wk(I_Pgfrz  )   ! [loss] r->g

          fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
          fac    = (          fac_sw ) &
                 + ( 1.0_RP - fac_sw ) * min( -wk(I_dqr_dt)/(net-fac_sw), 1.0_RP ) ! loss limiter

          wk(I_Praut  ) = wk(I_Praut  ) * fac
          wk(I_Pracw  ) = wk(I_Pracw  ) * fac
          wk(I_Psmlt  ) = wk(I_Psmlt  ) * fac
          wk(I_Pgmlt  ) = wk(I_Pgmlt  ) * fac
          wk(I_Prevp  ) = wk(I_Prevp  ) * fac
          wk(I_Piacr_s) = wk(I_Piacr_s) * fac
          wk(I_Psacr_s) = wk(I_Psacr_s) * fac
          wk(I_Piacr_g) = wk(I_Piacr_g) * fac
          wk(I_Psacr_g) = wk(I_Psacr_g) * fac
          wk(I_Pgacr  ) = wk(I_Pgacr  ) * fac
          wk(I_Pgfrz  ) = wk(I_Pgfrz  ) * fac

          ! [QV]
          net = &
              + wk(I_Prevp  ) & ! [prod] r->v
              + wk(I_Pisub  ) & ! [prod] i->v
              + wk(I_Pssub  ) & ! [prod] s->v
              + wk(I_Pgsub  ) & ! [prod] g->v
              - wk(I_Pigen  ) & ! [loss] v->i
              - wk(I_Pidep  ) & ! [loss] v->i
              - wk(I_Psdep  ) & ! [loss] v->s
              - wk(I_Pgdep  )   ! [loss] v->g

          fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
          fac    = (          fac_sw ) &
                 + ( 1.0_RP - fac_sw ) * min( -wk(I_dqv_dt)/(net-fac_sw), 1.0_RP ) ! loss limiter

          wk(I_Prevp  ) = wk(I_Prevp  ) * fac
          wk(I_Pisub  ) = wk(I_Pisub  ) * fac
          wk(I_Pssub  ) = wk(I_Pssub  ) * fac
          wk(I_Pgsub  ) = wk(I_Pgsub  ) * fac
          wk(I_Pigen  ) = wk(I_Pigen  ) * fac
          wk(I_Pidep  ) = wk(I_Pidep  ) * fac
          wk(I_Psdep  ) = wk(I_Psdep  ) * fac
          wk(I_Pgdep  ) = wk(I_Pgdep  ) * fac

          ! [QS]
          net = &
              + wk(I_Psdep  ) & ! [prod] v->s
              + wk(I_Psacw  ) & ! [prod] c->s
              + wk(I_Psfw   ) & ! [prod] c->s
              + wk(I_Piacr_s) & ! [prod] r->s
              + wk(I_Psacr_s) & ! [prod] r->s
              + wk(I_Psaut  ) & ! [prod] i->s
              + wk(I_Praci_s) & ! [prod] i->s
              + wk(I_Psaci  ) & ! [prod] i->s
              + wk(I_Psfi   ) & ! [prod] i->s
              - wk(I_Pssub  ) & ! [loss] s->v
              - wk(I_Psmlt  ) & ! [loss] s->r
              - wk(I_Pgaut  ) & ! [loss] s->g
              - wk(I_Pracs  ) & ! [loss] s->g
              - wk(I_Pgacs  )   ! [loss] s->g

          fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
          fac    = (          fac_sw ) &
                 + ( 1.0_RP - fac_sw ) * min( -wk(I_dqs_dt)/(net-fac_sw), 1.0_RP ) ! loss limiter

          wk(I_Psdep  ) = wk(I_Psdep  ) * fac
          wk(I_Psacw  ) = wk(I_Psacw  ) * fac
          wk(I_Psfw   ) = wk(I_Psfw   ) * fac
          wk(I_Piacr_s) = wk(I_Piacr_s) * fac
          wk(I_Psacr_s) = wk(I_Psacr_s) * fac
          wk(I_Psaut  ) = wk(I_Psaut  ) * fac
          wk(I_Praci_s) = wk(I_Praci_s) * fac
          wk(I_Psaci  ) = wk(I_Psaci  ) * fac
          wk(I_Psfi   ) = wk(I_Psfi   ) * fac
          wk(I_Pssub  ) = wk(I_Pssub  ) * fac
          wk(I_Psmlt  ) = wk(I_Psmlt  ) * fac
          wk(I_Pgaut  ) = wk(I_Pgaut  ) * fac
          wk(I_Pracs  ) = wk(I_Pracs  ) * fac
          wk(I_Pgacs  ) = wk(I_Pgacs  ) * fac

          ! [QG]
          net = &
              + wk(I_Pgdep  ) & ! [prod] v->g
              + wk(I_Pgacw  ) & ! [prod] c->g
              + wk(I_Piacr_g) & ! [prod] r->g
              + wk(I_Psacr_g) & ! [prod] r->g
              + wk(I_Pgacr  ) & ! [prod] r->g
              + wk(I_Pgfrz  ) & ! [prod] r->g
              + wk(I_Praci_g) & ! [prod] i->g
              + wk(I_Pgaci  ) & ! [prod] i->g
              + wk(I_Pgaut  ) & ! [prod] s->g
              + wk(I_Pracs  ) & ! [prod] s->g
              + wk(I_Pgacs  ) & ! [prod] s->g
              - wk(I_Pgsub  ) & ! [loss] g->v
              - wk(I_Pgmlt  )   ! [loss] g->r

          fac_sw = 0.5_RP + sign( 0.5_RP, net+EPS ) ! if production > loss , fac_sw=1
          fac    = (          fac_sw ) &
                 + ( 1.0_RP - fac_sw ) * min( -wk(I_dqg_dt)/(net-fac_sw), 1.0_RP ) ! loss limiter

          wk(I_Pgdep  ) = wk(I_Pgdep  ) * fac
          wk(I_Pgacw  ) = wk(I_Pgacw  ) * fac
          wk(I_Piacr_g) = wk(I_Piacr_g) * fac
          wk(I_Psacr_g) = wk(I_Psacr_g) * fac
          wk(I_Pgacr  ) = wk(I_Pgacr  ) * fac
          wk(I_Pgfrz  ) = wk(I_Pgfrz  ) * fac
          wk(I_Praci_g) = wk(I_Praci_g) * fac
          wk(I_Pgaci  ) = wk(I_Pgaci  ) * fac
          wk(I_Pgaut  ) = wk(I_Pgaut  ) * fac
          wk(I_Pracs  ) = wk(I_Pracs  ) * fac
          wk(I_Pgacs  ) = wk(I_Pgacs  ) * fac
          wk(I_Pgsub  ) = wk(I_Pgsub  ) * fac
          wk(I_Pgmlt  ) = wk(I_Pgmlt  ) * fac

          qc_t = + wk(I_Pimlt  ) & ! [prod] i->c
                 - wk(I_Praut  ) & ! [loss] c->r
                 - wk(I_Pracw  ) & ! [loss] c->r
                 - wk(I_Pihom  ) & ! [loss] c->i
                 - wk(I_Pihtr  ) & ! [loss] c->i
                 - wk(I_Psacw  ) & ! [loss] c->s
                 - wk(I_Psfw   ) & ! [loss] c->s
                 - wk(I_Pgacw  )   ! [loss] c->g

          qr_t = + wk(I_Praut  ) & ! [prod] c->r
                 + wk(I_Pracw  ) & ! [prod] c->r
                 + wk(I_Psmlt  ) & ! [prod] s->r
                 + wk(I_Pgmlt  ) & ! [prod] g->r
                 - wk(I_Prevp  ) & ! [loss] r->v
                 - wk(I_Piacr_s) & ! [loss] r->s
                 - wk(I_Psacr_s) & ! [loss] r->s
                 - wk(I_Piacr_g) & ! [loss] r->g
                 - wk(I_Psacr_g) & ! [loss] r->g
                 - wk(I_Pgacr  ) & ! [loss] r->g
                 - wk(I_Pgfrz  )   ! [loss] r->g

          qi_t = + wk(I_Pigen  ) & ! [prod] v->i
                 + wk(I_Pidep  ) & ! [prod] v->i
                 + wk(I_Pihom  ) & ! [prod] c->i
                 + wk(I_Pihtr  ) & ! [prod] c->i
                 - wk(I_Pisub  ) & ! [loss] i->v
                 - wk(I_Pimlt  ) & ! [loss] i->c
                 - wk(I_Psaut  ) & ! [loss] i->s
                 - wk(I_Praci_s) & ! [loss] i->s
                 - wk(I_Psaci  ) & ! [loss] i->s
                 - wk(I_Psfi   ) & ! [loss] i->s
                 - wk(I_Praci_g) & ! [loss] i->g
                 - wk(I_Pgaci  )   ! [loss] i->g

          qs_t = + wk(I_Psdep  ) & ! [prod] v->s
                 + wk(I_Psacw  ) & ! [prod] c->s
                 + wk(I_Psfw   ) & ! [prod] c->s
                 + wk(I_Piacr_s) & ! [prod] r->s
                 + wk(I_Psacr_s) & ! [prod] r->s
                 + wk(I_Psaut  ) & ! [prod] i->s
                 + wk(I_Praci_s) & ! [prod] i->s
                 + wk(I_Psaci  ) & ! [prod] i->s
                 + wk(I_Psfi   ) & ! [prod] i->s
                 - wk(I_Pssub  ) & ! [loss] s->v
                 - wk(I_Psmlt  ) & ! [loss] s->r
                 - wk(I_Pgaut  ) & ! [loss] s->g
                 - wk(I_Pracs  ) & ! [loss] s->g
                 - wk(I_Pgacs  )   ! [loss] s->g

          qg_t = + wk(I_Pgdep  ) & ! [prod] v->g
                 + wk(I_Pgacw  ) & ! [prod] c->g
                 + wk(I_Piacr_g) & ! [prod] r->g
                 + wk(I_Psacr_g) & ! [prod] r->g
                 + wk(I_Pgacr  ) & ! [prod] r->g
                 + wk(I_Pgfrz  ) & ! [prod] r->g
                 + wk(I_Praci_g) & ! [prod] i->g
                 + wk(I_Pgaci  ) & ! [prod] i->g
                 + wk(I_Pgaut  ) & ! [prod] s->g
                 + wk(I_Pracs  ) & ! [prod] s->g
                 + wk(I_Pgacs  ) & ! [prod] s->g
                 - wk(I_Pgsub  ) & ! [loss] g->v
                 - wk(I_Pgmlt  )   ! [loss] g->r

          qc_t = max( qc_t, -wk(I_dqc_dt) )
          qr_t = max( qr_t, -wk(I_dqr_dt) )
          qi_t = max( qi_t, -wk(I_dqi_dt) )
          qs_t = max( qs_t, -wk(I_dqs_dt) )
          qg_t = max( qg_t, -wk(I_dqg_dt) )

          qv_t = - ( qc_t &
                   + qr_t &
                   + qi_t &
                   + qs_t &
                   + qg_t )

          drhogqv(ij,k) = rhog(ij,k) * qv_t * dt
          drhogqc(ij,k) = rhog(ij,k) * qc_t * dt
          drhogqr(ij,k) = rhog(ij,k) * qr_t * dt
          drhogqi(ij,k) = rhog(ij,k) * qi_t * dt
          drhogqs(ij,k) = rhog(ij,k) * qs_t * dt
          drhogqg(ij,k) = rhog(ij,k) * qg_t * dt

          Vt(ij,k,I_QR) = Vtr
          Vt(ij,k,I_QI) = Vti
          Vt(ij,k,I_QS) = Vts
          Vt(ij,k,I_QG) = Vtg

          !--- Individual tendency for each process
          !        do ip = 1, wk_nmax
          !           ml_wk(ij,k,ip) = wk(ip)
          !        enddo

          !--- Effective Radius of Liquid Water
          xf_qc = rhoqc / Nc(ij,k) * 1.E-9_RP                    ! mean mass   of qc [kg]
          rf_qc = ( coef_xf * xf_qc + EPS )**0.33333333_RP       ! mean radius of qc [m]
          r2_qc = coef_dgam * rf_qc**2 * ( Nc(ij,k) * 1.E+6_RP ) ! r^2 moment  of qc
          r2_qr = 0.25_RP * N0r * GAM_3 * RLMDr_3                ! r^2 moment  of qr
          r3_qc = coef_xf * dens * qc                            ! r^3 moment  of qc
          r3_qr = coef_xf * dens * qr                            ! r^3 moment  of qr

          zerosw = 0.5_RP - sign(0.5_RP, (r2_qc+r2_qr)-r2_min )
          rceff(ij,k) = ( r3_qc + r3_qr ) / ( r2_qc + r2_qr + zerosw ) * ( 1.0_RP - zerosw )

          sw = ( 0.5_RP + sign(0.5_RP, (qc+qr)-q_min ) ) * zerosw ! if qc+qr > 0.01[g/kg], sw=1

          rctop(ij,1) = (        sw ) * rceff(ij,k) &
                      + ( 1.0_RP-sw ) * UNDEF
          tctop(ij,1) = (        sw ) * temp &
                      + ( 1.0_RP-sw ) * UNDEF

          zerosw = 0.5_RP - sign(0.5_RP, r2_qc-r2_min )
          rceff_cld(ij,k) = r3_qc / ( r2_qc + zerosw ) * ( 1.0_RP - zerosw )

          sw = ( 0.5_RP + sign(0.5_RP, qc-q_min ) ) * zerosw ! if qc > 0.01[g/kg], sw=1

          rctop_cld(ij,1) = (        sw ) * rceff_cld(ij,k) &
                          + ( 1.0_RP-sw ) * UNDEF
          tctop_cld(ij,1) = (        sw ) * temp &
                          + ( 1.0_RP-sw ) * UNDEF

#if defined _OPENMP && _OPENMP >= 201511
       enddo
#endif
    enddo
    enddo
    !$omp end parallel do

if ( TIME_CSTEP == EX_STEP ) then
   call PROF_rapend  ('_kernel_partM12')
endif

    !--- update rhogq
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kmin,kmax,rhogq,drhogqv,drhogqc,drhogqr,drhogqi,drhogqs,drhogqg, &
    !$omp        I_QV,I_QC,I_QR,I_QI,I_QS,I_QG ) &
    !$omp collapse(2)
    do k  = kmin, kmax
    do ij = 1, ijdim
       rhogq(ij,k,I_QV) = rhogq(ij,k,I_QV) + drhogqv(ij,k)
       rhogq(ij,k,I_QC) = rhogq(ij,k,I_QC) + drhogqc(ij,k)
       rhogq(ij,k,I_QR) = rhogq(ij,k,I_QR) + drhogqr(ij,k)
       rhogq(ij,k,I_QI) = rhogq(ij,k,I_QI) + drhogqi(ij,k)
       rhogq(ij,k,I_QS) = rhogq(ij,k,I_QS) + drhogqs(ij,k)
       rhogq(ij,k,I_QG) = rhogq(ij,k,I_QG) + drhogqg(ij,k)
    enddo
    enddo
    !$omp end parallel do

    !--- update rhoge
    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kmin,kmax,rhoge,drhogqv,drhogqi,drhogqs,drhogqg,LHV,LHF) &
    !$omp collapse(2)
    do k  = kmin, kmax
    do ij = 1, ijdim
       rhoge(ij,k) = rhoge(ij,k) - LHV * drhogqv(ij,k) &
                                 + LHF * drhogqi(ij,k) &
                                 + LHF * drhogqs(ij,k) &
                                 + LHF * drhogqg(ij,k)
    enddo
    enddo
    !$omp end parallel do

    !$omp parallel do default(none),private(ij,sw), &
    !$omp shared(ijdim,kmin,kmax,rceff,tctop,rctop,rwtop,rceff_cld,tctop_cld,rctop_cld,rwtop_cld,UNDEF)
    do ij = 1, ijdim
       rceff    (ij,kmin-1) = 0.0_RP
       rceff_cld(ij,kmin-1) = 0.0_RP
       rceff    (ij,kmax+1) = 0.0_RP
       rceff_cld(ij,kmax+1) = 0.0_RP

       sw = ( 0.5_RP + sign(0.5_RP, tctop(ij,1)-TEM00 ) ) ! if T > 0C, sw=1

       rwtop(ij,1) = (        sw ) * rctop(ij,1) &
                   + ( 1.0_RP-sw ) * UNDEF

       sw = ( 0.5_RP + sign(0.5_RP, rctop_cld(ij,1)-TEM00 ) ) ! if T > 0C, sw=1

       rwtop_cld(ij,1) = (        sw ) * rctop_cld(ij,1) &
                       + ( 1.0_RP-sw ) * UNDEF
    enddo
    !$omp end parallel do

if ( TIME_CSTEP == EX_STEP ) then
   call PROF_rapstart('_kernel_partM3')
endif

    !---< preciptation >---

    !$omp parallel do default(none),private(ij,nq), &
    !$omp shared(nqmax,ijdim,kmin,kmax,Vt) &
    !$omp collapse(2)
    do nq = 1, nqmax
    do ij = 1, ijdim
       Vt(ij,kmin-1,nq) = 0.0_RP
       Vt(ij,kmax+1,nq) = 0.0_RP
    enddo
    enddo
    !$omp end parallel do

    !--- update mass concentration
    !$omp parallel do default(none),private(ij,k,nq), &
    !$omp shared(ijdim,kmin,kmax,NQW_STR,NQW_END,q,rhogq,rhog) &
    !$omp collapse(3)
    do nq = NQW_STR, NQW_END
    do k  = kmin, kmax
    do ij = 1, ijdim
       q(ij,k,nq) = rhogq(ij,k,nq) / rhog(ij,k)
    enddo
    enddo
    enddo
    !$omp end parallel do

    call thrmdyn_qd( ijdim,     & ! [IN]
                     kdim,      & ! [IN]
                     q (:,:,:), & ! [IN]
                     qd(:,:)    ) ! [OUT]

    call thrmdyn_cv( ijdim,      & ! [IN]
                     kdim,       & ! [IN]
                     qd (:,:),   & ! [IN]
                     q  (:,:,:), & ! [IN]
                     cva(:,:)    ) ! [OUT]

    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,tem,rhoge,rhog,cva,rgs,rgsh,gam2,gam2h,gsgam2,gsgam2h), &
    !$omp collapse(2)
    do k  = 1, kdim
    do ij = 1, ijdim
       tem (ij,k) = rhoge(ij,k) / ( rhog(ij,k) * cva(ij,k) )
       rgs (ij,k) = gam2 (ij,k) / gsgam2 (ij,k)
       rgsh(ij,k) = gam2h(ij,k) / gsgam2h(ij,k)
    enddo
    enddo
    !$omp end parallel do

    preciptation_flag(:)    = .false.
    preciptation_flag(I_QR) = .true.
    preciptation_flag(I_QI) = .true.
    preciptation_flag(I_QS) = .true.
    preciptation_flag(I_QG) = .true.

    if ( precip_transport_type == '3WATER' ) then
       preciptation_flag(I_QI) = .false.
    endif

!ESC!    if ( precip_scheme_type == 'Upwind-Euler' ) then
!ESC!
!ESC!       call precip_transport_euler( ijdim,             &
!ESC!                                    rhog,              &
!ESC!                                    rhogvx,            &
!ESC!                                    rhogvy,            &
!ESC!                                    rhogvz,            &
!ESC!                                    rhogw,             &
!ESC!                                    rhoge,             &
!ESC!                                    rhogq,             &
!ESC!                                    rho,               &
!ESC!                                    tem,               &
!ESC!                                    pre,               &
!ESC!                                    vx,                &
!ESC!                                    vy,                &
!ESC!                                    vz,                &
!ESC!                                    w,                 &
!ESC!                                    q,                 &
!ESC!                                    qd,                &
!ESC!                                    z,                 &
!ESC!                                    Vt,                &
!ESC!                                    preciptation_flag, &
!ESC!                                    precip,            &
!ESC!                                    precip_rhoe,       &
!ESC!                                    precip_lh_heat,    &
!ESC!                                    precip_rhophi,     &
!ESC!                                    precip_rhokin,     &
!ESC!                                    gprec,             &
!ESC!                                    gsgam2,            &
!ESC!                                    gsgam2h,           &
!ESC!                                    rgs,               &
!ESC!                                    rgsh,              &
!ESC!                                    ix,                &
!ESC!                                    iy,                &
!ESC!                                    iz,                &
!ESC!                                    jx,                &
!ESC!                                    jy,                &
!ESC!                                    jz,                &
!ESC!                                    dt                 )
!ESC!
!ESC!    elseif( precip_scheme_type == 'Flux-Semilag_new' ) then

       call precip_transport_new( ijdim,             &
                                  rhog,              &
                                  rhogvx,            &
                                  rhogvy,            &
                                  rhogvz,            &
                                  rhogw,             &
                                  rhoge,             &
                                  rhogq,             &
                                  rho,               &
                                  tem,               &
                                  pre,               &
                                  vx,                &
                                  vy,                &
                                  vz,                &
                                  w,                 &
                                  q,                 &
                                  qd,                &
                                  z,                 &
                                  Vt,                &
                                  preciptation_flag, &
                                  precip,            &
                                  precip_rhoe,       &
                                  precip_lh_heat,    &
                                  precip_rhophi,     &
                                  precip_rhokin,     &
                                  gprec,             &
                                  gsgam2,            &
                                  gsgam2h,           &
                                  rgs,               &
                                  rgsh,              &
                                  ix,                &
                                  iy,                &
                                  iz,                &
                                  jx,                &
                                  jy,                &
                                  jz,                &
                                  dt                 )

!ESC!    else
!ESC!
!ESC!       call precip_transport_nwater( ijdim,             &
!ESC!                                     rhog,              &
!ESC!                                     rhogvx,            &
!ESC!                                     rhogvy,            &
!ESC!                                     rhogvz,            &
!ESC!                                     rhogw,             &
!ESC!                                     rhoge,             &
!ESC!                                     rhogq,             &
!ESC!                                     rho,               &
!ESC!                                     tem,               &
!ESC!                                     pre,               &
!ESC!                                     vx,                &
!ESC!                                     vy,                &
!ESC!                                     vz,                &
!ESC!                                     w,                 &
!ESC!                                     q,                 &
!ESC!                                     qd,                &
!ESC!                                     z,                 &
!ESC!                                     Vt,                &
!ESC!                                     preciptation_flag, &
!ESC!                                     precip,            &
!ESC!                                     precip_rhoe,       &
!ESC!                                     precip_lh_heat,    &
!ESC!                                     precip_rhophi,     &
!ESC!                                     precip_rhokin,     &
!ESC!                                     gprec,             &
!ESC!                                     gsgam2,            &
!ESC!                                     gsgam2h,           &
!ESC!                                     rgs,               &
!ESC!                                     rgsh,              &
!ESC!                                     ix,                &
!ESC!                                     iy,                &
!ESC!                                     iz,                &
!ESC!                                     jx,                &
!ESC!                                     jy,                &
!ESC!                                     jz,                &
!ESC!                                     dt                 )
!ESC!
!ESC!    endif

if ( TIME_CSTEP == EX_STEP ) then
   call PROF_rapend  ('_kernel_partM3')
endif

    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,ml_Pconv,ml_Pconw,ml_Pconi,q,I_QV,I_QC,I_QI) &
    !$omp collapse(2)
    do k  = 1, kdim
    do ij = 1, ijdim
       ml_Pconv(ij,k) = q(ij,k,I_QV)
       ml_Pconw(ij,k) = q(ij,k,I_QC)
       ml_Pconi(ij,k) = q(ij,k,I_QI)
    enddo
    enddo
    !$omp end parallel do

if ( TIME_CSTEP == EX_STEP ) then
   call PROF_rapstart('_kernel_partM4')
endif

    if ( OPT_EXPLICIT_ICEGEN ) then

       call SATURATION_setrange( 100.0_RP, 90.0_RP )

       call SATURATION_adjustment( ijdim,             & ! [IN]
                                   kdim,              & ! [IN]
                                   rhog  (:,:),       & ! [IN]
                                   rhoge (:,:),       & ! [INOUT]
                                   rhogq (:,:,:),     & ! [INOUT]
                                   tem   (:,:),       & ! [INOUT]
                                   q     (:,:,:),     & ! [INOUT]
                                   qd    (:,:),       & ! [IN]
                                   gsgam2(:,:),       & ! [IN]
                                   ice_adjust=.false. ) ! [IN]

    else

       call SATURATION_setrange( 273.16_RP, 233.16_RP )

       call SATURATION_adjustment( ijdim,             & ! [IN]
                                   kdim,              & ! [IN]
                                   rhog  (:,:),       & ! [IN]
                                   rhoge (:,:),       & ! [INOUT]
                                   rhogq (:,:,:),     & ! [INOUT]
                                   tem   (:,:),       & ! [INOUT]
                                   q     (:,:,:),     & ! [INOUT]
                                   qd    (:,:),       & ! [IN]
                                   gsgam2(:,:),       & ! [IN]
                                   ice_adjust=.true.  ) ! [IN]

    endif

    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,ml_Pconv,ml_Pconw,ml_Pconi,q,I_QV,I_QC,I_QI) &
    !$omp collapse(2)
    do k  = 1, kdim
    do ij = 1, ijdim
       ml_Pconv(ij,k) = q(ij,k,I_QV) - ml_Pconv(ij,k)
       ml_Pconw(ij,k) = q(ij,k,I_QC) - ml_Pconw(ij,k)
       ml_Pconi(ij,k) = q(ij,k,I_QI) - ml_Pconi(ij,k)
    enddo
    enddo
    !$omp end parallel do

if ( TIME_CSTEP == EX_STEP ) then
   call PROF_rapend  ('_kernel_partM4')
endif

    ! Individual tendency output
!     do ip = 1, wk_nmax
!        call history_in( 'ml_'//trim(w_name(ip)), ml_wk(:,:,ip) )
!     enddo
!ESC!    call history_in( 'ml_Pconv', ml_Pconv(:,:) )
!ESC!    call history_in( 'ml_Pconw', ml_Pconw(:,:) )
!ESC!    call history_in( 'ml_Pconi', ml_Pconi(:,:) )

    call negative_filter( ijdim,         & ! [IN]
                          rhog  (:,:),   & ! [INOUT]
                          rhoge (:,:),   & ! [INOUT]
                          rhogq (:,:,:), & ! [INOUT]
                          rho   (:,:),   & ! [INOUT]
                          tem   (:,:),   & ! [IN]
                          pre   (:,:),   & ! [INOUT]
                          q     (:,:,:), & ! [INOUT]
                          gsgam2(:,:)    ) ! [IN]

    call PROF_rapend  ('____MP_NSW6')

    return
  end subroutine mp_nsw6

  !-----------------------------------------------------------------------
  subroutine negative_filter( &
       ijdim, &
       rhog,  &
       rhoge, &
       rhogq, &
       rho,   &
       tem,   &
       pre,   &
       q,     &
       gsgam2 )
!ESC!    use mod_adm, only: &
!ESC!       kdim => ADM_kall, &
!ESC!       kmin => ADM_kmin, &
!ESC!       kmax => ADM_kmax
!ESC!    use mod_const, only: &
!ESC!       CONST_Rdry, &
!ESC!       CONST_Rvap
!ESC!    use mod_runconf, only: &
!ESC!       nqmax => TRC_VMAX, &
!ESC!       NQW_STR,           &
!ESC!       NQW_END,           &
!ESC!       I_QV
    use mod_thrmdyn, only: &
       thrmdyn_cv, &
       thrmdyn_qd
    implicit none

    integer,  intent(in)    :: ijdim
    real(RP), intent(inout) :: rhog  (ijdim,kdim)
    real(RP), intent(inout) :: rhoge (ijdim,kdim)
    real(RP), intent(inout) :: rhogq (ijdim,kdim,nqmax)
    real(RP), intent(inout) :: rho   (ijdim,kdim)
    real(RP), intent(in)    :: tem   (ijdim,kdim)
    real(RP), intent(inout) :: pre   (ijdim,kdim)
    real(RP), intent(inout) :: q     (ijdim,kdim,nqmax)
    real(RP), intent(in)    :: gsgam2(ijdim,kdim)

    real(RP) :: qd (ijdim,kdim)
    real(RP) :: cva(ijdim,kdim)
    real(RP) :: diffq
    real(RP) :: Rdry, Rvap

    integer  :: ij, k, nq
    !---------------------------------------------------------------------------

    Rdry = CONST_Rdry
    Rvap = CONST_Rvap

    !$omp parallel default(none),private(ij,k,nq,diffq), &
    !$omp shared(ijdim,kmin,kmax,NQW_STR,NQW_END,rhogq,rhog,q,rho,gsgam2,I_QV)

    !$omp do
    do k  = kmin, kmax
    do ij = 1, ijdim
       diffq = 0.0_RP
       do nq = NQW_STR+1, NQW_END
          ! total hydrometeor (before correction)
          diffq = diffq + rhogq(ij,k,nq)
          ! remove negative value of hydrometeors (mass)
          rhogq(ij,k,nq) = max( rhogq(ij,k,nq), 0.0_RP )
       enddo

       do nq = NQW_STR+1, NQW_END
          ! difference between before and after correction
          diffq = diffq - rhogq(ij,k,nq)
       enddo

       ! Compensate for the lack of hydrometeors by the water vapor
       rhogq(ij,k,I_QV) = rhogq(ij,k,I_QV) + diffq
    enddo
    enddo
    !$omp end do

    !$omp do
    do k  = kmin, kmax
    do ij = 1, ijdim
       diffq = rhogq(ij,k,I_QV)
       ! remove negative value of water vapor (mass)
       rhogq(ij,k,I_QV) = max( rhogq(ij,k,I_QV), 0.0_RP )

       diffq = diffq - rhogq(ij,k,I_QV)

       ! Apply correction to total density
       rhog(ij,k) = rhog(ij,k) * ( 1.0_RP - diffq ) ! diffq is negative
       rho (ij,k) = rhog(ij,k) / gsgam2(ij,k)
    enddo
    enddo
    !$omp end do

    !$omp do
    do nq = NQW_STR, NQW_END
    do k  = kmin, kmax
    do ij = 1, ijdim
       !--- update mass concentration
       q(ij,k,nq) = rhogq(ij,k,nq) / rhog(ij,k)
    enddo
    enddo
    enddo
    !$omp end do

    !$omp end parallel

    call thrmdyn_qd( ijdim,     & ! [IN]
                     kdim,      & ! [IN]
                     q (:,:,:), & ! [IN]
                     qd(:,:)    ) ! [OUT]

    call thrmdyn_cv( ijdim,      & ! [IN]
                     kdim,       & ! [IN]
                     qd (:,:),   & ! [IN]
                     q  (:,:,:), & ! [IN]
                     cva(:,:)    ) ! [OUT]

    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kmin,kmax,rhoge,pre,rho,tem,rhog,qd,q,cva,I_QV,Rdry,Rvap)
    do k  = kmin, kmax
    do ij = 1, ijdim
       rhoge(ij,k) = tem(ij,k) * rhog(ij,k) * cva(ij,k)
       pre  (ij,k) = rho(ij,k) * ( qd(ij,k) * Rdry + q(ij,k,I_QV) * Rvap ) * tem(ij,k)
    enddo
    enddo
    !$omp end parallel do

    return
  end subroutine negative_filter

  !-----------------------------------------------------------------------
  subroutine Bergeron_param( &
       ijdim, &
       tem,   &
       a1,    &
       a2,    &
       ma2    )
!ESC!    use mod_adm, only: &
!ESC!       kdim => ADM_kall, &
!ESC!       kmin => ADM_kmin, &
!ESC!       kmax => ADM_kmax
!ESC!    use mod_const, only: &
!ESC!       TEM00 => CONST_TEM00
    implicit none

    integer,  intent(in)  :: ijdim
    real(RP), intent(in)  :: tem(ijdim,kdim)
    real(RP), intent(out) :: a1 (ijdim,kdim)
    real(RP), intent(out) :: a2 (ijdim,kdim)
    real(RP), intent(out) :: ma2(ijdim,kdim)

    real(RP) :: a1_tab(32)
    real(RP) :: a2_tab(32)

    data a1_tab / 0.0001E-7_RP, 0.7939E-7_RP, 0.7841E-6_RP, 0.3369E-5_RP, 0.4336E-5_RP, &
                  0.5285E-5_RP, 0.3728E-5_RP, 0.1852E-5_RP, 0.2991E-6_RP, 0.4248E-6_RP, &
                  0.7434E-6_RP, 0.1812E-5_RP, 0.4394E-5_RP, 0.9145E-5_RP, 0.1725E-4_RP, &
                  0.3348E-4_RP, 0.1725E-4_RP, 0.9175E-5_RP, 0.4412E-5_RP, 0.2252E-5_RP, &
                  0.9115E-6_RP, 0.4876E-6_RP, 0.3473E-6_RP, 0.4758E-6_RP, 0.6306E-6_RP, &
                  0.8573E-6_RP, 0.7868E-6_RP, 0.7192E-6_RP, 0.6513E-6_RP, 0.5956E-6_RP, &
                  0.5333E-6_RP, 0.4834E-6_RP  /

    data a2_tab / 0.0100_RP, 0.4006_RP, 0.4831_RP, 0.5320_RP, 0.5307_RP, &
                  0.5319_RP, 0.5249_RP, 0.4888_RP, 0.3849_RP, 0.4047_RP, &
                  0.4318_RP, 0.4771_RP, 0.5183_RP, 0.5463_RP, 0.5651_RP, &
                  0.5813_RP, 0.5655_RP, 0.5478_RP, 0.5203_RP, 0.4906_RP, &
                  0.4447_RP, 0.4126_RP, 0.3960_RP, 0.4149_RP, 0.4320_RP, &
                  0.4506_RP, 0.4483_RP, 0.4460_RP, 0.4433_RP, 0.4413_RP, &
                  0.4382_RP, 0.4361_RP  /

    real(RP) :: temc
    integer  :: itemc
    real(RP) :: fact

    integer  :: ij, k
    !---------------------------------------------------------------------------

    !$omp parallel do default(none),private(ij,k,temc,itemc,fact), &
    !$omp shared(ijdim,kmin,kmax,a1,a2,ma2,tem,a1_tab,a2_tab)
    do k  = kmin, kmax
    do ij = 1, ijdim
       temc  = min( max( tem(ij,k)-TEM00, -30.99_RP ), 0.0_RP ) ! 0C <= T  <  31C
       itemc = int( -temc ) + 1                                 ! 1  <= iT <= 31
       fact  = - ( temc + real(itemc-1,kind=8) )

       a1 (ij,k) = ( 1.0_RP-fact ) * a1_tab(itemc  ) &
                 + (        fact ) * a1_tab(itemc+1)

       a2 (ij,k) = ( 1.0_RP-fact ) * a2_tab(itemc  ) &
                 + (        fact ) * a2_tab(itemc+1)

       ma2(ij,k) = 1.0_RP - a2(ij,k)

       a1 (ij,k) = a1(ij,k) * 1.E-3_RP**ma2(ij,k) ! [g->kg]
    enddo
    enddo
    !$omp end parallel do

    return
  end subroutine Bergeron_param

!ESC!  !-----------------------------------------------------------------------------
!ESC!  subroutine mp_nsw6_terminal_velocity( &
!ESC!       ijdim, &
!ESC!       kdim,  &
!ESC!       rho,   &
!ESC!       tem,   &
!ESC!       w,     &
!ESC!       qc,    &
!ESC!       qr,    &
!ESC!       qi,    &
!ESC!       qs,    &
!ESC!       qg,    &
!ESC!       vt_qc, &
!ESC!       vt_qr, &
!ESC!       vt_qi, &
!ESC!       vt_qs, &
!ESC!       vt_qg  )
!ESC!    use mod_adm, only: &
!ESC!       kmin => adm_kmin, &
!ESC!       kmax => adm_kmax
!ESC!    use mod_const, only: &
!ESC!       CONST_EPS,            &
!ESC!       TEM00 => CONST_TEM00
!ESC!    implicit none
!ESC!
!ESC!    integer,  intent(in)  :: ijdim
!ESC!    integer,  intent(in)  :: kdim
!ESC!    real(RP), intent(in)  :: rho  (ijdim,kdim)   !< density
!ESC!    real(RP), intent(in)  :: tem  (ijdim,kdim)   !< temperature
!ESC!    real(RP), intent(in)  :: w    (ijdim,kdim)   !< vertical velocity [m/s]    ![Add] 2016/06/04 W.Roh
!ESC!    real(RP), intent(in)  :: qc   (ijdim,kdim)   !< mixing ratio of cloud
!ESC!    real(RP), intent(in)  :: qr   (ijdim,kdim)   !< mixing ratio of rain
!ESC!    real(RP), intent(in)  :: qi   (ijdim,kdim)   !< mixing ratio of ice
!ESC!    real(RP), intent(in)  :: qs   (ijdim,kdim)   !< mixing ratio of snow
!ESC!    real(RP), intent(in)  :: qg   (ijdim,kdim)   !< mixing ratio of graupel
!ESC!    real(RP), intent(out) :: vt_qc(ijdim,kdim)   !< terminal velocity of cloud
!ESC!    real(RP), intent(out) :: vt_qr(ijdim,kdim)   !< terminal velocity of rain
!ESC!    real(RP), intent(out) :: vt_qi(ijdim,kdim)   !< terminal velocity of ice
!ESC!    real(RP), intent(out) :: vt_qs(ijdim,kdim)   !< terminal velocity of snow
!ESC!    real(RP), intent(out) :: vt_qg(ijdim,kdim)   !< terminal velocity of graupel
!ESC!
!ESC!    real(RP) :: dens                             !< density
!ESC!    real(RP) :: rho_fact                         !< density factor
!ESC!
!ESC!    real(RP) :: RLMDr, RLMDr_dr
!ESC!    real(RP) :: RLMDs, RLMDs_ds, RLMDs_1bs
!ESC!    real(RP) :: RLMDg, RLMDg_dg
!ESC!
!ESC!    !---< Roh and Satoh (2014) >---
!ESC!    real(RP) :: tems, Xs2
!ESC!    real(RP) :: MOMs_0bs, RMOMs_Vt
!ESC!    real(RP) :: coef_at(4), coef_bt(4)
!ESC!    real(RP) :: loga_, b_, nm
!ESC!
!ESC!    real(RP) :: zerosw
!ESC!    real(RP) :: EPS
!ESC!
!ESC!    integer  :: ij,k
!ESC!    !---------------------------------------------------------------------------
!ESC!
!ESC!    EPS = CONST_EPS
!ESC!
!ESC!    !$omp parallel do default(none), &
!ESC!    !$omp private(ij,k,dens,rho_fact,RLMDr,RLMDr_dr,RLMDs,RLMDs_ds,RLMDs_1bs,RLMDg,RLMDg_dg, &
!ESC!    !$omp         tems,Xs2,MOMs_0bs,RMOMs_Vt,coef_at,coef_bt,loga_,b_,nm,zerosw),      &
!ESC!    !$omp shared(ijdim,kmin,kmax,vt_qc,vt_qr,vt_qi,vt_qs,vt_qg,rho,tem,qr,qi,qs,qg,          &
!ESC!    !$omp        GAM_1br,GAM_1brdr,GAM_1bs,GAM_1bsds,GAM_1bg,GAM_1bgdg,coef_a,coef_b,        &
!ESC!    !$omp        N0r,N0s,N0g,Ar,As,Ag,Cr,Cs,Cg,sw_roh2014,EPS)
!ESC!    do k  = kmin, kmax
!ESC!    do ij = 1, ijdim
!ESC!
!ESC!       dens     = rho(ij,k)
!ESC!       rho_fact = sqrt( dens00 / dens )
!ESC!
!ESC!       ! slope parameter lambda (Rain)
!ESC!       zerosw = 0.5_RP - sign(0.5_RP, qr(ij,k) - 1.E-12_RP )
!ESC!       RLMDr  = sqrt(sqrt( dens * qr(ij,k) / ( Ar * N0r * GAM_1br ) + zerosw )) * ( 1.0_RP - zerosw )
!ESC!       RLMDr_dr  = sqrt( RLMDr )       ! **Dr
!ESC!
!ESC!       ! slope parameter lambda (Snow)
!ESC!       zerosw = 0.5_RP - sign(0.5_RP, qs(ij,k) - 1.E-12_RP )
!ESC!       RLMDs  = sqrt(sqrt( dens * qs(ij,k) / ( As * N0s * GAM_1bs ) + zerosw )) * ( 1.0_RP - zerosw )
!ESC!
!ESC!       RLMDs_ds  = sqrt( sqrt(RLMDs) ) ! **Ds
!ESC!       RLMDs_1bs = RLMDs**4 ! (1+Bs)
!ESC!
!ESC!       MOMs_0bs   = N0s * GAM_1bs   * RLMDs_1bs       ! Ns * 0+bs
!ESC!       RMOMs_Vt   = GAM_1bsds / GAM_1bs * RLMDs_ds
!ESC!
!ESC!       ! bimodal size distribution of snow
!ESC!       zerosw = 0.5_RP - sign(0.5_RP, dens * qs(ij,k) - 1.E-12_RP )
!ESC!       Xs2    = dens * qs(ij,k) / As * ( 1.0_RP - zerosw )
!ESC!
!ESC!       tems       = min( -0.1_RP, tem(ij,k)-TEM00 )
!ESC!       coef_at(1) = coef_a( 1) + tems * ( coef_a( 2) + tems * ( coef_a( 5) + tems * coef_a( 9) ) )
!ESC!       coef_at(2) = coef_a( 3) + tems * ( coef_a( 4) + tems *   coef_a( 7) )
!ESC!       coef_at(3) = coef_a( 6) + tems *   coef_a( 8)
!ESC!       coef_at(4) = coef_a(10)
!ESC!       coef_bt(1) = coef_b( 1) + tems * ( coef_b( 2) + tems * ( coef_b( 5) + tems * coef_b( 9) ) )
!ESC!       coef_bt(2) = coef_b( 3) + tems * ( coef_b( 4) + tems *   coef_b( 7) )
!ESC!       coef_bt(3) = coef_b( 6) + tems *   coef_b( 8)
!ESC!       coef_bt(4) = coef_b(10)
!ESC!       ! 0 + Bs(=2) moment
!ESC!       nm = 2.0_RP
!ESC!       loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
!ESC!          b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
!ESC!       MOMs_0bs = (        sw_roh2014 ) * 10.0_RP**loga_ * Xs2**b_ &
!ESC!                + ( 1.0_RP-sw_roh2014 ) * MOMs_0bs
!ESC!       ! Bs(=2) + Ds(=0.25) moment
!ESC!       nm = 2.25_RP
!ESC!       loga_  = coef_at(1) + nm * ( coef_at(2) + nm * ( coef_at(3) + nm * coef_at(4) ) )
!ESC!          b_  = coef_bt(1) + nm * ( coef_bt(2) + nm * ( coef_bt(3) + nm * coef_bt(4) ) )
!ESC!
!ESC!       zerosw = 0.5_RP - sign(0.5_RP, abs(MOMs_0bs) - EPS )
!ESC!       RMOMs_Vt = (        sw_roh2014 ) * 10.0_RP**loga_ * Xs2**b_ / ( MOMs_0bs + zerosw ) &
!ESC!                + ( 1.0_RP-sw_roh2014 ) * RMOMs_Vt
!ESC!
!ESC!       ! slope parameter lambda (Graupel)
!ESC!       zerosw = 0.5_RP - sign(0.5_RP, qg(ij,k) - 1.E-12_RP )
!ESC!       RLMDg  = sqrt(sqrt( dens * qg(ij,k) / ( Ag * N0g * GAM_1bg ) + zerosw )) * ( 1.0_RP - zerosw )
!ESC!       RLMDg_dg  = sqrt( RLMDg )       ! **Dg
!ESC!
!ESC!       !---< terminal velocity >---
!ESC!
!ESC!       vt_qc(ij,k) = 0.0_RP
!ESC!       zerosw = 0.5_RP + sign(0.5_RP, qi(ij,k) - 1.E-8_RP )
!ESC!       vt_qi(ij,k) = -3.29_RP * ( dens * qi(ij,k) * zerosw )**0.16_RP
!ESC!       vt_qr(ij,k) = -Cr * rho_fact * GAM_1brdr / GAM_1br * RLMDr_dr
!ESC!       vt_qs(ij,k) = -Cs * rho_fact * RMOMs_Vt
!ESC!       vt_qg(ij,k) = -Cg * rho_fact * GAM_1bgdg / GAM_1bg * RLMDg_dg
!ESC!
!ESC!    enddo
!ESC!    enddo
!ESC!
!ESC!    !$omp parallel do default(none),private(ij), &
!ESC!    !$omp shared(ijdim,kmin,kmax,vt_qc,vt_qr,vt_qi,vt_qs,vt_qg)
!ESC!    do ij = 1, ijdim
!ESC!       vt_qc(ij,kmin-1) = 0.0_RP
!ESC!       vt_qr(ij,kmin-1) = 0.0_RP
!ESC!       vt_qi(ij,kmin-1) = 0.0_RP
!ESC!       vt_qs(ij,kmin-1) = 0.0_RP
!ESC!       vt_qg(ij,kmin-1) = 0.0_RP
!ESC!
!ESC!       vt_qc(ij,kmax+1) = 0.0_RP
!ESC!       vt_qr(ij,kmax+1) = 0.0_RP
!ESC!       vt_qi(ij,kmax+1) = 0.0_RP
!ESC!       vt_qs(ij,kmax+1) = 0.0_RP
!ESC!       vt_qg(ij,kmax+1) = 0.0_RP
!ESC!    enddo
!ESC!    !$omp end parallel do
!ESC!
!ESC!    return
!ESC!  end subroutine mp_nsw6_terminal_velocity
!ESC!
!ESC!  !-----------------------------------------------------------------------------
!ESC!  subroutine mp_nsw6_effective_radius( &
!ESC!       ijdim, kdim, kmin, kmax, & ! in
!ESC!       rho, tem, pre, ccn,      & ! in
!ESC!       w,                       & ! in ! [Add] 2016/03/29 T.Seiki
!ESC!       qc,  qr,  qi,  qs,  qg,  & ! in
!ESC!       rec, rer, rei, res, reg  )
!ESC!    use mod_const, only: &
!ESC!       EPS   => CONST_EPS,  &
!ESC!       PI    => CONST_PI,   &
!ESC!       rho_w => CONST_DWATR
!ESC!    implicit none
!ESC!
!ESC!    integer,  intent(in)  :: ijdim
!ESC!    integer,  intent(in)  :: kdim
!ESC!    integer,  intent(in)  :: kmin
!ESC!    integer,  intent(in)  :: kmax
!ESC!    real(RP), intent(in)  :: rho(ijdim,kdim) ! air density [kg/m3]
!ESC!    real(RP), intent(in)  :: tem(ijdim,kdim) ! temperature [K]
!ESC!    real(RP), intent(in)  :: pre(ijdim,kdim) ! pressure [Pa]
!ESC!    real(RP), intent(in)  :: ccn(ijdim,kdim) ! cloud condensation nuclei [1/m3]
!ESC!    real(RP), intent(in)  :: w  (ijdim,kdim) ! vertical velocity [m/s]  [Add] 2016/03/29 T.Seiki
!ESC!    real(RP), intent(in)  :: qc (ijdim,kdim)
!ESC!    real(RP), intent(in)  :: qr (ijdim,kdim)
!ESC!    real(RP), intent(in)  :: qi (ijdim,kdim)
!ESC!    real(RP), intent(in)  :: qs (ijdim,kdim)
!ESC!    real(RP), intent(in)  :: qg (ijdim,kdim)
!ESC!    real(RP), intent(out) :: rec(ijdim,kdim)
!ESC!    real(RP), intent(out) :: rer(ijdim,kdim)
!ESC!    real(RP), intent(out) :: rei(ijdim,kdim)
!ESC!    real(RP), intent(out) :: res(ijdim,kdim)
!ESC!    real(RP), intent(out) :: reg(ijdim,kdim)
!ESC!
!ESC!    ! [Add] 2016/03/29 T.Seiki
!ESC!    real(RP), parameter :: rho_i = 916.7_RP
!ESC!    real(RP) :: nc
!ESC!    real(RP) :: Xni
!ESC!    real(RP) :: RLMDr, LMDs, LMDg
!ESC!    real(RP) :: factor, OneThird, zerosw
!ESC!
!ESC!    integer  :: ij, k
!ESC!    !---------------------------------------------------------------------------
!ESC!
!ESC!    OneThird = 1.0_RP / 3.0_RP
!ESC!
!ESC!    ! [Add] 2016/03/29 T.Seiki
!ESC!    factor = 3.0_RP / ( 4.0_RP * PI * rho_w )
!ESC!
!ESC!    if ( OPT_INDIR ) then
!ESC!       do k  = 1, kdim
!ESC!       do ij = 1, ijdim
!ESC!          nc = ccn(ij,k) ! [/m3]
!ESC!
!ESC!          rec(ij,k) = 1.1_RP * ( factor * rho(ij,k) * max( qc(ij,k), EPS ) / nc )**OneThird
!ESC!          rec(ij,k) = min( 1.E-3_RP, max( 1.E-6_RP, rec(ij,k) ) ) ! limiter
!ESC!       enddo
!ESC!       enddo
!ESC!    else
!ESC!       nc = Nc_def * 1.E+6_RP ! [/m3]
!ESC!       do k  = 1, kdim
!ESC!       do ij = 1, ijdim
!ESC!          rec(ij,k) = 1.1_RP * ( factor * rho(ij,k) * max( qc(ij,k), EPS ) / nc )**OneThird
!ESC!          rec(ij,k) = min( 1.E-3_RP, max( 1.E-6_RP, rec(ij,k) ) ) ! limiter
!ESC!       enddo
!ESC!       enddo
!ESC!    endif
!ESC!
!ESC!    factor = 3.0_RP / ( 4.0_RP * PI * rho_i )
!ESC!
!ESC!    if ( OPT_EXPLICIT_ICEGEN ) then
!ESC!       do k  = 1, kdim
!ESC!       do ij = 1, ijdim
!ESC!          Xni = 5.38E+7_RP * ( rho(ij,k) * max( qi(ij,k), EPS ) )**0.75_RP
!ESC!          Xni = min( 1.E+6_RP, max( 1.E+3_RP, Xni ) ) ! limiter
!ESC!
!ESC!          rei(ij,k) = ( factor * rho(ij,k) * max( qi(ij,k), EPS ) / Xni )**OneThird
!ESC!          rei(ij,k) = min( 1.E-3_RP, max( 1.E-6_RP, rei(ij,k) ) ) ! limiter
!ESC!       enddo
!ESC!       enddo
!ESC!    else
!ESC!       do k  = 1, kdim
!ESC!       do ij = 1, ijdim
!ESC!          rei(ij,k) = 40.E-6_RP
!ESC!       enddo
!ESC!       enddo
!ESC!    endif
!ESC!
!ESC!    factor = 3.0_RP / ( 2.0_RP * PI * rho_i )
!ESC!
!ESC!    do k  = 1, kdim
!ESC!    do ij = 1, ijdim
!ESC!       ! slope parameter lambda (Rain)
!ESC!       zerosw = 0.5_RP - sign(0.5_RP, qr(ij,k) - 1.E-12_RP )
!ESC!       RLMDr  = sqrt(sqrt( rho(ij,k) * qr(ij,k) / ( Ar * N0r * GAM_1br ) + zerosw )) * ( 1.0_RP - zerosw )
!ESC!       ! slope parameter lambda (Graupel)
!ESC!       zerosw = 0.5_RP - sign(0.5_RP, qg(ij,k) - 1.E-12_RP )
!ESC!       LMDg   = sqrt(sqrt( Ag * N0g * GAM_1bg / ( rho(ij,k) * qg(ij,k) + zerosw ) )) * ( 1.0_RP - zerosw )
!ESC!
!ESC!       rer(ij,k) = 1.5_RP * RLMDr
!ESC!       reg(ij,k) = factor * rho(ij,k) * qg(ij,k) / N0g * LMDg**3
!ESC!    enddo
!ESC!    enddo
!ESC!
!ESC!    ! default value
!ESC!    if ( Roh_flag ) then ! [add] WS.Roh 20141011
!ESC!
!ESC!       do k  = 1, kdim
!ESC!       do ij = 1, ijdim
!ESC!          ! ratio of (pi/4) to 0.45 is a correction to circle with maximum dimension.
!ESC!          ! 0.45 is derived from A=0.45*D**2 (approximateion of 0.2285*D**2 from Mitchell, 1996, JAS)
!ESC!          res(ij,k) = ( 3.0_RP * As / ( PI * rho_i ) ) * ( PI * 0.25_RP / 0.45_RP ) ! = 125um when as=0.069
!ESC!       enddo
!ESC!       enddo
!ESC!
!ESC!    else
!ESC!       factor = 3.0_RP / ( 2.0_RP * PI * rho_i )
!ESC!
!ESC!       do k  = 1, kdim
!ESC!       do ij = 1, ijdim
!ESC!          ! slope parameter lambda (Graupel)
!ESC!          zerosw = 0.5_RP - sign(0.5_RP, qs(ij,k) - 1.E-12_RP )
!ESC!          LMDs   = sqrt(sqrt( As * N0s * GAM_1bs / ( rho(ij,k) * qs(ij,k) + zerosw ) )) * ( 1.0_RP - zerosw )
!ESC!
!ESC!          res(ij,k) = factor * rho(ij,k) * qs(ij,k) / N0s * LMDs**3
!ESC!       enddo
!ESC!       enddo
!ESC!    endif
!ESC!
!ESC!    return
!ESC!  end subroutine mp_nsw6_effective_radius
!ESC!
!ESC!  !-----------------------------------------------------------------------------
!ESC!  subroutine ml2sl( &
!ESC!       ijdim, &
!ESC!       l,     &
!ESC!       rho,   &
!ESC!       v3d,   &
!ESC!       v2d    )
!ESC!    use mod_adm, only: &
!ESC!       k0   => ADM_KNONE, &
!ESC!       kdim => ADM_kall,  &
!ESC!       ADM_kmin,          &
!ESC!       ADM_kmax
!ESC!    use mod_gmtr, only: &
!ESC!       GMTR_area
!ESC!    use mod_vmtr, only: &
!ESC!       VMTR_VOLUME
!ESC!    implicit none
!ESC!
!ESC!    integer, intent(in)  :: ijdim
!ESC!    integer, intent(in)  :: l
!ESC!    real(RP), intent(in)  :: rho(ijdim,kdim)
!ESC!    real(RP), intent(in)  :: v3d(ijdim,kdim)
!ESC!    real(RP), intent(out) :: v2d(ijdim,k0)
!ESC!
!ESC!    integer :: k
!ESC!    !---------------------------------------------------------------------------
!ESC!
!ESC!    v2d(:,k0) = 0.0_RP
!ESC!    do k = ADM_kmin, ADM_kmax
!ESC!       v2d(:,k0) = v2d(:,k0) + rho(:,k) * v3d(:,k) * VMTR_VOLUME(:,k,l) / GMTR_area(:,l)
!ESC!    enddo
!ESC!
!ESC!  end subroutine ml2sl

end module mod_mp_nsw6
