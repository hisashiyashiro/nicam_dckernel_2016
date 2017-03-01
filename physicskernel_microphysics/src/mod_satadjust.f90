!-------------------------------------------------------------------------------
!> Module saturation process
!!
!! @par Description
!!         This module is for saturation processes
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_satadjust
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
  !++ Public parameters
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: SATURATION_setup
  public :: SATURATION_setrange

  public :: SATURATION_alpha
  public :: SATURATION_dalphadT

  public :: SATURATION_psat_all
  public :: SATURATION_psat_liq
  public :: SATURATION_psat_ice

  public :: SATURATION_qsat_liq
  public :: SATURATION_qsat_ice

  interface SATURATION_alpha
     module procedure SATURATION_alpha_0D
     module procedure SATURATION_alpha_1D
     module procedure SATURATION_alpha_2D
  end interface SATURATION_alpha

  interface SATURATION_dalphadT
     module procedure SATURATION_dalphadT_0D
     module procedure SATURATION_dalphadT_1D
     module procedure SATURATION_dalphadT_2D
  end interface SATURATION_dalphadT

  interface SATURATION_psat_all
     module procedure SATURATION_psat_all_0D
     module procedure SATURATION_psat_all_1D
     module procedure SATURATION_psat_all_2D
  end interface SATURATION_psat_all

  interface SATURATION_psat_liq
     module procedure SATURATION_psat_liq_0D
     module procedure SATURATION_psat_liq_1D
     module procedure SATURATION_psat_liq_2D
  end interface SATURATION_psat_liq

  interface SATURATION_psat_ice
     module procedure SATURATION_psat_ice_0D
     module procedure SATURATION_psat_ice_1D
     module procedure SATURATION_psat_ice_2D
  end interface SATURATION_psat_ice

  interface SATURATION_qsat_liq
     module procedure SATURATION_qsat_liq_0D
     module procedure SATURATION_qsat_liq_1D
     module procedure SATURATION_qsat_liq_2D
  end interface SATURATION_qsat_liq

  interface SATURATION_qsat_ice
     module procedure SATURATION_qsat_ice_0D
     module procedure SATURATION_qsat_ice_1D
     module procedure SATURATION_qsat_ice_2D
  end interface SATURATION_qsat_ice

  public :: SATURATION_rh
  public :: SATURATION_dewtem

  public :: SATURATION_adjustment

  public :: SATURATION_enthalpy
  public :: moist_dqsw_dtem_rho
  public :: moist_dqsi_dtem_rho
  public :: moist_dqsw_dtem_dpre
  public :: moist_dqsi_dtem_dpre

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  real(RP), public :: CPovR_liq
  real(RP), public :: CPovR_ice
  real(RP), public :: CVovR_liq
  real(RP), public :: CVovR_ice
  real(RP), public :: LovR_liq
  real(RP), public :: LovR_ice

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  private :: satadjust_all
  private :: satadjust_liq

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  real(RP), private, parameter :: TEM_MIN   = 10.0_RP   !< minimum temperature [K]
  real(RP), private, parameter :: DTEM_EPS0 = 1.0E-8_RP ! temperature convergence criterion

  character(len=H_SHORT), private :: ALPHA_TYPE             = 'LINEAR'

  real(RP),               private :: SATURATION_ULIMIT_TEMP = 273.15_RP     !< upper limit temperature
  real(RP),               private :: SATURATION_LLIMIT_TEMP = 233.15_RP     !< lower limit temperature

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> Setup
  subroutine SATURATION_setup
!ESC!    use mod_adm, only: &
!ESC!       ADM_proc_stop
!ESC!    use mod_const, only: &
!ESC!       CONST_Rvap,  &
!ESC!       CONST_CPvap, &
!ESC!       CONST_CVvap, &
!ESC!       CONST_CL,    &
!ESC!       CONST_CI
!ESC!    use mod_runconf, only: &
!ESC!       EIN_TYPE, &
!ESC!       LHV,      &
!ESC!       LHS
    implicit none

    NAMELIST / SATURATIONPARAM / &
       ALPHA_TYPE

    integer :: ierr
    !---------------------------------------------------------------------------

!ESC!    !--- read parameters
!ESC!    write(IO_FID_LOG,*)
!ESC!    write(IO_FID_LOG,*) '+++ Module[saturation]/Category[nhm share]'
!ESC!    rewind(IO_FID_CONF)
!ESC!    read(IO_FID_CONF,nml=SATURATIONPARAM,iostat=ierr)
!ESC!    if ( ierr < 0 ) then
!ESC!       write(IO_FID_LOG,*) '*** SATURATIONPARAM is not specified. use default.'
!ESC!    elseif( ierr > 0 ) then
!ESC!       write(*         ,*) 'xxx Not appropriate names in namelist SATURATIONPARAM. STOP.'
!ESC!       write(IO_FID_LOG,*) 'xxx Not appropriate names in namelist SATURATIONPARAM. STOP.'
!ESC!       call ADM_proc_stop
!ESC!    endif
!ESC!    write(IO_FID_LOG,nml=SATURATIONPARAM)

    if ( EIN_TYPE == 'EXACT' ) then
       CPovR_liq = ( CONST_CPvap - CONST_CL ) / CONST_Rvap
       CPovR_ice = ( CONST_CPvap - CONST_CI ) / CONST_Rvap
       CVovR_liq = ( CONST_CVvap - CONST_CL ) / CONST_Rvap
       CVovR_ice = ( CONST_CVvap - CONST_CI ) / CONST_Rvap
    elseif(      EIN_TYPE == 'SIMPLE'  &
            .OR. EIN_TYPE == 'SIMPLE2' ) then
       CPovR_liq = 0.0_RP
       CPovR_ice = 0.0_RP
       CVovR_liq = 0.0_RP
       CVovR_ice = 0.0_RP
    endif
    LovR_liq = LHV / CONST_Rvap
    LovR_ice = LHS / CONST_Rvap

    return
  end subroutine SATURATION_setup

  !-----------------------------------------------------------------------------
  subroutine SATURATION_setrange( Tw, Ti )
    implicit none

    real(RP), intent(in) :: Tw
    real(RP), intent(in) :: Ti
    !---------------------------------------------------------------------------

    SATURATION_ULIMIT_TEMP = Tw
    SATURATION_LLIMIT_TEMP = Ti

    return
  end subroutine SATURATION_setrange

  !-----------------------------------------------------------------------------
  !> calc liquid/ice separation factor (0D)
  subroutine SATURATION_alpha_0D( &
       tem,   &
       alpha  )
!ESC!    use mod_const, only: &
!ESC!       CONST_PI
    implicit none

    real(RP), intent(in)  :: tem   !< temperature [K]
    real(RP), intent(out) :: alpha !< liquid/ice separation factor [0-1]

    real(RP) :: fact
    real(RP) :: PI
    !---------------------------------------------------------------------------

    if    ( ALPHA_TYPE == 'COS' ) then

       PI = CONST_PI

       if    ( tem > SATURATION_ULIMIT_TEMP ) then
          alpha = 1.0_RP
       elseif( tem < SATURATION_LLIMIT_TEMP ) then
          alpha = 0.0_RP
       else
          fact = ( SATURATION_ULIMIT_TEMP - tem                    ) &
               / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )

          alpha = 0.5_RP * ( 1.0_RP + cos( fact * PI ) )
       endif

    elseif( ALPHA_TYPE == 'LINEAR' ) then

       fact = ( tem                    - SATURATION_LLIMIT_TEMP ) &
            / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )

       alpha = min( max( fact, 0.0_RP ), 1.0_RP )

    endif

    return
  end subroutine SATURATION_alpha_0D

  !-----------------------------------------------------------------------------
  !> calc liquid/ice separation factor (1D)
  subroutine SATURATION_alpha_1D( &
       ijdim, &
       tem,   &
       alpha  )
!ESC!    use mod_const, only: &
!ESC!       CONST_PI
    implicit none

    integer,  intent(in)  :: ijdim
    real(RP), intent(in)  :: tem  (ijdim) !< temperature [K]
    real(RP), intent(out) :: alpha(ijdim) !< liquid/ice separation factor [0-1]

    real(RP) :: fact
    real(RP) :: PI

    integer  :: ij
    !---------------------------------------------------------------------------

    if    ( ALPHA_TYPE == 'COS' ) then

       PI = CONST_PI

       !$omp parallel default(none),private(ij,fact), &
       !$omp shared(ijdim,alpha,tem,SATURATION_ULIMIT_TEMP,SATURATION_LLIMIT_TEMP,PI)

       !$omp do
       do ij = 1, ijdim
          fact = ( SATURATION_ULIMIT_TEMP - tem(ij)                ) &
               / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )

          alpha(ij) = 0.5_RP * ( 1.0_RP + cos( fact * PI ) )
       enddo
       !$omp end do

       !$omp do
       do ij = 1, ijdim
          if    ( tem(ij) > SATURATION_ULIMIT_TEMP ) then
             alpha(ij) = 1.0_RP
          elseif( tem(ij) < SATURATION_LLIMIT_TEMP ) then
             alpha(ij) = 0.0_RP
          endif
       enddo
       !$omp end do

       !$omp end parallel

    elseif( ALPHA_TYPE == 'LINEAR' ) then

       !$omp parallel do default(none),private(ij,fact), &
       !$omp shared(ijdim,alpha,tem,SATURATION_ULIMIT_TEMP,SATURATION_LLIMIT_TEMP)
       do ij = 1, ijdim
          fact = ( tem(ij)                - SATURATION_LLIMIT_TEMP ) &
               / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )

          alpha(ij) = min( max( fact, 0.0_RP ), 1.0_RP )
       enddo
       !$omp end parallel do

    endif

    return
  end subroutine SATURATION_alpha_1D

  !-----------------------------------------------------------------------------
  !> calc liquid/ice separation factor (2D)
  subroutine SATURATION_alpha_2D( &
       ijdim, &
       kdim,  &
       tem,   &
       alpha  )
!ESC!    use mod_const, only: &
!ESC!       CONST_PI
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(RP), intent(in)  :: tem  (ijdim,kdim) !< temperature [K]
    real(RP), intent(out) :: alpha(ijdim,kdim) !< liquid/ice separation factor [0-1]

    real(RP) :: fact
    real(RP) :: PI

    integer  :: ij, k
    !---------------------------------------------------------------------------

    if    ( ALPHA_TYPE == 'COS' ) then

       PI = CONST_PI

       !$omp parallel default(none),private(ij,k,fact), &
       !$omp shared(ijdim,kdim,alpha,tem,SATURATION_ULIMIT_TEMP,SATURATION_LLIMIT_TEMP,PI)
       do k  = 1, kdim

          !$omp do
          do ij = 1, ijdim
             fact = ( SATURATION_ULIMIT_TEMP - tem(ij,k)              ) &
                  / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )

             alpha(ij,k) = 0.5_RP * ( 1.0_RP + cos( fact * PI ) )
          enddo
          !$omp end do

          !$omp do
          do ij = 1, ijdim
             if    ( tem(ij,k) > SATURATION_ULIMIT_TEMP ) then
                alpha(ij,k) = 1.0_RP
             elseif( tem(ij,k) < SATURATION_LLIMIT_TEMP ) then
                alpha(ij,k) = 0.0_RP
             endif
          enddo
          !$omp end do

       enddo
       !$omp end parallel

    elseif( ALPHA_TYPE == 'LINEAR' ) then

       !$omp parallel do default(none),private(ij,k,fact), &
       !$omp shared(ijdim,kdim,alpha,tem,SATURATION_ULIMIT_TEMP,SATURATION_LLIMIT_TEMP)
       do k  = 1, kdim
       do ij = 1, ijdim
          fact = ( tem(ij,k)              - SATURATION_LLIMIT_TEMP ) &
               / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )

          alpha(ij,k) = min( max( fact, 0.0_RP ), 1.0_RP )
       enddo
       enddo
       !$omp end parallel do

    endif

    return
  end subroutine SATURATION_alpha_2D

  !-----------------------------------------------------------------------------
  !> calc d(alpha)/d(tem), 0D
  subroutine SATURATION_dalphadT_0D( &
       tem,      &
       dalpha_dT )
!ESC!    use mod_const, only: &
!ESC!       CONST_PI
    implicit none

    real(RP), intent(in)  :: tem       !< temperature [K]
    real(RP), intent(out) :: dalpha_dT !< d(alpha)/d(T)

    real(RP) :: lim1, lim2, fact
    real(RP) :: PI
    !---------------------------------------------------------------------------

    if    ( ALPHA_TYPE == 'COS' ) then

       PI = CONST_PI

       ! if Tup < tem, dalpha/dT = 0 (no slope)
       lim1 = 0.5_RP + sign( 0.5_RP, SATURATION_ULIMIT_TEMP - tem )
       ! if Tdn > tem, dalpha/dT = 0 (no slope)
       lim2 = 0.5_RP + sign( 0.5_RP, tem - SATURATION_LLIMIT_TEMP )

       fact = ( SATURATION_ULIMIT_TEMP - tem                    ) &
            / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )

       dalpha_dT = lim1 * lim2 / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP ) &
                 * 0.5_RP * PI * sin( fact * PI )

    elseif( ALPHA_TYPE == 'LINEAR' ) then

       ! if Tup < tem, dalpha/dT = 0 (no slope)
       lim1 = 0.5_RP + sign( 0.5_RP, SATURATION_ULIMIT_TEMP - tem )
       ! if Tdn > tem, dalpha/dT = 0 (no slope)
       lim2 = 0.5_RP + sign( 0.5_RP, tem - SATURATION_LLIMIT_TEMP )

       dalpha_dT = lim1 * lim2 / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )

    endif

    return
  end subroutine SATURATION_dalphadT_0D

  !-----------------------------------------------------------------------------
  !> calc d(alpha)/d(tem), 1D
  subroutine SATURATION_dalphadT_1D( &
       ijdim,    &
       tem,      &
       dalpha_dT )
!ESC!    use mod_const, only: &
!ESC!       CONST_PI
    implicit none

    integer,  intent(in)  :: ijdim
    real(RP), intent(in)  :: tem      (ijdim) !< temperature [K]
    real(RP), intent(out) :: dalpha_dT(ijdim) !< d(alpha)/d(T)

    real(RP) :: lim1, lim2, fact
    real(RP) :: PI

    integer  :: ij
    !---------------------------------------------------------------------------

    if    ( ALPHA_TYPE == 'COS' ) then

       PI = CONST_PI

       !$omp parallel do default(none),private(ij,fact,lim1,lim2), &
       !$omp shared(ijdim,dalpha_dT,tem,SATURATION_ULIMIT_TEMP,SATURATION_LLIMIT_TEMP,PI)
       do ij = 1, ijdim
          ! if Tup < tem, dalpha/dT = 0 (no slope)
          lim1 = 0.5_RP + sign( 0.5_RP, SATURATION_ULIMIT_TEMP - tem(ij) )
          ! if Tdn > tem, dalpha/dT = 0 (no slope)
          lim2 = 0.5_RP + sign( 0.5_RP, tem(ij) - SATURATION_LLIMIT_TEMP )

          fact = ( SATURATION_ULIMIT_TEMP - tem(ij)                ) &
               / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )

          dalpha_dT(ij) = lim1 * lim2 / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP ) &
                        * 0.5_RP * PI * sin( fact * PI )
       enddo
       !$omp end parallel do

    elseif( ALPHA_TYPE == 'LINEAR' ) then

       !$omp parallel do default(none),private(ij,lim1,lim2), &
       !$omp shared(ijdim,dalpha_dT,tem,SATURATION_ULIMIT_TEMP,SATURATION_LLIMIT_TEMP)
       do ij = 1, ijdim
          ! if Tup < tem, dalpha/dT = 0 (no slope)
          lim1 = 0.5_RP + sign( 0.5_RP, SATURATION_ULIMIT_TEMP - tem(ij) )
          ! if Tdn > tem, dalpha/dT = 0 (no slope)
          lim2 = 0.5_RP + sign( 0.5_RP, tem(ij) - SATURATION_LLIMIT_TEMP )

          dalpha_dT(ij) = lim1 * lim2 / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )
       enddo
       !$omp end parallel do

    endif

    return
  end subroutine SATURATION_dalphadT_1D

  !-----------------------------------------------------------------------------
  !> calc d(alpha)/d(temp), 2D
  subroutine SATURATION_dalphadT_2D( &
       ijdim,    &
       kdim,     &
       tem,      &
       dalpha_dT )
!ESC!    use mod_const, only: &
!ESC!       CONST_PI
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(RP), intent(in)  :: tem      (ijdim,kdim) !< temperature [K]
    real(RP), intent(out) :: dalpha_dT(ijdim,kdim) !< d(alpha)/d(T)

    real(RP) :: lim1, lim2, fact
    real(RP) :: PI

    integer  :: ij, k
    !---------------------------------------------------------------------------

    if    ( ALPHA_TYPE == 'COS' ) then

       PI = CONST_PI

       !$omp parallel do default(none),private(ij,k,fact,lim1,lim2), &
       !$omp shared(ijdim,kdim,dalpha_dT,tem,SATURATION_ULIMIT_TEMP,SATURATION_LLIMIT_TEMP,PI)
       do k  = 1, kdim
       do ij = 1, ijdim
          ! if Tup < tem, dalpha/dT = 0 (no slope)
          lim1 = 0.5_RP + sign( 0.5_RP, SATURATION_ULIMIT_TEMP - tem(ij,k) )
          ! if Tdn > tem, dalpha/dT = 0 (no slope)
          lim2 = 0.5_RP + sign( 0.5_RP, tem(ij,k) - SATURATION_LLIMIT_TEMP )

          fact = ( SATURATION_ULIMIT_TEMP - tem(ij,k)              ) &
               / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )

          dalpha_dT(ij,k) = lim1 * lim2 / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP ) &
                            * 0.5_RP * PI * sin( fact * PI )
       enddo
       enddo
       !$omp end parallel do

    elseif( ALPHA_TYPE == 'LINEAR' ) then

       !$omp parallel do default(none),private(ij,k,lim1,lim2), &
       !$omp shared(ijdim,kdim,dalpha_dT,tem,SATURATION_ULIMIT_TEMP,SATURATION_LLIMIT_TEMP)
       do k  = 1, kdim
       do ij = 1, ijdim
          ! if Tup < tem, dalpha/dT = 0 (no slope)
          lim1 = 0.5_RP + sign( 0.5_RP, SATURATION_ULIMIT_TEMP - tem(ij,k) )
          ! if Tdn > tem, dalpha/dT = 0 (no slope)
          lim2 = 0.5_RP + sign( 0.5_RP, tem(ij,k) - SATURATION_LLIMIT_TEMP )

          dalpha_dT(ij,k) = lim1 * lim2 / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )
       enddo
       enddo
       !$omp end parallel do

    endif

    return
  end subroutine SATURATION_dalphadT_2D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure (liquid/ice mixture) (0D)
  subroutine SATURATION_psat_all_0D( &
       tem,   &
       psat   )
    implicit none

    real(RP), intent(in)  :: tem  !< temperature               [K]
    real(RP), intent(out) :: psat !< saturation vapor pressure [Pa]

    real(RP) :: alpha, psatl, psati
    !---------------------------------------------------------------------------

    call SATURATION_alpha   ( tem, alpha )
    call SATURATION_psat_liq( tem, psatl )
    call SATURATION_psat_ice( tem, psati )

    psat = psatl * (          alpha ) &
         + psati * ( 1.0_RP - alpha )

    return
  end subroutine SATURATION_psat_all_0D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure (liquid/ice mixture) (1D)
  subroutine SATURATION_psat_all_1D( &
       ijdim, &
       tem,   &
       psat   )
    implicit none

    integer,  intent(in)  :: ijdim
    real(RP), intent(in)  :: tem (ijdim) !< temperature               [K]
    real(RP), intent(out) :: psat(ijdim) !< saturation vapor pressure [Pa]

    real(RP) :: alpha(ijdim), psatl(ijdim), psati(ijdim)

    integer  :: ij
    !---------------------------------------------------------------------------

    call SATURATION_alpha   ( ijdim, tem(:), alpha(:) )
    call SATURATION_psat_liq( ijdim, tem(:), psatl(:) )
    call SATURATION_psat_ice( ijdim, tem(:), psati(:) )

    do ij = 1, ijdim
       psat(ij) = psatl(ij) * (          alpha(ij) ) &
                + psati(ij) * ( 1.0_RP - alpha(ij) )
    enddo

    return
  end subroutine SATURATION_psat_all_1D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure (liquid/ice mixture) (2D)
  subroutine SATURATION_psat_all_2D( &
       ijdim, &
       kdim,  &
       tem,   &
       psat   )
!ESC!    use mod_const, only: &
!ESC!       CONST_PI,    &
!ESC!       CONST_TEM00, &
!ESC!       CONST_PSAT0
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(RP), intent(in)  :: tem (ijdim,kdim) !< temperature               [K]
    real(RP), intent(out) :: psat(ijdim,kdim) !< saturation vapor pressure [Pa]

    real(RP) :: alpha, psatl, psati
    real(RP) :: rtem
    real(RP) :: PI, RTEM00, PSAT0

    integer  :: ij, k
    !---------------------------------------------------------------------------

    RTEM00 = 1.0_RP / CONST_TEM00
    PSAT0  = CONST_PSAT0

    if    ( ALPHA_TYPE == 'COS' ) then

       PI = CONST_PI

       !$omp parallel do default(none),private(ij,k,alpha,rtem,psatl,psati), &
       !$omp shared(ijdim,kdim,psat,tem,SATURATION_ULIMIT_TEMP,SATURATION_LLIMIT_TEMP,PI, &
       !$omp        CPovR_liq,LovR_liq,CPovR_ice,LovR_ice,RTEM00,PSAT0)
       do k  = 1, kdim
       do ij = 1, ijdim

          if    ( tem(ij,k) > SATURATION_ULIMIT_TEMP ) then
             alpha = 1.0_RP
          elseif( tem(ij,k) < SATURATION_LLIMIT_TEMP ) then
             alpha = 0.0_RP
          else
             alpha = ( SATURATION_ULIMIT_TEMP - tem(ij,k)              ) &
                   / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )

             alpha = 0.5_RP * ( 1.0_RP + cos( alpha * PI ) )
          endif

          rtem = 1.0_RP / max( tem(ij,k), TEM_MIN )

          psatl = PSAT0 * ( tem(ij,k) * RTEM00 )**CPovR_liq &
                        * exp( LovR_liq * ( RTEM00 - rtem ) )

          psati = PSAT0 * ( tem(ij,k) * RTEM00 )**CPovR_ice &
                        * exp( LovR_ice * ( RTEM00 - rtem ) )

          psat(ij,k) = psatl * (          alpha ) &
                     + psati * ( 1.0_RP - alpha )
       enddo
       enddo
       !$omp end parallel do

    elseif( ALPHA_TYPE == 'LINEAR' ) then

       !$omp parallel do default(none),private(ij,k,alpha,rtem,psatl,psati), &
       !$omp shared(ijdim,kdim,psat,tem,SATURATION_ULIMIT_TEMP,SATURATION_LLIMIT_TEMP, &
       !$omp        CPovR_liq,LovR_liq,CPovR_ice,LovR_ice,RTEM00,PSAT0)
       do k  = 1, kdim
       do ij = 1, ijdim
          alpha = ( tem(ij,k)              - SATURATION_LLIMIT_TEMP ) &
                / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )

          alpha = min( max( alpha, 0.0_RP ), 1.0_RP )

          rtem = 1.0_RP / max( tem(ij,k), TEM_MIN )

          psatl = PSAT0 * ( tem(ij,k) * RTEM00 )**CPovR_liq &
                        * exp( LovR_liq * ( RTEM00 - rtem ) )

          psati = PSAT0 * ( tem(ij,k) * RTEM00 )**CPovR_ice &
                        * exp( LovR_ice * ( RTEM00 - rtem ) )

          psat(ij,k) = psatl * (          alpha ) &
                     + psati * ( 1.0_RP - alpha )
       enddo
       enddo
       !$omp end parallel do

    endif

    return
  end subroutine SATURATION_psat_all_2D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (0D)
  subroutine SATURATION_psat_liq_0D( &
       tem,   &
       psat   )
!ESC!    use mod_const, only: &
!ESC!       CONST_TEM00, &
!ESC!       CONST_PSAT0
    implicit none

    real(RP), intent(in)  :: tem  !< temperature               [K]
    real(RP), intent(out) :: psat !< saturation vapor pressure [Pa]

    real(RP) :: rtem
    real(RP) :: RTEM00, PSAT0
    !---------------------------------------------------------------------------

    RTEM00 = 1.0_RP / CONST_TEM00
    PSAT0  = CONST_PSAT0

    rtem = 1.0_RP / max( tem, TEM_MIN )

    psat = PSAT0 * ( tem * RTEM00 )**CPovR_liq &
                 * exp( LovR_liq * ( RTEM00 - rtem ) )

    return
  end subroutine SATURATION_psat_liq_0D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (1D)
  subroutine SATURATION_psat_liq_1D( &
       ijdim, &
       tem,   &
       psat   )
!ESC!    use mod_const, only: &
!ESC!       CONST_TEM00, &
!ESC!       CONST_PSAT0
    implicit none

    integer,  intent(in)  :: ijdim
    real(RP), intent(in)  :: tem (ijdim) !< temperature               [K]
    real(RP), intent(out) :: psat(ijdim) !< saturation vapor pressure [Pa]

    real(RP) :: rtem
    real(RP) :: RTEM00, PSAT0

    integer  :: ij
    !---------------------------------------------------------------------------

    RTEM00 = 1.0_RP / CONST_TEM00
    PSAT0  = CONST_PSAT0

    !$omp parallel do default(none),private(ij,rtem), &
    !$omp shared(ijdim,psat,tem,CPovR_liq,LovR_liq,RTEM00,PSAT0)
    do ij = 1, ijdim
       rtem = 1.0_RP / max( tem(ij), TEM_MIN )

       psat(ij) = PSAT0 * ( tem(ij) * RTEM00 )**CPovR_liq &
                        * exp( LovR_liq * ( RTEM00 - rtem ) )
    enddo
    !$omp end parallel do

    return
  end subroutine SATURATION_psat_liq_1D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (2D)
  subroutine SATURATION_psat_liq_2D( &
       ijdim, &
       kdim,  &
       tem,   &
       psat   )
!ESC!    use mod_const, only: &
!ESC!       CONST_TEM00, &
!ESC!       CONST_PSAT0
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(RP), intent(in)  :: tem (ijdim,kdim) !< temperature               [K]
    real(RP), intent(out) :: psat(ijdim,kdim) !< saturation vapor pressure [Pa]

    real(RP) :: rtem
    real(RP) :: RTEM00, PSAT0

    integer  :: ij, k
    !---------------------------------------------------------------------------

    RTEM00 = 1.0_RP / CONST_TEM00
    PSAT0  = CONST_PSAT0

    !$omp parallel do default(none),private(ij,k,rtem), &
    !$omp shared(ijdim,kdim,psat,tem,CPovR_liq,LovR_liq,RTEM00,PSAT0)
    do k  = 1, kdim
    do ij = 1, ijdim
       rtem = 1.0_RP / max( tem(ij,k), TEM_MIN )

       psat(ij,k) = PSAT0 * ( tem(ij,k) * RTEM00 )**CPovR_liq &
                          * exp( LovR_liq * ( RTEM00 - rtem ) )
    enddo
    enddo
    !$omp end parallel do

    return
  end subroutine SATURATION_psat_liq_2D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (0D)
  subroutine SATURATION_psat_ice_0D( &
       tem,   &
       psat   )
!ESC!    use mod_const, only: &
!ESC!       CONST_TEM00, &
!ESC!       CONST_PSAT0
    implicit none

    real(RP), intent(in)  :: tem  !< temperature               [K]
    real(RP), intent(out) :: psat !< saturation vapor pressure [Pa]

    real(RP) :: rtem
    real(RP) :: RTEM00, PSAT0
    !---------------------------------------------------------------------------

    RTEM00 = 1.0_RP / CONST_TEM00
    PSAT0  = CONST_PSAT0

    rtem = 1.0_RP / max( tem, TEM_MIN )

    psat = PSAT0 * ( tem * RTEM00 )**CPovR_ice &
                 * exp( LovR_ice * ( RTEM00 - rtem ) )

    return
  end subroutine SATURATION_psat_ice_0D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (1D)
  subroutine SATURATION_psat_ice_1D( &
       ijdim, &
       tem,   &
       psat   )
!ESC!    use mod_const, only: &
!ESC!       CONST_TEM00, &
!ESC!       CONST_PSAT0
    implicit none

    integer,  intent(in)  :: ijdim
    real(RP), intent(in)  :: tem (ijdim) !< temperature               [K]
    real(RP), intent(out) :: psat(ijdim) !< saturation vapor pressure [Pa]

    real(RP) :: rtem
    real(RP) :: RTEM00, PSAT0

    integer  :: ij
    !---------------------------------------------------------------------------

    RTEM00 = 1.0_RP / CONST_TEM00
    PSAT0  = CONST_PSAT0

    !$omp parallel do default(none),private(ij,rtem), &
    !$omp shared(ijdim,psat,tem,CPovR_ice,LovR_ice,RTEM00,PSAT0)
    do ij = 1, ijdim
       rtem = 1.0_RP / max( tem(ij), TEM_MIN )

       psat(ij) = PSAT0 * ( tem(ij) * RTEM00 )**CPovR_ice &
                        * exp( LovR_ice * ( RTEM00 - rtem ) )
    enddo
    !$omp end parallel do

    return
  end subroutine SATURATION_psat_ice_1D

  !-----------------------------------------------------------------------------
  !> calc saturation vapor pressure from Clausius-Clapeyron equation (2D)
  subroutine SATURATION_psat_ice_2D( &
       ijdim, &
       kdim,  &
       tem,   &
       psat   )
!ESC!    use mod_const, only: &
!ESC!       CONST_TEM00, &
!ESC!       CONST_PSAT0
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(RP), intent(in)  :: tem (ijdim,kdim) !< temperature               [K]
    real(RP), intent(out) :: psat(ijdim,kdim) !< saturation vapor pressure [Pa]

    real(RP) :: rtem
    real(RP) :: RTEM00, PSAT0

    integer  :: ij, k
    !---------------------------------------------------------------------------

    RTEM00 = 1.0_RP / CONST_TEM00
    PSAT0  = CONST_PSAT0

    !$omp parallel do default(none),private(ij,k,rtem), &
    !$omp shared(ijdim,kdim,psat,tem,CPovR_ice,LovR_ice,RTEM00,PSAT0)
    do k  = 1, kdim
    do ij = 1, ijdim
       rtem = 1.0_RP / max( tem(ij,k), TEM_MIN )

       psat(ij,k) = PSAT0 * ( tem(ij,k) * RTEM00 )**CPovR_ice &
                          * exp( LovR_ice * ( RTEM00 - rtem ) )
    enddo
    enddo
    !$omp end parallel do

    return
  end subroutine SATURATION_psat_ice_2D

  !-----------------------------------------------------------------------------
  subroutine SATURATION_qsat_liq_0D( &
       tem,   &
       pre,   &
       qsat   )
!ESC!    use mod_const, only: &
!ESC!       CONST_EPSvap
    implicit none

    real(RP), intent(in)  :: pre
    real(RP), intent(in)  :: tem
    real(RP), intent(out) :: qsat

    real(RP) :: psat
    real(RP) :: EPSvap
    !---------------------------------------------------------------------------

    EPSvap = CONST_EPSvap

    call SATURATION_psat_liq_0D( tem, psat )

    qsat = EPSvap * psat / ( pre - ( 1.0_RP-EPSvap ) * psat )

    return
  end subroutine SATURATION_qsat_liq_0D

  !-----------------------------------------------------------------------------
  subroutine SATURATION_qsat_liq_1D( &
       ijdim, &
       tem,   &
       pre,   &
       qsat   )
!ESC!    use mod_const, only: &
!ESC!       CONST_EPSvap
    implicit none

    integer,  intent(in)  :: ijdim
    real(RP), intent(in)  :: pre (ijdim)
    real(RP), intent(in)  :: tem (ijdim)
    real(RP), intent(out) :: qsat(ijdim)

    real(RP) :: psat(ijdim)
    real(RP) :: EPSvap

    integer  :: ij
    !---------------------------------------------------------------------------

    EPSvap = CONST_EPSvap

    call SATURATION_psat_liq_1D( ijdim, tem(:), psat(:) )

    !$omp parallel do default(none),private(ij), &
    !$omp shared(ijdim,qsat,pre,psat,EPSvap)
    do ij = 1, ijdim
       qsat(ij) = EPSvap * psat(ij) / ( pre(ij) - ( 1.0_RP-EPSvap ) * psat(ij) )
    enddo
    !$omp end parallel do

    return
  end subroutine SATURATION_qsat_liq_1D

  !-----------------------------------------------------------------------------
  subroutine SATURATION_qsat_liq_2D( &
       ijdim, &
       kdim,  &
       tem,   &
       pre,   &
       qsat   )
!ESC!    use mod_const, only: &
!ESC!       CONST_EPSvap
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(RP), intent(in)  :: pre (ijdim,kdim)
    real(RP), intent(in)  :: tem (ijdim,kdim)
    real(RP), intent(out) :: qsat(ijdim,kdim)

    real(RP) :: psat(ijdim,kdim)
    real(RP) :: EPSvap

    integer  :: ij, k
    !---------------------------------------------------------------------------

    EPSvap = CONST_EPSvap

    call SATURATION_psat_liq_2D( ijdim, kdim, tem(:,:), psat(:,:) )

    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,qsat,pre,psat,EPSvap)
    do k  = 1, kdim
    do ij = 1, ijdim
       qsat(ij,k) = EPSvap * psat(ij,k) / ( pre(ij,k) - ( 1.0_RP-EPSvap ) * psat(ij,k) )
    enddo
    enddo
    !$omp end parallel do

    return
  end subroutine SATURATION_qsat_liq_2D

  !-----------------------------------------------------------------------------
  subroutine SATURATION_qsat_ice_0D( &
       tem,   &
       pre,   &
       qsat   )
!ESC!    use mod_const, only: &
!ESC!       CONST_EPSvap
    implicit none

    real(RP), intent(in)  :: pre
    real(RP), intent(in)  :: tem
    real(RP), intent(out) :: qsat

    real(RP) :: psat
    real(RP) :: EPSvap
    !---------------------------------------------------------------------------

    EPSvap = CONST_EPSvap

    call SATURATION_psat_ice_0D( tem, psat )

    qsat = EPSvap * psat / ( pre - ( 1.0_RP-EPSvap ) * psat )

    return
  end subroutine SATURATION_qsat_ice_0D

  !-----------------------------------------------------------------------------
  subroutine SATURATION_qsat_ice_1D( &
       ijdim, &
       tem,   &
       pre,   &
       qsat   )
!ESC!    use mod_const, only: &
!ESC!       CONST_EPSvap
    implicit none

    integer,  intent(in)  :: ijdim
    real(RP), intent(in)  :: pre (ijdim)
    real(RP), intent(in)  :: tem (ijdim)
    real(RP), intent(out) :: qsat(ijdim)

    real(RP) :: psat(ijdim)
    real(RP) :: EPSvap

    integer  :: ij
    !---------------------------------------------------------------------------

    EPSvap = CONST_EPSvap

    call SATURATION_psat_ice_1D( ijdim, tem(:), psat(:) )

    !$omp parallel do default(none),private(ij), &
    !$omp shared(ijdim,qsat,pre,psat,EPSvap)
    do ij = 1, ijdim
       qsat(ij) = EPSvap * psat(ij) / ( pre(ij) - ( 1.0_RP-EPSvap ) * psat(ij) )
    enddo
    !$omp end parallel do

    return
  end subroutine SATURATION_qsat_ice_1D

  !-----------------------------------------------------------------------------
  subroutine SATURATION_qsat_ice_2D( &
       ijdim, &
       kdim,  &
       tem,   &
       pre,   &
       qsat   )
!ESC!    use mod_const, only: &
!ESC!       CONST_EPSvap
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(RP), intent(in)  :: pre (ijdim,kdim)
    real(RP), intent(in)  :: tem (ijdim,kdim)
    real(RP), intent(out) :: qsat(ijdim,kdim)

    real(RP) :: psat(ijdim,kdim)
    real(RP) :: EPSvap

    integer  :: ij, k
    !---------------------------------------------------------------------------

    EPSvap = CONST_EPSvap

    call SATURATION_psat_ice_2D( ijdim, kdim, tem(:,:), psat(:,:) )

    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,qsat,pre,psat,EPSvap)
    do k  = 1, kdim
    do ij = 1, ijdim
       qsat(ij,k) = EPSvap * psat(ij,k) / ( pre(ij,k) - ( 1.0_RP-EPSvap ) * psat(ij,k) )
    enddo
    enddo
    !$omp end parallel do

    return
  end subroutine SATURATION_qsat_ice_2D

  !-----------------------------------------------------------------------------
  subroutine SATURATION_rh( &
       ijdim, &
       kdim,  &
       rho,   &
       tem,   &
       qv,    &
       rh     )
!ESC!    use mod_const, only: &
!ESC!       CONST_Rvap, &
!ESC!       CONST_TMELT
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(RP), intent(in)  :: rho(ijdim,kdim)
    real(RP), intent(in)  :: tem(ijdim,kdim)
    real(RP), intent(in)  :: qv (ijdim,kdim)
    real(RP), intent(out) :: rh (ijdim,kdim)

    real(RP) :: psat_liq(ijdim,kdim)
    real(RP) :: psat_ice(ijdim,kdim)
    real(RP) :: rh_liq, rh_ice, alpha
    real(RP) :: Rvap, TMELT

    integer :: ij, k
    !---------------------------------------------------------------------------

    Rvap  = CONST_Rvap
    TMELT = CONST_TMELT

    call SATURATION_psat_liq( ijdim, kdim, tem(:,:), psat_liq(:,:) )
    call SATURATION_psat_ice( ijdim, kdim, tem(:,:), psat_ice(:,:) )

    !$omp parallel do default(none),private(ij,k,rh_liq,rh_ice,alpha), &
    !$omp shared(ijdim,kdim,rh,rho,tem,qv,psat_liq,psat_ice,Rvap,TMELT)
    do k  = 1, kdim
    do ij = 1, ijdim
       rh_liq = qv(ij,k) / psat_liq(ij,k) * ( rho(ij,k) * Rvap * tem(ij,k) )
       rh_ice = qv(ij,k) / psat_ice(ij,k) * ( rho(ij,k) * Rvap * tem(ij,k) )
       alpha  = 0.5_RP + sign(0.5_RP,tem(ij,k)-TMELT)

       rh(ij,k) = (        alpha ) * rh_liq &
                + ( 1.0_RP-alpha ) * rh_ice
    enddo
    enddo
    !$omp end parallel do

    return
  end subroutine SATURATION_rh

  !-----------------------------------------------------------------------------
  subroutine SATURATION_dewtem( &
       ijdim, &
       kdim,  &
       tem,   &
       pre,   &
       qd,    &
       qv,    &
       wtem   )
!ESC!    use mod_const, only: &
!ESC!       CONST_EPSvap, &
!ESC!       CONST_PSAT0,  &
!ESC!       CONST_TEM00
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(RP), intent(in)  :: tem (ijdim,kdim)
    real(RP), intent(in)  :: pre (ijdim,kdim)
    real(RP), intent(in)  :: qd  (ijdim,kdim)
    real(RP), intent(in)  :: qv  (ijdim,kdim)
    real(RP), intent(out) :: wtem(ijdim,kdim)       ! dew point temperature

    real(RP) :: prev(ijdim,kdim)
    logical  :: flag(ijdim,kdim)
    real(RP) :: rtem, lv, psatl, f, dfdtem
    real(RP) :: EPSvap, RTEM00, PSAT0

    real(RP), parameter :: criteria = 1.E-8_RP
    integer,  parameter :: itelim   = 20

    integer  :: ite
    integer  :: ij, k
    !---------------------------------------------------------------------------

    EPSvap = CONST_EPSvap
    RTEM00 = 1.0_RP / CONST_TEM00
    PSAT0  = CONST_PSAT0

    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,prev,wtem,flag,tem,pre,qd,qv,EPSvap)
    do k  = 1, kdim
    do ij = 1, ijdim
       prev(ij,k) = pre(ij,k) * qv(ij,k) / ( EPSvap * qd(ij,k) + qv(ij,k) )

       wtem(ij,k) = tem(ij,k) * 0.98_RP

       flag(ij,k) = .false.
    enddo
    enddo
    !$omp end parallel do

    !$omp parallel default(none),private(ij,k,rtem,lv,psatl,f,dfdtem), &
    !$omp shared(ite,ijdim,kdim,wtem,flag,prev,CPovR_liq,LovR_liq,RTEM00,PSAT0)
    do ite = 1, itelim
       !$omp do
       do k  = 1, kdim
       do ij = 1, ijdim
          if ( .NOT. flag(ij,k) ) then
             rtem   = 1.0_RP / max( wtem(ij,k), TEM_MIN )

             lv     = LovR_liq + wtem(ij,k) * CPovR_liq

             psatl  = PSAT0 * ( wtem(ij,k) * RTEM00 )**CPovR_liq &
                            * exp( LovR_liq * ( RTEM00 - rtem ) )

             f      = psatl - prev(ij,k)
             dfdtem = lv * psatl * rtem**2

             wtem(ij,k) = wtem(ij,k) - f / dfdtem

             if ( abs(f) < prev(ij,k) * criteria ) then
                flag(ij,k) = .true.
             endif
          endif
       enddo
       enddo
       !$omp end do
    enddo
    !$omp end parallel

    return
  end subroutine SATURATION_dewtem

  !-----------------------------------------------------------------------------
  subroutine SATURATION_adjustment( &
       ijdim,     &
       kdim,      &
       rhog,      &
       rhoge,     &
       rhogq,     &
       tem,       &
       q,         &
       qd,        &
       gsgam2,    &
       ice_adjust )
!ESC!    use mod_runconf, only: &
!ESC!       nqmax => TRC_VMAX, &
!ESC!       I_QV,              &
!ESC!       I_QC,              &
!ESC!       I_QI,              &
!ESC!       LHV,               &
!ESC!       LHF
    use mod_thrmdyn, only: &
       thrmdyn_cv
    implicit none

    integer,  intent(in)    :: ijdim
    integer,  intent(in)    :: kdim
    real(RP), intent(in)    :: rhog  (ijdim,kdim)
    real(RP), intent(inout) :: rhoge (ijdim,kdim)
    real(RP), intent(inout) :: rhogq (ijdim,kdim,nqmax)
    real(RP), intent(inout) :: tem   (ijdim,kdim)
    real(RP), intent(inout) :: q     (ijdim,kdim,nqmax)
    real(RP), intent(in)    :: qd    (ijdim,kdim)
    real(RP), intent(in)    :: gsgam2(ijdim,kdim)
    logical,  intent(in)    :: ice_adjust

    real(RP) :: ein_moist(ijdim,kdim)
    real(RP) :: qsum     (ijdim,kdim)
    real(RP) :: CVtot    (ijdim,kdim)
    real(RP) :: rho      (ijdim,kdim)

    integer  :: ij, k
    !---------------------------------------------------------------------------

    call PROF_rapstart('____Saturation_Adjustment')

    ! ein_moist = U1(rho,qsum,T1) : "unsaturated temperature"
    if ( I_QI > 0 .AND. ice_adjust ) then
       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kdim,ein_moist,qsum,rhoge,rhog,q,LHV,LHF,I_QV,I_QC,I_QI)
       do k  = 1, kdim
       do ij = 1, ijdim
          ein_moist(ij,k) = rhoge(ij,k) / rhog(ij,k) &
                          + q(ij,k,I_QV) * LHV &
                          - q(ij,k,I_QI) * LHF

          qsum     (ij,k) = q(ij,k,I_QV) &
                          + q(ij,k,I_QC) &
                          + q(ij,k,I_QI)

          q(ij,k,I_QV) = qsum(ij,k)
          q(ij,k,I_QC) = 0.0_RP
          q(ij,k,I_QI) = 0.0_RP
       enddo
       enddo
       !$omp end parallel do
    else
       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kdim,ein_moist,qsum,rhoge,rhog,q,LHV,I_QV,I_QC)
       do k  = 1, kdim
       do ij = 1, ijdim
          ein_moist(ij,k) = rhoge(ij,k) / rhog(ij,k) &
                          + q(ij,k,I_QV) * LHV

          qsum     (ij,k) = q(ij,k,I_QV) &
                          + q(ij,k,I_QC)

          q(ij,k,I_QV) = qsum(ij,k)
          q(ij,k,I_QC) = 0.0_RP
       enddo
       enddo
       !$omp end parallel do
    endif

    call THRMDYN_cv( ijdim,        & ! [IN]
                     kdim,         & ! [IN]
                     qd   (:,:),   & ! [IN]
                     q    (:,:,:), & ! [IN]
                     CVtot(:,:)    ) ! [OUT]

    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,rho,tem,ein_moist,q,rhog,gsgam2,CVtot,I_QV,LHV)
    do k  = 1, kdim
    do ij = 1, ijdim
       rho(ij,k) = rhog(ij,k) / gsgam2(ij,k)
       tem(ij,k) = ( ein_moist(ij,k) - q(ij,k,I_QV) * LHV ) / CVtot(ij,k)
    enddo
    enddo
    !$omp end parallel do

    if ( I_QI > 0 .AND. ice_adjust ) then
       call satadjust_all( ijdim,           & ! [IN]
                           kdim,            & ! [IN]
                           rho      (:,:),  & ! [IN]
                           ein_moist(:,:),  & ! [IN]
                           qsum     (:,:),  & ! [IN]
                           tem      (:,:),  & ! [INOUT]
                           q        (:,:,:) ) ! [INOUT]
    else
       call satadjust_liq( ijdim,           & ! [IN]
                           kdim,            & ! [IN]
                           rho      (:,:),  & ! [IN]
                           ein_moist(:,:),  & ! [IN]
                           qsum     (:,:),  & ! [IN]
                           tem      (:,:),  & ! [INOUT]
                           q        (:,:,:) ) ! [INOUT]
    endif

    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,rhogq,rhog,q,I_QV,I_QC)
    do k  = 1, kdim
    do ij = 1, ijdim
       rhogq(ij,k,I_QV) = rhog(ij,k) * q(ij,k,I_QV)
       rhogq(ij,k,I_QC) = rhog(ij,k) * q(ij,k,I_QC)
    enddo
    enddo
    !$omp end parallel do

    if ( I_QI > 0 .AND. ice_adjust ) then
       !$omp parallel do default(none),private(ij,k), &
       !$omp shared(ijdim,kdim,rhogq,rhog,q,I_QI)
       do k  = 1, kdim
       do ij = 1, ijdim
          rhogq(ij,k,I_QI) = rhog(ij,k) * q(ij,k,I_QI)
       enddo
       enddo
       !$omp end parallel do
    endif

    call THRMDYN_cv( ijdim,        & ! [IN]
                     kdim,         & ! [IN]
                     qd   (:,:),   & ! [IN]
                     q    (:,:,:), & ! [IN]
                     CVtot(:,:)    ) ! [OUT]

    !$omp parallel do default(none),private(ij,k), &
    !$omp shared(ijdim,kdim,rhoge,rhog,tem,CVtot)
    do k  = 1, kdim
    do ij = 1, ijdim
       rhoge(ij,k) = rhog(ij,k) * tem(ij,k) * CVtot(ij,k)
    enddo
    enddo
    !$omp end parallel do

    call PROF_rapend  ('____Saturation_Adjustment')

    return
  end subroutine SATURATION_adjustment

  !-----------------------------------------------------------------------------
  subroutine satadjust_all( &
       ijdim,  &
       kdim,   &
       rho,    &
       Emoist, &
       qsum,   &
       tem,    &
       q       )
!ESC!    use mod_const, only: &
!ESC!       CONST_EPS,   &
!ESC!       CONST_CVdry, &
!ESC!       CONST_Rvap,  &
!ESC!       CONST_PSAT0, &
!ESC!       CONST_TEM00
!ESC!    use mod_adm, only: &
!ESC!       ADM_proc_stop,    &
!ESC!       kmin => ADM_kmin, &
!ESC!       kmax => ADM_kmax
!ESC!    use mod_runconf, only: &
!ESC!       nqmax => TRC_VMAX, &
!ESC!       NQW_STR,           &
!ESC!       NQW_END,           &
!ESC!       I_QV,              &
!ESC!       I_QC,              &
!ESC!       I_QI,              &
!ESC!       CVW,               &
!ESC!       LHV,               &
!ESC!       LHF
    use mod_thrmdyn, only: &
       thrmdyn_qd
    implicit none

    integer,  intent(in)    :: ijdim
    integer,  intent(in)    :: kdim
    real(RP), intent(in)    :: rho   (ijdim,kdim)
    real(RP), intent(in)    :: Emoist(ijdim,kdim)
    real(RP), intent(in)    :: qsum  (ijdim,kdim)
    real(RP), intent(inout) :: tem   (ijdim,kdim)
    real(RP), intent(inout) :: q     (ijdim,kdim,nqmax)

    real(RP) :: qd(ijdim,kdim)

    real(RP) :: rtem, alpha, psatl, psati, psat, qsatl, qsati, qsat
    real(RP) :: CVtot, Emoist_new, dtemp, lim1, lim2
    real(RP) :: dalpha_dT, dqsatl_dT, dqsati_dT, dqsat_dT, dqc_dT, dqi_dT, dCVtot_dT, dEmoist_dT
    real(RP) :: RTEM00, PSAT0, Rvap, CVdry, EPS

    real(RP) :: dtemp_criteria
    integer,  parameter :: itelim = 100

    logical  :: converged
    integer  :: ij, k, nq, ite
    !---------------------------------------------------------------------------

    dtemp_criteria = 10.0_RP**(-(RP_PREC+1)/2)

    EPS    = CONST_EPS
    CVdry  = CONST_CVdry
    Rvap   = CONST_Rvap
    PSAT0  = CONST_PSAT0
    RTEM00 = 1.0_RP / CONST_TEM00

    call THRMDYN_qd( ijdim, kdim, q(:,:,:), qd(:,:) )

    !$omp parallel do default(none), &
    !$omp private(ij,k,nq,ite,converged,rtem,lim1,lim2,                                    &
    !$omp         alpha,psatl,psati,psat,qsatl,qsati,qsat,CVtot,Emoist_new,dtemp,          &
    !$omp         dalpha_dT,dqsatl_dT,dqsati_dT,dqsat_dT,dqc_dT,dqi_dT,dCVtot_dT,dEmoist_dT), &
    !$omp shared (ijdim,kmin,kmax,tem,q,qd,rho,Emoist,qsum,                               &
    !$omp         CPovR_liq,CVovR_liq,LovR_liq,CPovR_ice,CVovR_ice,LovR_ice,              &
    !$omp         NQW_STR,NQW_END,I_QV,I_QC,I_QI,CVW,LHV,LHF,RTEM00,PSAT0,Rvap,CVdry,EPS, &
    !$omp         SATURATION_ULIMIT_TEMP,SATURATION_LLIMIT_TEMP,dtemp_criteria)
    do k  = kmin, kmax
    do ij = 1, ijdim

       alpha = ( tem(ij,k)              - SATURATION_LLIMIT_TEMP ) &
             / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )
       alpha = min( max( alpha, 0.0_RP ), 1.0_RP )

       rtem = 1.0_RP / max( tem(ij,k), TEM_MIN )

       psatl = PSAT0 * ( tem(ij,k) * RTEM00 )**CPovR_liq * exp( LovR_liq * ( RTEM00 - rtem ) )
       psati = PSAT0 * ( tem(ij,k) * RTEM00 )**CPovR_ice * exp( LovR_ice * ( RTEM00 - rtem ) )
       psat  = psatl * (          alpha ) &
             + psati * ( 1.0_RP - alpha )

       qsat  = psat  / ( rho(ij,k) * Rvap * tem(ij,k) )

       if ( qsum(ij,k)-qsat > EPS ) then
          converged = .false.

          do ite = 1, itelim

             alpha = ( tem(ij,k)              - SATURATION_LLIMIT_TEMP ) &
                   / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )
             alpha = min( max( alpha, 0.0_RP ), 1.0_RP )

             rtem = 1.0_RP / max( tem(ij,k), TEM_MIN )

             psatl = PSAT0 * ( tem(ij,k) * RTEM00 )**CPovR_liq * exp( LovR_liq * ( RTEM00 - rtem ) )
             psati = PSAT0 * ( tem(ij,k) * RTEM00 )**CPovR_ice * exp( LovR_ice * ( RTEM00 - rtem ) )
             psat  = psatl * (          alpha ) &
                   + psati * ( 1.0_RP - alpha )

             qsatl = psatl / ( rho(ij,k) * Rvap * tem(ij,k) )
             qsati = psati / ( rho(ij,k) * Rvap * tem(ij,k) )
             qsat  = psat  / ( rho(ij,k) * Rvap * tem(ij,k) )

             ! Separation
             q(ij,k,I_QV) = qsat
             q(ij,k,I_QC) = ( qsum(ij,k)-qsat ) * (        alpha )
             q(ij,k,I_QI) = ( qsum(ij,k)-qsat ) * ( 1.0_RP-alpha )

             CVtot = qd(ij,k) * CVdry
             do nq = NQW_STR, NQW_END
                CVtot = CVtot + q(ij,k,nq) * CVW(nq)
             enddo

             Emoist_new = tem(ij,k) * CVtot + q(ij,k,I_QV) * LHV - q(ij,k,I_QI) * LHF

             ! dX/dT
             lim1       = 0.5_RP + sign( 0.5_RP, SATURATION_ULIMIT_TEMP - tem(ij,k) )
             lim2       = 0.5_RP + sign( 0.5_RP, tem(ij,k) - SATURATION_LLIMIT_TEMP )
             dalpha_dT  = lim1 * lim2 / ( SATURATION_ULIMIT_TEMP - SATURATION_LLIMIT_TEMP )

             dqsatl_dT  = ( LovR_liq / tem(ij,k)**2 + CVovR_liq / tem(ij,k) ) * qsatl
             dqsati_dT  = ( LovR_ice / tem(ij,k)**2 + CVovR_ice / tem(ij,k) ) * qsati
             dqsat_dT   = qsatl * dalpha_dT + dqsatl_dT * (        alpha ) &
                        - qsati * dalpha_dT + dqsati_dT * ( 1.0_RP-alpha )

             dqc_dT     =  ( qsum(ij,k)-qsat ) * dalpha_dT - dqsat_dT * (        alpha )
             dqi_dT     = -( qsum(ij,k)-qsat ) * dalpha_dT - dqsat_dT * ( 1.0_RP-alpha )

             dCVtot_dT  = dqsat_dT * CVW(I_QV) &
                        + dqc_dT   * CVW(I_QC) &
                        + dqi_dT   * CVW(I_QI)

             dEmoist_dT = tem(ij,k) * dCVtot_dT &
                        + CVtot                 &
                        + dqsat_dT  * LHV       &
                        - dqi_dT    * LHF

             dtemp      = ( Emoist_new - Emoist(ij,k) ) / dEmoist_dT

             tem(ij,k)  = tem(ij,k) - dtemp

             if ( abs(dtemp) < dtemp_criteria ) then
                converged = .true.
                exit
             endif

             if( tem(ij,k)*0.0_RP /= 0.0_RP) exit
          enddo

          if ( .NOT. converged ) then
             write(*,*) rho(ij,k),tem(ij,k),q(ij,k,I_QV),q(ij,k,I_QC),q(ij,k,I_QI)
             write(*,*) 'xxx [satadjust_all] not converged! dtemp=', dtemp,ij,k,ite
             call ADM_proc_stop
          endif

       endif

    enddo
    enddo
    !$omp end parallel do

    return
  end subroutine satadjust_all

  !-----------------------------------------------------------------------------
  subroutine satadjust_liq( &
       ijdim,  &
       kdim,   &
       rho,    &
       Emoist, &
       qsum,   &
       tem,    &
       q       )
!ESC!    use mod_const, only: &
!ESC!       CONST_EPS,   &
!ESC!       CONST_CVdry, &
!ESC!       CONST_Rvap,  &
!ESC!       CONST_PSAT0, &
!ESC!       CONST_TEM00
!ESC!    use mod_adm, only: &
!ESC!       ADM_proc_stop,    &
!ESC!       kmin => ADM_kmin, &
!ESC!       kmax => ADM_kmax
!ESC!    use mod_runconf, only: &
!ESC!       nqmax => TRC_VMAX, &
!ESC!       NQW_STR,           &
!ESC!       NQW_END,           &
!ESC!       I_QV,              &
!ESC!       I_QC,              &
!ESC!       CVW,               &
!ESC!       LHV
    use mod_thrmdyn, only: &
       thrmdyn_qd
    implicit none

    integer,  intent(in)    :: ijdim
    integer,  intent(in)    :: kdim
    real(RP), intent(in)    :: rho   (ijdim,kdim)
    real(RP), intent(in)    :: Emoist(ijdim,kdim)
    real(RP), intent(in)    :: qsum  (ijdim,kdim)
    real(RP), intent(inout) :: tem   (ijdim,kdim)
    real(RP), intent(inout) :: q     (ijdim,kdim,nqmax)

    real(RP) :: qd(ijdim,kdim)

    real(RP) :: rtem, psat, qsat
    real(RP) :: CVtot, Emoist_new, dtemp
    real(RP) :: dqsat_dT, dCVtot_dT, dEmoist_dT
    real(RP) :: RTEM00, PSAT0, Rvap, CVdry, EPS

    real(RP) :: dtemp_criteria
    integer,  parameter :: itelim = 100

    logical  :: converged
    integer  :: ij, k, nq, ite
    !---------------------------------------------------------------------------

    dtemp_criteria = 10.0_RP**(-(RP_PREC+1)/2)

    EPS    = CONST_EPS
    CVdry  = CONST_CVdry
    Rvap   = CONST_Rvap
    PSAT0  = CONST_PSAT0
    RTEM00 = 1.0_RP / CONST_TEM00

    call THRMDYN_qd( ijdim, kdim, q(:,:,:), qd(:,:) )

    !$omp parallel do default(none), &
    !$omp private(ij,k,nq,ite,converged,rtem,                                   &
    !$omp         psat,qsat,CVtot,Emoist_new,dtemp,dqsat_dT,dCVtot_dT,dEmoist_dT), &
    !$omp shared (ijdim,kmin,kmax,tem,q,qd,rho,Emoist,qsum,                      &
    !$omp         CPovR_liq,CVovR_liq,LovR_liq,                                  &
    !$omp         NQW_STR,NQW_END,I_QV,I_QC,CVW,LHV,RTEM00,PSAT0,Rvap,CVdry,EPS, &
    !$omp         SATURATION_ULIMIT_TEMP,SATURATION_LLIMIT_TEMP,dtemp_criteria)
    do k  = kmin, kmax
    do ij = 1, ijdim

       rtem = 1.0_RP / max( tem(ij,k), TEM_MIN )

       psat = PSAT0 * ( tem(ij,k) * RTEM00 )**CPovR_liq * exp( LovR_liq * ( RTEM00 - rtem ) )

       qsat = psat / ( rho(ij,k) * Rvap * tem(ij,k) )

       if ( qsum(ij,k)-qsat > EPS ) then
          converged = .false.

          do ite = 1, itelim

             rtem = 1.0_RP / max( tem(ij,k), TEM_MIN )

             psat = PSAT0 * ( tem(ij,k) * RTEM00 )**CPovR_liq * exp( LovR_liq * ( RTEM00 - rtem ) )

             qsat = psat / ( rho(ij,k) * Rvap * tem(ij,k) )

             ! Separation
             q(ij,k,I_QV) = qsat
             q(ij,k,I_QC) = qsum(ij,k)-qsat

             CVtot = qd(ij,k) * CVdry
             do nq = NQW_STR, NQW_END
                CVtot = CVtot + q(ij,k,nq) * CVW(nq)
             enddo

             Emoist_new = tem(ij,k) * CVtot + q(ij,k,I_QV) * LHV

             ! dX/dT
             dqsat_dT   = ( LovR_liq / tem(ij,k)**2 + CVovR_liq / tem(ij,k) ) * qsat

             dCVtot_dT  = dqsat_dT * ( CVW(I_QV) - CVW(I_QC) )

             dEmoist_dT = tem(ij,k) * dCVtot_dT &
                        + CVtot                 &
                        + dqsat_dT  * LHV

             dtemp      = ( Emoist_new - Emoist(ij,k) ) / dEmoist_dT

             tem(ij,k) = tem(ij,k) - dtemp

             if ( abs(dtemp) < dtemp_criteria ) then
                converged = .true.
                exit
             endif

             if( tem(ij,k)*0.0_RP /= 0.0_RP) exit
          enddo

          if ( .NOT. converged ) then
             write(*,*) rho(ij,k),tem(ij,k),q(ij,k,I_QV),q(ij,k,I_QC)
             write(*,*) 'xxx [satadjust_liq] not converged! dtemp=', dtemp,ij,k,ite
             call ADM_proc_stop
          endif

       endif

    enddo
    enddo

    return
  end subroutine satadjust_liq

  !-----------------------------------------------------------------------------
  subroutine SATURATION_enthalpy( &
       ijdim,       &
       kdim,        &
       tem,         &
       ent,         &
       pre,         &
       qw,          &
       adiabat_type )
!ESC!    use mod_const, only: &
!ESC!       CONST_Rdry,   &
!ESC!       CONST_Rvap,   &
!ESC!       CONST_CPdry,  &
!ESC!       CONST_CPvap,  &
!ESC!       CONST_LHV0,   &
!ESC!       CONST_PSAT0,  &
!ESC!       CONST_EPSvap, &
!ESC!       CONST_PRE00,  &
!ESC!       CONST_TEM00
!ESC!    use mod_runconf, only: &
!ESC!       nqmax => TRC_VMAX, &
!ESC!       I_QV,              &
!ESC!       I_QC
    use mod_thrmdyn, only: &
       thrmdyn_ent
    implicit none

    integer,          intent(in)    :: ijdim
    integer,          intent(in)    :: kdim
    real(RP),         intent(inout) :: tem(ijdim,kdim)
    real(RP),         intent(in)    :: ent(ijdim,kdim)
    real(RP),         intent(in)    :: pre(ijdim,kdim)
    real(RP),         intent(in)    :: qw (ijdim,kdim)
    character(len=*), intent(in)    :: adiabat_type
    !                                = 'RMA' reversible: qv = qv*, qc = qw - qv*
    !                                = 'PMA' pseudo    : if qw > qv* : qv = qv*, qc=0
    !                                = 'DA'  qc = 0

    real(RP) :: qd     (ijdim,kdim)
    real(RP) :: qsat   (ijdim,kdim)
    real(RP) :: q      (ijdim,kdim,nqmax)
    real(RP) :: pred   (ijdim,kdim)
    real(RP) :: prev   (ijdim,kdim)
    real(RP) :: tem_old(ijdim,kdim)
    real(RP) :: ents   (ijdim,kdim)
    real(RP) :: dents  (ijdim,kdim)
    real(RP) :: dtem   (ijdim,kdim)
    logical  :: flag   (ijdim,kdim)

    real(RP) :: PREMIN     = 1.E-10_RP
    real(RP) :: DENTS_FACT = 1.0_RP
    real(RP) :: DENTS_MIN  = 1.E-10_RP
    real(RP) :: TEMMIN     = 1.0_RP

    integer, parameter :: itelim = 100

    integer  :: ite
    integer  :: ij, k
    !---------------------------------------------------------------------------

    flag(:,:) = .false.

    q(:,:,:) = 0.0_RP

    if ( adiabat_type == 'DA' ) then

       tem (:,:)      = exp( ( ent(:,:) + CONST_Rdry * log( pre(:,:) / CONST_PRE00 ) ) / CONST_CPdry ) * CONST_TEM00
       qd  (:,:)      = 1.0_RP-qw(:,:)
       q   (:,:,I_QV) =        qw(:,:)
       pred(:,:)      = pre(:,:) * CONST_EPSvap * qd(:,:) / ( CONST_EPSvap * qd(:,:) + q(:,:,I_QV) )
       pred(:,:)      = max( pred(:,:), PREMIN )
       prev(:,:)      = pre(:,:) * q(:,:,I_QV)            / ( CONST_EPSvap * qd(:,:) + q(:,:,I_QV) )
       prev(:,:)      = max( prev(:,:), PREMIN )

       tem(:,:) = exp( ( ent(:,:) &
                       + qd(:,:)      * ( CONST_Rdry * log( pred(:,:) / CONST_PRE00 )                              )&
                       + q (:,:,I_QV) * ( CONST_Rvap * log( prev(:,:) / CONST_PSAT0 ) - CONST_LHV0 / CONST_TEM00 ) ) &
                     / ( qd(:,:) * CONST_CPdry + q(:,:,I_QV) * CONST_CPvap ) ) * CONST_TEM00

    elseif( adiabat_type == 'RMA' .OR. adiabat_type == 'PMA' ) then

       tem(:,:) = exp( ( ent(:,:) + CONST_Rdry * log( pre(:,:) / CONST_PRE00 ) ) / CONST_CPdry ) * CONST_TEM00

       do ite = 1, itelim

          tem_old(:,:) = tem(:,:)

          call SATURATION_qsat_liq( ijdim, kdim, tem(:,:), pre(:,:), qsat(:,:) )

          q(:,:,I_QV) = min ( qw(:,:), qsat(:,:) )

          if    ( adiabat_type == 'RMA' ) then
             q (:,:,I_QC) = qw(:,:) - q(:,:,I_QV)
             qd(:,:)      = 1.0_RP - qw(:,:)
          elseif( adiabat_type == 'PMA' ) then
             q (:,:,I_QC) = 0.0_RP
             qd(:,:)      = 1.0_RP - q(:,:,I_QV)
          endif

          call thrmdyn_ent( ijdim, & ! [IN]
                            kdim,  & ! [IN]
                            ents,  & ! [OUT]
                            tem,   & ! [IN]
                            pre,   & ! [IN]
                            q,     & ! [IN]
                            qd     ) ! [IN]

          do k  = 1, kdim
          do ij = 1, ijdim
             dents(:,:) = ( CONST_CPdry * qd(:,:) + CONST_CPvap * q(:,:,I_QV) ) / tem(:,:)    &
                        + (                   CONST_LHV0  * q(:,:,I_QV) ) / tem(:,:)**2 &
                        * ( CONST_LHV0 / ( tem(:,:) * CONST_RVAP ) - 1.0_RP )
             dents(:,:) = max ( dents(:,:) * DENTS_FACT, DENTS_MIN )
          enddo
          enddo

          do k  = 1, kdim
          do ij = 1, ijdim
             if ( .NOT. flag(ij,k) ) then
                tem (ij,k) = tem(ij,k) - ( ents(ij,k) - ent(ij,k) ) / dents(ij,k)
                tem (ij,k) = max ( tem(ij,k), TEMMIN )
                dtem(ij,k) = abs( tem(ij,k) - tem_old(ij,k) )
             endif

             if ( dtem(ij,k) < DTEM_EPS0 ) then
                dtem(ij,k) = 0.0_RP
                flag(ij,k) = .true.
             endif
          enddo
          enddo

          if( maxval( dtem(:,:) ) < DTEM_EPS0 ) exit

       enddo ! iteration
    else
       write(IO_FID_LOG,*) '### SATURATION_enthalpy: invalid adiabat_type=', adiabat_type
    endif

    return
  end subroutine SATURATION_enthalpy

  !-----------------------------------------------------------------------------
  ! (d qsw/d T)_{rho}: partial difference of qsat_water
  subroutine moist_dqsw_dtem_rho( &
       ijdim,  &
       kdim,   &
       tem,    &
       rho,    &
       dqsdtem )
!ESC!    use mod_const, only: &
!ESC!       CONST_Rvap, &
!ESC!       CONST_TEM00
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(RP), intent(in)  :: tem    (ijdim,kdim)
    real(RP), intent(in)  :: rho    (ijdim,kdim)
    real(RP), intent(out) :: dqsdtem(ijdim,kdim)

    real(RP) :: psat(ijdim,kdim) ! saturation vapor pressure
    real(RP) :: LovR(ijdim,kdim) ! latent heat for condensation
    !---------------------------------------------------------------------------

    call SATURATION_psat_liq_2D( ijdim, kdim, tem(:,:), psat(:,:) )

    LovR   (:,:) = LovR_liq + CPovR_liq * ( tem(:,:) - CONST_TEM00 )
    dqsdtem(:,:) = psat(:,:) / ( rho (:,:) * CONST_Rvap * tem(:,:)**2 ) &
                             * ( LovR(:,:) / tem(:,:) - 1.0_RP )

    return
  end subroutine moist_dqsw_dtem_rho

  !-----------------------------------------------------------------------------
  ! (d qsi/d T)_{rho}: partial difference of qsat_ice
  subroutine moist_dqsi_dtem_rho( &
       ijdim,  &
       kdim,   &
       tem,    &
       rho,    &
       dqsdtem )
!ESC!    use mod_const, only: &
!ESC!       CONST_Rvap, &
!ESC!       CONST_TEM00
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(RP), intent(in)  :: tem    (ijdim,kdim)
    real(RP), intent(in)  :: rho    (ijdim,kdim)
    real(RP), intent(out) :: dqsdtem(ijdim,kdim)

    real(RP) :: psat(ijdim,kdim) ! saturation vapor pressure
    real(RP) :: LovR(ijdim,kdim) ! latent heat for condensation
    !---------------------------------------------------------------------------

    call SATURATION_psat_ice_2D( ijdim, kdim, tem(:,:), psat(:,:) )

    LovR   (:,:) = LovR_ice + CPovR_ice * ( tem(:,:) - CONST_TEM00 )
    dqsdtem(:,:) = psat(:,:) / ( rho (:,:) * CONST_Rvap * tem(:,:)**2 ) &
                             * ( LovR(:,:) / tem(:,:) - 1.0_RP )

    return
  end subroutine moist_dqsi_dtem_rho

  !-----------------------------------------------------------------------------
  ! (d qs/d T)_{p} and (d qs/d p)_{T}
  subroutine moist_dqsw_dtem_dpre( &
       ijdim,   &
       kdim,    &
       tem,     &
       pre,     &
       dqsdtem, &
       dqsdpre  )
!ESC!    use mod_const, only: &
!ESC!       CONST_EPSvap, &
!ESC!       CONST_TEM00
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(RP), intent(in)  :: tem    (ijdim,kdim)
    real(RP), intent(in)  :: pre    (ijdim,kdim)
    real(RP), intent(out) :: dqsdtem(ijdim,kdim)
    real(RP), intent(out) :: dqsdpre(ijdim,kdim)

    real(RP) :: psat(ijdim,kdim) ! saturation vapor pressure
    real(RP) :: LovR(ijdim,kdim) ! latent heat for condensation
    real(RP) :: den1(ijdim,kdim) ! denominator
    !---------------------------------------------------------------------------

    call SATURATION_psat_liq_2D( ijdim, kdim, tem(:,:), psat(:,:) )

    LovR   (:,:) = LovR_liq + CPovR_liq * ( tem(:,:) - CONST_TEM00 )
    den1   (:,:) = ( pre(:,:) - ( 1.0_RP-CONST_EPSvap ) * psat(:,:) )**2
    dqsdpre(:,:) = -CONST_EPSvap * psat(:,:) / den1(:,:)
    dqsdtem(:,:) =  CONST_EPSvap * psat(:,:) / den1(:,:) * pre(:,:) * LovR(:,:) / tem(:,:)**2

    return
  end subroutine moist_dqsw_dtem_dpre

  !-----------------------------------------------------------------------------
  ! (d qsi/d T)_{p} and (d qs/d p)_{T}
  subroutine moist_dqsi_dtem_dpre( &
       ijdim,   &
       kdim,    &
       tem,     &
       pre,     &
       dqsdtem, &
       dqsdpre  )
!ESC!    use mod_const, only: &
!ESC!       CONST_EPSvap, &
!ESC!       CONST_TEM00
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(RP), intent(in)  :: tem    (ijdim,kdim)
    real(RP), intent(in)  :: pre    (ijdim,kdim)
    real(RP), intent(out) :: dqsdtem(ijdim,kdim)
    real(RP), intent(out) :: dqsdpre(ijdim,kdim)

    real(RP) :: psat(ijdim,kdim) ! saturation vapor pressure
    real(RP) :: LovR(ijdim,kdim) ! latent heat for condensation
    real(RP) :: den1(ijdim,kdim) ! denominator
    !---------------------------------------------------------------------------

    call SATURATION_psat_ice_2D( ijdim, kdim, tem(:,:), psat(:,:) )

    LovR   (:,:) = LovR_ice + CPovR_ice * ( tem(:,:) - CONST_TEM00 )
    den1   (:,:) = ( pre(:,:) - ( 1.0_RP-CONST_EPSvap ) * psat(:,:) )**2
    dqsdpre(:,:) = -CONST_EPSvap * psat(:,:) / den1(:,:)
    dqsdtem(:,:) =  CONST_EPSvap * psat(:,:) / den1(:,:) * pre(:,:) * LovR(:,:) / tem(:,:)**2

    return
  end subroutine moist_dqsi_dtem_dpre

end module mod_satadjust
