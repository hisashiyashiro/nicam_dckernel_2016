!-------------------------------------------------------------------------------
module mod_misc
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public, including parameters
  !
  include 'problem_size.inc'

  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  public :: DEBUG_rapstart
  public :: DEBUG_rapend
  public :: DEBUG_rapreport

  public :: MISC_make_idstr        !--- make file name with a number
  public :: MISC_get_available_fid !--- get an available file ID

  public :: ADM_proc_stop

  public :: GRD_setup
  public :: GRD_input_vgrid

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer,             public :: EX_CSTEP_diffusion = 0
  integer,             public :: EX_TSTEP_diffusion = 60
  integer,             public :: EX_CSTEP_divdamp3d = 0
  integer,             public :: EX_TSTEP_divdamp3d = 140
  integer,             public :: EX_fid
  integer,             public :: EX_err
  character(len=1024), public :: EX_fname
  character(len=16),   public :: EX_item
  real(RP),            public :: EX_max
  real(RP),            public :: EX_min
  real(RP),            public :: EX_sum

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: DEBUG_rapid

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,                 private, parameter :: DEBUG_rapnlimit = 100
  integer,                 private            :: DEBUG_rapnmax   = 0
  character(len=ADM_NSYS), private            :: DEBUG_rapname(DEBUG_rapnlimit)
  real(8),                 private            :: DEBUG_raptstr(DEBUG_rapnlimit)
  real(8),                 private            :: DEBUG_rapttot(DEBUG_rapnlimit)
  integer,                 private            :: DEBUG_rapnstr(DEBUG_rapnlimit)
  integer,                 private            :: DEBUG_rapnend(DEBUG_rapnlimit)

#ifdef _FIXEDINDEX_
  real(DP), public              :: GRD_x    (ADM_gall   ,ADM_KNONE,ADM_lall   ,              ADM_nxyz)
  real(DP), public              :: GRD_x_pl (ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              ADM_nxyz)
  real(DP), public              :: GRD_xt   (ADM_gall   ,ADM_KNONE,ADM_lall   ,ADM_TI:ADM_TJ,ADM_nxyz)
  real(DP), public              :: GRD_xt_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              ADM_nxyz)
  real(DP), public              :: GRD_xr   (ADM_gall   ,ADM_KNONE,ADM_lall   ,ADM_AI:ADM_AJ,ADM_nxyz)
  real(DP), public              :: GRD_xr_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              ADM_nxyz)
#else
  real(DP), public, allocatable :: GRD_x    (:,:,:,:)
  real(DP), public, allocatable :: GRD_x_pl (:,:,:,:)
  real(DP), public, allocatable :: GRD_xt   (:,:,:,:,:)
  real(DP), public, allocatable :: GRD_xt_pl(:,:,:,:)
  real(DP), public, allocatable :: GRD_xr   (:,:,:,:,:)
  real(DP), public, allocatable :: GRD_xr_pl(:,:,:,:)
#endif

#ifdef _FIXEDINDEX_
  real(DP), public              :: GRD_gz   (ADM_kall)
  real(DP), public              :: GRD_gzh  (ADM_kall)
  real(DP), public              :: GRD_dgz  (ADM_kall)
  real(DP), public              :: GRD_dgzh (ADM_kall)
  real(DP), public              :: GRD_rdgz (ADM_kall)
  real(DP), public              :: GRD_rdgzh(ADM_kall)

  real(DP), public              :: GRD_afac(ADM_kall)
  real(DP), public              :: GRD_bfac(ADM_kall)
  real(DP), public              :: GRD_cfac(ADM_kall)
  real(DP), public              :: GRD_dfac(ADM_kall)
#else
  real(DP), public, allocatable :: GRD_gz   (:) ! gsi (z-star) coordinate
  real(DP), public, allocatable :: GRD_gzh  (:) ! gsi (z-star) coordinate at the half point
  real(DP), public, allocatable :: GRD_dgz  (:) ! d(gsi)
  real(DP), public, allocatable :: GRD_dgzh (:) ! d(gsi) at the half point
  real(DP), public, allocatable :: GRD_rdgz (:)
  real(DP), public, allocatable :: GRD_rdgzh(:)

  real(DP), public, allocatable :: GRD_afac (:) ! From the cell center value to the cell wall value
  real(DP), public, allocatable :: GRD_bfac (:) !    A(k-1/2) = ( afac(k) A(k) + bfac(k) * A(k-1) ) / 2
  real(DP), public, allocatable :: GRD_cfac (:) ! From the cell wall value to the cell center value
  real(DP), public, allocatable :: GRD_dfac (:) !    A(k) = ( cfac(k) A(k+1/2) + dfac(k) * A(k-1/2) ) / 2
#endif

#ifdef _FIXEDINDEX_
  real(DP), public              :: GMTR_P_var   (ADM_gall   ,ADM_KNONE,ADM_lall   ,              GMTR_P_nmax_var   )
  real(DP), public              :: GMTR_P_var_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              GMTR_P_nmax_var   )
  real(DP), public              :: GMTR_T_var   (ADM_gall   ,ADM_KNONE,ADM_lall   ,ADM_TI:ADM_TJ,GMTR_T_nmax_var   )
  real(DP), public              :: GMTR_T_var_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              GMTR_T_nmax_var   )
  real(DP), public              :: GMTR_A_var   (ADM_gall   ,ADM_KNONE,ADM_lall   ,ADM_AI:ADM_AJ,GMTR_A_nmax_var   )
  real(DP), public              :: GMTR_A_var_pl(ADM_gall_pl,ADM_KNONE,ADM_lall_pl,              GMTR_A_nmax_var_pl)

  real(DP), public              :: GMTR_area    (ADM_gall   ,ADM_lall   )
  real(DP), public              :: GMTR_area_pl (ADM_gall_pl,ADM_lall_pl)
  real(DP), public              :: GMTR_lat     (ADM_gall   ,ADM_lall   )
  real(DP), public              :: GMTR_lat_pl  (ADM_gall_pl,ADM_lall_pl)
  real(DP), public              :: GMTR_lon     (ADM_gall   ,ADM_lall   )
  real(DP), public              :: GMTR_lon_pl  (ADM_gall_pl,ADM_lall_pl)
#else
  real(DP), public, allocatable :: GMTR_P_var   (:,:,:,:)   ! geometrics for the cell point
  real(DP), public, allocatable :: GMTR_P_var_pl(:,:,:,:)
  real(DP), public, allocatable :: GMTR_T_var   (:,:,:,:,:) ! geometrics for the cell vertex
  real(DP), public, allocatable :: GMTR_T_var_pl(:,:,:,:)
  real(DP), public, allocatable :: GMTR_A_var   (:,:,:,:,:) ! geometrics for the cell arc
  real(DP), public, allocatable :: GMTR_A_var_pl(:,:,:,:)

  real(DP), public, allocatable :: GMTR_area    (:,:)       ! control area of the cell
  real(DP), public, allocatable :: GMTR_area_pl (:,:)
  real(DP), public, allocatable :: GMTR_lat     (:,:)       ! latitude  of the cell point
  real(DP), public, allocatable :: GMTR_lat_pl  (:,:)
  real(DP), public, allocatable :: GMTR_lon     (:,:)       ! longitude of the cell point
  real(DP), public, allocatable :: GMTR_lon_pl  (:,:)
#endif

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  function DEBUG_rapid( rapname ) result(id)
    implicit none

    character(len=*), intent(in) :: rapname

    integer :: id
    !---------------------------------------------------------------------------

    if ( DEBUG_rapnmax >= 1 ) then
       do id = 1, DEBUG_rapnmax
          if( trim(rapname) == trim(DEBUG_rapname(id)) ) return
       enddo
    endif

    DEBUG_rapnmax     = DEBUG_rapnmax + 1
    id                = DEBUG_rapnmax
    DEBUG_rapname(id) = trim(rapname)
    DEBUG_raptstr(id) = 0.D0
    DEBUG_rapttot(id) = 0.D0
    DEBUG_rapnstr(id) = 0
    DEBUG_rapnend(id) = 0

  end function DEBUG_rapid

  !-----------------------------------------------------------------------------
  subroutine DEBUG_rapstart( rapname )
    implicit none

    character(len=*), intent(in) :: rapname

    real(8) :: time
    integer :: id
    !---------------------------------------------------------------------------

    id = DEBUG_rapid( rapname )

    call CPU_TIME(time)

    DEBUG_raptstr(id) = time
    DEBUG_rapnstr(id) = DEBUG_rapnstr(id) + 1

    !write(ADM_LOG_FID,*) rapname, DEBUG_rapnstr(id)

#ifdef _FAPP_
    call fapp_start( rapname, id, 1 )
#endif

    return
  end subroutine DEBUG_rapstart

  !-----------------------------------------------------------------------------
  subroutine DEBUG_rapend( rapname )
    implicit none

    character(len=*), intent(in) :: rapname

    real(8) :: time
    integer :: id
    !---------------------------------------------------------------------------

    id = DEBUG_rapid( rapname )

    call CPU_TIME(time)

    DEBUG_rapttot(id) = DEBUG_rapttot(id) + ( time-DEBUG_raptstr(id) )
    DEBUG_rapnend(id) = DEBUG_rapnend(id) + 1

#ifdef _FAPP_
    call fapp_stop( rapname, id, 1 )
#endif

    return
  end subroutine DEBUG_rapend

  !-----------------------------------------------------------------------------
  subroutine DEBUG_rapreport
    implicit none

    integer :: id
    !---------------------------------------------------------------------------

    if ( DEBUG_rapnmax >= 1 ) then

       do id = 1, DEBUG_rapnmax
          if ( DEBUG_rapnstr(id) /= DEBUG_rapnend(id) ) then
              write(*,*) '*** Mismatch Report',id,DEBUG_rapname(id),DEBUG_rapnstr(id),DEBUG_rapnend(id)
          endif
       enddo

       write(ADM_LOG_FID,*)
       write(ADM_LOG_FID,*) '*** Computational Time Report'

       do id = 1, DEBUG_rapnmax
          write(ADM_LOG_FID,'(1x,A,I3.3,A,A,A,F10.3,A,I7)') &
          '*** ID=',id,' : ',DEBUG_rapname(id),' T=',DEBUG_rapttot(id),' N=',DEBUG_rapnstr(id)
       enddo

    else
       write(ADM_LOG_FID,*)
       write(ADM_LOG_FID,*) '*** Computational Time Report: NO item.'
    endif

    return
  end subroutine DEBUG_rapreport
  !-----------------------------------------------------------------------------
  !> make extention with process number
  subroutine MISC_make_idstr( &
       str,    &
       prefix, &
       ext,    &
       numID,  &
       digit   )
    implicit none

    character(len=*),  intent(out) :: str    !< combined extention string
    character(len=*),  intent(in)  :: prefix !< prefix
    character(len=*),  intent(in)  :: ext    !< extention ( e.g. .rgn )
    integer,           intent(in)  :: numID  !< number
    integer, optional, intent(in)  :: digit  !< digit

    logical, parameter            :: NSTR_ZERO_START = .true. ! number of separated file starts from 0 ?
    integer, parameter            :: NSTR_MAX_DIGIT  = 6      ! digit of separated file

    character(len=128) :: rankstr
    integer            :: setdigit
    !---------------------------------------------------------------------------

    if ( NSTR_ZERO_START ) then
       write(rankstr,'(I128.128)') numID-1
    else
       write(rankstr,'(I128.128)') numID
    endif

    if ( present(digit) ) then
       setdigit = digit
    else
       setdigit = NSTR_MAX_DIGIT
    endif

    rankstr(1:setdigit) = rankstr(128-(setdigit-1):128)
    rankstr(setdigit+1:128) = ' '

    str = trim(prefix)//'.'//trim(ext)//trim(rankstr) ! -> prefix.ext00000

    return
  end subroutine MISC_make_idstr

  !-----------------------------------------------------------------------------
  !> Search and get available machine id
  !> @return fid
  function MISC_get_available_fid() result(fid)
    implicit none

    integer :: fid

    integer, parameter :: min_fid =  7 !< minimum available fid
    integer, parameter :: max_fid = 99 !< maximum available fid

    logical :: i_opened
    !---------------------------------------------------------------------------

    do fid = min_fid, max_fid
       inquire(fid,opened=i_opened)
       if( .NOT. i_opened ) return
    enddo

  end function MISC_get_available_fid

  !-----------------------------------------------------------------------------
  subroutine ADM_proc_stop

    stop

  end subroutine ADM_proc_stop

  !-----------------------------------------------------------------------------
  subroutine GRD_setup
    implicit none

    integer :: k
    !---------------------------------------------------------------------------

    !--- < setting the vertical coordinate > ---
    allocate( GRD_gz   (ADM_kall) )
    allocate( GRD_gzh  (ADM_kall) )
    allocate( GRD_dgz  (ADM_kall) )
    allocate( GRD_dgzh (ADM_kall) )
    allocate( GRD_rdgz (ADM_kall) )
    allocate( GRD_rdgzh(ADM_kall) )

    call GRD_input_vgrid(vgrid_fname)

    ! calculation of grid intervals ( cell center )
    do k = ADM_kmin-1, ADM_kmax
       GRD_dgz(k) = GRD_gzh(k+1) - GRD_gzh(k)
    enddo
    GRD_dgz(ADM_kmax+1) = GRD_dgz(ADM_kmax)

    ! calculation of grid intervals ( cell wall )
    do k = ADM_kmin, ADM_kmax+1
       GRD_dgzh(k) = GRD_gz(k) - GRD_gz(k-1)
    enddo
    GRD_dgzh(ADM_kmin-1) = GRD_dgzh(ADM_kmin)

    ! calculation of 1/dgz and 1/dgzh
    do k = 1, ADM_kall
       GRD_rdgz (k) = 1.D0 / grd_dgz(k)
       GRD_rdgzh(k) = 1.D0 / grd_dgzh(k)
    enddo

    !---< vertical interpolation factor >---
    allocate( GRD_afac (ADM_kall) )
    allocate( GRD_bfac (ADM_kall) )
    allocate( GRD_cfac (ADM_kall) )
    allocate( GRD_dfac (ADM_kall) )

    ! From the cell center value to the cell wall value
    ! A(k-1/2) = ( afac(k) A(k) + bfac(k) * A(k-1) ) / 2
    do k = ADM_kmin, ADM_kmax+1
       GRD_afac(k) = 2.D0 * ( GRD_gzh(k) - GRD_gz(k-1) ) / ( GRD_gz(k) - GRD_gz(k-1) )
    enddo
    GRD_afac(ADM_kmin-1) = 2.D0

    GRD_bfac(:) = 2.D0 - GRD_afac(:)

    ! From the cell wall value to the cell center value
    ! A(k) = ( cfac(k) A(k+1/2) + dfac(k) * A(k-1/2) ) / 2
    do k = ADM_kmin, ADM_kmax
       GRD_cfac(k) = 2.D0 * ( GRD_gz(k) - GRD_gzh(k) ) / ( GRD_gzh(k+1) - GRD_gzh(k) )
    enddo
    GRD_cfac(ADM_kmin-1) = 2.D0
    GRD_cfac(ADM_kmax+1) = 0.D0

    GRD_dfac(:) = 2.D0 - GRD_cfac(:)

    return
  end subroutine GRD_setup

  !-----------------------------------------------------------------------------
  subroutine GRD_input_vgrid( &
       fname  )
!ESC!    use mod_misc,  only :&
!ESC!       MISC_get_available_fid
!ESC!    use mod_adm, only: &
!ESC!       ADM_proc_stop, &
!ESC!       ADM_vlayer
    implicit none

    character(len=*), intent(in) :: fname

    integer :: num_of_layer
    integer :: fid, ierr
    !---------------------------------------------------------------------------

    fid = MISC_get_available_fid()
    open( unit   = fid,           &
          file   = trim(fname),   &
          status = 'old',         &
          form   = 'unformatted', &
          access = 'sequential',  &
          iostat = ierr           )

       if ( ierr /= 0 ) then
          write(ADM_LOG_FID,*) 'xxx No vertical grid file : ', trim(fname)
          call ADM_proc_stop
       endif

       read(fid) num_of_layer
       if ( num_of_layer /= ADM_vlayer ) then
          write(ADM_LOG_FID,*)          'Msg : Sub[GRD_input_vgrid]/Mod[grid]'
          write(ADM_LOG_FID,*)          '   *** inconsistency in number of vertical layers.'
          call ADM_proc_stop
       endif

       read(fid) GRD_gz
       read(fid) GRD_gzh

       close(fid)

    return
  end subroutine GRD_input_vgrid

end module mod_misc
