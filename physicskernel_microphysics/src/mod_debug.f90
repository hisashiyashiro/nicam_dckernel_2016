!-------------------------------------------------------------------------------
!>
!! Debug utility module
!!
!! @par Description
!!         This module is for dubug.
!!
!! @author  H.Tomita
!!
!! @par History
!! @li      2012-06-29 (H.Yashiro)  [NEW]
!<
module mod_debug
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mpi
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
  public :: PROF_setup
  public :: PROF_setprefx
  public :: PROF_rapstart
  public :: PROF_rapend
  public :: PROF_rapreport

  public :: PROF_valcheck

  interface PROF_valcheck
     module procedure PROF_valcheck_SP_1D
     module procedure PROF_valcheck_SP_2D
     module procedure PROF_valcheck_SP_3D
     module procedure PROF_valcheck_SP_4D
     module procedure PROF_valcheck_SP_5D
     module procedure PROF_valcheck_SP_6D
     module procedure PROF_valcheck_DP_1D
     module procedure PROF_valcheck_DP_2D
     module procedure PROF_valcheck_DP_3D
     module procedure PROF_valcheck_DP_4D
     module procedure PROF_valcheck_DP_5D
     module procedure PROF_valcheck_DP_6D
  end interface PROF_valcheck

  public :: MISC_make_idstr        !--- make file name with a number
  public :: IO_get_available_fid !--- get an available file ID
  public :: MISC_gammafunc         !--- Gamma function

  public :: ADM_proc_stop
  public :: ADM_MPItime
  public :: GRD_setup
  public :: GRD_input_vgrid
  public :: cnvvar_rhogkin_in

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  integer,             public :: EX_STEP = 49
  integer,             public :: EX_rgnid
  integer,             public :: EX_fid
  integer,             public :: EX_err
  character(len=1024), public :: EX_fname

  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  private :: get_rapid

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer,                  private, parameter :: PROF_rapnlimit = 300
  character(len=H_SHORT),   private            :: PROF_prefix    = ''
  integer,                  private            :: PROF_rapnmax   = 0
  character(len=H_SHORT*2), private            :: PROF_rapname (PROF_rapnlimit)
  integer,                  private            :: PROF_grpnmax   = 0
  character(len=H_SHORT),   private            :: PROF_grpname (PROF_rapnlimit)
  integer,                  private            :: PROF_grpid   (PROF_rapnlimit)
  real(DP),                 private            :: PROF_raptstr (PROF_rapnlimit)
  real(DP),                 private            :: PROF_rapttot (PROF_rapnlimit)
  integer,                  private            :: PROF_rapnstr (PROF_rapnlimit)
  integer,                  private            :: PROF_rapnend (PROF_rapnlimit)
  integer,                  private            :: PROF_raplevel(PROF_rapnlimit)

  integer,                  private, parameter :: PROF_default_rap_level = 2
  integer,                  private            :: PROF_rap_level         = 2
  logical,                  private            :: PROF_mpi_barrier       = .false.

  character(len=7),         private            :: PROF_header
  character(len=16),        private            :: PROF_item
  real(DP),                 private            :: PROF_max
  real(DP),                 private            :: PROF_min
  real(DP),                 private            :: PROF_sum

  integer, private, parameter :: min_fid = 7
  integer, private, parameter :: max_fid = 99
  logical, private, parameter :: NSTR_ZERO_START = .true.
  integer, private            :: NSTR_MAX_DIGIT  = 5

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> setup
  subroutine PROF_setup
!ESC!    use mod_adm, only: &
!ESC!       ADM_proc_stop
    implicit none

    namelist / PARAM_PROF / &
       PROF_rap_level, &
       PROF_mpi_barrier

    integer :: ierr
    !---------------------------------------------------------------------------

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '++++++ Module[PROF] / Categ[COMMON] / Origin[SCALElib]'

!ESC!    !--- read namelist
!ESC!    rewind(IO_FID_CONF)
!ESC!    read(IO_FID_CONF,nml=PARAM_PROF,iostat=ierr)
!ESC!    if( ierr < 0 ) then !--- missing
!ESC!       write(IO_FID_LOG,*) '*** Not found namelist. Default used.'
!ESC!    elseif( ierr > 0 ) then !--- fatal error
!ESC!       write(*,*) 'xxx Not appropriate names in namelist PARAM_PROF. Check!'
!ESC!       call ADM_proc_stop
!ESC!    endif
!ESC!    write(IO_FID_LOG,nml=PARAM_PROF)
!ESC!
!ESC!    write(IO_FID_LOG,*) '*** Rap output level              = ', PROF_rap_level
!ESC!    write(IO_FID_LOG,*) '*** Add MPI_barrier in every rap? = ', PROF_mpi_barrier

    PROF_prefix = ''

    return
  end subroutine PROF_setup

  !-----------------------------------------------------------------------------
  subroutine PROF_setprefx( &
       prefxname )
    implicit none

    character(len=*), intent(in) :: prefxname !< prefix
    !---------------------------------------------------------------------------

    if ( prefxname == '' ) then !--- no prefix
       PROF_prefix = ''
    else
       PROF_prefix = trim(prefxname)//'_'
    endif

    return
  end subroutine PROF_setprefx

  !-----------------------------------------------------------------------------
  !> Start raptime
  subroutine PROF_rapstart( rapname_base, level )
!ESC!    use mod_adm, only: &
!ESC!       ADM_MPIbarrier, &
!ESC!       ADM_MPItime
    implicit none

    character(len=*), intent(in) :: rapname_base    !< name of item

    integer,          intent(in), optional :: level !< level of item (default is 2)

    character(len=H_SHORT*2) :: rapname             !< name of item with prefix
    integer                  :: id, level_
    !---------------------------------------------------------------------------

    if ( present(level) ) then
       level_ = level
    else
       level_ = PROF_default_rap_level
    endif

    if( level_ > PROF_rap_level ) return

    rapname = trim(PROF_prefix)//trim(rapname_base)

    id = get_rapid( rapname, level_ )

!ESC!    if(PROF_mpi_barrier) call ADM_MPIbarrier

    PROF_raptstr(id) = ADM_MPItime()
    PROF_rapnstr(id) = PROF_rapnstr(id) + 1

    !write(IO_FID_LOG,*) rapname, PROF_rapnstr(id)

#ifdef _FAPP_
    call FAPP_START( trim(rapname), id, level_ )
#endif
#ifdef _FINEPA_
    call START_COLLECTION( trim(rapname) )
#endif

    return
  end subroutine PROF_rapstart

  !-----------------------------------------------------------------------------
  !> Save raptime
  subroutine PROF_rapend( rapname_base, level )
!ESC!    use mod_adm, only: &
!ESC!       ADM_MPIbarrier, &
!ESC!       ADM_MPItime
    implicit none

    character(len=*), intent(in) :: rapname_base    !< name of item

    integer,          intent(in), optional :: level !< level of item

    character(len=H_SHORT*2) :: rapname             !< name of item with prefix
    integer                  :: id, level_
    !---------------------------------------------------------------------------

    if ( present(level) ) then
       if( level > PROF_rap_level ) return
    endif

    rapname = trim(PROF_prefix)//trim(rapname_base)

    id = get_rapid( rapname, level_ )

    if( level_ > PROF_rap_level ) return

!ESC!    if(PROF_mpi_barrier) call ADM_MPIbarrier

    PROF_rapttot(id) = PROF_rapttot(id) + ( ADM_MPItime()-PROF_raptstr(id) )
    PROF_rapnend(id) = PROF_rapnend(id) + 1

#ifdef _FINEPA_
    call STOP_COLLECTION( trim(rapname) )
#endif
#ifdef _FAPP_
    call FAPP_STOP( trim(rapname), id, level_ )
#endif

    return
  end subroutine PROF_rapend

  !-----------------------------------------------------------------------------
  !> Report raptime
  subroutine PROF_rapreport
!ESC!    use mod_adm, only: &
!ESC!       ADM_MPItimestat, &
!ESC!       ADM_IsMaster
    implicit none

    real(DP) :: avgvar(PROF_rapnlimit)
    real(DP) :: maxvar(PROF_rapnlimit)
    real(DP) :: minvar(PROF_rapnlimit)
    integer  :: maxidx(PROF_rapnlimit)
    integer  :: minidx(PROF_rapnlimit)

    integer :: id, gid
    integer :: fid
    !---------------------------------------------------------------------------

    do id = 1, PROF_rapnmax
       if ( PROF_rapnstr(id) /= PROF_rapnend(id) ) then
           write(*,*) '*** Mismatch Report',id,PROF_rapname(id),PROF_rapnstr(id),PROF_rapnend(id)
       endif
    enddo

    write(IO_FID_LOG,*)
    write(IO_FID_LOG,*) '*** Computational Time Report'
    write(IO_FID_LOG,*) '*** Rap level is ', PROF_rap_level

!ESC!    if ( .false. ) then ! report for each node

       do gid = 1, PROF_rapnmax
       do id  = 1, PROF_rapnmax
          if (       PROF_raplevel(id) <= PROF_rap_level &
               .AND. PROF_grpid(id)    == gid            ) then
             write(IO_FID_LOG,'(1x,A,I3.3,A,A,A,F10.3,A,I9)') &
                  '*** ID=',id,' : ',PROF_rapname(id),' T=',PROF_rapttot(id),' N=',PROF_rapnstr(id)
          endif
       enddo
       enddo

!ESC!    else
!ESC!
!ESC!       call ADM_MPItimestat( avgvar      (1:PROF_rapnmax), & ! [OUT]
!ESC!                             maxvar      (1:PROF_rapnmax), & ! [OUT]
!ESC!                             minvar      (1:PROF_rapnmax), & ! [OUT]
!ESC!                             maxidx      (1:PROF_rapnmax), & ! [OUT]
!ESC!                             minidx      (1:PROF_rapnmax), & ! [OUT]
!ESC!                             PROF_rapttot(1:PROF_rapnmax)  ) ! [IN]
!ESC!
!ESC!       fid = -1
!ESC!       if ( .false. ) then ! report to STDOUT
!ESC!          if ( ADM_IsMaster ) then
!ESC!             write(*,*) '*** Computational Time Report'
!ESC!             fid = 6 ! master node
!ESC!          endif
!ESC!       else
!ESC!          fid = IO_FID_LOG
!ESC!       endif
!ESC!
!ESC!       do gid = 1, PROF_rapnmax
!ESC!       do id  = 1, PROF_rapnmax
!ESC!          if (       PROF_raplevel(id) <= PROF_rap_level &
!ESC!               .AND. PROF_grpid(id)    == gid            &
!ESC!               .AND. fid > 0                             ) then
!ESC!             write(IO_FID_LOG,'(1x,A,I3.3,A,A,A,F10.3,A,F10.3,A,I5,A,A,F10.3,A,I5,A,A,I9)') &
!ESC!                  '*** ID=',id,' : ',PROF_rapname(id), &
!ESC!                  ' T(avg)=',avgvar(id), &
!ESC!                  ', T(max)=',maxvar(id),'[',maxidx(id),']', &
!ESC!                  ', T(min)=',minvar(id),'[',minidx(id),']', &
!ESC!                  ' N=',PROF_rapnstr(id)
!ESC!          endif
!ESC!       enddo
!ESC!       enddo
!ESC!
!ESC!    endif

    return
  end subroutine PROF_rapreport

  !-----------------------------------------------------------------------------
  !> Get item ID or register item
  function get_rapid( rapname, level ) result(id)
    implicit none

    character(len=*), intent(in)    :: rapname !< name of item
    integer,          intent(inout) :: level   !< level of item
    integer                         :: id

    character(len=H_SHORT*2) :: trapname
    character(len=H_SHORT)   :: trapname2
    !---------------------------------------------------------------------------

    trapname  = trim(rapname)
    trapname2 = trim(rapname)

    do id = 1, PROF_rapnmax
       if ( trapname == PROF_rapname(id) ) then
          level = PROF_raplevel(id)
          return
       endif
    enddo

    PROF_rapnmax     = PROF_rapnmax + 1
    id               = PROF_rapnmax
    PROF_rapname(id) = trapname

    PROF_rapnstr(id) = 0
    PROF_rapnend(id) = 0
    PROF_rapttot(id) = 0.0_DP

    PROF_grpid   (id) = get_grpid(trapname2)
    PROF_raplevel(id) = level

    return
  end function get_rapid

  !-----------------------------------------------------------------------------
  !> Get group ID
  function get_grpid( rapname ) result(gid)
    implicit none

    character(len=*), intent(in) :: rapname !< name of item
    integer                      :: gid

    character(len=H_SHORT) :: grpname
    integer                :: idx
    !---------------------------------------------------------------------------

    idx = index(rapname," ")
    if ( idx > 1 ) then
       grpname = rapname(1:idx-1)
    else
       grpname = rapname
    endif

    do gid = 1, PROF_grpnmax
       if( grpname == PROF_grpname(gid) ) return
    enddo

    PROF_grpnmax      = PROF_grpnmax + 1
    gid               = PROF_grpnmax
    PROF_grpname(gid) = grpname

    return
  end function get_grpid

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_SP_1D( &
       header,  &
       varname, &
       var      )
    implicit none

    character(len=*),  intent(in)  :: header
    character(len=*),  intent(in)  :: varname
    real(SP),          intent(in)  :: var(:)
    !---------------------------------------------------------------------------

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    PROF_max    = real(maxval(var),kind=DP)
    PROF_min    = real(minval(var),kind=DP)
    PROF_sum    = real(sum   (var),kind=DP)
    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,',min=',PROF_min,',sum=',PROF_sum

    return
  end subroutine PROF_valcheck_SP_1D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_SP_2D( &
       header,  &
       varname, &
       var      )
    implicit none

    character(len=*),  intent(in)  :: header
    character(len=*),  intent(in)  :: varname
    real(SP),          intent(in)  :: var(:,:)
    !---------------------------------------------------------------------------

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    PROF_max    = real(maxval(var),kind=DP)
    PROF_min    = real(minval(var),kind=DP)
    PROF_sum    = real(sum   (var),kind=DP)
    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,',min=',PROF_min,',sum=',PROF_sum

    return
  end subroutine PROF_valcheck_SP_2D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_SP_3D( &
       header,  &
       varname, &
       var      )
    implicit none

    character(len=*),  intent(in)  :: header
    character(len=*),  intent(in)  :: varname
    real(SP),          intent(in)  :: var(:,:,:)
    !---------------------------------------------------------------------------

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    PROF_max    = real(maxval(var),kind=DP)
    PROF_min    = real(minval(var),kind=DP)
    PROF_sum    = real(sum   (var),kind=DP)
    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,',min=',PROF_min,',sum=',PROF_sum

    return
  end subroutine PROF_valcheck_SP_3D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_SP_4D( &
       header,  &
       varname, &
       var      )
    implicit none

    character(len=*),  intent(in)  :: header
    character(len=*),  intent(in)  :: varname
    real(SP),          intent(in)  :: var(:,:,:,:)
    !---------------------------------------------------------------------------

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    PROF_max    = real(maxval(var),kind=DP)
    PROF_min    = real(minval(var),kind=DP)
    PROF_sum    = real(sum   (var),kind=DP)
    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,',min=',PROF_min,',sum=',PROF_sum

    return
  end subroutine PROF_valcheck_SP_4D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_SP_5D( &
       header,  &
       varname, &
       var      )
    implicit none

    character(len=*),  intent(in)  :: header
    character(len=*),  intent(in)  :: varname
    real(SP),          intent(in)  :: var(:,:,:,:,:)
    !---------------------------------------------------------------------------

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    PROF_max    = real(maxval(var),kind=DP)
    PROF_min    = real(minval(var),kind=DP)
    PROF_sum    = real(sum   (var),kind=DP)
    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,',min=',PROF_min,',sum=',PROF_sum

    return
  end subroutine PROF_valcheck_SP_5D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_SP_6D( &
       header,  &
       varname, &
       var      )
    implicit none

    character(len=*),  intent(in)  :: header
    character(len=*),  intent(in)  :: varname
    real(SP),          intent(in)  :: var(:,:,:,:,:,:)
    !---------------------------------------------------------------------------

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    PROF_max    = real(maxval(var),kind=DP)
    PROF_min    = real(minval(var),kind=DP)
    PROF_sum    = real(sum   (var),kind=DP)
    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,',min=',PROF_min,',sum=',PROF_sum

    return
  end subroutine PROF_valcheck_SP_6D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_DP_1D( &
       header,  &
       varname, &
       var      )
    implicit none

    character(len=*),  intent(in)  :: header
    character(len=*),  intent(in)  :: varname
    real(DP),          intent(in)  :: var(:)
    !---------------------------------------------------------------------------

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    PROF_max    = real(maxval(var),kind=DP)
    PROF_min    = real(minval(var),kind=DP)
    PROF_sum    = real(sum   (var),kind=DP)
    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,',min=',PROF_min,',sum=',PROF_sum

    return
  end subroutine PROF_valcheck_DP_1D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_DP_2D( &
       header,  &
       varname, &
       var      )
    implicit none

    character(len=*),  intent(in)  :: header
    character(len=*),  intent(in)  :: varname
    real(DP),          intent(in)  :: var(:,:)
    !---------------------------------------------------------------------------

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    PROF_max    = real(maxval(var),kind=DP)
    PROF_min    = real(minval(var),kind=DP)
    PROF_sum    = real(sum   (var),kind=DP)
    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,',min=',PROF_min,',sum=',PROF_sum

    return
  end subroutine PROF_valcheck_DP_2D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_DP_3D( &
       header,  &
       varname, &
       var      )
    implicit none

    character(len=*),  intent(in)  :: header
    character(len=*),  intent(in)  :: varname
    real(DP),          intent(in)  :: var(:,:,:)
    !---------------------------------------------------------------------------

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    PROF_max    = real(maxval(var),kind=DP)
    PROF_min    = real(minval(var),kind=DP)
    PROF_sum    = real(sum   (var),kind=DP)
    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,',min=',PROF_min,',sum=',PROF_sum

    return
  end subroutine PROF_valcheck_DP_3D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_DP_4D( &
       header,  &
       varname, &
       var      )
    implicit none

    character(len=*),  intent(in)  :: header
    character(len=*),  intent(in)  :: varname
    real(DP),          intent(in)  :: var(:,:,:,:)
    !---------------------------------------------------------------------------

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    PROF_max    = real(maxval(var),kind=DP)
    PROF_min    = real(minval(var),kind=DP)
    PROF_sum    = real(sum   (var),kind=DP)
    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,',min=',PROF_min,',sum=',PROF_sum

    return
  end subroutine PROF_valcheck_DP_4D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_DP_5D( &
       header,  &
       varname, &
       var      )
    implicit none

    character(len=*),  intent(in)  :: header
    character(len=*),  intent(in)  :: varname
    real(DP),          intent(in)  :: var(:,:,:,:,:)
    !---------------------------------------------------------------------------

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    PROF_max    = real(maxval(var),kind=DP)
    PROF_min    = real(minval(var),kind=DP)
    PROF_sum    = real(sum   (var),kind=DP)
    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,',min=',PROF_min,',sum=',PROF_sum

    return
  end subroutine PROF_valcheck_DP_5D

  !-----------------------------------------------------------------------------
  subroutine PROF_valcheck_DP_6D( &
       header,  &
       varname, &
       var      )
    implicit none

    character(len=*),  intent(in)  :: header
    character(len=*),  intent(in)  :: varname
    real(DP),          intent(in)  :: var(:,:,:,:,:,:)
    !---------------------------------------------------------------------------

    PROF_header = trim(header)
    PROF_item   = trim(varname)
    PROF_max    = real(maxval(var),kind=DP)
    PROF_min    = real(minval(var),kind=DP)
    PROF_sum    = real(sum   (var),kind=DP)
    write(IO_FID_LOG,'(1x,A,A7,A,A16,3(A,ES24.16))') &
    '+',PROF_header,'[',PROF_item,'] max=',PROF_max,',min=',PROF_min,',sum=',PROF_sum

    return
  end subroutine PROF_valcheck_DP_6D

  !-----------------------------------------------------------------------------
  subroutine MISC_make_idstr( &
       str,     & !--- [OUT]
       prefix,  & !--- [IN]
       ext,     & !--- [IN]
       numID,   & !--- [IN]
       digit    ) !--- [IN]
    implicit none

    character(len=*), intent(out) :: str    ! combined string (file name)
    character(len=*), intent(in)  :: prefix ! prefix
    character(len=*), intent(in)  :: ext    ! extention( e.g. .rgn )
    integer,          intent(in)  :: numID  ! number

    integer, optional, intent(in) :: digit  ! digit

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
  function IO_get_available_fid()  &
       result(fid)                     !--- file id
    implicit none

    integer :: fid
    logical :: i_opened

    do fid = min_fid,max_fid
       INQUIRE (fid, OPENED=I_OPENED)
       if(.not.I_OPENED) return
    enddo

  end function IO_get_available_fid

  !-----------------------------------------------------------------------------
  function MISC_gammafunc( xx ) result(f)
    implicit none
    real(RP), intent(in) :: xx
    real(RP) :: f
    real(RP) :: coef(6)=(/&
         +76.18009172947146_RP,&
         -86.50532032941677_RP,&
         +24.01409824083091_RP,&
         -1.231739572450155_RP,&
         +0.1208650973866179E-2_RP,&
         -0.5395239384953E-5_RP&
         /)
    integer :: j
    real(RP) :: x,y,tmp,ser

    x=xx
    y=x
    tmp=x+5.5_RP
    tmp = tmp - (x+0.5)*log(tmp)
    ser=1.000000000190015_RP
    do j=1,6
       y=y+1
       ser = ser+coef(j)/y
    enddo
    f = exp(-tmp+log(2.5066282746310005_RP*ser/x))
  end function MISC_gammafunc

  !-----------------------------------------------------------------------------
  subroutine ADM_proc_stop
    stop
  end subroutine ADM_proc_stop

  !-----------------------------------------------------------------------------
  !> Get MPI time
  !> @return time
  function ADM_MPItime() result(time)
    implicit none

    real(DP) :: time
    !---------------------------------------------------------------------------

!ESC!    if ( ADM_myprc_is_run ) then
!ESC!       time = real(MPI_WTIME(), kind=8)
!ESC!    else
       call cpu_time(time)
!ESC!    endif

  end function ADM_MPItime

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
       GRD_rdgz (k) = 1.0_RP / grd_dgz(k)
       GRD_rdgzh(k) = 1.0_RP / grd_dgzh(k)
    enddo

    !---< vertical interpolation factor >---
    allocate( GRD_afact(ADM_kall) )
    allocate( GRD_bfact(ADM_kall) )
    allocate( GRD_cfact(ADM_kall) )
    allocate( GRD_dfact(ADM_kall) )

    ! vertical interpolation factor
    do k = ADM_kmin, ADM_kmax+1
       GRD_afact(k) = ( GRD_gzh(k) - GRD_gz(k-1) ) &
                    / ( GRD_gz (k) - GRD_gz(k-1) )
    enddo
    GRD_afact(ADM_kmin-1) = 1.0_RP

    GRD_bfact(:) = 1.0_RP - GRD_afact(:)

    do k = ADM_kmin, ADM_kmax
       GRD_cfact(k) = ( GRD_gz (k  ) - GRD_gzh(k) ) &
                    / ( GRD_gzh(k+1) - GRD_gzh(k) )
    enddo
    GRD_cfact(ADM_kmin-1) = 1.0_RP
    GRD_cfact(ADM_kmax+1) = 0.0_RP

    GRD_dfact(:) = 1.0_RP - GRD_cfact(:)

    return
  end subroutine GRD_setup

  !-----------------------------------------------------------------------------
  subroutine GRD_input_vgrid( &
       fname  )
!ESC!    use mod_misc,  only :&
!ESC!       IO_get_available_fid
!ESC!    use mod_adm, only: &
!ESC!       ADM_proc_stop, &
!ESC!       ADM_vlayer
    implicit none

    character(len=*), intent(in) :: fname

    real(DP) :: GRD_gz_DP (ADM_kall)
    real(DP) :: GRD_gzh_DP(ADM_kall)

    integer :: num_of_layer
    integer :: fid, ierr
    !---------------------------------------------------------------------------

    fid = IO_get_available_fid()
    open( unit   = fid,           &
          file   = trim(fname),   &
          status = 'old',         &
          form   = 'unformatted', &
          access = 'sequential',  &
          iostat = ierr           )

       if ( ierr /= 0 ) then
          write(IO_FID_LOG,*) 'xxx No vertical grid file : ', trim(fname)
          call ADM_proc_stop
       endif

       read(fid) num_of_layer
       if ( num_of_layer /= ADM_vlayer ) then
          write(IO_FID_LOG,*)          'Msg : Sub[GRD_input_vgrid]/Mod[grid]'
          write(IO_FID_LOG,*)          '   *** inconsistency in number of vertical layers.'
          call ADM_proc_stop
       endif

       read(fid) GRD_gz_DP
       read(fid) GRD_gzh_DP

    close(fid)

    GRD_gz  = real(GRD_gz_DP ,kind=RP)
    GRD_gzh = real(GRD_gzh_DP,kind=RP)

    return
  end subroutine GRD_input_vgrid

  !-----------------------------------------------------------------------------
  subroutine cnvvar_rhogkin_in( &
       ijdim,     &
       kdim,      &
       rhog,      &
       rhogvx,    &
       rhogvy,    &
       rhogvz,    &
       rhogw,     &
       C2Wfact,   &
       W2Cfact,   &
       rhogkin,   &
       rhogkin_h, &
       rhogkin_v  )
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    real(RP), intent(in)  :: rhog     (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 )
    real(RP), intent(in)  :: rhogvx   (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X vx
    real(RP), intent(in)  :: rhogvy   (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X vy
    real(RP), intent(in)  :: rhogvz   (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X vz
    real(RP), intent(in)  :: rhogw    (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X w
    real(RP), intent(in)  :: C2Wfact  (ijdim,kdim,2) ! rho X ( G^1/2 X gamma2 ) X w
    real(RP), intent(in)  :: W2Cfact  (ijdim,kdim,2) ! rho X ( G^1/2 X gamma2 ) X w
    real(RP), intent(out) :: rhogkin  (ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X kin
    real(RP), intent(out) :: rhogkin_h(ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X kin (horizontal)
    real(RP), intent(out) :: rhogkin_v(ijdim,kdim)   ! rho X ( G^1/2 X gamma2 ) X kin (vertical)

    integer  :: gall, kmin, kmax

    integer  :: g, k
    !---------------------------------------------------------------------------

    call PROF_rapstart('CNV_rhogkin',2)

    gall = ijdim
    kmin = 2
    kmax = kdim-1

    !$omp parallel default(none),private(g,k),                           &
    !$omp shared(gall,kmin,kmax,rhog,rhogvx,rhogvy,rhogvz,rhogw,rhogkin, &
    !$omp        rhogkin_h,rhogkin_v,C2Wfact,W2Cfact)

    !--- horizontal kinetic energy
    !$omp do
    do k = kmin, kmax
    do g = 1, gall
       rhogkin_h(g,k) = 0.5_RP * ( rhogvx(g,k) * rhogvx(g,k) &
                                 + rhogvy(g,k) * rhogvy(g,k) &
                                 + rhogvz(g,k) * rhogvz(g,k) ) / rhog(g,k)
    enddo
    enddo
    !$omp end do

    !$omp workshare
    rhogkin_h(:,kmin-1) = 0.0_RP
    rhogkin_h(:,kmax+1) = 0.0_RP
    !$omp end workshare

    !--- vertical kinetic energy
    !$omp do
    do k = kmin+1, kmax
    do g = 1, gall
       rhogkin_v(g,k) = 0.5_RP * ( rhogw(g,k) * rhogw(g,k) ) &
                      / ( C2Wfact(g,k,1) * rhog(g,k  ) &
                        + C2Wfact(g,k,2) * rhog(g,k-1) )
    enddo
    enddo
    !$omp end do

    !$omp workshare
    rhogkin_v(:,kmin-1) = 0.0_RP
    rhogkin_v(:,kmin  ) = 0.0_RP
    rhogkin_v(:,kmax+1) = 0.0_RP
    !$omp end workshare

    !--- total kinetic energy
    !$omp do
    do k = kmin, kmax
    do g = 1, gall
       rhogkin(g,k) = rhogkin_h(g,k)                      & ! horizontal
                    + ( W2Cfact(g,k,1) * rhogkin_v(g,k+1) & ! vertical
                      + W2Cfact(g,k,2) * rhogkin_v(g,k  ) )
    enddo
    enddo
    !$omp end do

    !$omp workshare
    rhogkin(:,kmin-1) = 0.0_RP
    rhogkin(:,kmax+1) = 0.0_RP
    !$omp end workshare

    !$omp end parallel

    call PROF_rapend('CNV_rhogkin',2)

    return
  end subroutine cnvvar_rhogkin_in

end module mod_debug
