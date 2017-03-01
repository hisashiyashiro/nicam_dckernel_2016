!-------------------------------------------------------------------------------
!>
!! Calendar module
!!
!! @par Description
!!         This module provides the subroutines for calendar.
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_calendar
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
!ESC!  use mod_stdio
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedures
  !
  public :: calendar_setup
  public :: calendar_ointvl
  public :: calendar_ss2cc
  public :: calendar_ssaft
  public :: calendar_secdy
  public :: calendar_ss2ds
  public :: calendar_ss2yd
  public :: calendar_ym2dd
  public :: calendar_ds2ss
  public :: calendar_dayyr
  public :: calendar_ss2yh
  public :: calendar_ss2ym
  public :: calendar_xx2ss
  public :: calendar_cc2yh
  public :: calendar_yh2ss
  public :: calendar_PERPR
  public :: calendar_dd2ym
  public :: calendar_daymo

  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  integer, private :: isecdy, isechr, jy, jy4, jcent, jcent4
  integer, private :: idays0, idy, ileap, id, m, idayyr, jyear, jmonth
  !
  !--- flag of automatic or not
  logical, private, save ::  oauto = .true.
  !<---                      yr=0-999     : 360day
  !<---                      yr=1000-1899 : 365day
  !<---                      yr=1900-     : gregorian
  !
  !--- flag of the gregorian calendar or not
  logical, private, save :: ogrego = .true.
  !
  !--- number of days in a month
  integer, private, save :: monday ( 12,2 )
  data   monday /                            &
       31,28,31,30,31,30,31,31,30,31,30,31,  &
       31,29,31,30,31,30,31,31,30,31,30,31 /
  !
  !--- flag of ideal calender (n day per month)
  logical, private, save :: oideal = .false.
  !------ 1 month = x days in the ideal case
  integer, private, save :: idaymo = 30
  !------ 1 year = x months in the ideal case
  integer, private, save :: imonyr = 12
  !
  !--- flag of perpetual or not
  logical, private, save :: operpt = .false.
  !------ perpetual date(year)
  integer, private, save :: iyrpp  = 0
  !------ perpetual date(month)
  integer, private, save :: imonpp = 3
  !------ perpetual date(day)
  integer, private, save :: idaypp = 21
  !
  !--- 1 minute = x sec.
  integer, private, save :: isecmn = 60
  !--- 1 hour = x minutes
  integer, private, save :: iminhr = 60
  !--- 1 day = x hours
  integer, private, save :: ihrday = 24
  !
  logical, private, save :: ooperz = .false.
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> setup
  subroutine calendar_setup
    implicit none

    namelist  /nm_calendar/  &
         oauto ,             &
         ogrego,             &
         oideal,             &
         operpt,             &
         idaymo,             &
         imonyr,             &
         iyrpp,              &
         imonpp,             &
         idaypp,             &
         isecmn,             &
         iminhr,             &
         ihrday

    integer :: ierr

!ESC!    rewind(IO_FID_CONF)
!ESC!    read(IO_FID_CONF, nm_calendar, iostat=ierr)
!ESC!    if ( ierr < 0 ) then
!ESC!       write(IO_FID_LOG,*) 'Msg : Sub[calendar_setup]/Mod[calendar]'
!ESC!       write(IO_FID_LOG,*) ' *** Not found namelist.'
!ESC!       write(IO_FID_LOG,*) ' *** Use default values.'
!ESC!    elseif( ierr > 0 ) then
!ESC!       write(*,*) 'Msg : Sub[calendar_setup]/Mod[calendar]'
!ESC!       write(*,*) ' *** WARNING : Not appropriate names in namelist!! CHECK!!'
!ESC!    endif

    return
  end subroutine calendar_setup

  !-----------------------------------------------------------------------------
  subroutine calendar_PERPR( &    !! calendar, refer to fixed date
       IYEAR ,               & !--- OUT
       IMONTH,               & !--- OUT
       IDAY,                 & !--- OUT
       OOPERP                & !--- OUT
       )
    implicit none

    integer, intent(out) :: IYEAR
    integer, intent(out) :: IMONTH
    integer, intent(out) :: IDAY
    logical, intent(out) :: OOPERP

    OOPERP = OPERPT
    IYEAR  = IYRPP
    IMONTH = IMONPP
    IDAY   = IDAYPP

    RETURN
  end subroutine calendar_PERPR

  !-----------------------------------------------------------------------------
  subroutine calendar_daymo(&
       ndaymo,              & !--- OUT : day
       iyear,               & !--- IN : year
       imonth               & !--- IN : month
       )
    !--- calendar, no.of day in a month
    implicit none

    integer, intent(out) :: ndaymo
    integer, intent(in) :: iyear
    integer, intent(in) :: imonth

    if ( oauto  ) then
       if ( iyear >= 1900 ) then
          ogrego = .true.
       else
          ogrego = .false.
          if ( iyear >= 1000 ) then
             oideal = .false.
          else
             oideal = .true.
             idaymo = 30
             imonyr = 12
          endif
       endif
    endif

    if ( ogrego ) then
       if ( ocleap( iyear ) ) then
          ndaymo = monday( imonth,2 )
       else
          ndaymo = monday( imonth,1 )
       endif
    elseif( .not. oideal ) then
       ndaymo = monday( imonth,1 )
    else
       ndaymo = idaymo
    endif

    return
  end subroutine calendar_daymo

  !-----------------------------------------------------------------------------
  subroutine calendar_dayyr(&
       ndayyr,              & !-- OUT : day
       iyear                & !---IN : year
       )
    !--- calendar, no.of day in an year
    implicit none

    integer, intent(out) :: ndayyr
    integer, intent(in) :: iyear

    if ( oauto  ) then
       if ( iyear >= 1900 ) then
          ogrego = .true.
       else
          ogrego = .false.
          if ( iyear >= 1000 ) then
             oideal = .false.
          else
             oideal = .true.
             idaymo = 30
             imonyr = 12
          endif
       endif
    endif

    if ( ogrego ) then
       if ( ocleap( iyear ) ) then
          ndayyr = 366
       else
          ndayyr = 365
       endif
    elseif( .not. oideal ) then
       ndayyr = 365
    else
       ndayyr = idaymo*imonyr
    endif

    return
  end subroutine calendar_dayyr

  !-----------------------------------------------------------------------------
  subroutine calendar_monyr( &
       nmonyr,               & !--- OUT : month
       iyear                 & !--- IN : year
       )
    !--- calendar, no.of month in an year
    implicit none

    integer, intent(out) :: nmonyr
    integer, intent(in) :: iyear

    if ( oauto  ) then
       nmonyr = 12
    else
       nmonyr = imonyr
    endif

    return
  end subroutine calendar_monyr

  !-----------------------------------------------------------------------------
  subroutine calendar_secdy( &
       nsecdy                & !--- OUT : sec
       )
    !--- calendar, no.of sec. in a day
    implicit none

    integer, intent(out) :: nsecdy
    nsecdy = isecmn*iminhr*ihrday

    return
  end subroutine calendar_secdy

  !-----------------------------------------------------------------------------
  subroutine calendar_secmi( &
       nsecmi                & !--- OUT
       )
    !--- calendar, no of sec. in a minute
    implicit none

    integer, intent(out) :: nsecmi
    nsecmi = isecmn

    return
  end subroutine calendar_secmi

  !-----------------------------------------------------------------------------
  subroutine calendar_sechr( &
       nsechr                &
       )
    !--- calendar, no.of sec. in an hour
    implicit none

    integer, intent(out) :: nsechr
    nsechr = isecmn*iminhr

    return
  end subroutine calendar_sechr

  !-----------------------------------------------------------------------------
  subroutine calendar_ss2ds( &
       idays, &
       rsec,  &
       dsec   )
    implicit none

    integer,  intent(out) :: idays
    real(DP), intent(out) :: rsec
    real(DP), intent(in)  :: dsec
    !---------------------------------------------------------------------------

    isecdy = isecmn * iminhr * ihrday
    idays  = int( dsec / real(isecdy,kind=DP) ) + 1
    rsec   = dsec - real((idays-1),kind=DP) * real(isecdy,kind=DP)

    if ( nint(rsec) >= isecdy ) then
       idays = idays + 1
       rsec  = rsec - real(isecdy,kind=DP)
    endif

    return
  end subroutine calendar_ss2ds

  !-----------------------------------------------------------------------------
  subroutine calendar_ds2ss( &
       dsec,  &
       idays, &
       rsec   )
    implicit none

    real(DP), intent(out) :: dsec
    integer,  intent(in)  :: idays
    real(DP), intent(in)  :: rsec
    !---------------------------------------------------------------------------

    isecdy = isecmn * iminhr * ihrday
    dsec   = real(idays-1,kind=DP) * real(isecdy,kind=DP) + real(rsec,kind=DP)

    return
  end subroutine calendar_ds2ss

  !-----------------------------------------------------------------------------
  subroutine calendar_rs2hm( &
       ihour, &
       imin,  &
       isec,  &
       rsec   )
    implicit none

    integer,  intent(out) :: ihour
    integer,  intent(out) :: imin
    integer,  intent(out) :: isec
    real(DP), intent(in)  :: rsec
    !---------------------------------------------------------------------------

    isechr = isecmn * iminhr
    ihour  = int ( rsec / real(isechr,kind=DP) )
    imin   = int ( ( rsec - real(ihour*isechr,kind=DP) ) / real(isecmn,kind=DP) )
    isec   = nint( rsec - real(ihour*isechr,kind=DP) - real(imin*isecmn,kind=DP) )

    if ( isec >= isecmn ) then
       imin  = imin + 1
       isec  = isec - isecmn
    endif

    if ( imin == iminhr ) then
       ihour = ihour + 1
       imin  = imin  - iminhr
    endif

    return
  end subroutine calendar_rs2hm

  !-----------------------------------------------------------------------------
  subroutine calendar_hm2rs( &
       rsec  ,               & !--- OUT
       ihour ,               & !--- IN
       imin  ,               & !--- IN
       isec                  & !--- IN
       )
    !--- calendar, hhmmss -> sec.
    implicit none

    real(DP), intent(out) :: rsec
    integer, intent(in) :: ihour
    integer, intent(in) :: imin
    integer, intent(in) :: isec

    rsec = real( ihour*isecmn*iminhr + imin*isecmn + isec, kind=DP )

    return
  end subroutine calendar_hm2rs

  !-----------------------------------------------------------------------------
  subroutine calendar_dd2ym( &
       iyear,  &
       imonth, &
       iday,   &
       idays   )
    implicit none

    integer, intent(out) :: iyear
    integer, intent(out) :: imonth
    integer, intent(out) :: iday
    integer, intent(in)  :: idays
    !---------------------------------------------------------------------------

    if ( oauto ) then
       if ( idays >= 693961 ) then       !" 1900*365+1900/4-19+5
          ogrego = .true.
       else
          ogrego = .false.
          if ( idays >= 1000*365 ) then
             oideal = .false.
          else
             oideal = .true.
             idaymo = 30
             imonyr = 12
          endif
       endif
    endif

    if ( operpt .AND. ooperz ) then
       iyear  = iyrpp
       imonth = imonpp
       iday   = idaypp
       return
    endif

    if ( ogrego ) then
       jy     = int(real(idays,kind=DP)/365.24_DP)
1100   continue
       jy4    = (jy+3)/4
       jcent  = (jy+99)/100
       jcent4 = (jy+399)/400
       idays0 = jy*365+jy4-jcent+jcent4
       if ( idays <= idays0 ) then
          jy = jy -1
          if ( jy >= 0 ) goto 1100
       endif
       iyear = jy
       idy   = idays - idays0
       if ( ocleap( iyear ) ) then
          ileap  = 2
       else
          ileap  = 1
       endif
    elseif( .not. oideal ) then
       iyear = idays/365
       idy   = idays - iyear*365
       ileap = 1
    endif
    if ( ogrego .OR. .not. oideal ) then
       id = 0
       do m = 1, 12
          id = id + monday(m,ileap)
          if ( idy <= id ) then
             imonth = m
             iday   = idy-id+monday(m,ileap)
             exit
          endif
       enddo
    else
       idayyr = idaymo*imonyr
       iyear  = ( idays-1 ) / idayyr
       imonth = ( idays-1 - iyear*idayyr )/idaymo+1
       iday   = idays - iyear*idayyr - (imonth-1)*idaymo
    endif

    return
  end subroutine calendar_dd2ym

  !-----------------------------------------------------------------------------
  subroutine calendar_ym2dd( &
       idays ,               & !--- OUT
       iyear ,               & !--- IN
       imonth,               & !--- IN
       iday                  & !--- IN
       )
    !--- calendar, yymmdd -> day
    implicit none
    integer, intent(out) :: idays
    integer, intent(in) :: iyear
    integer, intent(in) :: imonth
    integer, intent(in) :: iday
    if ( oauto  ) then
       if ( iyear >= 1900 ) then
          ogrego = .true.
       else
          ogrego = .false.
          if ( iyear >= 1000 ) then
             oideal = .false.
          else
             oideal = .true.
             idaymo = 30
             imonyr = 12
          endif
       endif
    endif
    if ( ogrego .OR. .not. oideal ) then
       if ( imonth > 0 ) then
          jyear  = iyear + (imonth-1)/12
          jmonth = mod(imonth-1,12)+1
       else
          jyear  = iyear - (-imonth)/12 - 1
          jmonth = 12-mod(-imonth,12)
       endif
    endif

    if ( ogrego ) then
       jy4    = (jyear+3)/4
       jcent  = (jyear+99)/100
       jcent4 = (jyear+399)/400
       idays0 = jyear*365+jy4-jcent+jcent4
       if ( ocleap( jyear ) ) then
          ileap = 2
       else
          ileap = 1
       endif
    elseif( .not. oideal ) then
       idays0 = jyear*365
       ileap  = 1
    endif

    if ( ogrego .OR. .not. oideal ) then
       id = 0
       do m = 1, jmonth-1
          id = id + monday(m,ileap)
       enddo
    else
       idays0 = iyear*idaymo*imonyr
       id     = (imonth-1)*idaymo
    endif
    idays = idays0 + id + iday

    return
  end subroutine calendar_ym2dd
  !-----------------------------------------------------------------------------
  subroutine calendar_ym2yd( &
       idaysy,               & !--- OUT
       iyear ,               & !--- IN
       imonth,               & !--- IN
       iday                  & !--- IN
       )
    !--- calendar, yymmdd -> yydd
    implicit none
    integer, intent(out) :: idaysy
    integer, intent(in) :: iyear
    integer, intent(in) :: imonth
    integer, intent(in) :: iday
    if ( oauto  ) then
       if ( iyear >= 1900 ) then
          ogrego = .true.
       else
          ogrego = .false.
          if ( iyear >= 1000 ) then
             oideal = .false.
          else
             oideal = .true.
             idaymo = 30
             imonyr = 12
          endif
       endif
    endif
    if ( ogrego .OR. .not. oideal ) then
       if ( imonth > 0 ) then
          jyear  = iyear + (imonth-1) / 12
          jmonth = mod(imonth-1,12) + 1
       else
          jyear  = iyear - (-imonth) / 12 - 1
          jmonth = 12 - mod(-imonth,12)
       endif
    endif

    if ( ogrego ) then
       if ( ocleap( jyear ) ) then
          ileap = 2
       else
          ileap = 1
       endif
    elseif( .not. oideal ) then
       ileap  = 1
    endif

    if ( ogrego .OR. .not. oideal ) then
       id = 0
       do m = 1, jmonth-1
          id = id + monday(m,ileap)
       enddo
    else
       id     = (imonth-1)*idaymo
    endif
    idaysy = id + iday

    return
  end subroutine calendar_ym2yd

  !-----------------------------------------------------------------------------
  subroutine calendar_ss2yh( &
       idate, &
       dsec   )
    implicit none

    integer,  intent(out) :: idate(6) ! yymmddhhmmss
    real(DP), intent(in)  :: dsec     ! time(second)

    integer  :: idays ! serial number of day
    real(DP) :: rsec  ! seconds in a day
    !---------------------------------------------------------------------------

    call calendar_ss2ds( idays, rsec, dsec )
    call calendar_dd2ym( idate(1), idate(2), idate(3), idays )
    call calendar_rs2hm( idate(4), idate(5), idate(6), rsec  )

    return
  end subroutine calendar_ss2yh

  !-----------------------------------------------------------------------------
  subroutine calendar_yh2ss( &
       dsec  ,               & !--- OUT
       idate                 & !--- IN
       )
    !--- calendar, date -> sec.
    implicit none

    real(DP), intent(out) :: dsec    !" time
    integer, intent(in) :: idate(6) !" yymmddhhmmss

    integer    idays             !" serial no.of day
    real(DP)    rsec              !" no. of sec. in a day

    call calendar_ym2dd &
         ( idays   , &
         idate(1), idate(2), idate(3) )

    call calendar_hm2rs &
         ( rsec    , &
         idate(4), idate(5), idate(6) )

    call calendar_ds2ss &
         ( dsec  , &
         idays , rsec   )

    return
  end subroutine calendar_yh2ss

  !-----------------------------------------------------------------------------
  subroutine calendar_ss2yd( &
       iyear ,               & !--- OUT
       idaysy,               & !--- OUT
       dsec                  & !--- IN
       )
    !--- calendar, sec. -> yydd
    implicit none

    integer, intent(out) :: iyear
    integer, intent(out) :: idaysy
    real(DP), intent(in) :: dsec

    integer :: imonth, iday
    call calendar_ss2ym &
         ( iyear , imonth, iday  , &
         dsec                    )
    call calendar_ym2yd &
         ( idaysy, &
         iyear , imonth, iday    )

    return
  end subroutine calendar_ss2yd

  !-----------------------------------------------------------------------------
  subroutine calendar_ss2ym( &
       iyear ,               & !--- OUT
       imonth,               & !--- OUT
       iday  ,               & !--- OUT
       dsec                  & !--- IN
       )
    !--- calendar, sec. -> yymmdd
    implicit none

    integer, intent(out) :: iyear
    integer, intent(out) :: imonth
    integer, intent(out) :: iday
    real(DP), intent(in) :: dsec

    integer :: idays
    real(DP) :: rsec

    call calendar_ss2ds &
         ( idays , rsec  , &
         dsec            )
    call calendar_dd2ym &
         ( iyear , imonth, iday  , &
         idays                  )

    return
  end subroutine calendar_ss2ym

  !-----------------------------------------------------------------------------
  subroutine calendar_xx2ss( &
       ddsec, &
       rtdur, &
       hunit, &
       dsec   )
    implicit none

    real(DP),         intent(out) :: ddsec
    real(DP),         intent(in)  :: rtdur
    character(len=*), intent(in)  :: hunit
    real(DP),         intent(in)  :: dsec

    character(len=10) :: hunitx
    integer :: isecmi, isechr, isecdy
    integer :: iyear, imonth, iday, ndaymo, ndayyr
    !---------------------------------------------------------------------------

    hunitx = hunit

    if    ( hunitx(1:1) == 's' .OR. hunitx(1:1) == 'S' ) then

       ddsec = real(rtdur,kind=DP)

    elseif( hunitx(1:2) == 'mi' .OR. hunitx(1:2) == 'MI' ) then

       call calendar_secmi( isecmi )
       ddsec = real(rtdur,kind=DP)*real(isecmi,kind=DP)

    elseif( hunitx(1:1) == 'h' .OR. hunitx(1:1) == 'H' ) then

       call calendar_sechr( isechr )
       ddsec = real(rtdur,kind=DP)*real(isechr,kind=DP)

    elseif( hunitx(1:1) == 'd' .OR. hunitx(1:1) == 'D' ) then

       call calendar_secdy( isecdy )
       ddsec = real(rtdur,kind=DP)*real(isecdy,kind=DP)

    elseif( hunitx(1:2) == 'mo' .OR. hunitx(1:2) == 'MO' ) then

       call calendar_ss2ym( iyear, imonth, iday, dsec )
       call calendar_daymo( ndaymo, iyear, imonth )
       call calendar_secdy( isecdy )
       ddsec = real(rtdur,kind=DP) * real(ndaymo,kind=DP) * real(isecdy,kind=DP)

    elseif( hunitx(1:1) == 'y' .OR. hunitx(1:1) == 'Y' ) then

       call calendar_ss2ym( iyear, imonth, iday, dsec )
       call calendar_dayyr( ndayyr, iyear )
       call calendar_secdy( isecdy )
       ddsec = real(rtdur,kind=DP)*real(ndayyr,kind=DP)*real(isecdy,kind=DP)

    else
       write(*,*) 'xxx cxx2ss: invalid unit : ', hunit, ' [sec] assumed'
       ddsec = rtdur
    endif

    return
  end subroutine calendar_xx2ss

  !-----------------------------------------------------------------------------
  subroutine calendar_cc2yh( &
       itime ,               & !--- OUT
       htime                 & !--- IN
       )
    !--- calendar, character -> date
    implicit none

    integer, intent(out) :: itime ( 6 )
    character(len=*), intent(in) :: htime

    integer :: i
    read ( htime, 2600 ) (itime(i),i=1,6)
2600 format( i4.4,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2 )

    return
  end subroutine calendar_cc2yh

  !-----------------------------------------------------------------------------
  !> second to character
  subroutine calendar_ss2cc( &
       htime,      &
       dsec,       &
       number_only )
    implicit none

    character(len=*), intent(out) :: htime
    real(DP),         intent(in)  :: dsec
    logical,          intent(in), optional :: number_only

    integer :: itime(6)
    !---------------------------------------------------------------------------

    call calendar_ss2yh( itime, dsec )

    if ( present(number_only) ) then
       if ( number_only ) then
          write(htime,'(I4.4,5(I2.2))') itime(1), itime(2), itime(3), &
                                        itime(4), itime(5), itime(6)
          return
       endif
    endif

    write(htime,'(I4.4,5(A,I2.2))') itime(1), '/', itime(2), '/', itime(3), &
                               '-', itime(4), ':', itime(5), ':', itime(6)

    return
  end subroutine calendar_ss2cc

  !-----------------------------------------------------------------------------
  function ocleap( &
       iyear       & !--- IN
       )
    !--- calendar :bissextile or not
    implicit none

    logical :: ocleap
    integer, intent(in) :: iyear

    integer :: iy, iycen, icent

    iy     = mod(iyear,4)
    iycen  = mod(iyear,100)
    icent  = mod(iyear/100,4)
    if ( iy == 0 .AND. ( iycen /= 0 .OR. icent == 0 ) ) then
       ocleap = .true.
    else
       ocleap = .false.
    endif

    return
  end function ocleap

  !-----------------------------------------------------------------------------
  subroutine calendar_ssaft( &
       dseca ,               & !--- OUT
       dsec  ,               & !--- IN
       raftr ,               & !--- IN
       hunit                 & !--- IN
       )
    !--- calendar, time advancing
    implicit none

    real(DP), intent(out) :: dseca
    real(DP), intent(in) :: dsec
    real(DP), intent(in) :: raftr
    character(len=*), intent(in) :: hunit

    integer :: idays, iyear, imonth, iday
    real(DP) :: rsec
    real(DP) :: ddtime

    if ( hunit(1:1) == 'y' .OR. hunit(1:1) == 'Y' &
         .OR. hunit(1:1) == 'mo' .OR. hunit(1:2) == 'MO' ) then
       call calendar_ss2ds &
            ( idays , rsec  , &
            dsec            )
       call calendar_dd2ym &
            ( iyear , imonth, iday  , &
            idays                  )
       if ( hunit(1:1) == 'y' .OR. hunit(1:1) == 'Y' ) then
          iyear  = iyear  + int(raftr)
       elseif( hunit(1:2) == 'mo' .OR. hunit(1:2) == 'MO' ) then
          imonth = imonth + int(raftr)
       endif
       call calendar_ym2dd &
            ( idays , &
            iyear , imonth, iday   )
       call calendar_ds2ss &
            ( dseca , &
            idays , rsec    )
    else
       call calendar_xx2ss( ddtime, raftr, hunit, dsec )
       dseca = dsec + ddtime
    endif

    return
  end subroutine calendar_ssaft

  !-----------------------------------------------------------------------------
  function   calendar_ointvl( &
       dtime ,                & !--- IN
       dtprev ,               & !--- IN
       dtorgn,                & !--- IN
       rintv ,                & !--- IN
       htunit                 & !--- IN
       )
    !--- time step passed ?
    implicit none

    logical :: calendar_ointvl
    real(DP), intent(in) :: dtime
    real(DP), intent(in) :: dtprev
    real(DP), intent(in) :: dtorgn
    real(DP), intent(in) :: rintv
    character(len=*), intent(in) :: htunit

    real(DP) :: ddtime
    character(len=5) :: hunit
    integer :: iyear, imon, iday, iyearp, imonp, idayp
    integer :: iy, imo
    integer :: nmonyr, ndayyr, ndaymo
    real(DP) :: ry, rmo

    hunit = htunit

    if ( dtime == dtprev ) then
       calendar_ointvl = .true.
       return
    endif

    calendar_ointvl = .false.

    call calendar_ss2ym ( iyear, imon, iday, dtime  )
    call calendar_ss2ym ( iyearp, imonp, idayp, dtprev )
    call calendar_xx2ss ( ddtime, rintv, hunit, dtime  )

    if ( dtime >= dtorgn ) then
       if      ( hunit(1:1) == 'y' .OR. hunit(1:1) == 'Y' ) then
          call calendar_monyr( nmonyr, iyear )
          call calendar_dayyr( ndayyr, iyear )
          ry = real( iyear-iyearp, kind=DP )&
             + real( imon-imonp, kind=DP ) / real( nmonyr, kind=DP ) &
             + real( iday-idayp, kind=DP ) / real( ndayyr, kind=DP )
          if ( ry >= rintv ) then
             calendar_ointvl = .true.
          endif
       elseif( hunit(1:2) == 'mo' .OR. hunit(1:2) == 'MO' ) then
          imo = 0
          do iy = iyearp, iyear-1
             call calendar_monyr( nmonyr, iy )
             imo = imo + nmonyr
          enddo
          call calendar_daymo( ndaymo, iyear, imon )
          rmo = real( imon-imonp+imo, kind=DP ) &
              + real( iday-idayp, kind=DP ) / real( ndaymo, kind=DP )
          if ( rmo >= rintv ) then
             calendar_ointvl = .true.
          endif
       elseif(      calendar_dgaus((dtime -dtorgn)/ddtime)  &
            > calendar_dgaus((dtprev-dtorgn)/ddtime) ) then
          calendar_ointvl = .true.
       endif
    endif

    return
  end function calendar_ointvl

  !-----------------------------------------------------------------------------
  function calendar_dgaus( &
        dx                 & !--- IN
        )
    !--- dicard gaussian
    implicit none
    real(DP) :: calendar_dgaus
    real(DP), intent(in) :: dx

    calendar_dgaus = aint(dx) + aint(dx - aint(dx) + 1.0_DP) - 1.0_DP

  end function calendar_dgaus

end module mod_calendar
