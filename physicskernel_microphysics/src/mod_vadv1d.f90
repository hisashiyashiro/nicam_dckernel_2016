!-------------------------------------------------------------------------------
!> Module vertical advection
!!
!! @par Description
!!         This module is for 1-d vertical advection
!!
!! @author NICAM developers
!<
!-------------------------------------------------------------------------------
module mod_vadv1d
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
  public :: vadv1d_prep
  public :: vadv1d_getflux_new

  !-----------------------------------------------------------------------------
  !
  !++ Private procedures
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private variables
  !

  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  subroutine vadv1d_prep( &
       ijdim,     &
       kdim,      &
       kmin,      &
       kmax,      &
       dz,        &
       zh,        &
       wp,        &
       zdis,      &
       kcell,     &
       kcell_max, &
       kcell_min, &
       dt         )
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: kmin
    integer,  intent(in)  :: kmax
    real(RP), intent(in)  :: dz       (kdim)
    real(RP), intent(in)  :: zh       (kdim)
    real(RP), intent(in)  :: wp       (ijdim,kdim)
    real(RP), intent(out) :: zdis     (ijdim,kdim) ! [bugfix] H.Yashiro 20120606
    integer,  intent(out) :: kcell    (ijdim,kdim)
    integer,  intent(out) :: kcell_max(kdim)
    integer,  intent(out) :: kcell_min(kdim)
    real(RP), intent(in)  :: dt

    real(RP) :: wh(ijdim,kdim)
    real(RP) :: zzmax, zzmin

    integer  :: ij, k, k2
    !---------------------------------------------------------------------------

    ! vetical velocity at the half level
    do k  = kmin+1, kmax
    do ij = 1, ijdim
       wh(ij,k) = 0.5_RP * ( wp(ij,k-1)+wp(ij,k) )
    enddo
    enddo

    ! bottom boundary for wh
    ! top    boundary for wh : same as inner region
    do ij = 1, ijdim
      wh(ij,kmin  ) = wp(ij,kmin  )
      wh(ij,kmin-1) = wp(ij,kmin-1)
      wh(ij,kmax+1) = wp(ij,kmax  )
    enddo

    ! calculation of distance of cell wall during dt
    do k = kmin+1, kmax
    do ij = 1, ijdim
       zdis(ij,k) = dt    * wh(ij,k)                                                              &
                  - dt**2 * wh(ij,k) *     ( wh(ij,k+1)-wh(ij,k-1) ) / ( dz(k-1)+dz(k) ) / 2.0_RP &
                  + dt**3 * wh(ij,k) * ( ( ( wh(ij,k+1)-wh(ij,k-1) ) / ( dz(k-1)+dz(k) ) )**2     &
                                       + wh(ij,k) * ( ( ( wh(ij,k+1)-wh(ij,k) ) / dz(k)           &
                                                      - ( wh(ij,k)-wh(ij,k-1) ) / dz(k-1) )       &
                                                    / ( dz(k-1)+dz(k) ) * 2.0_RP )                &
                                       ) / 6.0_RP
    enddo
    enddo

    ! bottom and top boundary for zdis
    do ij = 1, ijdim
       zdis(ij,kmin-1) = 0.0_RP
       zdis(ij,kmin  ) = dt    * wh(ij,kmin  ) &
                       - dt**2 * wh(ij,kmin  ) * ( wh(ij,kmin+1)-wh(ij,kmin) ) / dz(kmin) / 2.0_RP
       zdis(ij,kmax+1) = dt    * wh(ij,kmax+1) &
                       - dt**2 * wh(ij,kmax+1) * ( wh(ij,kmax+1)-wh(ij,kmax) ) / dz(kmax) / 2.0_RP
    enddo

    ! calculation of kcell
    ! top boundary: rigid [kcell(:,kmax+1) = kmax+1]
    do k  = 1, kdim
    do ij = 1, ijdim
       kcell(ij,k)  = k
       kcell_min(k) = k
       kcell_max(k) = k
    enddo
    enddo

    ! setup limiter of max and min of kcell
    do k = kmin, kmax
       zzmax = maxval( zdis(:,k) )
       zzmin = minval( zdis(:,k) )

       if ( zzmax > 0.0_RP ) then
          do k2 = k, kmin, -1
             if (       zh(k2)   <= zh(k)-zzmax &
                  .AND. zh(k2+1) >  zh(k)-zzmax ) then

                kcell_min(k) = k2
                exit

             endif
          enddo
       endif

       if ( zzmin < 0.0_RP ) then
          do k2 = k, kmax
             if (       zh(k2)   <= zh(k)-zzmin &
                  .AND. zh(k2+1) >  zh(k)-zzmin ) then

                kcell_max(k) = k2
                exit

             endif
          enddo
       endif
    enddo

    ! determine the kcell at each point.
    do k  = kmin, kmax
    do ij = 1, ijdim
       if (       kcell_min(k) == k &
            .AND. kcell_max(k) == k ) then

          kcell(ij,k) = k

       else

          kcell(ij,k) = 0

          do k2 = kcell_min(k), kcell_max(k)
             kcell(ij,k) = max ( kcell(ij,k), &
                                 int( k2  * sign ( 1.0_RP, (zh(k)-zdis(ij,k))-zh(k2)   ) &
                                          * sign ( 1.0_RP, zh(k2+1)-(zh(k)-zdis(ij,k)) ) ) &
                               )
          enddo

       endif
    enddo
    enddo

    do k  = 1, kdim
    do ij = 1, ijdim
       if ( kcell(ij,k) == 0 ) then
          kcell(ij,k) = kmin
       endif

       if ( kcell(ij,k) == kmax+1 ) then
          kcell(ij,k) = kmax
       endif
    enddo
    enddo

    return
  end subroutine vadv1d_prep

  !-----------------------------------------------------------------------------
  subroutine vadv1d_getflux_new( &
       ijdim,     &
       kdim,      &
       kmin,      &
       kmax,      &
       dz,        &
       rhof,      &
       zdis0,     &
       kcell,     &
       kcell_max, &
       kcell_min, &
       frhof      )
!ESC!    use mod_const, only: &
!ESC!       CONST_EPS
    implicit none

    integer,  intent(in)  :: ijdim
    integer,  intent(in)  :: kdim
    integer,  intent(in)  :: kmin
    integer,  intent(in)  :: kmax
    real(RP), intent(in)  :: dz       (kdim)
    real(RP), intent(in)  :: rhof     (ijdim,kdim)
    real(RP), intent(in)  :: zdis0    (ijdim,kdim) ! [bugfix] H.Yashiro 20120606
    integer,  intent(in)  :: kcell    (ijdim,kdim)
    integer,  intent(in)  :: kcell_max(kdim)
    integer,  intent(in)  :: kcell_min(kdim)
    real(RP), intent(out) :: frhof    (ijdim,kdim)

    real(RP) :: zdis(ijdim,kdim)
    real(RP) :: fact

    integer  :: ij, k, k2
    !---------------------------------------------------------------------------

    !------ integration in the integer cells
    !$omp parallel &
    !$omp default(none) &
    !$omp shared(frhof,kcell_min,kcell_max,zdis,zdis0,dz,rhof,kcell) &
    !$omp shared(kdim,kmin,kmax,ijdim,CONST_EPS) &
    !$omp private(k,ij,k2,fact)

    do k  = 1, kdim
       !$omp do
       do ij = 1, ijdim
          frhof(ij,k) = 0.0_RP
       enddo
       !$omp end do nowait
    enddo


    do k = kmin, kmax
       if (       kcell_min(k) == k &
            .AND. kcell_max(k) == k ) then

          !$omp do
          do ij = 1, ijdim
             zdis(ij,k) = zdis0(ij,k)
          enddo
          !$omp end do nowait
       else
          do k2 = kcell_min(k), kcell_max(k)

             ! sum up over k2 = kcell(ij,k)+1, k-1 : if w > 0
             ! or          k2 = k, kcell(ij,k)-1   : if w < 0
             !$omp do
             do ij = 1, ijdim
                fact = dz(k2) * 0.25_RP &
                     * ( ( sign(1,k2-(kcell(ij,k)+1)) + 1.0_RP ) * ( sign(1,(k          -1)-k2) + 1.0_RP ) &
                       - ( sign(1,k2-k              ) + 1.0_RP ) * ( sign(1,(kcell(ij,k)-1)-k2) + 1.0_RP ) )


                frhof(ij,k) = frhof(ij,k) + rhof(ij,k2) * fact
                zdis (ij,k) = zdis0(ij,k) - fact
             enddo
             !$omp end do nowait
          enddo
       endif

       !ocl XFILL
       !$omp do
       do ij = 1, ijdim
          frhof(ij,k) = frhof(ij,k) + rhof(ij,kcell(ij,k)) * zdis(ij,k)
       enddo
       !$omp end do nowait
    enddo

    do k  = 1, kdim
       !$omp do
       do ij = 1, ijdim
          frhof(ij,k) = frhof(ij,k) * ( 0.5_RP + sign(0.5_RP, abs(frhof(ij,k)) - CONST_EPS ) ) ! small negative filter
       enddo
       !$omp end do nowait
    enddo

    !$omp end parallel

    return
  end subroutine vadv1d_getflux_new

end module mod_vadv1d
