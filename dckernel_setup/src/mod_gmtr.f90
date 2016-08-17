!-------------------------------------------------------------------------------
!>
!! Geometrics module
!!
!! @par Description
!!          In this module, the geometrics of the icosahedral grid such as
!!          area are calculated.
!!
!! @author  H.Tomita
!!
!! @par History
!! @li      2004-02-17 (H.Tomita)  Imported from igdc-4.33
!! @li      2010-06-08 (S.Iga)     a new grid is implemented
!! @li      2011-07-22 (T.Ohno)    a new grid is implemented
!! @li      2011-08-18 (T.Ohno)    bugfix
!!
!<
module mod_gmtr
  !-----------------------------------------------------------------------------
  !
  !++ Used modules
  !
  use mod_precision
  use mod_misc
!ESC!  use mod_stdio
!ESC!  use mod_adm, only: &
!ESC!     TI,      &
!ESC!     TJ,      &
!ESC!     AI,      &
!ESC!     AIJ,     &
!ESC!     AJ,      &
!ESC!     K0,   &
!ESC!     ADM_lall,    &
!ESC!     ADM_lall_pl, &
!ESC!     ADM_gall,    &
!ESC!     ADM_gall_pl
  !-----------------------------------------------------------------------------
  implicit none
  private
  !-----------------------------------------------------------------------------
  !
  !++ Public procedure
  !
  !-----------------------------------------------------------------------------
  !
  !++ Public parameters & variables
  !
  !-----------------------------------------------------------------------------
  !
  !++ Private procedure
  !
  public :: GMTR_p_setup
  public :: GMTR_t_setup
  public :: GMTR_a_setup
  private :: GMTR_TNvec

  !-----------------------------------------------------------------------------
  !
  !++ Private parameters & variables
  !
  !-----------------------------------------------------------------------------
contains
  !-----------------------------------------------------------------------------
  !> calc geometrical information for cell point
  subroutine GMTR_p_setup( &
       GRD_x,  GRD_x_pl,  &
       GRD_xt, GRD_xt_pl, &
       GRD_s,  GRD_s_pl,  &
       GMTR_p, GMTR_p_pl, &
       GRD_rscale         )
!ESC!    use mod_adm, only: &
!ESC!       ADM_nxyz,     &
!ESC!       ADM_have_pl,  &
!ESC!       ADM_have_sgp, &
!ESC!       ADM_vlink,    &
!ESC!       ADM_gmin,     &
!ESC!       ADM_gmax,     &
!ESC!       ADM_gslf_pl
!ESC!    use mod_grd, only: &
!ESC!       GRD_XDIR,     &
!ESC!       GRD_YDIR,     &
!ESC!       GRD_ZDIR,     &
!ESC!       GRD_LON,      &
!ESC!       GRD_grid_type
    use mod_vector, only: &
       VECTR_triangle,      &
       VECTR_triangle_plane
    implicit none

    real(RP), intent(in)  :: GRD_x    (ADM_gall   ,K0,ADM_lall   ,      ADM_nxyz)
    real(RP), intent(in)  :: GRD_x_pl (ADM_gall_pl,K0,ADM_lall_pl,      ADM_nxyz)
    real(RP), intent(in)  :: GRD_xt   (ADM_gall   ,K0,ADM_lall   ,TI:TJ,ADM_nxyz)
    real(RP), intent(in)  :: GRD_xt_pl(ADM_gall_pl,K0,ADM_lall_pl,      ADM_nxyz)
    real(RP), intent(in)  :: GRD_s    (ADM_gall   ,K0,ADM_lall   ,      2)
    real(RP), intent(in)  :: GRD_s_pl (ADM_gall_pl,K0,ADM_lall_pl,      2)
    real(RP), intent(out) :: GMTR_p   (ADM_gall   ,K0,ADM_lall   ,GMTR_p_nmax)
    real(RP), intent(out) :: GMTR_p_pl(ADM_gall_pl,K0,ADM_lall_pl,GMTR_p_nmax)
    real(RP), intent(in)  :: GRD_rscale

    real(RP) :: wk   (ADM_nxyz,0:7,ADM_gall)
    real(RP) :: wk_pl(ADM_nxyz,0:ADM_vlink+1)

    real(RP) :: area
    real(RP) :: cos_lambda, sin_lambda

    integer  :: ij
    integer  :: ip1j, ijp1, ip1jp1
    integer  :: im1j, ijm1, im1jm1

    integer  :: i, j, l, d, v, n
    !---------------------------------------------------------------------------

    write(ADM_LOG_FID,*) '*** setup metrics for hexagonal/pentagonal mesh'

    GMTR_p   (:,:,:,:)   = 0.0_RP
    GMTR_p_pl(:,:,:,:)   = 0.0_RP

    do l = 1, ADM_lall
       do j = ADM_gmin, ADM_gmax
       do i = ADM_gmin, ADM_gmax
          ij     = suf(i  ,j  )
          ip1j   = suf(i+1,j  )
          ip1jp1 = suf(i+1,j+1)
          ijp1   = suf(i  ,j+1)
          im1j   = suf(i-1,j  )
          im1jm1 = suf(i-1,j-1)
          ijm1   = suf(i  ,j-1)

          !--- prepare 1 center and 6 vertices
          do d = 1, ADM_nxyz
             wk(d,0,ij) = GRD_x(ij,k0,l,d)

             wk(d,1,ij) = GRD_xt(ijm1  ,k0,l,TJ,d)
             wk(d,2,ij) = GRD_xt(ij    ,k0,l,TI,d)
             wk(d,3,ij) = GRD_xt(ij    ,k0,l,TJ,d)
             wk(d,4,ij) = GRD_xt(im1j  ,k0,l,TI,d)
             wk(d,5,ij) = GRD_xt(im1jm1,k0,l,TJ,d)
             wk(d,6,ij) = GRD_xt(im1jm1,k0,l,TI,d)
             wk(d,7,ij) = wk(d,1,ij)
          enddo
       enddo ! i loop
       enddo ! j loop

       if ( ADM_have_sgp(l) ) then ! pentagon
          wk(:,6,suf(ADM_gmin,ADM_gmin)) = wk(:,1,suf(ADM_gmin,ADM_gmin))
          wk(:,7,suf(ADM_gmin,ADM_gmin)) = wk(:,1,suf(ADM_gmin,ADM_gmin))
       endif

       !--- calc control area
       if ( GRD_grid_type == 'ON_PLANE' ) then
          do j = ADM_gmin, ADM_gmax
          do i = ADM_gmin, ADM_gmax
             ij = suf(i,j)

             area = 0.0_RP
             do v = 1, 6
                area = area + VECTR_triangle_plane( wk(:,0,ij), wk(:,v,ij), wk(:,v+1,ij) )
             enddo

             GMTR_p(ij,k0,l,P_AREA)  = area
             GMTR_p(ij,k0,l,P_RAREA) = 1.0_RP / GMTR_p(ij,k0,l,P_AREA)

          enddo ! i loop
          enddo ! j loop
       else
          do j = ADM_gmin, ADM_gmax
          do i = ADM_gmin, ADM_gmax
             ij = suf(i,j)

             wk(:,:,ij) = wk(:,:,ij) / GRD_rscale

             area = 0.0_RP
             do v = 1, 6
                area = area + VECTR_triangle( wk(:,0,ij), wk(:,v,ij), wk(:,v+1,ij), GMTR_polygon_type, GRD_rscale )
             enddo

             GMTR_p(ij,k0,l,P_AREA)  = area
             GMTR_p(ij,k0,l,P_RAREA) = 1.0_RP / GMTR_p(ij,k0,l,P_AREA)

          enddo ! i loop
          enddo ! j loop
       endif

       !--- calc coefficient between xyz <-> latlon
       if ( GRD_grid_type == 'ON_PLANE' ) then
          GMTR_p(:,k0,l,P_IX) = 1.0_RP
          GMTR_p(:,k0,l,P_IY) = 0.0_RP
          GMTR_p(:,k0,l,P_IZ) = 0.0_RP
          GMTR_p(:,k0,l,P_JX) = 0.0_RP
          GMTR_p(:,k0,l,P_JY) = 1.0_RP
          GMTR_p(:,k0,l,P_JZ) = 0.0_RP
       else
          do j = ADM_gmin, ADM_gmax
          do i = ADM_gmin, ADM_gmax
             ij = suf(i,j)

             sin_lambda = sin( GRD_s(ij,k0,l,GRD_LON) )
             cos_lambda = cos( GRD_s(ij,k0,l,GRD_LON) )

             GMTR_p(ij,k0,l,P_IX) = -sin_lambda
             GMTR_p(ij,k0,l,P_IY) =  cos_lambda
             GMTR_p(ij,k0,l,P_IZ) = 0.0_RP
             GMTR_p(ij,k0,l,P_JX) = -( GRD_x(ij,k0,l,ZDIR) * cos_lambda ) / GRD_rscale
             GMTR_p(ij,k0,l,P_JY) = -( GRD_x(ij,k0,l,ZDIR) * sin_lambda ) / GRD_rscale
             GMTR_p(ij,k0,l,P_JZ) =  ( GRD_x(ij,k0,l,XDIR) * cos_lambda &
                                     + GRD_x(ij,k0,l,YDIR) * sin_lambda ) / GRD_rscale
          enddo ! i loop
          enddo ! j loop
       endif
    enddo ! l loop

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl
          !--- prepare 1 center and * vertices
          do d = 1, ADM_nxyz
             wk_pl(d,0) = GRD_x_pl(n,k0,l,d)
             do v = 1, ADM_vlink ! (ICO=5)
                wk_pl(d,v) = GRD_xt_pl(v+1,k0,l,d)
             enddo
             wk_pl(d,ADM_vlink+1) = wk_pl(d,1)
          enddo

          wk_pl(:,:) = wk_pl(:,:) / GRD_rscale

          !--- calc control area
          area = 0.0_RP
          do v = 1, ADM_vlink ! (ICO=5)
             area = area + VECTR_triangle( wk_pl(:,0), wk_pl(:,v), wk_pl(:,v+1), GMTR_polygon_type, GRD_rscale )
          enddo

          GMTR_p_pl(n,k0,l,P_AREA)  = area
          GMTR_p_pl(n,k0,l,P_RAREA) = 1.0_RP / GMTR_p_pl(n,k0,l,P_AREA)

          !--- calc coefficient between xyz <-> latlon
          sin_lambda = sin( GRD_s_pl(n,k0,l,GRD_LON) )
          cos_lambda = cos( GRD_s_pl(n,k0,l,GRD_LON) )

          GMTR_p_pl(n,k0,l,P_IX) = -sin_lambda
          GMTR_p_pl(n,k0,l,P_IY) =  cos_lambda
          GMTR_p_pl(n,k0,l,P_IZ) = 0.0_RP
          GMTR_p_pl(n,k0,l,P_JX) = -( GRD_x_pl(n,k0,l,ZDIR) * cos_lambda ) / GRD_rscale
          GMTR_p_pl(n,k0,l,P_JY) = -( GRD_x_pl(n,k0,l,ZDIR) * sin_lambda ) / GRD_rscale
          GMTR_p_pl(n,k0,l,P_JZ) =  ( GRD_x_pl(n,k0,l,XDIR) * cos_lambda &
                                    + GRD_x_pl(n,k0,l,YDIR) * sin_lambda ) / GRD_rscale
       enddo ! l loop
    endif

    return
  end subroutine GMTR_p_setup

  !-----------------------------------------------------------------------------
  !> calc geometrical information for cell vertex (triangle)
  subroutine GMTR_t_setup( &
       GRD_x,  GRD_x_pl,  &
       GRD_xt, GRD_xt_pl, &
       GMTR_t, GMTR_t_pl, &
       GRD_rscale         )
!ESC!    use mod_adm, only: &
!ESC!       ADM_nxyz,     &
!ESC!       ADM_have_pl,  &
!ESC!       ADM_have_sgp, &
!ESC!       ADM_gmin,     &
!ESC!       ADM_gmax,     &
!ESC!       ADM_gslf_pl,  &
!ESC!       ADM_gmin_pl,  &
!ESC!       ADM_gmax_pl
!ESC!    use mod_grd, only: &
!ESC!       GRD_grid_type
    use mod_vector, only: &
       VECTR_triangle,      &
       VECTR_triangle_plane
    implicit none

    real(RP), intent(in)  :: GRD_x    (ADM_gall   ,K0,ADM_lall   ,      ADM_nxyz)
    real(RP), intent(in)  :: GRD_x_pl (ADM_gall_pl,K0,ADM_lall_pl,      ADM_nxyz)
    real(RP), intent(in)  :: GRD_xt   (ADM_gall   ,K0,ADM_lall   ,TI:TJ,ADM_nxyz)
    real(RP), intent(in)  :: GRD_xt_pl(ADM_gall_pl,K0,ADM_lall_pl,      ADM_nxyz)
    real(RP), intent(out) :: GMTR_t   (ADM_gall   ,K0,ADM_lall   ,TI:TJ,GMTR_t_nmax)
    real(RP), intent(out) :: GMTR_t_pl(ADM_gall_pl,K0,ADM_lall_pl,      GMTR_t_nmax)
    real(RP), intent(in)  :: GRD_rscale

    real(RP) :: wk   (ADM_nxyz,0:3,ADM_gall,TI:TJ)
    real(RP) :: wk_pl(ADM_nxyz,0:3)

    real(RP) :: area, area1, area2, area3

    integer  :: ij
    integer  :: ip1j, ijp1, ip1jp1

    integer  :: i, j, l, d, v, n, t
    !---------------------------------------------------------------------------

    write(ADM_LOG_FID,*) '*** setup metrics for triangle mesh'

    GMTR_t   (:,:,:,:,:) = 0.0_RP
    GMTR_t_pl(:,:,:,:)   = 0.0_RP

    do l = 1,ADM_lall
       do j = ADM_gmin-1, ADM_gmax
       do i = ADM_gmin-1, ADM_gmax
          ij     = suf(i  ,j  )
          ip1j   = suf(i+1,j  )
          ip1jp1 = suf(i+1,j+1)
          ijp1   = suf(i  ,j+1)

          !--- prepare 1 center and 3 vertices for 2 triangles
          do d = 1, ADM_nxyz
             wk(d,0,ij,TI) = GRD_xt(ij,k0,l,TI,d)

             wk(d,1,ij,TI) = GRD_x(ij    ,k0,l,d)
             wk(d,2,ij,TI) = GRD_x(ip1j  ,k0,l,d)
             wk(d,3,ij,TI) = GRD_x(ip1jp1,k0,l,d)

             wk(d,0,ij,TJ) = GRD_xt(ij,k0,l,TJ,d)

             wk(d,1,ij,TJ) = GRD_x(ij    ,k0,l,d)
             wk(d,2,ij,TJ) = GRD_x(ip1jp1,k0,l,d)
             wk(d,3,ij,TJ) = GRD_x(ijp1  ,k0,l,d)
          enddo
       enddo
       enddo

       !--- treat unused triangle
       wk(:,:,suf(ADM_gmax,ADM_gmin-1),TI) = wk(:,:,suf(ADM_gmax,ADM_gmin-1),TJ)
       wk(:,:,suf(ADM_gmin-1,ADM_gmax),TJ) = wk(:,:,suf(ADM_gmin-1,ADM_gmax),TI)

       if ( ADM_have_sgp(l) ) then ! pentagon
          wk(:,:,suf(ADM_gmin-1,ADM_gmin-1),TI) = wk(:,:,suf(ADM_gmin,ADM_gmin-1),TJ)
       endif

       if ( GRD_grid_type == 'ON_PLANE' ) then
          do t = TI,TJ
          do j = ADM_gmin-1, ADM_gmax
          do i = ADM_gmin-1, ADM_gmax
             ij = suf(i,j)

             area1 = VECTR_triangle_plane( wk(:,0,ij,t), wk(:,2,ij,t), wk(:,3,ij,t) )
             area2 = VECTR_triangle_plane( wk(:,0,ij,t), wk(:,3,ij,t), wk(:,1,ij,t) )
             area3 = VECTR_triangle_plane( wk(:,0,ij,t), wk(:,1,ij,t), wk(:,2,ij,t) )

             area = area1 + area2 + area3

             GMTR_t(ij,k0,l,t,T_AREA)  = area
             GMTR_t(ij,k0,l,t,T_RAREA) = 1.0_RP / area

             GMTR_t(ij,k0,l,t,W1)    = area1 / area
             GMTR_t(ij,k0,l,t,W2)    = area2 / area
             GMTR_t(ij,k0,l,t,W3)    = area3 / area
          enddo
          enddo
          enddo
       else
          do t = TI,TJ
          do j = ADM_gmin-1, ADM_gmax
          do i = ADM_gmin-1, ADM_gmax
             ij = suf(i,j)

             wk(:,:,ij,t) = wk(:,:,ij,t) / GRD_rscale

             area1 = VECTR_triangle( wk(:,0,ij,t), wk(:,2,ij,t), wk(:,3,ij,t), GMTR_polygon_type, GRD_rscale )
             area2 = VECTR_triangle( wk(:,0,ij,t), wk(:,3,ij,t), wk(:,1,ij,t), GMTR_polygon_type, GRD_rscale )
             area3 = VECTR_triangle( wk(:,0,ij,t), wk(:,1,ij,t), wk(:,2,ij,t), GMTR_polygon_type, GRD_rscale )

             area = area1 + area2 + area3

             GMTR_t(ij,k0,l,t,T_AREA)  = area
             GMTR_t(ij,k0,l,t,T_RAREA) = 1.0_RP / area

             GMTR_t(ij,k0,l,t,W1)    = area1 / area
             GMTR_t(ij,k0,l,t,W2)    = area2 / area
             GMTR_t(ij,k0,l,t,W3)    = area3 / area
          enddo
          enddo
          enddo
       endif

    enddo

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1,ADM_lall_pl
          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v + 1
             if( ijp1 == ADM_gmax_pl+1 ) ijp1 = ADM_gmin_pl

             do d = 1, ADM_nxyz
                wk_pl(d,0) = GRD_xt_pl(ij,k0,l,d)

                wk_pl(d,1) = GRD_x_pl(n   ,k0,l,d)
                wk_pl(d,2) = GRD_x_pl(ij  ,k0,l,d)
                wk_pl(d,3) = GRD_x_pl(ijp1,k0,l,d)
             enddo

             wk_pl(:,:) = wk_pl(:,:) / GRD_rscale

             area1 = VECTR_triangle( wk_pl(:,0), wk_pl(:,2), wk_pl(:,3), GMTR_polygon_type, GRD_rscale )
             area2 = VECTR_triangle( wk_pl(:,0), wk_pl(:,3), wk_pl(:,1), GMTR_polygon_type, GRD_rscale )
             area3 = VECTR_triangle( wk_pl(:,0), wk_pl(:,1), wk_pl(:,2), GMTR_polygon_type, GRD_rscale )

             area = area1 + area2 + area3

             GMTR_t_pl(ij,k0,l,T_AREA)  = area
             GMTR_t_pl(ij,k0,l,T_RAREA) = 1.0_RP / area

             GMTR_t_pl(ij,k0,l,W1)    = area1 / area
             GMTR_t_pl(ij,k0,l,W2)    = area2 / area
             GMTR_t_pl(ij,k0,l,W3)    = area3 / area
          enddo
       enddo
    endif

    return
  end subroutine GMTR_t_setup

  !-----------------------------------------------------------------------------
  !> calc geometrical information for cell arc
  subroutine GMTR_a_setup( &
       GRD_x,  GRD_x_pl,  &
       GRD_xt, GRD_xt_pl, &
       GMTR_a, GMTR_a_pl, &
       GRD_rscale         )
!ESC!    use mod_adm, only: &
!ESC!       ADM_nxyz,     &
!ESC!       ADM_have_pl,  &
!ESC!       ADM_have_sgp, &
!ESC!       ADM_gmin,     &
!ESC!       ADM_gmax,     &
!ESC!       ADM_gslf_pl,  &
!ESC!       ADM_gmin_pl,  &
!ESC!       ADM_gmax_pl
!ESC!    use mod_grd, only: &
!ESC!       GRD_grid_type
    implicit none

    real(RP), intent(in)  :: GRD_x    (ADM_gall   ,K0,ADM_lall   ,      ADM_nxyz)
    real(RP), intent(in)  :: GRD_x_pl (ADM_gall_pl,K0,ADM_lall_pl,      ADM_nxyz)
    real(RP), intent(in)  :: GRD_xt   (ADM_gall   ,K0,ADM_lall   ,TI:TJ,ADM_nxyz)
    real(RP), intent(in)  :: GRD_xt_pl(ADM_gall_pl,K0,ADM_lall_pl,      ADM_nxyz)
    real(RP), intent(out) :: GMTR_a   (ADM_gall   ,K0,ADM_lall   ,AI:AJ,GMTR_a_nmax   )
    real(RP), intent(out) :: GMTR_a_pl(ADM_gall_pl,K0,ADM_lall_pl,      GMTR_a_nmax_pl)
    real(RP), intent(in)  :: GRD_rscale

    real(RP) :: wk   (ADM_nxyz,2,ADM_gall)
    real(RP) :: wk_pl(ADM_nxyz,2)

    real(RP) :: Tvec(3), Nvec(3)

    integer  :: ij
    integer  :: ip1j, ijp1, ip1jp1
    integer  :: im1j, ijm1

    integer  :: i, j, l, d, v, n
    !---------------------------------------------------------------------------

    write(ADM_LOG_FID,*) '*** setup metrics for cell arcs'

    GMTR_a   (:,:,:,:,:) = 0.0_RP
    GMTR_a_pl(:,:,:,:)   = 0.0_RP

    !--- Triangle
    do l = 1, ADM_lall

       !--- AI
       do j = ADM_gmin-1, ADM_gmax+1
       do i = ADM_gmin-1, ADM_gmax
          ij   = suf(i  ,j  )
          ip1j = suf(i+1,j  )

          do d = 1, ADM_nxyz
             wk(d,1,ij) = GRD_x(ij  ,k0,l,d)
             wk(d,2,ij) = GRD_x(ip1j,k0,l,d)
          enddo
       enddo
       enddo

       ! treat arc of unused triangle
       wk(:,1,suf(ADM_gmax  ,ADM_gmin-1)) = GRD_x(suf(ADM_gmax  ,ADM_gmin-1),k0,l,:)
       wk(:,2,suf(ADM_gmax  ,ADM_gmin-1)) = GRD_x(suf(ADM_gmax  ,ADM_gmin  ),k0,l,:)
       wk(:,1,suf(ADM_gmin-1,ADM_gmax+1)) = GRD_x(suf(ADM_gmin  ,ADM_gmax+1),k0,l,:)
       wk(:,2,suf(ADM_gmin-1,ADM_gmax+1)) = GRD_x(suf(ADM_gmin  ,ADM_gmax  ),k0,l,:)

       if ( ADM_have_sgp(l) ) then ! pentagon
          wk(:,1,suf(ADM_gmin-1,ADM_gmin-1)) = GRD_x(suf(ADM_gmin  ,ADM_gmin-1),k0,l,:)
          wk(:,2,suf(ADM_gmin-1,ADM_gmin-1)) = GRD_x(suf(ADM_gmin+1,ADM_gmin  ),k0,l,:)
       endif

       do j = ADM_gmin-1, ADM_gmax+1
       do i = ADM_gmin-1, ADM_gmax
          ij = suf(i,j)

          call GMTR_TNvec( Tvec(:), Nvec(:),                            & ! [OUT]
                           wk(:,1,ij), wk(:,2,ij),                      & ! [IN]
                           GRD_grid_type, GMTR_polygon_type, GRD_rscale ) ! [IN]

          GMTR_a(ij,k0,l,AI,TNX) = Nvec(1)
          GMTR_a(ij,k0,l,AI,TNY) = Nvec(2)
          GMTR_a(ij,k0,l,AI,TNZ) = Nvec(3)
          GMTR_a(ij,k0,l,AI,TTX) = Tvec(1)
          GMTR_a(ij,k0,l,AI,TTY) = Tvec(2)
          GMTR_a(ij,k0,l,AI,TTZ) = Tvec(3)
       enddo
       enddo

       !--- AIJ
       do j = ADM_gmin-1, ADM_gmax
       do i = ADM_gmin-1, ADM_gmax
          ij     = suf(i  ,j  )
          ip1jp1 = suf(i+1,j+1)

          do d = 1, ADM_nxyz
             wk(d,1,ij) = GRD_x(ij    ,k0,l,d)
             wk(d,2,ij) = GRD_x(ip1jp1,k0,l,d)
          enddo
       enddo
       enddo

       do j = ADM_gmin-1, ADM_gmax
       do i = ADM_gmin-1, ADM_gmax
          ij = suf(i,j)

          call GMTR_TNvec( Tvec(:), Nvec(:),                            & ! [OUT]
                           wk(:,1,ij), wk(:,2,ij),                      & ! [IN]
                           GRD_grid_type, GMTR_polygon_type, GRD_rscale ) ! [IN]

          GMTR_a(ij,k0,l,AIJ,TNX) = Nvec(1)
          GMTR_a(ij,k0,l,AIJ,TNY) = Nvec(2)
          GMTR_a(ij,k0,l,AIJ,TNZ) = Nvec(3)
          GMTR_a(ij,k0,l,AIJ,TTX) = Tvec(1)
          GMTR_a(ij,k0,l,AIJ,TTY) = Tvec(2)
          GMTR_a(ij,k0,l,AIJ,TTZ) = Tvec(3)
       enddo
       enddo

       !--- AJ
       do j = ADM_gmin-1, ADM_gmax
       do i = ADM_gmin-1, ADM_gmax+1
          ij   = suf(i  ,j  )
          ijp1 = suf(i  ,j+1)

          do d = 1, ADM_nxyz
             wk(d,1,ij) = GRD_x(ij  ,k0,l,d)
             wk(d,2,ij) = GRD_x(ijp1,k0,l,d)
          enddo
       enddo
       enddo

       ! treat arc of unused triangle
       wk(:,1,suf(ADM_gmax+1,ADM_gmin-1)) = GRD_x(suf(ADM_gmax+1,ADM_gmin),k0,l,:)
       wk(:,2,suf(ADM_gmax+1,ADM_gmin-1)) = GRD_x(suf(ADM_gmax  ,ADM_gmin),k0,l,:)
       wk(:,1,suf(ADM_gmin-1,ADM_gmax  )) = GRD_x(suf(ADM_gmin-1,ADM_gmax),k0,l,:)
       wk(:,2,suf(ADM_gmin-1,ADM_gmax  )) = GRD_x(suf(ADM_gmin  ,ADM_gmax),k0,l,:)

       do j = ADM_gmin-1, ADM_gmax
       do i = ADM_gmin-1, ADM_gmax+1
          ij = suf(i,j)

          call GMTR_TNvec( Tvec(:), Nvec(:),                            & ! [OUT]
                           wk(:,1,ij), wk(:,2,ij),                      & ! [IN]
                           GRD_grid_type, GMTR_polygon_type, GRD_rscale ) ! [IN]

          GMTR_a(ij,k0,l,AJ,TNX) = Nvec(1)
          GMTR_a(ij,k0,l,AJ,TNY) = Nvec(2)
          GMTR_a(ij,k0,l,AJ,TNZ) = Nvec(3)
          GMTR_a(ij,k0,l,AJ,TTX) = Tvec(1)
          GMTR_a(ij,k0,l,AJ,TTY) = Tvec(2)
          GMTR_a(ij,k0,l,AJ,TTZ) = Tvec(3)
       enddo
       enddo

    enddo ! l loop

    !--- Hexagon
    do l = 1, ADM_lall

       !--- AI
       do j = ADM_gmin,   ADM_gmax
       do i = ADM_gmin-1, ADM_gmax
          ij   = suf(i  ,j  )
          ijm1 = suf(i  ,j-1)

          do d = 1, ADM_nxyz
             wk(d,1,ij) = GRD_xt(ij  ,k0,l,TI,d)
             wk(d,2,ij) = GRD_xt(ijm1,k0,l,TJ,d)
          enddo
       enddo
       enddo

       do j = ADM_gmin,   ADM_gmax
       do i = ADM_gmin-1, ADM_gmax
          ij = suf(i,j)

          call GMTR_TNvec( Tvec(:), Nvec(:),                            & ! [OUT]
                           wk(:,1,ij), wk(:,2,ij),                      & ! [IN]
                           GRD_grid_type, GMTR_polygon_type, GRD_rscale ) ! [IN]

          GMTR_a(ij,k0,l,AI,HNX) = Nvec(1)
          GMTR_a(ij,k0,l,AI,HNY) = Nvec(2)
          GMTR_a(ij,k0,l,AI,HNZ) = Nvec(3)
          GMTR_a(ij,k0,l,AI,HTX) = Tvec(1)
          GMTR_a(ij,k0,l,AI,HTY) = Tvec(2)
          GMTR_a(ij,k0,l,AI,HTZ) = Tvec(3)
       enddo
       enddo

       !--- AIJ
       do j = ADM_gmin-1, ADM_gmax
       do i = ADM_gmin-1, ADM_gmax
          ij   = suf(i  ,j  )

          do d = 1, ADM_nxyz
             wk(d,1,ij) = GRD_xt(ij  ,k0,l,TJ,d)
             wk(d,2,ij) = GRD_xt(ij  ,k0,l,TI,d)
          enddo
       enddo
       enddo

       ! treat arc of unused hexagon
       wk(:,1,suf(ADM_gmax  ,ADM_gmin-1)) = GRD_xt(suf(ADM_gmax  ,ADM_gmin-1),k0,l,TJ,:)
       wk(:,2,suf(ADM_gmax  ,ADM_gmin-1)) = GRD_xt(suf(ADM_gmax  ,ADM_gmin  ),k0,l,TI,:)
       wk(:,1,suf(ADM_gmin-1,ADM_gmax  )) = GRD_xt(suf(ADM_gmin  ,ADM_gmax  ),k0,l,TJ,:)
       wk(:,2,suf(ADM_gmin-1,ADM_gmax  )) = GRD_xt(suf(ADM_gmin-1,ADM_gmax  ),k0,l,TI,:)

       do j = ADM_gmin-1, ADM_gmax
       do i = ADM_gmin-1, ADM_gmax
          ij = suf(i,j)

          call GMTR_TNvec( Tvec(:), Nvec(:),                            & ! [OUT]
                           wk(:,1,ij), wk(:,2,ij),                      & ! [IN]
                           GRD_grid_type, GMTR_polygon_type, GRD_rscale ) ! [IN]

          GMTR_a(ij,k0,l,AIJ,HNX) = Nvec(1)
          GMTR_a(ij,k0,l,AIJ,HNY) = Nvec(2)
          GMTR_a(ij,k0,l,AIJ,HNZ) = Nvec(3)
          GMTR_a(ij,k0,l,AIJ,HTX) = Tvec(1)
          GMTR_a(ij,k0,l,AIJ,HTY) = Tvec(2)
          GMTR_a(ij,k0,l,AIJ,HTZ) = Tvec(3)
       enddo
       enddo

       !--- AJ
       do j = ADM_gmin-1, ADM_gmax
       do i = ADM_gmin,   ADM_gmax
          ij   = suf(i  ,j  )
          im1j = suf(i-1,j  )

          do d = 1, ADM_nxyz
             wk(d,1,ij) = GRD_xt(im1j,k0,l,TI,d)
             wk(d,2,ij) = GRD_xt(ij  ,k0,l,TJ,d)
          enddo
       enddo
       enddo

       if ( ADM_have_sgp(l) ) then ! pentagon
          wk(:,1,suf(ADM_gmin  ,ADM_gmin-1)) = GRD_xt(suf(ADM_gmin  ,ADM_gmin  ),k0,l,TI,:)
          wk(:,2,suf(ADM_gmin  ,ADM_gmin-1)) = GRD_xt(suf(ADM_gmin  ,ADM_gmin-1),k0,l,TJ,:)
       endif

       do j = ADM_gmin-1, ADM_gmax
       do i = ADM_gmin,   ADM_gmax
          ij = suf(i,j)

          call GMTR_TNvec( Tvec(:), Nvec(:),                            & ! [OUT]
                           wk(:,1,ij), wk(:,2,ij),                      & ! [IN]
                           GRD_grid_type, GMTR_polygon_type, GRD_rscale ) ! [IN]

          GMTR_a(ij,k0,l,AJ,HNX) = Nvec(1)
          GMTR_a(ij,k0,l,AJ,HNY) = Nvec(2)
          GMTR_a(ij,k0,l,AJ,HNZ) = Nvec(3)
          GMTR_a(ij,k0,l,AJ,HTX) = Tvec(1)
          GMTR_a(ij,k0,l,AJ,HTY) = Tvec(2)
          GMTR_a(ij,k0,l,AJ,HTZ) = Tvec(3)
       enddo
       enddo

    enddo ! l loop

    if ( ADM_have_pl ) then
       n = ADM_gslf_pl

       do l = 1, ADM_lall_pl

          !--- Triangle (arc 1)
          do v = ADM_gmin_pl, ADM_gmax_pl
             ij = v

             do d = 1, ADM_nxyz
                wk_pl(d,1) = GRD_x_pl(n ,k0,l,d)
                wk_pl(d,2) = GRD_x_pl(ij,k0,l,d)
             enddo

             call GMTR_TNvec( Tvec(:), Nvec(:),                            & ! [OUT]
                              wk_pl(:,1), wk_pl(:,2),                      & ! [IN]
                              GRD_grid_type, GMTR_polygon_type, GRD_rscale ) ! [IN]

             GMTR_a_pl(ij,k0,l,TNX) = Nvec(1)
             GMTR_a_pl(ij,k0,l,TNY) = Nvec(2)
             GMTR_a_pl(ij,k0,l,TNZ) = Nvec(3)
             GMTR_a_pl(ij,k0,l,TTX) = Tvec(1)
             GMTR_a_pl(ij,k0,l,TTY) = Tvec(2)
             GMTR_a_pl(ij,k0,l,TTZ) = Tvec(3)
          enddo

          !--- Triangle (arc 2)
          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijp1 = v+1
             if ( ijp1 == ADM_gmax_pl+1 ) ijp1 = ADM_gmin_pl

             do d = 1, ADM_nxyz
                wk_pl(d,1) = GRD_x_pl(ij  ,k0,l,d)
                wk_pl(d,2) = GRD_x_pl(ijp1,k0,l,d)
             enddo

             call GMTR_TNvec( Tvec(:), Nvec(:),                            & ! [OUT]
                              wk_pl(:,1), wk_pl(:,2),                      & ! [IN]
                              GRD_grid_type, GMTR_polygon_type, GRD_rscale ) ! [IN]

             GMTR_a_pl(ij,k0,l,TN2X) = Nvec(1)
             GMTR_a_pl(ij,k0,l,TN2Y) = Nvec(2)
             GMTR_a_pl(ij,k0,l,TN2Z) = Nvec(3)
             GMTR_a_pl(ij,k0,l,TT2X) = Tvec(1)
             GMTR_a_pl(ij,k0,l,TT2Y) = Tvec(2)
             GMTR_a_pl(ij,k0,l,TT2Z) = Tvec(3)
          enddo

          !--- Hexagon
          do v = ADM_gmin_pl, ADM_gmax_pl
             ij   = v
             ijm1 = v-1
             if ( ijm1 == ADM_gmin_pl-1 ) ijp1 = ADM_gmax_pl

             do d = 1, ADM_nxyz
                wk_pl(d,1) = GRD_xt_pl(ijm1,k0,l,d)
                wk_pl(d,2) = GRD_xt_pl(ij  ,k0,l,d)
             enddo

             call GMTR_TNvec( Tvec(:), Nvec(:),                            & ! [OUT]
                              wk_pl(:,1), wk_pl(:,2),                      & ! [IN]
                              GRD_grid_type, GMTR_polygon_type, GRD_rscale ) ! [IN]

             GMTR_a_pl(ij,k0,l,HNX) = Nvec(1)
             GMTR_a_pl(ij,k0,l,HNY) = Nvec(2)
             GMTR_a_pl(ij,k0,l,HNZ) = Nvec(3)
             GMTR_a_pl(ij,k0,l,HTX) = Tvec(1)
             GMTR_a_pl(ij,k0,l,HTY) = Tvec(2)
             GMTR_a_pl(ij,k0,l,HTZ) = Tvec(3)
          enddo

       enddo
    endif

    return
  end subroutine GMTR_a_setup

  !-----------------------------------------------------------------------------
  subroutine GMTR_TNvec( &
       vT,           &
       vN,           &
       vFrom,        &
       vTo,          &
       grid_type,    &
       polygon_type, &
       radius        )
    use mod_vector, only: &
       VECTR_dot,   &
       VECTR_cross, &
       VECTR_abs,   &
       VECTR_angle
    implicit none

    real(RP),         intent(out) :: vT   (3)     ! tangential vector
    real(RP),         intent(out) :: vN   (3)     ! normal     vector
    real(RP),         intent(in)  :: vFrom(3)
    real(RP),         intent(in)  :: vTo  (3)
    character(len=*), intent(in)  :: grid_type    ! ON_SPHERE or ON_PLANE
    character(len=*), intent(in)  :: polygon_type ! ON_SPHERE or ON_PLANE
    real(RP),         intent(in)  :: radius

    real(RP), parameter :: o(3) = 0.0_RP

    real(RP) :: angle, length
    real(RP) :: distance
    !---------------------------------------------------------------------------

    ! calculate tangential vector
    vT(:) = vTo(:) - vFrom(:)

    if ( grid_type == 'ON_PLANE' ) then ! treat as point on the plane

       ! calculate normal vector
       vN(1) = -vT(2)
       vN(2) =  vT(1)
       vN(3) = 0.0_RP

    elseif( grid_type == 'ON_SPHERE' ) then ! treat as point on the sphere

       if ( polygon_type == 'ON_PLANE' ) then ! length of a line

          call VECTR_dot( distance, vFrom(:), vTo(:), vFrom(:), vTo(:) )
          distance = sqrt( distance )

       elseif( polygon_type == 'ON_SPHERE' ) then ! length of a geodesic line ( angle * radius )

          call VECTR_angle( angle, vFrom(:), o(:), vTo(:) )
          distance = angle * radius

       endif

       call VECTR_abs( length, vT(:) )
       vT(:) = vT(:) * distance / length

       ! calculate normal vector
       call VECTR_cross( vN(:), o(:), vFrom(:), o(:), vTo(:) )

       call VECTR_abs( length, vN(:) )
       vN(:) = vN(:) * distance / length

    endif

    return
  end subroutine GMTR_TNvec

  !-----------------------------------------------------------------------------
  integer function suf(i,j)
!ESC!    use mod_adm, only: &
!ESC!       ADM_gall_1d
    implicit none

    integer :: i, j
    !---------------------------------------------------------------------------

    suf = ADM_gall_1d * (j-1) + i

  end function suf

end module mod_gmtr
!-------------------------------------------------------------------------------
