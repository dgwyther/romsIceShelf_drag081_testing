      MODULE analytical_mod
!
!svn $Id$
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!   PACKAGE:                                                 !
!                                                                      !
!  This package is used to provide various analytical fields to the    !
!  model when appropriate.                                             !
!                                                                      !
!=======================================================================
!
      implicit none
!
      CONTAINS
!
      SUBROUTINE ana_btflux (ng, tile, model, itrc)
!
!=======================================================================
!                                                                      !
!  This routine sets kinematic bottom flux of tracer type variables    !
!  (tracer units m/s).                                                 !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, itrc
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
!
!  Set horizontal starting and ending indices for automatic private
!  storage arrays.
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
!  Determine array lower and upper bounds in the I- and J-directions.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
!  Set array lower and upper bounds for MIN(I,J) directions and
!  MAX(I,J) directions.
!
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
!
      CALL ana_btflux_tile (ng, tile, model, itrc,                      &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      FORCES(ng) % btflx)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        ANANAME( 3)="ROMS/Functionals/ana_btflux.h"
      END IF
      RETURN
      END SUBROUTINE ana_btflux
!
!***********************************************************************
      SUBROUTINE ana_btflux_tile (ng, tile, model, itrc,                &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            btflx)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, itrc
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
      real(r8), intent(inout) :: btflx(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: i, j
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrB, IstrP, IstrR, IstrT, IstrM, IstrU
      integer :: Iend, IendB, IendP, IendR, IendT
      integer :: Jstr, JstrB, JstrP, JstrR, JstrT, JstrM, JstrV
      integer :: Jend, JendB, JendP, JendR, JendT
      integer :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer :: Jendp1, Jendp2, Jendp2i, Jendp3
!
      Istr   =BOUNDS(ng) % Istr   (tile)
      IstrB  =BOUNDS(ng) % IstrB  (tile)
      IstrM  =BOUNDS(ng) % IstrM  (tile)
      IstrP  =BOUNDS(ng) % IstrP  (tile)
      IstrR  =BOUNDS(ng) % IstrR  (tile)
      IstrT  =BOUNDS(ng) % IstrT  (tile)
      IstrU  =BOUNDS(ng) % IstrU  (tile)
      Iend   =BOUNDS(ng) % Iend   (tile)
      IendB  =BOUNDS(ng) % IendB  (tile)
      IendP  =BOUNDS(ng) % IendP  (tile)
      IendR  =BOUNDS(ng) % IendR  (tile)
      IendT  =BOUNDS(ng) % IendT  (tile)
      Jstr   =BOUNDS(ng) % Jstr   (tile)
      JstrB  =BOUNDS(ng) % JstrB  (tile)
      JstrM  =BOUNDS(ng) % JstrM  (tile)
      JstrP  =BOUNDS(ng) % JstrP  (tile)
      JstrR  =BOUNDS(ng) % JstrR  (tile)
      JstrT  =BOUNDS(ng) % JstrT  (tile)
      JstrV  =BOUNDS(ng) % JstrV  (tile)
      Jend   =BOUNDS(ng) % Jend   (tile)
      JendB  =BOUNDS(ng) % JendB  (tile)
      JendP  =BOUNDS(ng) % JendP  (tile)
      JendR  =BOUNDS(ng) % JendR  (tile)
      JendT  =BOUNDS(ng) % JendT  (tile)
!
      Istrm3 =BOUNDS(ng) % Istrm3 (tile)            ! Istr-3
      Istrm2 =BOUNDS(ng) % Istrm2 (tile)            ! Istr-2
      Istrm1 =BOUNDS(ng) % Istrm1 (tile)            ! Istr-1
      IstrUm2=BOUNDS(ng) % IstrUm2(tile)            ! IstrU-2
      IstrUm1=BOUNDS(ng) % IstrUm1(tile)            ! IstrU-1
      Iendp1 =BOUNDS(ng) % Iendp1 (tile)            ! Iend+1
      Iendp2 =BOUNDS(ng) % Iendp2 (tile)            ! Iend+2
      Iendp2i=BOUNDS(ng) % Iendp2i(tile)            ! Iend+2 interior
      Iendp3 =BOUNDS(ng) % Iendp3 (tile)            ! Iend+3
      Jstrm3 =BOUNDS(ng) % Jstrm3 (tile)            ! Jstr-3
      Jstrm2 =BOUNDS(ng) % Jstrm2 (tile)            ! Jstr-2
      Jstrm1 =BOUNDS(ng) % Jstrm1 (tile)            ! Jstr-1
      JstrVm2=BOUNDS(ng) % JstrVm2(tile)            ! JstrV-2
      JstrVm1=BOUNDS(ng) % JstrVm1(tile)            ! JstrV-1
      Jendp1 =BOUNDS(ng) % Jendp1 (tile)            ! Jend+1
      Jendp2 =BOUNDS(ng) % Jendp2 (tile)            ! Jend+2
      Jendp2i=BOUNDS(ng) % Jendp2i(tile)            ! Jend+2 interior
      Jendp3 =BOUNDS(ng) % Jendp3 (tile)            ! Jend+3
!
!-----------------------------------------------------------------------
!  Set kinematic bottom heat flux (degC m/s) at horizontal RHO-points.
!-----------------------------------------------------------------------
!
      IF (itrc.eq.itemp) THEN
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            btflx(i,j,itrc)=0.0_r8
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Set kinematic bottom salt flux (m/s) at horizontal RHO-points,
!  scaling by bottom salinity is done elsewhere.
!-----------------------------------------------------------------------
!
      ELSE IF (itrc.eq.isalt) THEN
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            btflx(i,j,itrc)=0.0_r8
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Set kinematic bottom flux (T m/s) of passive tracers, if any.
!-----------------------------------------------------------------------
!
      ELSE
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            btflx(i,j,itrc)=0.0_r8
          END DO
        END DO
      END IF
      RETURN
      END SUBROUTINE ana_btflux_tile
      SUBROUTINE ana_icedrag (ng, tile, model)
!
!================================================== David E. Gwyther ===
!                                                                      !
!  This routine sets analytical, spatially varying ice bottom roughness!
!  length (m), or linear drag coefficients (m/s), or quadratic drag    !
!  coefficients (nondimensional) at RHO-points.                        !
!                                                                      !
! Currently, only  is implemented.                             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ncparam
      USE mod_iceshelfvar
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
!
!  Set horizontal starting and ending indices for automatic private
!  storage arrays.
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
!  Determine array lower and upper bounds in the I- and J-directions.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
!  Set array lower and upper bounds for MIN(I,J) directions and
!  MAX(I,J) directions.
!
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
!
      CALL ana_icedrag_tile (ng, tile, model,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    ICESHELFVAR(ng) % idrag2)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        ANANAME(38)="ROMS/Functionals/ana_icedrag.h"
      END IF
      RETURN
      END SUBROUTINE ana_icedrag
!
!***********************************************************************
      SUBROUTINE ana_icedrag_tile (ng, tile, model,                     &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          idrag2)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
      USE mp_exchange_mod, ONLY : mp_exchange2d
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
      real(r8), intent(out) :: idrag2(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: i, j
      real(r8) :: cff
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrB, IstrP, IstrR, IstrT, IstrM, IstrU
      integer :: Iend, IendB, IendP, IendR, IendT
      integer :: Jstr, JstrB, JstrP, JstrR, JstrT, JstrM, JstrV
      integer :: Jend, JendB, JendP, JendR, JendT
      integer :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer :: Jendp1, Jendp2, Jendp2i, Jendp3
!
      Istr   =BOUNDS(ng) % Istr   (tile)
      IstrB  =BOUNDS(ng) % IstrB  (tile)
      IstrM  =BOUNDS(ng) % IstrM  (tile)
      IstrP  =BOUNDS(ng) % IstrP  (tile)
      IstrR  =BOUNDS(ng) % IstrR  (tile)
      IstrT  =BOUNDS(ng) % IstrT  (tile)
      IstrU  =BOUNDS(ng) % IstrU  (tile)
      Iend   =BOUNDS(ng) % Iend   (tile)
      IendB  =BOUNDS(ng) % IendB  (tile)
      IendP  =BOUNDS(ng) % IendP  (tile)
      IendR  =BOUNDS(ng) % IendR  (tile)
      IendT  =BOUNDS(ng) % IendT  (tile)
      Jstr   =BOUNDS(ng) % Jstr   (tile)
      JstrB  =BOUNDS(ng) % JstrB  (tile)
      JstrM  =BOUNDS(ng) % JstrM  (tile)
      JstrP  =BOUNDS(ng) % JstrP  (tile)
      JstrR  =BOUNDS(ng) % JstrR  (tile)
      JstrT  =BOUNDS(ng) % JstrT  (tile)
      JstrV  =BOUNDS(ng) % JstrV  (tile)
      Jend   =BOUNDS(ng) % Jend   (tile)
      JendB  =BOUNDS(ng) % JendB  (tile)
      JendP  =BOUNDS(ng) % JendP  (tile)
      JendR  =BOUNDS(ng) % JendR  (tile)
      JendT  =BOUNDS(ng) % JendT  (tile)
!
      Istrm3 =BOUNDS(ng) % Istrm3 (tile)            ! Istr-3
      Istrm2 =BOUNDS(ng) % Istrm2 (tile)            ! Istr-2
      Istrm1 =BOUNDS(ng) % Istrm1 (tile)            ! Istr-1
      IstrUm2=BOUNDS(ng) % IstrUm2(tile)            ! IstrU-2
      IstrUm1=BOUNDS(ng) % IstrUm1(tile)            ! IstrU-1
      Iendp1 =BOUNDS(ng) % Iendp1 (tile)            ! Iend+1
      Iendp2 =BOUNDS(ng) % Iendp2 (tile)            ! Iend+2
      Iendp2i=BOUNDS(ng) % Iendp2i(tile)            ! Iend+2 interior
      Iendp3 =BOUNDS(ng) % Iendp3 (tile)            ! Iend+3
      Jstrm3 =BOUNDS(ng) % Jstrm3 (tile)            ! Jstr-3
      Jstrm2 =BOUNDS(ng) % Jstrm2 (tile)            ! Jstr-2
      Jstrm1 =BOUNDS(ng) % Jstrm1 (tile)            ! Jstr-1
      JstrVm2=BOUNDS(ng) % JstrVm2(tile)            ! JstrV-2
      JstrVm1=BOUNDS(ng) % JstrVm1(tile)            ! JstrV-1
      Jendp1 =BOUNDS(ng) % Jendp1 (tile)            ! Jend+1
      Jendp2 =BOUNDS(ng) % Jendp2 (tile)            ! Jend+2
      Jendp2i=BOUNDS(ng) % Jendp2i(tile)            ! Jend+2 interior
      Jendp3 =BOUNDS(ng) % Jendp3 (tile)            ! Jend+3
!
!-----------------------------------------------------------------------
!  Set spatially varying quadratic ice-ocean drag coefficient (nondimensional)
!-----------------------------------------------------------------------
!
      DO j=JstrT,JendT          ! based on Chezy coefficient (g/c^2)
        DO i=IstrT,IendT
          idrag2(i,j)=0.1_r8
        END DO
      END DO
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          idrag2)
      END IF
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    idrag2)
      RETURN
      END SUBROUTINE ana_icedrag_tile
      SUBROUTINE ana_grid (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine sets model grid using an analytical expressions.       !
!                                                                      !
!  On Output:  stored in common blocks:                                !
!                                                                      !
!                           "grid"    (file grid.h)                    !
!                           "scalars" (file scalar.h)                  !
!                                                                      !
!     el       Length (m) of domain box in the ETA-direction.          !
!     f        Coriolis parameter (1/seconds) at RHO-points.           !
!     h        Bathymetry (meters; positive) at RHO-points.            !
!     hmin     Minimum depth of bathymetry (m).                        !
!     hmax     Maximum depth of bathymetry (m).                        !
!     pm       Coordinate transformation metric "m" (1/meters)         !
!              associated with the differential distances in XI        !
!              at RHO-points.                                          !
!     pn       Coordinate transformation metric "n" (1/meters)         !
!              associated with the differential distances in ETA.      !
!              at RHO-points.                                          !
!     xl       Length (m) of domain box in the XI-direction.           !
!     xp       XI-coordinates (m) at PSI-points.                       !
!     xr       XI-coordinates (m) at RHO-points.                       !
!     yp       ETA-coordinates (m) at PSI-points.                      !
!     yr       ETA-coordinates (m) at RHO-points.                      !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
!
!  Set horizontal starting and ending indices for automatic private
!  storage arrays.
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
!  Determine array lower and upper bounds in the I- and J-directions.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
!  Set array lower and upper bounds for MIN(I,J) directions and
!  MAX(I,J) directions.
!
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
!
      CALL ana_grid_tile (ng, tile, model,                              &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    GRID(ng) % angler,                            &
     &                    GRID(ng) % dmde,                              &
     &                    GRID(ng) % dndx,                              &
     &                    GRID(ng) % zice,                              &
     &                    GRID(ng) % lonp,                              &
     &                    GRID(ng) % lonr,                              &
     &                    GRID(ng) % lonu,                              &
     &                    GRID(ng) % lonv,                              &
     &                    GRID(ng) % latp,                              &
     &                    GRID(ng) % latr,                              &
     &                    GRID(ng) % latu,                              &
     &                    GRID(ng) % latv,                              &
     &                    GRID(ng) % pn,                                &
     &                    GRID(ng) % pm,                                &
     &                    GRID(ng) % f,                                 &
     &                    GRID(ng) % h)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        ANANAME( 7)="ROMS/Functionals/ana_grid.h"
      END IF
      RETURN
      END SUBROUTINE ana_grid
!
!***********************************************************************
      SUBROUTINE ana_grid_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          angler,                                 &
     &                          dmde, dndx,                             &
     &                          zice,                                   &
     &                          lonp, lonr, lonu, lonv,                 &
     &                          latp, latr, latu, latv,                 &
     &                          pn, pm, f, h)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_scalars
!
      USE distribute_mod, ONLY : mp_reduce
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
      USE mp_exchange_mod, ONLY : mp_exchange2d
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
      real(r8), intent(out) :: angler(LBi:,LBj:)
      real(r8), intent(out) :: dmde(LBi:,LBj:)
      real(r8), intent(out) :: dndx(LBi:,LBj:)
      real(r8), intent(out) :: zice(LBi:,LBj:)
      real(r8), intent(out) :: lonp(LBi:,LBj:)
      real(r8), intent(out) :: lonr(LBi:,LBj:)
      real(r8), intent(out) :: lonu(LBi:,LBj:)
      real(r8), intent(out) :: lonv(LBi:,LBj:)
      real(r8), intent(out) :: latp(LBi:,LBj:)
      real(r8), intent(out) :: latr(LBi:,LBj:)
      real(r8), intent(out) :: latu(LBi:,LBj:)
      real(r8), intent(out) :: latv(LBi:,LBj:)
      real(r8), intent(out) :: pn(LBi:,LBj:)
      real(r8), intent(out) :: pm(LBi:,LBj:)
      real(r8), intent(out) :: f(LBi:,LBj:)
      real(r8), intent(out) :: h(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: Imin, Imax, Jmin, Jmax
      integer :: NSUB, i, ival, j, k
      real(r8), parameter :: twopi = 2.0_r8*pi
      real(r8) :: Esize, Xsize, beta, cff, depth, dth
      real(r8) :: dx, dy, f0, my_min, my_max, r, theta, val1, val2
      real(r8), dimension(2) :: buffer
      character (len=3), dimension(2) :: op_handle
      real(r8) :: wrkX(IminS:ImaxS,JminS:JmaxS)
      real(r8) :: wrkY(IminS:ImaxS,JminS:JmaxS)
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrB, IstrP, IstrR, IstrT, IstrM, IstrU
      integer :: Iend, IendB, IendP, IendR, IendT
      integer :: Jstr, JstrB, JstrP, JstrR, JstrT, JstrM, JstrV
      integer :: Jend, JendB, JendP, JendR, JendT
      integer :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer :: Jendp1, Jendp2, Jendp2i, Jendp3
!
      Istr   =BOUNDS(ng) % Istr   (tile)
      IstrB  =BOUNDS(ng) % IstrB  (tile)
      IstrM  =BOUNDS(ng) % IstrM  (tile)
      IstrP  =BOUNDS(ng) % IstrP  (tile)
      IstrR  =BOUNDS(ng) % IstrR  (tile)
      IstrT  =BOUNDS(ng) % IstrT  (tile)
      IstrU  =BOUNDS(ng) % IstrU  (tile)
      Iend   =BOUNDS(ng) % Iend   (tile)
      IendB  =BOUNDS(ng) % IendB  (tile)
      IendP  =BOUNDS(ng) % IendP  (tile)
      IendR  =BOUNDS(ng) % IendR  (tile)
      IendT  =BOUNDS(ng) % IendT  (tile)
      Jstr   =BOUNDS(ng) % Jstr   (tile)
      JstrB  =BOUNDS(ng) % JstrB  (tile)
      JstrM  =BOUNDS(ng) % JstrM  (tile)
      JstrP  =BOUNDS(ng) % JstrP  (tile)
      JstrR  =BOUNDS(ng) % JstrR  (tile)
      JstrT  =BOUNDS(ng) % JstrT  (tile)
      JstrV  =BOUNDS(ng) % JstrV  (tile)
      Jend   =BOUNDS(ng) % Jend   (tile)
      JendB  =BOUNDS(ng) % JendB  (tile)
      JendP  =BOUNDS(ng) % JendP  (tile)
      JendR  =BOUNDS(ng) % JendR  (tile)
      JendT  =BOUNDS(ng) % JendT  (tile)
!
      Istrm3 =BOUNDS(ng) % Istrm3 (tile)            ! Istr-3
      Istrm2 =BOUNDS(ng) % Istrm2 (tile)            ! Istr-2
      Istrm1 =BOUNDS(ng) % Istrm1 (tile)            ! Istr-1
      IstrUm2=BOUNDS(ng) % IstrUm2(tile)            ! IstrU-2
      IstrUm1=BOUNDS(ng) % IstrUm1(tile)            ! IstrU-1
      Iendp1 =BOUNDS(ng) % Iendp1 (tile)            ! Iend+1
      Iendp2 =BOUNDS(ng) % Iendp2 (tile)            ! Iend+2
      Iendp2i=BOUNDS(ng) % Iendp2i(tile)            ! Iend+2 interior
      Iendp3 =BOUNDS(ng) % Iendp3 (tile)            ! Iend+3
      Jstrm3 =BOUNDS(ng) % Jstrm3 (tile)            ! Jstr-3
      Jstrm2 =BOUNDS(ng) % Jstrm2 (tile)            ! Jstr-2
      Jstrm1 =BOUNDS(ng) % Jstrm1 (tile)            ! Jstr-1
      JstrVm2=BOUNDS(ng) % JstrVm2(tile)            ! JstrV-2
      JstrVm1=BOUNDS(ng) % JstrVm1(tile)            ! JstrV-1
      Jendp1 =BOUNDS(ng) % Jendp1 (tile)            ! Jend+1
      Jendp2 =BOUNDS(ng) % Jendp2 (tile)            ! Jend+2
      Jendp2i=BOUNDS(ng) % Jendp2i(tile)            ! Jend+2 interior
      Jendp3 =BOUNDS(ng) % Jendp3 (tile)            ! Jend+3
!
!-----------------------------------------------------------------------
!  Set grid parameters:
!
!     Xsize    Length (m) of domain box in the XI-direction.
!     Esize    Length (m) of domain box in the ETA-direction.
!     depth    Maximum depth of bathymetry (m).
!     f0       Coriolis parameter, f-plane constant (1/s).
!     beta     Coriolis parameter, beta-plane constant (1/s/m).
!-----------------------------------------------------------------------
!
      Xsize=570.466E+03_r8
      Esize=1120.0E+03_r8
      depth=900.0_r8
      f0=0.0_r8
      beta=0.0_r8
!
!  Load grid parameters to global storage.
!
      IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
        xl(ng)=Xsize
        el(ng)=Esize
      END IF
!
!-----------------------------------------------------------------------
!  Compute the (XI,ETA) coordinates at PSI- and RHO-points.
!  Set grid spacing (m).
!-----------------------------------------------------------------------
!
!  Determine I- and J-ranges for computing grid data.  These ranges
!  are special in periodic boundary conditons since periodicity cannot
!  be imposed in the grid coordinates.
!
      IF (DOMAIN(ng)%Western_Edge(tile)) THEN
        Imin=Istr-1
      ELSE
        Imin=Istr
      END IF
      IF (DOMAIN(ng)%Eastern_Edge(tile)) THEN
        Imax=Iend+1
      ELSE
        Imax=Iend
      END IF
      IF (DOMAIN(ng)%Southern_Edge(tile)) THEN
        Jmin=Jstr-1
      ELSE
        Jmin=Jstr
      END IF
      IF (DOMAIN(ng)%Northern_Edge(tile)) THEN
        Jmax=Jend+1
      ELSE
        Jmax=Jend
      END IF
!
!  Spherical coordinates set-up.
!
      spherical=.TRUE.
      DO j=JstrR,JendR
        cff=-80.0_r8+0.1_r8*REAL(j-1,r8)
        DO i=IstrR,IendR
          lonr(i,j)=0.3_r8*REAL(i-1,r8)
          latr(i,j)=cff
          lonu(i,j)=0.3_r8*REAL(i-1,r8)+0.15_r8
          lonp(i,j)=lonu(i,j)
          latu(i,j)=latr(i,j)
          lonv(i,j)=lonr(i,j)
          latv(i,j)=cff+0.05_r8
          latp(i,j)=latv(i,j)
        END DO
      END DO
!
!  Exchange boundary data.
!
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, .FALSE., .FALSE.,               &
     &                    lonp, lonr, lonu, lonv)
      CALL mp_exchange2d (ng, tile, model, 4,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints, .FALSE., .FALSE.,               &
     &                    latp, latr, latu, latv)
!
!-----------------------------------------------------------------------
! Compute coordinate transformation metrics at RHO-points "pm" and
! "pn"  (1/m) associated with the differential distances in XI and
! ETA, respectively.
!-----------------------------------------------------------------------
!
!
!  Spherical coordinates set-up.
!
      Eradius=6371020.0_r8
      val1=360.0_r8/(0.3_r8*(2.0_r8*pi*Eradius))
      val2=360.0_r8/(0.1_r8*(2.0_r8*pi*Eradius))
      DO j=MIN(JstrT,Jstr-1),MAX(Jend+1,JendT)
        DO i=MIN(IstrT,Istr-1),MAX(Iend+1,IendT)
          wrkX(i,j)=val1/COS(latr(i,j)*deg2rad)
          wrkY(i,j)=val2
        END DO
      END DO
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          pm(i,j)=wrkX(i,j)
          pn(i,j)=wrkY(i,j)
        END DO
      END DO
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pm)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          pn)
      END IF
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    pm, pn)
!
!-----------------------------------------------------------------------
!  Compute d(1/n)/d(xi) and d(1/m)/d(eta) at RHO-points.
!-----------------------------------------------------------------------
!
      DO j=Jstr,Jend
        DO i=Istr,Iend
          dndx(i,j)=0.5_r8*((1.0_r8/wrkY(i+1,j  ))-                     &
     &                      (1.0_r8/wrkY(i-1,j  )))
          dmde(i,j)=0.5_r8*((1.0_r8/wrkX(i  ,j+1))-                     &
     &                      (1.0_r8/wrkX(i  ,j-1)))
        END DO
      END DO
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          dndx)
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          dmde)
      END IF
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    dndx, dmde)
!
!-----------------------------------------------------------------------
! Angle (radians) between XI-axis and true EAST at RHO-points.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          f(i,j)=(4.0_r8*pi/86164.1_r8)*                                &
     &           SIN(latr(i,j)*deg2rad)
        END DO
      END DO
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          angler)
      END IF
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    angler)
!
!-----------------------------------------------------------------------
!  Compute Coriolis parameter (1/s) at RHO-points.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          f(i,j)=(4.0_r8*pi/86164.1_r8)*                                &
     &           SIN(latr(i,j)*deg2rad)
        END DO
      END DO
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          f)
      END IF
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    f)
!
!-----------------------------------------------------------------------
!  Set bathymetry (meters; positive) at RHO-points.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          h(i,j)=900.0_r8
        END DO
      END DO
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          h)
      END IF
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    h)
!
! Determine minimum depth: first, determine minimum values of depth
! within each subdomain, then determine global minimum by comparing
! these subdomain minima.
!
      my_min=h(IstrT,JstrT)
      my_max=h(IstrT,JstrT)
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          my_min=MIN(my_min,h(i,j))
          my_max=MAX(my_max,h(i,j))
        END DO
      END DO
      NSUB=1                             ! distributed-memory
      IF (tile_count.eq.0) THEN
        hmin(ng)=my_min
        hmax(ng)=my_max
      ELSE
        hmin(ng)=MIN(hmin(ng),my_min)
        hmax(ng)=MAX(hmax(ng),my_max)
      END IF
      tile_count=tile_count+1
      IF (tile_count.eq.NSUB) THEN
        tile_count=0
        buffer(1)=hmin(ng)
        buffer(2)=hmax(ng)
        op_handle(1)='MIN'
        op_handle(2)='MAX'
        CALL mp_reduce (ng, model, 2, buffer, op_handle)
        hmin(ng)=buffer(1)
        hmax(ng)=buffer(2)
      END IF
!
!-----------------------------------------------------------------------
!  Set depth of ice shelf (meters; negative) at RHO-points.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=IstrR,IendR
          IF (j.eq.0) THEN
            zice(i,j)=-700.0_r8
          ELSE IF (j.le.41) THEN
            zice(i,j)=-700.0_r8+(500.0_r8/40.0_r8)*REAL(j-1,r8)
          ELSE
            zice(i,j)=0.0_r8
          END IF
!          IF (j.le.10) THEN
!            zice(i,j)=-700.0_r8
!          ELSE IF (j.lt.41) THEN
!            zice(i,j)=-700.0_r8+(500.0_r8/30.0_r8)*REAL(j-11,r8)
!          ELSE
!            zice(i,j)=-200.0_r8
!          END IF
        END DO
      END DO
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          zice)
      END IF
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    zice)
      RETURN
      END SUBROUTINE ana_grid_tile
      SUBROUTINE ana_initial (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This subroutine sets initial conditions for momentum and tracer     !
!  type variables using analytical expressions.                        !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
!
!  Set horizontal starting and ending indices for automatic private
!  storage arrays.
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
!  Determine array lower and upper bounds in the I- and J-directions.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
!  Set array lower and upper bounds for MIN(I,J) directions and
!  MAX(I,J) directions.
!
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
!
      IF (model.eq.iNLM) THEN
        CALL ana_NLMinitial_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            GRID(ng) % h,                         &
     &                            GRID(ng) % lonr,                      &
     &                            GRID(ng) % latr,                      &
     &                            GRID(ng) % z_r,                       &
     &                            OCEAN(ng) % u,                        &
     &                            OCEAN(ng) % v,                        &
     &                            OCEAN(ng) % t,                        &
     &                            OCEAN(ng) % ubar,                     &
     &                            OCEAN(ng) % vbar,                     &
     &                            OCEAN(ng) % zeta)
      END IF
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        ANANAME(10)="ROMS/Functionals/ana_initial.h"
      END IF
      RETURN
      END SUBROUTINE ana_initial
!
!***********************************************************************
      SUBROUTINE ana_NLMinitial_tile (ng, tile, model,                  &
     &                                LBi, UBi, LBj, UBj,               &
     &                                IminS, ImaxS, JminS, JmaxS,       &
     &                                h,                                &
     &                                lonr, latr,                       &
     &                                z_r,                              &
     &                                u, v, t,                          &
     &                                ubar, vbar, zeta)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_scalars
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
      real(r8), intent(in) :: h(LBi:,LBj:)
      real(r8), intent(in) :: lonr(LBi:,LBj:)
      real(r8), intent(in) :: latr(LBi:,LBj:)
      real(r8), intent(in) :: z_r(LBi:,LBj:,:)
      real(r8), intent(out) :: u(LBi:,LBj:,:,:)
      real(r8), intent(out) :: v(LBi:,LBj:,:,:)
      real(r8), intent(out) :: t(LBi:,LBj:,:,:,:)
      real(r8), intent(out) :: ubar(LBi:,LBj:,:)
      real(r8), intent(out) :: vbar(LBi:,LBj:,:)
      real(r8), intent(out) :: zeta(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: Iless, Iplus, i, itrc, j, k
      real(r8) :: depth, dx, val1, val2, val3, val4, x, x0, y, y0
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrB, IstrP, IstrR, IstrT, IstrM, IstrU
      integer :: Iend, IendB, IendP, IendR, IendT
      integer :: Jstr, JstrB, JstrP, JstrR, JstrT, JstrM, JstrV
      integer :: Jend, JendB, JendP, JendR, JendT
      integer :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer :: Jendp1, Jendp2, Jendp2i, Jendp3
!
      Istr   =BOUNDS(ng) % Istr   (tile)
      IstrB  =BOUNDS(ng) % IstrB  (tile)
      IstrM  =BOUNDS(ng) % IstrM  (tile)
      IstrP  =BOUNDS(ng) % IstrP  (tile)
      IstrR  =BOUNDS(ng) % IstrR  (tile)
      IstrT  =BOUNDS(ng) % IstrT  (tile)
      IstrU  =BOUNDS(ng) % IstrU  (tile)
      Iend   =BOUNDS(ng) % Iend   (tile)
      IendB  =BOUNDS(ng) % IendB  (tile)
      IendP  =BOUNDS(ng) % IendP  (tile)
      IendR  =BOUNDS(ng) % IendR  (tile)
      IendT  =BOUNDS(ng) % IendT  (tile)
      Jstr   =BOUNDS(ng) % Jstr   (tile)
      JstrB  =BOUNDS(ng) % JstrB  (tile)
      JstrM  =BOUNDS(ng) % JstrM  (tile)
      JstrP  =BOUNDS(ng) % JstrP  (tile)
      JstrR  =BOUNDS(ng) % JstrR  (tile)
      JstrT  =BOUNDS(ng) % JstrT  (tile)
      JstrV  =BOUNDS(ng) % JstrV  (tile)
      Jend   =BOUNDS(ng) % Jend   (tile)
      JendB  =BOUNDS(ng) % JendB  (tile)
      JendP  =BOUNDS(ng) % JendP  (tile)
      JendR  =BOUNDS(ng) % JendR  (tile)
      JendT  =BOUNDS(ng) % JendT  (tile)
!
      Istrm3 =BOUNDS(ng) % Istrm3 (tile)            ! Istr-3
      Istrm2 =BOUNDS(ng) % Istrm2 (tile)            ! Istr-2
      Istrm1 =BOUNDS(ng) % Istrm1 (tile)            ! Istr-1
      IstrUm2=BOUNDS(ng) % IstrUm2(tile)            ! IstrU-2
      IstrUm1=BOUNDS(ng) % IstrUm1(tile)            ! IstrU-1
      Iendp1 =BOUNDS(ng) % Iendp1 (tile)            ! Iend+1
      Iendp2 =BOUNDS(ng) % Iendp2 (tile)            ! Iend+2
      Iendp2i=BOUNDS(ng) % Iendp2i(tile)            ! Iend+2 interior
      Iendp3 =BOUNDS(ng) % Iendp3 (tile)            ! Iend+3
      Jstrm3 =BOUNDS(ng) % Jstrm3 (tile)            ! Jstr-3
      Jstrm2 =BOUNDS(ng) % Jstrm2 (tile)            ! Jstr-2
      Jstrm1 =BOUNDS(ng) % Jstrm1 (tile)            ! Jstr-1
      JstrVm2=BOUNDS(ng) % JstrVm2(tile)            ! JstrV-2
      JstrVm1=BOUNDS(ng) % JstrVm1(tile)            ! JstrV-1
      Jendp1 =BOUNDS(ng) % Jendp1 (tile)            ! Jend+1
      Jendp2 =BOUNDS(ng) % Jendp2 (tile)            ! Jend+2
      Jendp2i=BOUNDS(ng) % Jendp2i(tile)            ! Jend+2 interior
      Jendp3 =BOUNDS(ng) % Jendp3 (tile)            ! Jend+3
!
!-----------------------------------------------------------------------
!  Initial conditions for 2D momentum (m/s) components.
!-----------------------------------------------------------------------
!
      DO j=JstrT,JendT
        DO i=IstrP,IendT
          ubar(i,j,1)=0.0_r8
        END DO
      END DO
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          vbar(i,j,1)=0.0_r8
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Initial conditions for free-surface (m).
!-----------------------------------------------------------------------
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          zeta(i,j,1)=0.0_r8
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Initial conditions for 3D momentum components (m/s).
!-----------------------------------------------------------------------
!
      DO k=1,N(ng)
       DO j=JstrT,JendT
         DO i=IstrP,IendT
            u(i,j,k,1)=0.0_r8
          END DO
        END DO
        DO j=JstrP,JendT
          DO i=IstrT,IendT
            v(i,j,k,1)=0.0_r8
          END DO
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Initial conditions for tracer type variables.
!-----------------------------------------------------------------------
!
!  Set initial conditions for potential temperature (Celsius) and
!  salinity (PSU).
!
      DO k=1,N(ng)
        DO j=JstrR,JendR
          DO i=IstrR,IendR
            t(i,j,k,1,itemp)=-1.9_r8
            t(i,j,k,1,isalt)=34.4_r8
            t(i,j,k,2,itemp)=t(i,j,k,1,itemp)
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt)
          END DO
        END DO
      END DO
      RETURN
      END SUBROUTINE ana_NLMinitial_tile
      SUBROUTINE ana_smflux (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine sets kinematic surface momentum flux (wind stress)     !
!  "sustr" and "svstr" (m2/s2) using an analytical expression.         !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
!
!  Set horizontal starting and ending indices for automatic private
!  storage arrays.
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
!  Determine array lower and upper bounds in the I- and J-directions.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
!  Set array lower and upper bounds for MIN(I,J) directions and
!  MAX(I,J) directions.
!
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
!
      CALL ana_smflux_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      GRID(ng) % angler,                          &
     &                      GRID(ng) % lonr,                            &
     &                      GRID(ng) % latr,                            &
     &                      FORCES(ng) % sustr,                         &
     &                      FORCES(ng) % svstr)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        ANANAME(24)="ROMS/Functionals/ana_smflux.h"
      END IF
      RETURN
      END SUBROUTINE ana_smflux
!
!***********************************************************************
      SUBROUTINE ana_smflux_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            angler,                               &
     &                            lonr, latr,                           &
     &                            sustr, svstr)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod
      USE mp_exchange_mod, ONLY : mp_exchange2d
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
      real(r8), intent(in) :: angler(LBi:,LBj:)
      real(r8), intent(in) :: lonr(LBi:,LBj:)
      real(r8), intent(in) :: latr(LBi:,LBj:)
      real(r8), intent(out) :: sustr(LBi:,LBj:)
      real(r8), intent(out) :: svstr(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: i, j
      real(r8) :: Ewind, Nwind, cff, val1, val2, windamp, winddir
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrB, IstrP, IstrR, IstrT, IstrM, IstrU
      integer :: Iend, IendB, IendP, IendR, IendT
      integer :: Jstr, JstrB, JstrP, JstrR, JstrT, JstrM, JstrV
      integer :: Jend, JendB, JendP, JendR, JendT
      integer :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer :: Jendp1, Jendp2, Jendp2i, Jendp3
!
      Istr   =BOUNDS(ng) % Istr   (tile)
      IstrB  =BOUNDS(ng) % IstrB  (tile)
      IstrM  =BOUNDS(ng) % IstrM  (tile)
      IstrP  =BOUNDS(ng) % IstrP  (tile)
      IstrR  =BOUNDS(ng) % IstrR  (tile)
      IstrT  =BOUNDS(ng) % IstrT  (tile)
      IstrU  =BOUNDS(ng) % IstrU  (tile)
      Iend   =BOUNDS(ng) % Iend   (tile)
      IendB  =BOUNDS(ng) % IendB  (tile)
      IendP  =BOUNDS(ng) % IendP  (tile)
      IendR  =BOUNDS(ng) % IendR  (tile)
      IendT  =BOUNDS(ng) % IendT  (tile)
      Jstr   =BOUNDS(ng) % Jstr   (tile)
      JstrB  =BOUNDS(ng) % JstrB  (tile)
      JstrM  =BOUNDS(ng) % JstrM  (tile)
      JstrP  =BOUNDS(ng) % JstrP  (tile)
      JstrR  =BOUNDS(ng) % JstrR  (tile)
      JstrT  =BOUNDS(ng) % JstrT  (tile)
      JstrV  =BOUNDS(ng) % JstrV  (tile)
      Jend   =BOUNDS(ng) % Jend   (tile)
      JendB  =BOUNDS(ng) % JendB  (tile)
      JendP  =BOUNDS(ng) % JendP  (tile)
      JendR  =BOUNDS(ng) % JendR  (tile)
      JendT  =BOUNDS(ng) % JendT  (tile)
!
      Istrm3 =BOUNDS(ng) % Istrm3 (tile)            ! Istr-3
      Istrm2 =BOUNDS(ng) % Istrm2 (tile)            ! Istr-2
      Istrm1 =BOUNDS(ng) % Istrm1 (tile)            ! Istr-1
      IstrUm2=BOUNDS(ng) % IstrUm2(tile)            ! IstrU-2
      IstrUm1=BOUNDS(ng) % IstrUm1(tile)            ! IstrU-1
      Iendp1 =BOUNDS(ng) % Iendp1 (tile)            ! Iend+1
      Iendp2 =BOUNDS(ng) % Iendp2 (tile)            ! Iend+2
      Iendp2i=BOUNDS(ng) % Iendp2i(tile)            ! Iend+2 interior
      Iendp3 =BOUNDS(ng) % Iendp3 (tile)            ! Iend+3
      Jstrm3 =BOUNDS(ng) % Jstrm3 (tile)            ! Jstr-3
      Jstrm2 =BOUNDS(ng) % Jstrm2 (tile)            ! Jstr-2
      Jstrm1 =BOUNDS(ng) % Jstrm1 (tile)            ! Jstr-1
      JstrVm2=BOUNDS(ng) % JstrVm2(tile)            ! JstrV-2
      JstrVm1=BOUNDS(ng) % JstrVm1(tile)            ! JstrV-1
      Jendp1 =BOUNDS(ng) % Jendp1 (tile)            ! Jend+1
      Jendp2 =BOUNDS(ng) % Jendp2 (tile)            ! Jend+2
      Jendp2i=BOUNDS(ng) % Jendp2i(tile)            ! Jend+2 interior
      Jendp3 =BOUNDS(ng) % Jendp3 (tile)            ! Jend+3
!
!-----------------------------------------------------------------------
!  Set kinematic surface momentum flux (wind stress) component in the
!  XI-direction (m2/s2) at horizontal U-points.
!-----------------------------------------------------------------------
!
      DO j=JstrR,JendR
        DO i=Istr,IendR
          sustr(i,j)=0.0_r8
        END DO
      END DO
!
!-----------------------------------------------------------------------
!  Set kinematic surface momentum flux (wind stress) component in the
!  ETA-direction (m2/s2) at horizontal V-points.
!-----------------------------------------------------------------------
!
      DO j=JstrP,JendT
        DO i=IstrT,IendT
          svstr(i,j)=0.0_r8
        END DO
      END DO
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_u2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          sustr)
        CALL exchange_v2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          svstr)
      END IF
      CALL mp_exchange2d (ng, tile, model, 2,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    sustr, svstr)
      RETURN
      END SUBROUTINE ana_smflux_tile
      SUBROUTINE ana_srflux (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This subroutine sets kinematic surface solar shortwave radiation    !
!  flux "srflx" (degC m/s) using an analytical expression.             !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
!
!  Set horizontal starting and ending indices for automatic private
!  storage arrays.
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
!  Determine array lower and upper bounds in the I- and J-directions.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
!  Set array lower and upper bounds for MIN(I,J) directions and
!  MAX(I,J) directions.
!
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
!
      CALL ana_srflux_tile (ng, tile, model,                            &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      GRID(ng) % lonr,                            &
     &                      GRID(ng) % latr,                            &
     &                      FORCES(ng) % srflx)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        ANANAME(27)="ROMS/Functionals/ana_srflux.h"
      END IF
      RETURN
      END SUBROUTINE ana_srflux
!
!***********************************************************************
      SUBROUTINE ana_srflux_tile (ng, tile, model,                      &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            lonr, latr,                           &
     &                            srflx)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
      USE mp_exchange_mod, ONLY : mp_exchange2d
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
      real(r8), intent(in) :: lonr(LBi:,LBj:)
      real(r8), intent(in) :: latr(LBi:,LBj:)
      real(r8), intent(out) :: srflx(LBi:,LBj:)
!
!  Local variable declarations.
!
      integer :: i, j
      real(r8) :: cff
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrB, IstrP, IstrR, IstrT, IstrM, IstrU
      integer :: Iend, IendB, IendP, IendR, IendT
      integer :: Jstr, JstrB, JstrP, JstrR, JstrT, JstrM, JstrV
      integer :: Jend, JendB, JendP, JendR, JendT
      integer :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer :: Jendp1, Jendp2, Jendp2i, Jendp3
!
      Istr   =BOUNDS(ng) % Istr   (tile)
      IstrB  =BOUNDS(ng) % IstrB  (tile)
      IstrM  =BOUNDS(ng) % IstrM  (tile)
      IstrP  =BOUNDS(ng) % IstrP  (tile)
      IstrR  =BOUNDS(ng) % IstrR  (tile)
      IstrT  =BOUNDS(ng) % IstrT  (tile)
      IstrU  =BOUNDS(ng) % IstrU  (tile)
      Iend   =BOUNDS(ng) % Iend   (tile)
      IendB  =BOUNDS(ng) % IendB  (tile)
      IendP  =BOUNDS(ng) % IendP  (tile)
      IendR  =BOUNDS(ng) % IendR  (tile)
      IendT  =BOUNDS(ng) % IendT  (tile)
      Jstr   =BOUNDS(ng) % Jstr   (tile)
      JstrB  =BOUNDS(ng) % JstrB  (tile)
      JstrM  =BOUNDS(ng) % JstrM  (tile)
      JstrP  =BOUNDS(ng) % JstrP  (tile)
      JstrR  =BOUNDS(ng) % JstrR  (tile)
      JstrT  =BOUNDS(ng) % JstrT  (tile)
      JstrV  =BOUNDS(ng) % JstrV  (tile)
      Jend   =BOUNDS(ng) % Jend   (tile)
      JendB  =BOUNDS(ng) % JendB  (tile)
      JendP  =BOUNDS(ng) % JendP  (tile)
      JendR  =BOUNDS(ng) % JendR  (tile)
      JendT  =BOUNDS(ng) % JendT  (tile)
!
      Istrm3 =BOUNDS(ng) % Istrm3 (tile)            ! Istr-3
      Istrm2 =BOUNDS(ng) % Istrm2 (tile)            ! Istr-2
      Istrm1 =BOUNDS(ng) % Istrm1 (tile)            ! Istr-1
      IstrUm2=BOUNDS(ng) % IstrUm2(tile)            ! IstrU-2
      IstrUm1=BOUNDS(ng) % IstrUm1(tile)            ! IstrU-1
      Iendp1 =BOUNDS(ng) % Iendp1 (tile)            ! Iend+1
      Iendp2 =BOUNDS(ng) % Iendp2 (tile)            ! Iend+2
      Iendp2i=BOUNDS(ng) % Iendp2i(tile)            ! Iend+2 interior
      Iendp3 =BOUNDS(ng) % Iendp3 (tile)            ! Iend+3
      Jstrm3 =BOUNDS(ng) % Jstrm3 (tile)            ! Jstr-3
      Jstrm2 =BOUNDS(ng) % Jstrm2 (tile)            ! Jstr-2
      Jstrm1 =BOUNDS(ng) % Jstrm1 (tile)            ! Jstr-1
      JstrVm2=BOUNDS(ng) % JstrVm2(tile)            ! JstrV-2
      JstrVm1=BOUNDS(ng) % JstrVm1(tile)            ! JstrV-1
      Jendp1 =BOUNDS(ng) % Jendp1 (tile)            ! Jend+1
      Jendp2 =BOUNDS(ng) % Jendp2 (tile)            ! Jend+2
      Jendp2i=BOUNDS(ng) % Jendp2i(tile)            ! Jend+2 interior
      Jendp3 =BOUNDS(ng) % Jendp3 (tile)            ! Jend+3
!
!-----------------------------------------------------------------------
!  Set incoming solar shortwave radiation (degC m/s).  Usually, the
!  shortwave radiation from input files is Watts/m2 and then converted
!  to degC m/s by multiplying by conversion factor 1/(rho0*Cp) during
!  reading (Fscale). However, we are already inside ROMS kernel here
!  and all the fluxes are kinematic so shortwave radiation units need
!  to be degC m/s.
!-----------------------------------------------------------------------
!
      cff=1.0_r8/(rho0*cp)
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          srflx(i,j)=0.0_r8
        END DO
      END DO
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          srflx)
      END IF
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    srflx)
      RETURN
      END SUBROUTINE ana_srflux_tile
      SUBROUTINE ana_stflux (ng, tile, model, itrc)
!
!=======================================================================
!                                                                      !
!  This routine sets kinematic surface flux of tracer type variables   !
!  "stflx" (tracer units m/s) using analytical expressions.            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_forces
      USE mod_ncparam
!
! Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, itrc
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
!
!  Set horizontal starting and ending indices for automatic private
!  storage arrays.
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
!  Determine array lower and upper bounds in the I- and J-directions.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
!  Set array lower and upper bounds for MIN(I,J) directions and
!  MAX(I,J) directions.
!
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
!
      CALL ana_stflux_tile (ng, tile, model, itrc,                      &
     &                      LBi, UBi, LBj, UBj,                         &
     &                      IminS, ImaxS, JminS, JmaxS,                 &
     &                      FORCES(ng) % srflx,                         &
     &                      FORCES(ng) % stflx)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        ANANAME(31)="ROMS/Functionals/ana_stflux.h"
      END IF
      RETURN
      END SUBROUTINE ana_stflux
!
!***********************************************************************
      SUBROUTINE ana_stflux_tile (ng, tile, model, itrc,                &
     &                            LBi, UBi, LBj, UBj,                   &
     &                            IminS, ImaxS, JminS, JmaxS,           &
     &                            srflx,                                &
     &                            stflx)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
      USE mp_exchange_mod, ONLY : mp_exchange2d
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model, itrc
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
      real(r8), intent(in) :: srflx(LBi:,LBj:)
      real(r8), intent(inout) :: stflx(LBi:,LBj:,:)
!
!  Local variable declarations.
!
      integer :: i, j
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrB, IstrP, IstrR, IstrT, IstrM, IstrU
      integer :: Iend, IendB, IendP, IendR, IendT
      integer :: Jstr, JstrB, JstrP, JstrR, JstrT, JstrM, JstrV
      integer :: Jend, JendB, JendP, JendR, JendT
      integer :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer :: Jendp1, Jendp2, Jendp2i, Jendp3
!
      Istr   =BOUNDS(ng) % Istr   (tile)
      IstrB  =BOUNDS(ng) % IstrB  (tile)
      IstrM  =BOUNDS(ng) % IstrM  (tile)
      IstrP  =BOUNDS(ng) % IstrP  (tile)
      IstrR  =BOUNDS(ng) % IstrR  (tile)
      IstrT  =BOUNDS(ng) % IstrT  (tile)
      IstrU  =BOUNDS(ng) % IstrU  (tile)
      Iend   =BOUNDS(ng) % Iend   (tile)
      IendB  =BOUNDS(ng) % IendB  (tile)
      IendP  =BOUNDS(ng) % IendP  (tile)
      IendR  =BOUNDS(ng) % IendR  (tile)
      IendT  =BOUNDS(ng) % IendT  (tile)
      Jstr   =BOUNDS(ng) % Jstr   (tile)
      JstrB  =BOUNDS(ng) % JstrB  (tile)
      JstrM  =BOUNDS(ng) % JstrM  (tile)
      JstrP  =BOUNDS(ng) % JstrP  (tile)
      JstrR  =BOUNDS(ng) % JstrR  (tile)
      JstrT  =BOUNDS(ng) % JstrT  (tile)
      JstrV  =BOUNDS(ng) % JstrV  (tile)
      Jend   =BOUNDS(ng) % Jend   (tile)
      JendB  =BOUNDS(ng) % JendB  (tile)
      JendP  =BOUNDS(ng) % JendP  (tile)
      JendR  =BOUNDS(ng) % JendR  (tile)
      JendT  =BOUNDS(ng) % JendT  (tile)
!
      Istrm3 =BOUNDS(ng) % Istrm3 (tile)            ! Istr-3
      Istrm2 =BOUNDS(ng) % Istrm2 (tile)            ! Istr-2
      Istrm1 =BOUNDS(ng) % Istrm1 (tile)            ! Istr-1
      IstrUm2=BOUNDS(ng) % IstrUm2(tile)            ! IstrU-2
      IstrUm1=BOUNDS(ng) % IstrUm1(tile)            ! IstrU-1
      Iendp1 =BOUNDS(ng) % Iendp1 (tile)            ! Iend+1
      Iendp2 =BOUNDS(ng) % Iendp2 (tile)            ! Iend+2
      Iendp2i=BOUNDS(ng) % Iendp2i(tile)            ! Iend+2 interior
      Iendp3 =BOUNDS(ng) % Iendp3 (tile)            ! Iend+3
      Jstrm3 =BOUNDS(ng) % Jstrm3 (tile)            ! Jstr-3
      Jstrm2 =BOUNDS(ng) % Jstrm2 (tile)            ! Jstr-2
      Jstrm1 =BOUNDS(ng) % Jstrm1 (tile)            ! Jstr-1
      JstrVm2=BOUNDS(ng) % JstrVm2(tile)            ! JstrV-2
      JstrVm1=BOUNDS(ng) % JstrVm1(tile)            ! JstrV-1
      Jendp1 =BOUNDS(ng) % Jendp1 (tile)            ! Jend+1
      Jendp2 =BOUNDS(ng) % Jendp2 (tile)            ! Jend+2
      Jendp2i=BOUNDS(ng) % Jendp2i(tile)            ! Jend+2 interior
      Jendp3 =BOUNDS(ng) % Jendp3 (tile)            ! Jend+3
!
!-----------------------------------------------------------------------
!  Set kinematic surface heat flux (degC m/s) at horizontal
!  RHO-points.
!-----------------------------------------------------------------------
!
      IF (itrc.eq.itemp) THEN
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            stflx(i,j,itrc)=0.0_r8
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Set kinematic surface freshwater flux (m/s) at horizontal
!  RHO-points, scaling by surface salinity is done in STEP3D.
!-----------------------------------------------------------------------
!
      ELSE IF (itrc.eq.isalt) THEN
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            stflx(i,j,itrc)=0.0_r8
          END DO
        END DO
!
!-----------------------------------------------------------------------
!  Set kinematic surface flux (T m/s) of passive tracers, if any.
!-----------------------------------------------------------------------
!
      ELSE
        DO j=JstrT,JendT
          DO i=IstrT,IendT
            stflx(i,j,itrc)=0.0_r8
          END DO
        END DO
      END IF
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          stflx(:,:,itrc))
      END IF
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
     &                    stflx(:,:,itrc))
      RETURN
      END SUBROUTINE ana_stflux_tile
      SUBROUTINE ana_nudgcoef (ng, tile, model)
!
!=======================================================================
!                                                                      !
!  This routine set nudging coefficients time-scales (1/s).            !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_ncparam
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
!
!  Local variable declarations.
!
      integer :: IminS, ImaxS, JminS, JmaxS
      integer :: LBi, UBi, LBj, UBj, LBij, UBij
!
!  Set horizontal starting and ending indices for automatic private
!  storage arrays.
!
      IminS=BOUNDS(ng)%Istr(tile)-3
      ImaxS=BOUNDS(ng)%Iend(tile)+3
      JminS=BOUNDS(ng)%Jstr(tile)-3
      JmaxS=BOUNDS(ng)%Jend(tile)+3
!
!  Determine array lower and upper bounds in the I- and J-directions.
!
      LBi=BOUNDS(ng)%LBi(tile)
      UBi=BOUNDS(ng)%UBi(tile)
      LBj=BOUNDS(ng)%LBj(tile)
      UBj=BOUNDS(ng)%UBj(tile)
!
!  Set array lower and upper bounds for MIN(I,J) directions and
!  MAX(I,J) directions.
!
      LBij=BOUNDS(ng)%LBij
      UBij=BOUNDS(ng)%UBij
!
      CALL ana_nudgcoef_tile (ng, tile, model,                          &
     &                        LBi, UBi, LBj, UBj,                       &
     &                        IminS, ImaxS, JminS, JmaxS)
!
! Set analytical header file name used.
!
      IF (Lanafile) THEN
        ANANAME(16)="ROMS/Functionals/ana_nudgcoef.h"
      END IF
      RETURN
      END SUBROUTINE ana_nudgcoef
!
!***********************************************************************
      SUBROUTINE ana_nudgcoef_tile (ng, tile, model,                    &
     &                              LBi, UBi, LBj, UBj,                 &
     &                              IminS, ImaxS, JminS, JmaxS)
!***********************************************************************
!
      USE mod_param
      USE mod_parallel
      USE mod_boundary
      USE mod_grid
      USE mod_ncparam
      USE mod_scalars
!
      USE distribute_mod, ONLY : mp_collect
      USE mp_exchange_mod, ONLY : mp_exchange2d
      USE mp_exchange_mod, ONLY : mp_exchange3d
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
!  Local variable declarations.
!
      integer :: Iwrk, i, itrc, j
      real(r8) :: cff1, cff2, cff3
      real(r8), parameter :: IniVal = 0.0_r8
      real(r8), dimension(IminS:ImaxS,JminS:JmaxS) :: wrk
!
!-----------------------------------------------------------------------
!  Set lower and upper tile bounds and staggered variables bounds for
!  this horizontal domain partition.  Notice that if tile=-1, it will
!  set the values for the global grid.
!-----------------------------------------------------------------------
!
      integer :: Istr, IstrB, IstrP, IstrR, IstrT, IstrM, IstrU
      integer :: Iend, IendB, IendP, IendR, IendT
      integer :: Jstr, JstrB, JstrP, JstrR, JstrT, JstrM, JstrV
      integer :: Jend, JendB, JendP, JendR, JendT
      integer :: Istrm3, Istrm2, Istrm1, IstrUm2, IstrUm1
      integer :: Iendp1, Iendp2, Iendp2i, Iendp3
      integer :: Jstrm3, Jstrm2, Jstrm1, JstrVm2, JstrVm1
      integer :: Jendp1, Jendp2, Jendp2i, Jendp3
!
      Istr   =BOUNDS(ng) % Istr   (tile)
      IstrB  =BOUNDS(ng) % IstrB  (tile)
      IstrM  =BOUNDS(ng) % IstrM  (tile)
      IstrP  =BOUNDS(ng) % IstrP  (tile)
      IstrR  =BOUNDS(ng) % IstrR  (tile)
      IstrT  =BOUNDS(ng) % IstrT  (tile)
      IstrU  =BOUNDS(ng) % IstrU  (tile)
      Iend   =BOUNDS(ng) % Iend   (tile)
      IendB  =BOUNDS(ng) % IendB  (tile)
      IendP  =BOUNDS(ng) % IendP  (tile)
      IendR  =BOUNDS(ng) % IendR  (tile)
      IendT  =BOUNDS(ng) % IendT  (tile)
      Jstr   =BOUNDS(ng) % Jstr   (tile)
      JstrB  =BOUNDS(ng) % JstrB  (tile)
      JstrM  =BOUNDS(ng) % JstrM  (tile)
      JstrP  =BOUNDS(ng) % JstrP  (tile)
      JstrR  =BOUNDS(ng) % JstrR  (tile)
      JstrT  =BOUNDS(ng) % JstrT  (tile)
      JstrV  =BOUNDS(ng) % JstrV  (tile)
      Jend   =BOUNDS(ng) % Jend   (tile)
      JendB  =BOUNDS(ng) % JendB  (tile)
      JendP  =BOUNDS(ng) % JendP  (tile)
      JendR  =BOUNDS(ng) % JendR  (tile)
      JendT  =BOUNDS(ng) % JendT  (tile)
!
      Istrm3 =BOUNDS(ng) % Istrm3 (tile)            ! Istr-3
      Istrm2 =BOUNDS(ng) % Istrm2 (tile)            ! Istr-2
      Istrm1 =BOUNDS(ng) % Istrm1 (tile)            ! Istr-1
      IstrUm2=BOUNDS(ng) % IstrUm2(tile)            ! IstrU-2
      IstrUm1=BOUNDS(ng) % IstrUm1(tile)            ! IstrU-1
      Iendp1 =BOUNDS(ng) % Iendp1 (tile)            ! Iend+1
      Iendp2 =BOUNDS(ng) % Iendp2 (tile)            ! Iend+2
      Iendp2i=BOUNDS(ng) % Iendp2i(tile)            ! Iend+2 interior
      Iendp3 =BOUNDS(ng) % Iendp3 (tile)            ! Iend+3
      Jstrm3 =BOUNDS(ng) % Jstrm3 (tile)            ! Jstr-3
      Jstrm2 =BOUNDS(ng) % Jstrm2 (tile)            ! Jstr-2
      Jstrm1 =BOUNDS(ng) % Jstrm1 (tile)            ! Jstr-1
      JstrVm2=BOUNDS(ng) % JstrVm2(tile)            ! JstrV-2
      JstrVm1=BOUNDS(ng) % JstrVm1(tile)            ! JstrV-1
      Jendp1 =BOUNDS(ng) % Jendp1 (tile)            ! Jend+1
      Jendp2 =BOUNDS(ng) % Jendp2 (tile)            ! Jend+2
      Jendp2i=BOUNDS(ng) % Jendp2i(tile)            ! Jend+2 interior
      Jendp3 =BOUNDS(ng) % Jendp3 (tile)            ! Jend+3
!
!-----------------------------------------------------------------------
!  Set up nudging towards data time-scale coefficients (1/s).
!-----------------------------------------------------------------------
!
!  Initialize.
!
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          wrk(i,j)=0.0_r8
        END DO
      END DO
!
!  Default nudging coefficients.  Set nudging coefficients uniformly to
!  the values specified in the standard input file.
!
!
!-----------------------------------------------------------------------
!  Set nudging coefficients (1/s) for passive/active (outflow/inflow)
!  open boundary conditions.  Weak nudging is expected in passive
!  outflow conditions and strong nudging is expected in active inflow
!  conditions.  Notice that interior nudging coefficient defined
!  above are zero out when boundary condition nudging.  The USER needs
!  to adapt this to his/her application!
!-----------------------------------------------------------------------
!
      IF (NudgingCoeff(ng)) THEN
!
!  Free-surface nudging coefficients.
!
        IF (LBC(iwest,isFsur,ng)%nudging) THEN
          IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
            FSobc_out(ng,iwest)=Znudg(ng)
            FSobc_in (ng,iwest)=obcfac(ng)*Znudg(ng)
          END IF
        END IF
!
        IF (LBC(ieast,isFsur,ng)%nudging) THEN
          IF (DOMAIN(ng)%NorthEast_Test(tile)) THEN
            FSobc_out(ng,ieast)=Znudg(ng)
            FSobc_in (ng,ieast)=obcfac(ng)*Znudg(ng)
          END IF
        END IF
!
        IF (LBC(isouth,isFsur,ng)%nudging) THEN
          IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
            FSobc_out(ng,isouth)=Znudg(ng)
            FSobc_in (ng,isouth)=obcfac(ng)*Znudg(ng)
          END IF
        END IF
!
        IF (LBC(inorth,isFsur,ng)%nudging) THEN
          IF (DOMAIN(ng)%NorthEast_Test(tile)) THEN
            FSobc_out(ng,inorth)=Znudg(ng)
            FSobc_in (ng,inorth)=obcfac(ng)*Znudg(ng)
          END IF
        END IF
!
!  2D momentum nudging coefficients.
!
        IF (LBC(iwest,isUbar,ng)%nudging.or.                            &
     &      LBC(iwest,isVbar,ng)%nudging) THEN
          IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
            M2obc_out(ng,iwest)=M2nudg(ng)
            M2obc_in (ng,iwest)=obcfac(ng)*M2nudg(ng)
          END IF
        END IF
!
        IF (LBC(ieast,isUbar,ng)%nudging.or.                            &
     &      LBC(ieast,isVbar,ng)%nudging) THEN
          IF (DOMAIN(ng)%NorthEast_Test(tile)) THEN
            M2obc_out(ng,ieast)=M2nudg(ng)
            M2obc_in (ng,ieast)=obcfac(ng)*M2nudg(ng)
          END IF
        END IF
!
        IF (LBC(isouth,isUbar,ng)%nudging.or.                           &
     &      LBC(isouth,isVbar,ng)%nudging) THEN
          IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
            M2obc_out(ng,isouth)=M2nudg(ng)
            M2obc_in (ng,isouth)=obcfac(ng)*M2nudg(ng)
          END IF
        END IF
!
        IF (LBC(inorth,isUbar,ng)%nudging.or.                           &
     &      LBC(inorth,isVbar,ng)%nudging) THEN
          IF (DOMAIN(ng)%NorthEast_Test(tile)) THEN
            M2obc_out(ng,inorth)=M2nudg(ng)
            M2obc_in (ng,inorth)=obcfac(ng)*M2nudg(ng)
          END IF
        END IF
!
!  Tracers nudging coefficients.
!
        DO itrc=1,NT(ng)
          IF (LBC(iwest,isTvar(itrc),ng)%nudging) THEN
            IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
              Tobc_out(itrc,ng,iwest)=Tnudg(itrc,ng)
              Tobc_in (itrc,ng,iwest)=obcfac(ng)*Tnudg(itrc,ng)
            END IF
          END IF
        END DO
!
        DO itrc=1,NT(ng)
          IF (LBC(ieast,isTvar(itrc),ng)%nudging) THEN
            IF (DOMAIN(ng)%NorthEast_Test(tile)) THEN
              Tobc_out(itrc,ng,ieast)=Tnudg(itrc,ng)
              Tobc_in (itrc,ng,ieast)=obcfac(ng)*Tnudg(itrc,ng)
            END IF
          END IF
        END DO
!
        DO itrc=1,NT(ng)
          IF (LBC(isouth,isTvar(itrc),ng)%nudging) THEN
            IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
              Tobc_out(itrc,ng,isouth)=Tnudg(itrc,ng)
              Tobc_in (itrc,ng,isouth)=obcfac(ng)*Tnudg(itrc,ng)
            END IF
          END IF
        END DO
!
        DO itrc=1,NT(ng)
          IF (LBC(inorth,isTvar(itrc),ng)%nudging) THEN
            IF (DOMAIN(ng)%NorthEast_Test(tile)) THEN
              Tobc_out(itrc,ng,inorth)=Tnudg(itrc,ng)
              Tobc_in (itrc,ng,inorth)=obcfac(ng)*Tnudg(itrc,ng)
            END IF
          END IF
        END DO
!
!  3D momentum nudging coefficients.
!
        IF (LBC(iwest,isUvel,ng)%nudging.or.                            &
     &      LBC(iwest,isVvel,ng)%nudging) THEN
          IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
            M3obc_out(ng,iwest)=M3nudg(ng)
            M3obc_in (ng,iwest)=obcfac(ng)*M3nudg(ng)
          END IF
        END IF
!
        IF (LBC(ieast,isUvel,ng)%nudging.or.                            &
     &      LBC(ieast,isVvel,ng)%nudging) THEN
          IF (DOMAIN(ng)%NorthEast_Test(tile)) THEN
            M3obc_out(ng,ieast)=M3nudg(ng)
            M3obc_in (ng,ieast)=obcfac(ng)*M3nudg(ng)
          END IF
        END IF
!
        IF (LBC(isouth,isUvel,ng)%nudging.or.                           &
     &      LBC(isouth,isVvel,ng)%nudging) THEN
          IF (DOMAIN(ng)%SouthWest_Test(tile)) THEN
            M3obc_out(ng,isouth)=M3nudg(ng)
            M3obc_in (ng,isouth)=obcfac(ng)*M3nudg(ng)
          END IF
        END IF
!
        IF (LBC(inorth,isUvel,ng)%nudging.or.                           &
     &      LBC(inorth,isVvel,ng)%nudging) THEN
          IF (DOMAIN(ng)%NorthEast_Test(tile)) THEN
            M3obc_out(ng,inorth)=M3nudg(ng)
            M3obc_in (ng,inorth)=obcfac(ng)*M3nudg(ng)
          END IF
        END IF
      END IF
      RETURN
      END SUBROUTINE ana_nudgcoef_tile
      END MODULE analytical_mod