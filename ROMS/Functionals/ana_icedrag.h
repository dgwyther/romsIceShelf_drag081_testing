      SUBROUTINE ana_icedrag (ng, tile, model)
!
!! svn $Id$
!!======================================================================
!! Copyright (c) 2002-2013 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!================================================== David E. Gwyther ===
!                                                                      !
!  This routine sets analytical, spatially varying ice bottom roughness!
!  length (m), or linear drag coefficients (m/s), or quadratic drag    !
!  coefficients (nondimensional) at RHO-points.                        !
!                                                                      !
! Currently, only UV_QDRAG is implemented.                             !
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

#include "tile.h"
!
      CALL ana_icedrag_tile (ng, tile, model,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
#if defined UV_LOGDRAG
     &                    GRID(ng) % ZoBot)
#elif defined UV_LDRAG
     &                    GRID(ng) % rdrag)
#elif defined UV_QDRAG
     &                    ICESHELFVAR(ng) % idrag2)
#endif
!
! Set analytical header file name used.
!
#ifdef DISTRIBUTE
      IF (Lanafile) THEN
#else
      IF (Lanafile.and.(tile.eq.0)) THEN
#endif
        ANANAME(38)=__FILE__
      END IF

      RETURN
      END SUBROUTINE ana_icedrag
!
!***********************************************************************
      SUBROUTINE ana_icedrag_tile (ng, tile, model,                        &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
#if defined UV_LOGDRAG
     &                          ZoBot)
#elif defined UV_LDRAG
     &                          rdrag)
#elif defined UV_QDRAG
     &                          idrag2)
#endif
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_scalars
!
      USE exchange_2d_mod, ONLY : exchange_r2d_tile
#ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
#endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
!
#ifdef ASSUMED_SHAPE
# if defined UV_LOGDRAG
      real(r8), intent(out) :: ZoBot(LBi:,LBj:)
# elif defined UV_LDRAG
      real(r8), intent(out) :: rdrag(LBi:,LBj:)
# elif defined UV_QDRAG
      real(r8), intent(out) :: idrag2(LBi:,LBj:)
# endif

#else

# if defined UV_LOGDRAG
      real(r8), intent(out) :: ZoBot(LBi:UBi,LBj:UBj)
# elif defined UV_LDRAG
      real(r8), intent(out) :: rdrag(LBi:UBi,LBj:UBj)
# elif defined UV_QDRAG
      real(r8), intent(out) :: idrag2(LBi:UBi,LBj:UBj)
# endif
#endif
!
!  Local variable declarations.
!
      integer :: i, j

      real(r8) :: cff

#include "set_bounds.h"
!
!-----------------------------------------------------------------------
#if defined UV_LOGDRAG
!  Set spatially varying bottom roughness length (m).
#elif defined UV_LDRAG
!  Set spatially varying linear drag coefficient (m/s).
#elif defined UV_QDRAG
!  Set spatially varying quadratic ice-ocean drag coefficient (nondimensional)
#endif
!-----------------------------------------------------------------------
!
#if defined UPWELLING
# if defined UV_LOGDRAG
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ZoBot(i,j)=0.05_r8*(1.0_r8+TANH(GRID(ng)%h(i,j)/50.0_r8))
        END DO
      END DO
# elif defined UV_LDRAG
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          rdrag(i,j)=0.002_r8*(1.0_r8-TANH(GRID(ng)%h(i,j)/150.0_r8))
        END DO
      END DO
# elif defined UV_QDRAG
      DO j=JstrT,JendT          ! based on Chezy coefficient (g/c^2)
        DO i=IstrT,IendT
          cff=1.8_r8*GRID(ng)%h(i,j)*LOG(GRID(ng)%h(i,j))
          idrag2(i,j)=g/(cff*cff)
        END DO
      END DO
# endif
#elif defined ICETEST
# if defined UV_LOGDRAG
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ZoBot(i,j)=0.05_r8*(1.0_r8+TANH(GRID(ng)%h(i,j)/50.0_r8))
            write(6,*) "UV_LOGDRAG not implemented with UV_ICEDRAG_GRID yet."
        END DO
      END DO
# elif defined UV_LDRAG
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          rdrag(i,j)=0.002_r8*(1.0_r8-TANH(GRID(ng)%h(i,j)/150.0_r8))
            write(6,*) "UV_LDRAG not implemented with UV_ICEDRAG_GRID yet."
        END DO
      END DO
# elif defined UV_QDRAG
      DO j=JstrT,JendT          ! based on Chezy coefficient (g/c^2)
        DO i=IstrT,IendT
          idrag2(i,j)=0.1_r8
        END DO
      END DO
# endif

#else
# if defined UV_LOGDRAG
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          ZoBot(i,j)=???
        END DO
      END DO
# elif defined UV_LDRAG
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          rdrag(i,j)=???
        END DO
      END DO
# elif defined UV_QDRAG
      DO j=JstrT,JendT
        DO i=IstrT,IendT
          idrag2(i,j)=???
        END DO
      END DO
# endif
#endif
!
!  Exchange boundary data.
!
      IF (EWperiodic(ng).or.NSperiodic(ng)) THEN
        CALL exchange_r2d_tile (ng, tile,                               &
     &                          LBi, UBi, LBj, UBj,                     &
#if defined UV_LOGDRAG
     &                          ZoBot)
#elif defined UV_LDRAG
     &                          rdrag)
#elif defined UV_QDRAG
     &                          idrag2)
#endif
      END IF

#ifdef DISTRIBUTE
      CALL mp_exchange2d (ng, tile, model, 1,                           &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    NghostPoints,                                 &
     &                    EWperiodic(ng), NSperiodic(ng),               &
# if defined UV_LOGDRAG
     &                    ZoBot)
# elif defined UV_LDRAG
     &                    rdrag)
# elif defined UV_QDRAG
     &                    idrag2)
# endif
#endif

      RETURN
      END SUBROUTINE ana_icedrag_tile
