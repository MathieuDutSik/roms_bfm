#include "cppdefs.h"
      MODULE simple_growth_mod
#ifdef SIMPLE_GROWTH
      implicit none

      PRIVATE
      PUBLIC  :: simple_growth

      CONTAINS
!
!***********************************************************************
      SUBROUTINE simple_growth (ng, tile)
!***********************************************************************
!
      USE mod_param
      USE mod_forces
      USE mod_grid
      USE mod_mixing
      USE mod_ocean
      USE mod_stepping
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
# include "tile.h"
!
      CALL simple_growth_tile (ng, tile,                                &
     &                    LBi, UBi, LBj, UBj,                           &
     &                    IminS, ImaxS, JminS, JmaxS,                   &
     &                    nstp(ng),                                     &
# ifdef MASKING
     &                    GRID(ng) % rmask,                             &
# endif
     &                    GRID(ng) % expfact,                           &
     &                    OCEAN(ng) % t)
      RETURN
      END SUBROUTINE simple_growth
!
!***********************************************************************
      SUBROUTINE simple_growth_tile (ng, tile,                          &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nstp,                                   &
# ifdef MASKING
     &                          rmask,                                  &
# endif
     &                          expfact,                                &
     &                          t)
!***********************************************************************
!
      USE mod_param
      USE mod_scalars
!
      USE bc_2d_mod, ONLY : bc_r2d_tile
# ifdef DISTRIBUTE
      USE mp_exchange_mod, ONLY : mp_exchange2d
# endif
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp
!
# ifdef ASSUMED_SHAPE
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  endif
      real(r8), intent(in) :: expfact(LBi:,LBj:)
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: expfact(LBi:UBi,LBj:UBj)
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, j, k, itrc
      real(r8) :: increment
# include "set_bounds.h"
      DO itrc=1,NT(ng)
        IF (LtracerGrowth(itrc,ng) .gt.0) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              DO k=1,N(ng)
                increment = 1.0_r8 + LtracerGrowth(itrc,ng)             &
     &                    * expfact(i,j) * * dt(ng)
                t(i,j,k,nstp,itrc) = t(i,j,k,nstp,itrc) * increment
              END DO
            END DO
          END DO
        END IF
      END DO
      END SUBROUTINE simple_growth_tile
#endif
      END MODULE simple_growth_mod
