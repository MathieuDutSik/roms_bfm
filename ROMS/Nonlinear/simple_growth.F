#include "cppdefs.h"
      MODULE simple_growth
#ifdef SIMPLE_GROWTH
      implicit none

      PRIVATE
      PUBLIC  :: lmd_skpp

      CONTAINS
!
!***********************************************************************
      SUBROUTINE lmd_skpp (ng, tile)
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
     &                    OCEAN(ng) % t)
      RETURN
      END SUBROUTINE lmd_skpp
!
!***********************************************************************
      SUBROUTINE simple_growth_tile (ng, tile,                          &
     &                          LBi, UBi, LBj, UBj,                     &
     &                          IminS, ImaxS, JminS, JmaxS,             &
     &                          nstp,                                   &
# ifdef MASKING
     &                          rmask,                                  &
# endif
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
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
# else
#  ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  endif
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,N(ng),3,NT(ng))
# endif
!
!  Local variable declarations.
!
      integer :: i, itrc, j, k
# include "set_bounds.h"
      DO itrc=1,NT(ng)
        IF (LtracerGrowth(itrc,ng) .gt.0) THEN
          DO j=Jstr,Jend
            DO i=Istr,Iend
              t(i,j
            END DO
          END DO
        END IF
      END DO
      END SUBROUTINE simple_growth_tile
#endif
      END MODULE simple_growth_mod
