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
     &                    nstp(ng), nnew(ng), nrhs(ng),                 &
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
     &                          nstp, nnew, nrhs,                       &
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
      integer, intent(in) :: nstp, nnew, nrhs
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
      real(r8) e_val, sum_expfact
      real(r8) min_expfact, max_expfact
      integer iter, idxw
      integer, SAVE :: idx_iter = 0
      real(r8) sum_VAL, max_increment, egrow
      LOGICAL DO_DEBUG
# include "set_bounds.h"
      DO_DEBUG = .FALSE.
      IF (DO_DEBUG) THEN
         sum_expfact = 0
         min_expfact =  40000
         max_expfact = -40000
         DO j=Jstr,Jend
            DO i=Istr,Iend
               e_val = expfact(i,j)
               sum_expfact = sum_expfact + e_val
               if (min_expfact .gt. e_val) min_expfact = e_val
               if (max_expfact .lt. e_val) max_expfact = e_val
            END DO
         END DO
         Print *, 'sum_expfact = ', sum_expfact
         Print *, 'min_expfact = ', min_expfact
         Print *, 'max_expfact = ', max_expfact
         Print *, 'dt=', dt(ng)
         DO itrc=1,NT(ng)
            Print *,'itrc=',itrc,' LtracerGrowth=',LtracerGrowth(itrc,ng)
         END DO
         Print *, 'nstp=', nstp
         Print *, 'nrhs=', nrhs
         Print *, 'nnew=', nnew
         Print *, 'idx_iter=', idx_iter
         idx_iter = idx_iter + 1
         DO iter=1,3
            sum_VAL=0
            DO j=Jstr,Jend
               DO i=Istr,Iend
                  DO k=1,N(ng)
                     sum_VAL = sum_VAL + t(i,j,k,iter,3)
                  END DO
               END DO
            END DO
            Print *, 'iter=', iter, ' sum_VAL=', sum_VAL
         END DO
      END IF
      idxw = nnew
      DO itrc=1,NT(ng)
        egrow = LtracerGrowth(itrc,ng)
!        egrow = 0
        IF (LtracerGrowth(itrc,ng) .gt.0) THEN
          max_increment = 0
          DO j=Jstr,Jend
            DO i=Istr,Iend
!              e_val = 1
              e_val = expfact(i,j)
              DO k=1,N(ng)
                increment = 1.0_r8 + egrow * e_val * dt(ng)
                IF (increment .gt. max_increment) THEN
                   max_increment = increment
                END IF
                t(i,j,k,idxw,itrc) = t(i,j,k,idxw,itrc) * increment
              END DO
            END DO
          END DO
        END IF
      END DO
      END SUBROUTINE simple_growth_tile
#endif
      END MODULE simple_growth_mod
