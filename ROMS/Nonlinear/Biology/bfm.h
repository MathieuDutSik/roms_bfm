!
! Geometrical stuff
!
      SUBROUTINE SPHERICAL_COORDINATE_DISTANCE(lon1,lon2, lat1,lat2, dist)
      USE mod_kinds
      IMPLICIT NONE
      REAL(r8), intent(in) :: lon1, lon2, lat1, lat2
      REAL(r8), intent(out) :: dist
      !
      REAL(r8) :: lon1_r, lon2_r, lat1_r, lat2_r
      REAL(r8) :: X1, Y1, Z1, X2, Y2, Z2
      REAL(r8) :: COEFF, scalprod
      COEFF = 3.14159_r8 / 180.0_r8
      lon1_r = lon1 * COEFF
      lon2_r = lon2 * COEFF
      lat1_r = lat1 * COEFF
      lat2_r = lat2 * COEFF
      X1 = COS(lon1_r) * COS(lat1_r)
      Y1 = SIN(lon1_r) * COS(lat1_r)
      Z1 = SIN(lat1_r)
      X2 = COS(lon2_r) * COS(lat2_r)
      Y2 = SIN(lon2_r) * COS(lat2_r)
      Z2 = SIN(lat2_r)
      scalprod = X1 * X2 + Y1 * Y2 + Z1 * Z2;
      if (scalprod .GE. 1) THEN
        dist = 0;
      ELSE
        dist = ACOS(scalprod)
      END IF
      END SUBROUTINE
!
      SUBROUTINE COMPUTE_SPHERICAL_AREA(lon1, lon2, lon3,               &
     &     lat1, lat2, lat3, area)
      USE mod_kinds
      IMPLICIT NONE
      REAL(r8), intent(in) :: LON1, LON2, LON3, LAT1, LAT2, LAT3
      REAL(r8), intent(out) :: AREA
      REAL(r8) :: DistA, DistB, DistC, DistS
      REAL(r8) :: eTan1, eTan2, eTan3, eTan4
      REAL(r8) :: eProd, sqrtProd
      CALL SPHERICAL_COORDINATE_DISTANCE(LON1, LON2, LAT1, LAT2, DistA)
      CALL SPHERICAL_COORDINATE_DISTANCE(LON1, LON3, LAT1, LAT3, DistB)
      CALL SPHERICAL_COORDINATE_DISTANCE(LON2, LON3, LAT2, LAT3, DistC)
      DistS=(DistA + DistB + DistC)/2.0_r8
      eTan1=tan(DistS/2.0_r8)
      eTan2=tan((DistS - DistA)/2.0_r8)
      eTan3=tan((DistS - DistB)/2.0_r8)
      eTan4=tan((DistS - DistC)/2.0_r8)
      eProd=eTan1*eTan2*eTan3*eTan4
      sqrtProd=SQRT(eProd)
      AREA=4.0_r8*ATAN(sqrtProd)
      END SUBROUTINE
!
!
!
      SUBROUTINE GetTime_for_daylength(ng, yday)
      USE mod_kinds
      USE mod_scalars, only : tdays
      IMPLICIT NONE
      INTEGER, intent(in) :: ng
      REAL(r8), intent(out) :: yday
      REAL(r8) tday_2017_01_01
      ! 1968-05-23 is by convention 0
      ! 2017-01-01 is thus 17755
      tday_2017_01_01 = 17755
      yday = tdays(ng) - tday_2017_01_01
      END SUBROUTINE
!
! Now the model part.
!
      SUBROUTINE biology (ng,tile)
!
!svn $Id: hypoxia_srm.h 864 2017-08-10 04:11:10Z arango $
!************************************************** Hernan G. Arango ***
!  Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!***********************************************************************
!                                                                      !
!  This routine computes the biological sources and sinks for the      !
!  Hypoxia Simple Respiration Model for dissolved oxygen. Then, it     !
!  adds those terms to the global biological fields.                   !
!                                                                      !
!***********************************************************************
!
      USE mod_param
#ifdef DIAGNOSTICS_BIO
      USE mod_diags
#endif
      USE mod_forces
      USE mod_grid
      USE mod_ncparam
      USE mod_ocean
      USE mod_stepping
!
      implicit none
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
!
!  Local variable declarations.
!
#include "tile.h"
!
!  Set header file name.
!
#ifdef PROFILE
      CALL wclock_on (ng, iNLM, 15, __LINE__, __FILE__)
#endif
      CALL biology_tile (ng, tile,                                      &
     &                   LBi, UBi, LBj, UBj, N(ng), NT(ng),             &
     &                   IminS, ImaxS, JminS, JmaxS,                    &
     &                   nstp(ng), nnew(ng),                            &
#ifdef MASKING
     &                   GRID(ng) % rmask,                              &
# if defined WET_DRY && defined DIAGNOSTICS_BIO
     &                   GRID(ng) % rmask_full,                         &
# endif
#endif
     &                   GRID(ng) % Hz,                                 &
#ifdef BULK_FLUXES
     &                   FORCES(ng) % Uwind,                            &
     &                   FORCES(ng) % Vwind,                            &
#else
     &                   FORCES(ng) % sustr,                            &
     &                   FORCES(ng) % svstr,                            &
#endif
#ifdef DIAGNOSTICS_BIO
     &                   DIAGS(ng) % DiaBio2d,                          &
#endif
     &                   OCEAN(ng) % t)

#ifdef PROFILE
      CALL wclock_off (ng, iNLM, 15, __LINE__, __FILE__)
#endif

      RETURN
      END SUBROUTINE biology

      subroutine BFM_ComputeSourceTerms
      use global_mem, only: RLEN
      use envforcing
      IMPLICIT NONE
      ! Assign external data
      call external_data
      ! Assign external event data
      call event_data
#ifdef INCLUDE_SEAICE
      call external_seaice
#endif
      call CalcVerticalExtinction( ) !     Compute extinction coefficient
      call EcologyDynamics      !     Compute reaction terms
      END SUBROUTINE
!-----------------------------------------------------------------------
! Use Magnus formula from https://en.wikipedia.org/wiki/Dew_point
! in order to get the dew point temperature.
      SUBROUTINE SET_BFM_FIELDS_FROM_ROMS(ng, tile)
      USE mod_grid
      USE mod_param
      USE mod_biology
      USE mod_ocean
      USE mod_scalars, only : Cp, itemp, isalt, rho0
      USE mod_stepping
      USE mod_parallel
      USE api_bfm
      USE envforcing, only : botdep_c, botdep_n, botdep_p,              &
     &    botdep_si, daylength, density
      USE mod_forces
#ifdef INCLUDE_PELCO2
      use mem_CO2, only: AtmCO20, AtmCO2, AtmSLP, AtmTDP
#endif
      USE mem, only : Volume, SUNQ,                                     &
     &     jbotR6c, jbotR6n, jbotR6p, jbotR6s,                          &
     &     Area, Area2d, ETW, EIR, ESS, EWIND, ESW,                     &
     &      EICE, NO_BOXES_XY, xEPS, ERHO, Depth
      use mem_PAR,    only: ChlAttenFlag, P_PARRGB, P_PAR,              &
     &     R_EPS, B_EPS, G_EPS,                                         &
     &     EIRR, EIRB, EIRG
      use constants,  only: E2W
      IMPLICIT NONE
      integer, intent(in) :: ng, tile
      integer tileS
      integer idx, idxB, i, j
      REAL(r8) :: ONE = 1.0_r8
      REAL(R8) lat
      REAL(RLEN) b_dew, c_dew
      REAL(RLEN) eTempDewPoint
      REAL(r8) eTemp, eSalt, gamma, RH
      REAL(r8) wlight
      REAL(r8) ux, uy
      REAL(r8) sumWlight, avgWlight
      REAL(r8) sumSUNQ, avgSUNQ
      REAL(r8) Hscale
      REAL(r8) eHz, eVolume
      REAL(r8) area1, area2, TotArea
      REAL(r8) lon1, lon2, lon3, lon4, lat1, lat2, lat3, lat4
      REAL(r8) :: Rearth = 6371000
      REAL(r8) yday
      REAL(r8) z1, z2
      integer k
# include "set_bounds.h"
      Print *, 'SET_BFM_FIELDS_FROM_ROMS, ChlAttenFlag=', ChlAttenFlag
      b_dew=18.678_r8
      c_dew=157.14_r8
      tileS = tile - first_tile(ng) + 1
      sumWlight = 0
      sumSUNQ = 0
      Hscale = rho0 * Cp
      CALL GetTime_for_daylength(ng, yday)
      DO idx=1,NO_BOXES_XY
         i = ListArrayWet(ng) % TheArr(tileS) % ListI(idx)
         j = ListArrayWet(ng) % TheArr(tileS) % ListJ(idx)
         lat = GRID(ng) % latr(i,j)
         ! The sun
         DO k=1,N(ng)
            idxB = N(ng) * (idx - 1) + k
            SUNQ(idx) = daylength(yday, lat)
            sumSUNQ = sumSUNQ + SUNQ(idx)
         END DO
         ! Now computing the area
         lat1 = GRID(ng) % latp(i,j)
         lat2 = GRID(ng) % latp(i+1,j)
         lat3 = GRID(ng) % latp(i+1,j+1)
         lat4 = GRID(ng) % latp(i,j+1)
         lon1 = GRID(ng) % lonp(i,j)
         lon2 = GRID(ng) % lonp(i+1,j)
         lon3 = GRID(ng) % lonp(i+1,j+1)
         lon4 = GRID(ng) % lonp(i,j+1)
         CALL COMPUTE_SPHERICAL_AREA(lon1,lon2,lon3,                    &
     &        lat1, lat2, lat3, area1)
         CALL COMPUTE_SPHERICAL_AREA(lon4,lon2,lon3,                    &
     &        lat4, lat2, lat3, area2)
         TotArea = Rearth * Rearth * (area1 + area2)
         Area2d(idx) = TotArea
         DO k=1,N(ng)
            idxB = N(ng) * (idx - 1) + k
            Area(idxB) = TotArea
         END DO
         ! Now computing the volume
         DO k=1,N(ng)
            eHz = GRID(ng) % Hz(i,j,k)
            idxB = N(ng) * (idx - 1) + k
            eVolume = eHz * TotArea
            Volume(idxB) = eVolume
         END DO
         ! Setting up the depths
         DO k=1,N(ng)
            idxB = N(ng) * (idx-1) + k
            z1 = OCEAN(ng) % zeta(i,j,nrhs(ng))
            z2 = GRID(ng) % z_r(i,j,k)
            Depth(idxB) = z1 - z2
         END DO
         ! The temperature and salinity
         DO k=1,N(ng)
            idxB = N(ng) * (idx - 1) + k
            eTemp = OCEAN(ng) % t(i, j, k, nrhs(ng), itemp)
            eSalt = OCEAN(ng) % t(i, j, k, nrhs(ng), isalt)
            ETW(idxB) = eTemp
            ESW(idxB) = eSalt
         END DO
         ! The irradiance
         ! wlight in the BFM is in W m^{-2}
         ! srflx is in degC m^{-2} so we need to convert to the W m^{-2}
         ! See the bulk_flux code.
         wlight = FORCES(ng) % srflx(i,j) * Hscale
         sumWlight = sumWlight + wlight
         IF (ChlAttenFlag .eq. 1) THEN
           DO k=1,N(ng)
             idxB = N(ng) * (idx - 1) + k
             EIR(idxB) = wlight * p_PAR / E2W
           END DO
         END IF
         IF (ChlAttenFlag .eq. 2) THEN
           DO k=1,N(ng)
              idxB = N(ng) * (idx - 1) + k
              EIR(idxB) = p_PARRGB * wlight / E2W
              EIRR(idxB) = EIR(idxB) * exp ( -R_eps(idxB) )
              EIRG(idxB) = EIR(idxB) * exp ( -G_eps(idxB) )
              EIRB(idxB) = EIR(idxB) * exp ( -B_eps(idxB) )
              EIR(idxB) = EIRB(idxB) + EIRG(idxB) + EIRR(idxB)
              ! weighted broadband diffuse attenuation coefficient for diagnostics
              xEPS(idxB) = (EIRB(idxB)*B_eps(idxB) +                    &
     &             EIRG(idxB)*G_eps(idxB) +                             &
     &             EIRR(idxB)*R_eps(idxB))/EIR(idxB)
           END DO
         END IF
         ! The wind
         ux = FORCES(ng) % Uwind(i,j)
         uy = FORCES(ng) % Uwind(i,j)
         EWIND(idx) = sqrt(ux * ux + uy * uy)
         ! The ICE
         EICE(idx) = 0.0
#ifdef INCLUDE_PELCO2
         ! The CO2 in ppm
         AtmCO2 % fnow(idx) = 402 ! We need better values
         ! The sea surface pressure
         AtmSLP % fnow(idx) = FORCES(ng) % Pair(i,j)
         ! Temperature of dew point
         RH = FORCES(ng) % Hair(i,j) ! The relative humidity. Between 0 and 1.
         gamma = log(RH) + b_dew * eTemp / (c_dew - eTemp)
         eTempDewPoint = c_dew * gamma / (b_dew - gamma)
         AtmTDP % fnow(idx) = eTempDewPoint
#endif
         ! deposition
         if (botdep_c>ZERO) jbotR6c(idx) = 0.0
         if (botdep_n>ZERO) jbotR6n(idx) = 0.0
         if (botdep_p>ZERO) jbotR6p(idx) = 0.0
         if (botdep_si>ZERO) jbotR6s(idx) = 0.0
#ifdef INCLUDE_SEAICE
         ! Sea ice
         EVB(idx) = ONE
         ETB(idx) = ONE
         ESB(idx) = ONE
         ! convert from irradiance to PAR in uE/m2/s
         EIB(idx) = ONE/E2W
         EHB(idx) = ONE
         ESI(idx) = ONE
#endif
         ! Sediment, set it to zero
         DO k=1,N(ng)
            idxB = N(ng) * (idx - 1) + k
            ESS(idxB) = 0
         END DO
      END DO
      ERHO(:) = density(ETW(:),ESW(:),Depth(:)/2.0_RLEN)
      avgWlight = sumWlight / NO_BOXES_XY
      avgSUNQ = sumSUNQ / NO_BOXES_XY
#ifdef BFM_DEBUG
      Print *, 'avgWlight=', avgWlight, 'avgSUNQ=', avgSUNQ
#endif
      END SUBROUTINE
!
!-----------------------------------------------------------------------
      SUBROUTINE READ_BFM_DIAGNOSTICS(LBi, UBi, LBj, UBj, UBk, UBt, ng, tile)
      USE mod_param
      USE mod_grid
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
      USE mod_parallel
      USE mod_ocean
      USE mem
      USE api_bfm
      IMPLICIT NONE
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer, intent(in) :: ng, tile
      integer tileS, i, j, k, idx, iNode, NO_BOXES_XY_loc
      real(RLEN) ARR(NO_BOXES)
      REAL(8) sumABS_ruPTc
      integer :: idx_ruptc = 1
      integer :: idx_ruztc = 2
      REAL(RLEN) val, eavg
      REAL(RLEN) :: sumABS(NO_BOXES_Z)
      REAL(RLEN) :: maxABS(NO_BOXES_Z)
      REAL(RLEN) :: minABS(NO_BOXES_Z)

      tileS = tile - first_tile(ng) + 1
      NO_BOXES_XY_loc = ListArrayWet(ng) % TheArr(tileS) % Nwetpoint
!      Print *, 'NO_BOXES=', NO_BOXES, ' RLEN=', RLEN
!      Print *, 'Before correct_flux_output idx_ruptc=', idx_ruptc
!      Print *, 'NO_BOXES=', NO_BOXES
!      Print *, 'NO_BOXES_XY=', NO_BOXES_XY
!      Print *, 'NO_BOXES_Z=', NO_BOXES_Z


      CALL correct_flux_output(1,idx_ruptc,1,NO_BOXES,ARR)
      sumABS(:) = 0
      maxABS(:) = 0
      minABS(:) = 100000
      DO iNode=1,NO_BOXES_XY_loc
         i = ListArrayWet(ng) % TheArr(tileS) % ListI(iNode)
         j = ListArrayWet(ng) % TheArr(tileS) % ListJ(iNode)
         DO k=1,NO_BOXES_Z
            idx = k + NO_BOXES_Z * (iNode-1)
            OCEAN(ng) % ROMS_ruPTC(i,j,k) = ARR(idx)
            val = ABS(ARR(idx))
            sumABS(k) = sumABS(k) + val
            maxABS(k) = MAX(maxABS(k), val)
            minABS(k) = MIN(minABS(k), val)
         END DO
      END DO
#ifdef BFM_DEBUG
      Print *, 'Statistics for RUPTC'
      DO k=1,NO_BOXES_Z
        eavg = sumABS(k) / NO_BOXES_XY_loc
        Print *, 'k=',k,'avg/max/min=', eavg,maxABS(k),minABS(k)
      END DO
      Print *, 'Before correct_flux_output idx_ruztc=', idx_ruztc
#endif
      CALL correct_flux_output(1,idx_ruztc,1,NO_BOXES,ARR)
      DO iNode=1,NO_BOXES_XY_loc
         i = ListArrayWet(ng) % TheArr(tileS) % ListI(iNode)
         j = ListArrayWet(ng) % TheArr(tileS) % ListJ(iNode)
         DO k=1,NO_BOXES_Z
            idx = k + NO_BOXES_Z * (iNode-1)
            OCEAN(ng) % ROMS_ruZTC(i,j,k) = ARR(idx)
         END DO
      END DO

      DO iNode=1,NO_BOXES_XY_loc
         i = ListArrayWet(ng) % TheArr(tileS) % ListI(iNode)
         j = ListArrayWet(ng) % TheArr(tileS) % ListJ(iNode)
         DO k=1,NO_BOXES_Z
            idx = k + NO_BOXES_Z * (iNode-1)
            OCEAN(ng) % ROMS_EIR(i,j,k) = EIR(idx)
            OCEAN(ng) % ROMS_DIC(i,j,k) = DIC(idx)
            OCEAN(ng) % ROMS_Chlo(i,j,k)  = Chla(idx)
            OCEAN(ng) % ROMS_ixEPS(i,j,k) = xEPS(idx)
         END DO
      END DO
      END SUBROUTINE
!-----------------------------------------------------------------------
      SUBROUTINE PRINT_T_KEYS(LBi, UBi, LBj, UBj, UBk, UBt, ng, tile, t)
      USE mod_param
      USE mod_grid
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
      USE mod_parallel
      USE mod_ocean
      USE mem
      USE api_bfm
      IMPLICIT NONE
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer, intent(in) :: ng, tile
#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: t(LBi:,LBj:,:,:,:)
#else
      real(r8), intent(in) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
#endif
      REAL(r8), allocatable :: F(:)
      REAL(r8) emaxval, eminval, eavgval, eVal
      integer iVar, itrc, ibio, pos, iNode, i, j, k, idx
      integer NO_BOXES_XY_loc, tileS, posmax
      Print *, 'PRINT_T_KEYS ng=', ng, ' tile=', tile
      posmax = 1
      tileS = tile - first_tile(ng) + 1
      NO_BOXES_XY_loc = ListArrayWet(ng) % TheArr(tileS) % Nwetpoint
      Npt = NO_BOXES_XY_loc * NO_BOXES_Z
      allocate(F(Npt))
      DO iVar=stPelStateS, stPelStateE
         itrc = iVar - stPelStateS + 1
         ibio = idbio(itrc)
         DO pos=1,posmax
            idx = 0
            DO iNode=1,NO_BOXES_XY_loc
               i = ListArrayWet(ng) % TheArr(tileS) % ListI(iNode)
               j = ListArrayWet(ng) % TheArr(tileS) % ListJ(iNode)
               DO k=1,NO_BOXES_Z
                  eVal = t(i, j, k, pos, ibio)
                  idx = idx + 1
                  F(idx) = eVal
               END DO
            END DO
            emaxval = maxval(F)
            eminval = minval(F)
            eavgval = sum(F) / Npt
            Print *, 'iVar/pos=', iVar, pos, 'avg/min/max=',            &
     &      eavgval, eminval, emaxval
         END DO
      END DO
      deallocate(F)
      END SUBROUTINE
!-----------------------------------------------------------------------
      SUBROUTINE COPY_T_to_D3STATE(LBi, UBi, LBj, UBj, UBk, UBt, ng, tile, eTimeIdx, t)
      USE mod_param
      USE mod_grid
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
      USE mod_parallel
      USE mod_ocean
      USE mem
      USE api_bfm
      IMPLICIT NONE
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer, intent(in) :: ng, tile, eTimeIdx
#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#else
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
#endif
      integer iNode, i, j, k, idx
      integer iVar, itrc, ibio, NO_BOXES_XY_loc, tileS
      REAL(r8) eVal
!      Print *, 'CP_T_D3 : stPelStateS=', stPelStateS, ' stPelStateE=', stPelStateE
!      Print *, 'first_tile=', first_tile(ng)
!      Print *, 'stPelStateS=', stPelStateS, ' stPelStateE=', stPelStateE
      tileS = tile - first_tile(ng) + 1
      NO_BOXES_XY_loc = ListArrayWet(ng) % TheArr(tileS) % Nwetpoint
      DO iNode=1,NO_BOXES_XY_loc
         i = ListArrayWet(ng) % TheArr(tileS) % ListI(iNode)
         j = ListArrayWet(ng) % TheArr(tileS) % ListJ(iNode)
         DO k=1,NO_BOXES_Z
            idx = k + NO_BOXES_Z * (iNode-1)
!            Print *, '1: iZ=', iZ, ' iNode=', iNode, ' k=', k, ' idx=', idx
            DO iVar=stPelStateS, stPelStateE
               itrc = iVar - stPelStateS + 1
               ibio = idbio(itrc)
!               Print *, 'iVar=', iVar, ' itrc=', itrc, ' ibio=', ibio
!               Print *, 'i=', i, ' j=', j, ' k=', k, ' nstp=', nstp
               eVal = t(i, j, k, eTimeIdx, ibio)
!               Print *, ' eVal=', eVal
               D3STATE(itrc, idx) = eVal
            END DO
         END DO
      END DO
      END SUBROUTINE

      SUBROUTINE INCREASE_T_NNEW(LBi, UBi, LBj, UBj, UBk, UBt, ng, tile, nnew, t)
      USE mod_param
      USE mod_grid
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
      USE mod_parallel
      USE mem
      USE api_bfm
      IMPLICIT NONE
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer, intent(in) :: ng, tile, nnew
#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#else
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
#endif
      integer tileS, NO_BOXES_XY_loc, iNode, i, j, k
      integer idx, itrc, ibio, iVar
      real(r8) cff
      tileS = tile - first_tile(ng) + 1
      NO_BOXES_XY_loc = ListArrayWet(ng) % TheArr(tileS) % Nwetpoint
      DO iNode=1,NO_BOXES_XY_loc
         i = ListArrayWet(ng) % TheArr(tileS) % ListI(iNode)
         j = ListArrayWet(ng) % TheArr(tileS) % ListJ(iNode)
         DO k=1,NO_BOXES_Z
            idx = k + NO_BOXES_Z * (iNode-1)
            DO iVar=stPelStateS, stPelStateE
               itrc = iVar - stPelStateS + 1
               ibio = idbio(itrc)
               IF (D3STATETYPE(itrc).ge.0) THEN
#ifndef EXPLICIT_SINK
                  cff = D3SOURCE(itrc,idx)
#else
                  cff = 0
                  DO u=1,NO_D3_BOX_STATES
                     cff = cff + D3SOURCE(itrc,u,idx)                   &
     &                         - D3SINK(itrc,u,idx)
                  END DO
#endif
                  t(i,j,k,nnew,ibio) = t(i,j,k,nnew,ibio) +             &
     &               cff * delt_bfm * GRID(ng) % Hz(i,j,k)
               END IF
            END DO
         END DO
      END DO
      END SUBROUTINE
!
      SUBROUTINE BFM_STATE_UPDATE(ng, tile)
      USE mod_parallel
      USE mod_biology
      USE mem
      USE api_bfm
      IMPLICIT NONE
      integer, intent(in) :: ng, tile
      integer tileS, i, j, k
      integer NO_BOXES_XY_loc
      tileS = tile - first_tile(ng) + 1
      NO_BOXES_XY_loc = ListArrayWet(ng) % TheArr(tileS) % Nwetpoint
      DO j=1,NO_D3_BOX_STATES
         IF (D3STATETYPE(j).ge.0) THEN
            DO i=1,NO_BOXES
#ifndef EXPLICIT_SINK
               D3STATE(j,i) = D3STATE(j,i) + delt_bfm * D3SOURCE(j,i)
#else
               DO k=1,NO_D3_BOX_STATES
                  D3STATE(j,i) = D3STATE(j,i) +                         &
     &                 delt_bfm * (D3SOURCE(j,k,i) - D3SINK(j,k,i))
               END DO
#endif
            END DO
         END IF
      END DO
      END SUBROUTINE
!
      SUBROUTINE COPY_D3STATE_to_T(LBi, UBi, LBj, UBj, UBk, UBt, ng, tile, eTimeIdx, t)
      USE mod_param
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
      USE mod_parallel
      USE mem
      USE api_bfm
      IMPLICIT NONE
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer, intent(in) :: ng, tile, eTimeIdx
#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#else
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
#endif
      integer iNode, i, j, k, idx
      integer iVar, itrc, ibio, NO_BOXES_XY_loc, tileS
      REAL(r8) eVal
!      Print *, 'stPelStateS=', stPelStateS, ' stPelStateE=', stPelStateE
      tileS = tile - first_tile(ng) + 1
      NO_BOXES_XY_loc = ListArrayWet(ng) % TheArr(tileS) % Nwetpoint
      DO iNode=1,NO_BOXES_XY_loc
         i = ListArrayWet(ng) % TheArr(tileS) % ListI(iNode)
         j = ListArrayWet(ng) % TheArr(tileS) % ListJ(iNode)
         DO k=1,NO_BOXES_Z
            idx = k + NO_BOXES_Z * (iNode-1)
!            Print *, '2: iZ=', iZ, ' iNode=', iNode, ' k=', k, ' idx=', idx
            DO iVar=stPelStateS, stPelStateE
               itrc = iVar - stPelStateS + 1
               ibio = idbio(itrc)
!               Print *, 'iVar=', iVar, ' itrc=', itrc, ' ibio=', ibio
               t(i, j, k, eTimeIdx, ibio) = D3STATE(itrc, idx)
            END DO
         END DO
      END DO
      END SUBROUTINE

      SUBROUTINE PRINT_BFM_STATE_KEYS(ng, tile)
      USE mod_param
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
      USE mem
      USE api_bfm
      IMPLICIT NONE
      integer, intent(in) :: ng, tile
      integer iNode, i, j, k, idx
      integer iVar, itrc, ibio, siz, TotalNb
      REAL(r8), allocatable :: ArrSum(:), ArrMin(:), ArrMax(:)
      REAL(r8) eVal, eAvg, eMin, eMax
      REAL(r8) eminval, emaxval, eavgval
      REAL(r8), allocatable :: F(:)
!      Print *, 'PRINT_AVERAGE: stPelStateS=', stPelStateS, ' stPelStateE=', stPelStateE
!      Print *, 'stPelFluxS=', stPelFluxS, ' stPelFluxE=', stPelFluxE
      DO iVar=stPelStateS, stPelStateE
        itrc = iVar - stPelStateS + 1
        ibio = idbio(itrc)
!        Print *, 'iVar=', iVar, ' itrc=', itrc, ' ibio=', ibio
      END DO
!      Print *, 'size(D3STATE,1:2)=', size(D3STATE,1), size(D3STATE,2)
!      Print *, 'Computing the average NO_BOXES_XY/Z=', NO_BOXES_XY, NO_BOXES_Z
      siz = stPelStateE + 1 - stPelStateS
      Allocate(ArrSum(siz), ArrMin(siz), ArrMax(siz))
      DO i=1,siz
        ArrSum(i) = 0
        ArrMin(i) =  100000
        ArrMax(i) = -100000
      END DO
      DO iNode=1,NO_BOXES_XY
         DO k=1,NO_BOXES_Z
            idx = k + NO_BOXES_Z * (iNode-1)
            DO iVar=stPelStateS, stPelStateE
               itrc = iVar - stPelStateS + 1
               eVal = D3STATE(itrc, idx)
               ArrSum(itrc) = ArrSum(itrc) + eVal
               IF (eVal .lt. ArrMin(itrc)) THEN
                 ArrMin(itrc) = eVal
               END IF
               IF (eVal .gt. ArrMax(itrc)) THEN
                 ArrMax(itrc) = eVal
               END IF
            END DO
         END DO
      END DO
      TotalNb = NO_BOXES_Z * NO_BOXES_XY
      Print *, 'PRINT_BFM_STATE_KEYS : The D3STATE variables'
      DO i=1,siz
        eAvg = ArrSum(i) / TotalNb
        eMin = ArrMin(i)
        eMax = ArrMax(i)
        Print *, 'i=', i, ' min/max/avg=', eMin, eMax, eAvg
      END DO
      deallocate(ArrSum, ArrMin, ArrMax)
      Print *, 'PRINT_BFM_STATE_KEYS : The D3DIAGNOS variables'
      allocate(F(TotalNb))
      DO j=1,NO_D3_BOX_DIAGNOSS
         DO i=1,TotalNb
            F(i) = D3DIAGNOS(j,i)
         END DO
         eminval = minval(F)
         emaxval = maxval(F)
         eavgval = sum(F) / TotalNb
         Print *, 'j=', j, ' min/max/avg=', eminval, emaxval, eavgval
      END DO

      DO j=1,NO_D3_BOX_STATES
         DO i=1,TotalNb
            F(i) = D3STATE(j,i)
         END DO
         eminval = minval(F)
         emaxval = maxval(F)
         eavgval = sum(F) / TotalNb
!         Print *, 'j=', j, ' min/max/avg=', eminval, emaxval, eavgval
      END DO
      deallocate(F)
      END SUBROUTINE

      SUBROUTINE COMPUTE_LOCAL_NB_WET(ng, tile, Nwetpoint)
      USE mod_grid
      USE mod_param
      USE mod_biology
      USE api_bfm
      USE mem
      implicit none
      integer, intent(in) :: ng, tile
      integer, intent(out) :: Nwetpoint
      integer i, j
# include "set_bounds.h"
      Nwetpoint=0
      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
#ifdef MASKING
          IF (GRID(ng) % rmask(i,j) .eq. 1) THEN
#endif
             Nwetpoint=Nwetpoint + 1
#ifdef MASKING
          END IF
#endif
        END DO
      END DO
      END SUBROUTINE
!
!
!
      SUBROUTINE COMPUTE_NO_BOXES_XY_Z_ALL
      use mod_grid
      USE mod_biology
      use api_bfm
      USE mod_parallel
      USE mod_scalars
      USE mem
      implicit none
      REAL(8) delta
      integer tile, tileS
      integer nl, ig, ng, eProd
      integer, parameter :: file_id = 1453
      integer istat, i
      LOGICAL DoNestLayer
      integer Nwetpoint, tileLen
      integer NO_BOXES_Z_max, NO_BOXES_XY_max
      NAMELIST /SETTING_BFM_COUPL/ delt_bfm, AnalyticalInitD3STATE,     &
     &   CopyInitialToD3STATE, CopyD3STATEtoInitial,                    &
     &   SourceTermD3STATE, AdvectionD3STATE
      eFile = TRIM(BFM_Prefix_NML) // 'bfm_input.nml'
      OPEN(file_id, FILE = eFile)
      READ(file_id, NML = SETTING_BFM_COUPL)
      CLOSE(file_id)
      DO ng=1,Ngrids
        MULTIPLIER(ng) = INT(delt_bfm / dt(ng))
        delta = ABS(MULTIPLIER(ng) * dt(ng) - delt_bfm)
        IF (delta .gt. 1.0) THEN
          Print *, "delta=", delta, " which should be 0"
          STOP
          END IF
      END DO
      DoNestLayer = .TRUE.
      !
      nl=0
      allocate(NO_BOXES_Z_arr(Ngrids), stat=istat)
      allocate(NO_BOXES_XY_arr(Ngrids), stat=istat)
      allocate(NO_BOXES_arr(Ngrids), stat=istat)
      allocate(ListArrayWet(Ngrids), stat=istat)
      NO_BOXES_Z_max  = 0
      NO_BOXES_XY_max = 0
      NEST_LAYER : DO WHILE (DoNestLayer)
        nl = nl + 1
        IF ((nl.le.0).or.(nl.gt.NestLayers)) EXIT
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          tileLen = 1 + last_tile(ng) - first_tile(ng)
          Print *, 'nl=', nl, ' ig=', ig, ' ng=', ng, ' tileLen=', tileLen
          allocate(NO_BOXES_Z_arr(ng) % ArrInt(tileLen), stat=istat)
          allocate(NO_BOXES_XY_arr(ng) % ArrInt(tileLen), stat=istat)
          allocate(NO_BOXES_arr(ng) % ArrInt(tileLen), stat=istat)
          allocate(ListArrayWet(ng) % TheArr(tileLen), stat=istat)
          DO tile=first_tile(ng),last_tile(ng),+1
            tileS = tile - first_tile(ng) + 1
            CALL COMPUTE_LOCAL_NB_WET(ng, tile, Nwetpoint)
            eProd = Nwetpoint * N(ng)
            Print *, '  tile=', tile, ' Nwetpoint=', Nwetpoint
            !
            NO_BOXES_XY_arr(ng) % ArrInt(tileS) = Nwetpoint
            NO_BOXES_Z_arr(ng) % ArrInt(tileS) = N(ng)
            IF (Nwetpoint .gt. NO_BOXES_XY_max) NO_BOXES_XY_max = Nwetpoint
            IF (N(ng) .gt. NO_BOXES_Z_max) NO_BOXES_Z_max = N(ng)
            NO_BOXES_arr(ng) % ArrInt(tileS) = eProd
            ListArrayWet(ng) % TheArr(tileS) % Nwetpoint = Nwetpoint
          END DO
        END DO
      END DO NEST_LAYER
      NO_BOXES_Y  = 1
      NO_BOXES_Z  = NO_BOXES_Z_max
      NO_BOXES_XY = NO_BOXES_XY_max
      NO_BOXES_X  = NO_BOXES_XY
      NO_BOXES    = NO_BOXES_XY * NO_BOXES_Z
      NO_STATES   = NO_D3_BOX_STATES * NO_BOXES + NO_BOXES_XY
      !
      allocate(BOTindices(Nwetpoint), stat=istat)
      allocate(SRFindices(Nwetpoint), stat=istat)
      END SUBROUTINE

      SUBROUTINE SET_BOT_SURFINDICES(ng, tile)
      use mod_grid
      use api_bfm
      use mod_parallel
      implicit none
      integer, intent(in) :: ng, tile
      integer eN, iNode, idxBOT, idxSRF
      integer tileS, NO_BOXES_XY_loc
      tileS = tile - first_tile(ng) + 1
      eN = NO_BOXES_Z_arr(ng) % ArrInt(tileS)
      NO_BOXES_XY_loc = NO_BOXES_XY_arr(ng) % ArrInt(tileS)
      DO iNode=1,NO_BOXES_XY_loc
        idxBOT = 1  + eN * (iNode-1)
        idxSRF = eN + eN * (iNode-1)
        BOTindices(iNode) = idxBOT
        SRFindices(iNode) = idxSRF
      END DO
      END SUBROUTINE
!
!
!
      SUBROUTINE FULL_INIT_OF_BFM_MODEL
#ifdef INCLUDE_PELCO2
      use mem_CO2, only: AtmCO20, AtmCO2, AtmSLP, AtmTDP
      use SystemForcing, only : FieldAllocate
#endif
      implicit none
      CALL COMPUTE_NO_BOXES_XY_Z_ALL
#ifdef INCLUDE_PELCO2
      CALL FieldAllocate(AtmCO2)
      CALL FieldAllocate(AtmSLP)
      CALL FieldAllocate(AtmTDP)
#endif
      CALL Allocate_GRID_ARRAY
      CALL INIT_BFM_MODEL
      CALL INIT_BFM_SYSTEM_VARIABLE
      END SUBROUTINE



      SUBROUTINE Allocate_GRID_ARRAY
      use mod_grid
      use mod_param
      USE mod_parallel
      implicit none
      integer nl, ig, ng, tile
      LOGICAL DoNestLayer
      DoNestLayer = .TRUE.
      nl=0
      NEST_LAYER : DO WHILE (DoNestLayer)
        nl = nl + 1
        IF ((nl.le.0).or.(nl.gt.NestLayers)) EXIT
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL Allocate_GRID_ARRAY_LOCAL(ng, tile)
          END DO
        END DO
      END DO NEST_LAYER
      END SUBROUTINE


      SUBROUTINE Allocate_GRID_ARRAY_LOCAL(ng, tile)
      use mod_biology
      use mod_grid
      use mod_param
      USE mod_parallel
      implicit none
      integer, intent(in) :: ng, tile
      integer tileS
      integer idx, j, i, Nwetpoint, istat
#include "set_bounds.h"
      tileS = tile - first_tile(ng) + 1
      Print *, 'ng=', ng, ' tile=', tile, ' tileS=', tileS
      Print *, 'allocated(ListArrayWet)=', allocated(ListArrayWet)
      Print *, 'allocated(ListArrayWet(ng) % TheArr)=', allocated(ListArrayWet(ng) % TheArr)
      Nwetpoint = ListArrayWet(ng) % TheArr(tileS) % Nwetpoint
      allocate(ListArrayWet(ng) % TheArr(tileS) % ListI(Nwetpoint), stat=istat)
      allocate(ListArrayWet(ng) % TheArr(tileS) % ListJ(Nwetpoint), stat=istat)
! The WET arrays.
      idx=0
      DO j=Jstr-1,JendR
        DO i=Istr-1,IendR
#ifdef MASKING
          IF (GRID(ng) % rmask(i,j) .eq. 1) THEN
#endif
            idx = idx + 1
            ListArrayWet(ng) % TheArr(tileS) % ListI(idx) = i
            ListArrayWet(ng) % TheArr(tileS) % ListJ(idx) = j
#ifdef MASKING
          END IF
#endif
        END DO
      END DO
      END SUBROUTINE

      SUBROUTINE INIT_BFM_MODEL
      USE mod_grid
      USE mod_param
      USE mod_biology
      USE api_bfm
      USE mem
      USE init_var_bfm_local, only : init_organic_constituents
      IMPLICIT NONE
      integer bio_setup_loc
      integer, parameter :: namlst = 10
      ! Initialise the BFM with standalone settings
      call init_bfm(namlst)
      ! Initialise state variable names and diagnostics
      call set_var_info_bfm
      ! Allocate memory and give initial values
      ! to the pelagic system
      ! We need the variable well set for this to work.
      bio_setup_loc=1
      call init_var_bfm(bio_setup_loc, AnalyticalInitD3STATE)
      ! Initialize internal constitutents of functional groups
      call init_organic_constituents(AnalyticalInitD3STATE)
      ! Need to set up bfmtime adequately for the runs.
      ! Need also to set up the wind.
      call init_envforcing_bfm
      ! Read restart file (if flag)
      ! Overwrite previous initialization
      ! Initialise the diagnostic variables
      call CalcVerticalExtinction( )
      call CalcChlorophylla( )
      END SUBROUTINE



      SUBROUTINE INIT_BFM_SYSTEM_VARIABLE
      use mod_grid
      use mod_param
      USE mod_parallel
      implicit none
      integer nl, ig, ng, tile
      LOGICAL DoNestLayer
      DoNestLayer = .TRUE.
      !
      nl=0
      NEST_LAYER : DO WHILE (DoNestLayer)
        nl = nl + 1
        IF ((nl.le.0).or.(nl.gt.NestLayers)) EXIT
        DO ig=1,GridsInLayer(nl)
          ng=GridNumber(ig,nl)
          DO tile=first_tile(ng),last_tile(ng),+1
            CALL INIT_BFM_SYSTEM_VARIABLE_LOCAL(ng, tile)
          END DO
        END DO
      END DO NEST_LAYER
      END SUBROUTINE


      SUBROUTINE INIT_BFM_SYSTEM_VARIABLE_LOCAL(ng, tile)
      use mod_param
      use mod_ocean
      implicit none
      integer, intent(in) :: ng, tile
# include "tile.h"
      CALL INIT_BFM_SYSTEM_VARIABLE_LOCAL_tile(LBi, UBi, LBj, UBj, N(ng), NT(ng), ng, tile, OCEAN(ng) % t)
      CALL SET_BFM_FIELDS_FROM_ROMS(ng, tile)
      END SUBROUTINE

      SUBROUTINE INIT_BFM_SYSTEM_VARIABLE_LOCAL_tile(LBi, UBi, LBj, UBj, UBk, UBt, ng, tile, t)
      USE mod_grid
      USE mod_param
      USE mod_biology
      USE api_bfm
      USE mem
      USE init_var_bfm_local, only : init_organic_constituents
      IMPLICIT NONE
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer, intent(in) :: ng, tile
#ifdef ASSUMED_SHAPE
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#else
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
#endif
      integer idx, i, j
      integer idxBOT, idxSRF, eN, iNode
      integer bio_setup_loc
# include "set_bounds.h"
      Print *, 'LBi=', LBi, ' UBi=', UBi
      Print *, 'LBj=', LBj, ' UBj=', UBj
      Print *, 'UBk=', UBk, ' UBt=', UBt
      Print *, 'ng=', ng, ' tile=', tile
      !
      Print *, 'Printing average of T'
      CALL PRINT_T_KEYS(LBi, UBi, LBj, UBj, UBk, UBt, ng, tile, t)
      Print *, 'Printing average of D3STATE'
      Print *, 'CopyInitialToD3STATE=', CopyInitialToD3STATE
      Print *, 'CopyD3STATEtoInitial=', CopyD3STATEtoInitial
      IF (CopyD3STATEtoInitial) THEN
        Print *, 'Copying D3STATE to T'

!        Print *, 'Printing T avg in INIT_BFM_SYSTEM_VARIABLE, Before'
!        CALL Print_t_average(ng, tile, 1)
!        CALL Print_t_average(ng, tile, 2)

        CALL COPY_D3STATE_to_T(LBi, UBi, LBj, UBj, UBk, UBt, ng, tile, 1, t)
        CALL COPY_D3STATE_to_T(LBi, UBi, LBj, UBj, UBk, UBt, ng, tile, 2, t)
!        Print *, 'Printing T avg in INIT_BFM_SYSTEM_VARIABLE, After'
!        CALL Print_t_average(ng, tile, 1)
!        CALL Print_t_average(ng, tile, 2)
      END IF
      IF (CopyInitialToD3STATE) THEN
        Print *, 'Copying T to D3STATE'
        CALL COPY_T_to_D3STATE(LBi, UBi, LBj, UBj, UBk, UBt, ng, tile, 1, t)
      END IF
      CALL PRINT_BFM_STATE_KEYS(ng, tile)
      Print *, 'End of INIT_BFM_SYSTEM_VARIABLE_local'
      END SUBROUTINE

      SUBROUTINE PRINT_CRITICAL(pos)
      USE mod_biology
      use mod_ocean
      use mod_mixing
      implicit none
      integer, intent(in) :: pos
      integer ibio, ng
      real(r8) val1, val2, val3, maxakt
      ng = 1
      ibio = idbio(1)
      val1 = maxval(OCEAN(ng) % t(:,:,:, 1, ibio))
      val2 = maxval(OCEAN(ng) % t(:,:,:, 2, ibio))
      val3 = maxval(OCEAN(ng) % t(:,:,:, 3, ibio))
      maxakt = maxval(MIXING(ng) % Akt)
!      Print *, 'Critical pos=', pos, ' val123=', val1, val2, val3
!      Print *, 'Critical pos=', pos, ' val12=', val1, val2, ' m(akt)=', maxakt
      END SUBROUTINE

!
!-----------------------------------------------------------------------
      SUBROUTINE biology_tile (ng, tile,                                &
     &                         LBi, UBi, LBj, UBj, UBk, UBt,            &
     &                         IminS, ImaxS, JminS, JmaxS,              &
     &                         nstp, nnew,                              &
#ifdef MASKING
     &                         rmask,                                   &
# if defined WET_DRY && defined DIAGNOSTICS_BIO
     &                         rmask_full,                              &
# endif
#endif
     &                         Hz,                                      &
#ifdef BULK_FLUXES
     &                         Uwind, Vwind,                            &
#else
     &                         sustr, svstr,                            &
#endif
#ifdef DIAGNOSTICS_BIO
     &                         DiaBio2d,                                &
#endif
     &                         t)
!-----------------------------------------------------------------------
!
      USE mod_param
      USE mod_ocean
      USE mod_grid
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
      USE mod_parallel
      USE mem
      USE api_bfm
!
!  Imported variable declarations.
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj, UBk, UBt
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: nstp, nnew

#ifdef ASSUMED_SHAPE
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:,LBj:)
#  if defined WET_DRY && defined DIAGNOSTICS_BIO
      real(r8), intent(in) :: rmask_full(LBi:,LBj:)
#  endif
# endif
      real(r8), intent(in) :: Hz(LBi:,LBj:,:)
# ifdef BULK_FLUXES
      real(r8), intent(in) :: Uwind(LBi:,LBj:)
      real(r8), intent(in) :: Vwind(LBi:,LBj:)
# else
      real(r8), intent(in) :: sustr(LBi:,LBj:)
      real(r8), intent(in) :: svstr(LBi:,LBj:)
# endif
# ifdef DIAGNOSTICS_BIO
      real(r8), intent(inout) :: DiaBio2d(LBi:,LBj:,:)
# endif
      real(r8), intent(inout) :: t(LBi:,LBj:,:,:,:)
#else
# ifdef MASKING
      real(r8), intent(in) :: rmask(LBi:UBi,LBj:UBj)
#  if defined WET_DRY && defined DIAGNOSTICS_BIO
      real(r8), intent(in) :: rmask_full(LBi:UBi,LBj:UBj)
#  endif
# endif
      real(r8), intent(in) :: Hz(LBi:UBi,LBj:UBj,UBk)
# ifdef BULK_FLUXES
      real(r8), intent(in) :: Uwind(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: Vwind(LBi:UBi,LBj:UBj)
# else
      real(r8), intent(in) :: sustr(LBi:UBi,LBj:UBj)
      real(r8), intent(in) :: svstr(LBi:UBi,LBj:UBj)
# endif
      real(r8), intent(inout) :: repiration(LBi:UBi,LBj:UBj,UBk)
# ifdef DIAGNOSTICS_BIO
      real(r8), intent(inout) :: DiaBio2d(LBi:UBi,LBj:UBj,NDbio2d)
# endif
      real(r8), intent(inout) :: t(LBi:UBi,LBj:UBj,UBk,3,UBt)
#endif
      real(r8) eVal
!
! We have two entries, nstp and nnew for the timing.
! How to reconcile both.
! If we use only nstp, then we are overwritten after by the nnew.
! on the other hand, said increase must be correct.
!
! We face two different extremes
! ---- Update only nstp. But then we face te problem of nnew overwritting
!   the nstp and so not significant growth occurring.
! ---- Put the updated value in the nnew. But then it explodes.
! ---- Do the increase on both nstp and nnew.
!
! The field values on nnew and nstp are not the same.
!
! What are other doing? The model like fennel.h is doing something like
! ---First Bio = Bio_old
! ---Then doing some work on Bio
! ---Then put the contribution to the nnew with something like
!        cff=Bio(i,k,ibio)-Bio_old(i,k,ibio)
!        t(i,j,k,nnew,ibio)=t(i,j,k,nnew,ibio)+cff*Hz(i,j,k)
!
! ecosim works on Bio / Bio_old / Bio_new
! nemuro is like fennel (Bio / Bio_old)
! npzd_Franks is also the same
! npzd_Powell is also the same
! npzd_iron is also the same
! oyster_float is float base, not tracer based
! red_tide is also like that.
!



!
!  Assigning the STATE variables from the t array
!  ! We need to determine if the diagnostics need to be recomputed.
!

!      Print *, 'Printing T average in biology_tile'
!      CALL Print_t_average(ng, tile, nstp)
!      CALL Print_t_average(ng, tile, nnew)
#ifdef BFM_DEBUG
      Print *, "PosMultiplier=", PosMultiplier, " MULTIPLIER(ng)=", MULTIPLIER(ng)
      Print *, "SourceTermD3STATE=", SourceTermD3STATE
      Print *, "AdvectionD3STATE=", AdvectionD3STATE
#endif
      PosMultiplier = PosMultiplier + 1
      IF (PosMultiplier == MULTIPLIER(ng)) THEN
        PosMultiplier = 0
        CALL SET_BOT_SURFINDICES(ng, tile)
        call SET_BFM_FIELDS_FROM_ROMS(ng, tile)

#ifdef BFM_DEBUG
        CALL PRINT_BFM_STATE_KEYS(ng, tile)
#endif
        IF (AdvectionD3STATE) THEN
           CALL COPY_T_to_D3STATE(LBi, UBi, LBj, UBj, UBk, UBt, ng, tile, nstp, t)
        ENd IF
        IF (SourceTermD3STATE) THEN
!
!          Now the time stepping operations
!          Right now we do Euler forward algorithm
!          It is probably needed to subdivide further the time interval
!          or to use a better integration method
!
           CALL BFM_ComputeSourceTerms
#ifdef BFM_DEBUG
# ifndef EXPLICIT_SINK
           Print *, "        Case ifndef EXPLICIT_SINK"
# else
           Print *, "        Case ifdef EXPLICIT_SINK"
# endif
#endif
           CALL BFM_STATE_UPDATE(ng, tile)
        END IF
        IF (AdvectionD3STATE) THEN
!
!          Doing the time step increase. Here, we follow directly the
!          schemes of redtide, fennel, npzd_*, nemuro.
!          We do not need to copy to nstp in that case.
!
           CALL INCREASE_T_NNEW(LBi, UBi, LBj, UBj, UBk, UBt,           &
     &         ng, tile, nnew, t)
        ELSE
!
!          Now copying back the field values
!          This is needed even if in absence of advection because we use this for
!          for outputting results.
!
           CALL COPY_D3STATE_to_T(LBi, UBi, LBj, UBj, UBk, UBt,         &
     &         ng, tile, nstp, t)
        END IF
        CALL READ_BFM_DIAGNOSTICS(LBi, UBi, LBj, UBj, UBk, UBt,         &
     &         ng, tile)
        IF (SourceTermD3STATE) THEN
!
!  Now set up the D3SOURCE and D3SINK to zero.
!
          CALL ResetFluxes
        END IF
!       Need to put code for the diagnostics. We do not put yet the dlux. Maybe never.
      END IF
      END SUBROUTINE biology_tile
