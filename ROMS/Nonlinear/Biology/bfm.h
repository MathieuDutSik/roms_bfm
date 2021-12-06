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
!  This simple formulation is based on the work of Malcolm Scully. The !
!  dissolved oxygen (DO) is respired at a constant rate.  A constant   !
!  magnitude of DO is respired during each time step and the DO is not !
!  allowed to go negative.                                             !
!                                                                      !
!  Dissolved Oxygen can either be added to the surface as a flux or    !
!  the surface layer DO can be set to saturation based on the layer    !
!  temperature and salinity.                                           !
!                                                                      !
!  The surface oxygen saturation formulation is the same as in the     !
!  Fennel ecosystem model.  The surface oxygen flux is different from  !
!  the method used by the below citations.  Setting the surface DO to  !
!  saturation is very similar to the ChesROMS-CRM method in Irby et.   !
!  al. (2016).                                                         !
!                                                                      !
!    Scully, M.E., 2010: Wind Modulation of Dissolved Oxygen in        !
!      Chesapeake Bay, Estuaries and Coasts, 33, 1164-1175.            !
!                                                                      !
!    Scully, M.E., 2013: Physical control on hypoxia in Chesapeake     !
!      Bay: A numerical modeling study, J. Geophys. Res., 118,         !
!      1239-1256.                                                      !
!                                                                      !
!    Irby, I. et al., 2016: Challenges associated with modeling low-   !
!      oxygen waters in Chesapeake Bay: a multiple model comparison,   !
!      Biogeosciences, 13, 2011-2028                                   !
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
!
      subroutine envforcing_bfm(ng, tile, step)
      use global_mem, only: RLEN
      use envforcing
      IMPLICIT NONE
      integer, intent(in)  :: ng, tile
      integer, intent(inout)  :: step
      call SET_BFM_FIELDS_FROM_ROMS(ng, tile)
      ! Assign external data
      call external_data
      ! Assign external event data
      call event_data
#ifdef INCLUDE_SEAICE
      call external_seaice
#endif
      if (init_forcing_vars) init_forcing_vars=.false.
      end subroutine envforcing_bfm
!-----------------------------------------------------------------------
! Use Magnus formula from https://en.wikipedia.org/wiki/Dew_point
! in order to get the dew point temperature.
      SUBROUTINE SET_BFM_FIELDS_FROM_ROMS(ng, tile)
      USE mod_grid
      USE mod_param
      USE mod_biology
      USE mod_ocean
      USE mod_scalars
      USE mod_stepping
      USE mod_parallel
      USE api_bfm
      USE envforcing, only : botdep_c, botdep_n, botdep_p,              &
     &    botdep_si, daylength
      USE mod_forces
#ifdef INCLUDE_PELCO2
      use mem_CO2, only: AtmCO20, AtmCO2, AtmSLP, AtmTDP
#endif
      USE mem
      use mem_PAR,    only: ChlAttenFlag, P_PARRGB, P_PAR,              &
     &                      R_EPS, B_EPS, G_EPS,                        &
     &                      EIRR, EIRB, EIRG
      use constants,  only: E2W
      IMPLICIT NONE
      integer, intent(in) :: ng, tile
      integer tileS
      integer idx, idxB, i, j
      REAL(r8) :: ONE = 1.0_r8
      REAL(R8) lat
      REAL(RLEN) b_dew, c_dew
      REAL(RLEN) eTempDewPoint
      REAL(r8) eTemp, eSalt
      REAL(r8) wlight
      REAL(r8) ux, uy
      integer k
# include "set_bounds.h"
      b_dew=18.678_r8
      c_dew=157.14_r8
      tileS = tile - first_tile(ng) + 1
      DO idx=1,NO_BOXES_XY
         i = ListArrayWet(ng) % TheArr(tileS) % ListI(idx)
         j = ListArrayWet(ng) % TheArr(tileS) % ListJ(idx)
         lat = GRID(ng) % latp(i,j)
         ! The sun
         SUNQ(idx) = daylength(REAL(Rclock % yday, r8), lat)
         ! The temperature and salinity
         eTemp = OCEAN(ng) % t(i, j, N(ng), nrhs(ng), itemp)
         eSalt = OCEAN(ng) % t(i, j, N(ng), nrhs(ng), isalt)
         ETW(idx) = eTemp
         ESW(idx) = eSalt
         ! The irradiance
         wlight = FORCES(ng) % srflx(i,j)
         IF (ChlAttenFlag .eq. 1) THEN
           DO k=1,N(ng)
             idxB = N(ng) * (idx - 1) + k
             EIR(idxB) = wlight * p_PAR / E2W
           END DO
         END IF
         IF (ChlAttenFlag .eq. 2) THEN
           DO k=1,N(ng)
              idxB = N(ng) * (idx - 1) + k
              EIR(idxB) = p_PARRGB* wlight / E2W
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
      END DO
      END SUBROUTINE
!
!-----------------------------------------------------------------------
      SUBROUTINE COPY_T_to_D3STATE(LBi, UBi, LBj, UBj, UBk, UBt, ng, tile, eTimeIdx, t)
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
      integer iNode, i, j, k, iZ, idx
      integer iVar, itrc, ibio, NO_BOXES_XY_loc, tileS
      REAL(r8) eVal
      Print *, 'CP_T_D3 : stPelStateS=', stPelStateS, ' stPelStateE=', stPelStateE
      tileS = tile - first_tile(ng) + 1
      NO_BOXES_XY_loc = ListArrayWet(ng) % TheArr(tileS) % Nwetpoint
      DO iNode=1,NO_BOXES_XY_loc
         i = ListArrayWet(ng) % TheArr(tileS) % ListI(iNode)
         j = ListArrayWet(ng) % TheArr(tileS) % ListJ(iNode)
         DO k=1,NO_BOXES_Z
            iZ = k
            idx = iZ + NO_BOXES_Z * (iNode-1)
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
      integer iNode, i, j, k, iZ, idx
      integer iVar, itrc, ibio, NO_BOXES_XY_loc, tileS
      REAL(r8) eVal
      Print *, 'CP_D3_T : stPelStateS=', stPelStateS, ' stPelStateE=', stPelStateE
      tileS = tile - first_tile(ng) + 1
      NO_BOXES_XY_loc = ListArrayWet(ng) % TheArr(tileS) % Nwetpoint
      DO iNode=1,NO_BOXES_XY_loc
         i = ListArrayWet(ng) % TheArr(tileS) % ListI(iNode)
         j = ListArrayWet(ng) % TheArr(tileS) % ListJ(iNode)
         DO k=1,NO_BOXES_Z
            iZ = k
            idx = iZ + NO_BOXES_Z * (iNode-1)
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


      SUBROUTINE PRINT_AVERAGE_D3STATE(ng, tile)
      USE mod_param
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
      USE mem
      USE api_bfm
      IMPLICIT NONE
      integer, intent(in) :: ng, tile
      integer iNode, i, j, k, iZ, idx
      integer iVar, itrc, ibio, siz, TotalNb
      REAL(r8), allocatable :: ArrSum(:), ArrMin(:), ArrMax(:)
      REAL(r8) eVal, eAvg, eMin, eMax
!      Print *, 'PRINT_AVERAGE: stPelStateS=', stPelStateS, ' stPelStateE=', stPelStateE
!      Print *, 'stPelFluxS=', stPelFluxS, ' stPelFluxE=', stPelFluxE
      DO iVar=stPelStateS, stPelStateE
        itrc = iVar - stPelStateS + 1
        ibio = idbio(itrc)
        Print *, 'iVar=', iVar, ' itrc=', itrc, ' ibio=', ibio
      END DO
      Print *, 'Computing the average NO_BOXES_XY/Z=', NO_BOXES_XY, NO_BOXES_Z
      siz = stPelStateE + 1 - stPelStateS
      Allocate(ArrSum(siz), ArrMin(siz), ArrMax(siz))
      DO i=1,siz
        ArrSum(i) = 0
        ArrMin(i) =  100000
        ArrMax(i) = -100000
      END DO
      DO iNode=1,NO_BOXES_XY
         DO k=1,NO_BOXES_Z
            iZ = k
            idx = iZ + NO_BOXES_Z * (iNode-1)
            DO iVar=stPelStateS, stPelStateE
               itrc = iVar - stPelStateS + 1
               ibio = idbio(itrc)
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
      DO i=1,siz
        eAvg = ArrSum(i) / TotalNb
        eMin = ArrMin(i)
        eMax = ArrMax(i)
        Print *, 'i=', i, ' avg/min/max=', eAvg, eMin, eMax
      END DO
      deallocate(ArrSum, ArrMin, ArrMax)
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
      integer istat
      LOGICAL DoNestLayer
      integer Nwetpoint, tileLen
      integer NO_BOXES_Z_max, NO_BOXES_XY_max
      NAMELIST /SETTING_BFM_COUPL/ delt_bfm, AnalyticalInitD3STATE, CopyInitialToD3STATE, CopyD3STATEtoInitial, SourceTermD3STATE
      OPEN(file_id, FILE="bfm_input.nml")
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
      NO_BOXES_Z = NO_BOXES_Z_max
      NO_BOXES_XY = NO_BOXES_XY_max
      NO_BOXES_X = NO_BOXES_XY
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
      implicit none
      CALL COMPUTE_NO_BOXES_XY_Z_ALL
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
      !
      Print *, 'Printing average of D3STATE'
      CALL PRINT_AVERAGE_D3STATE(ng, tile)
      Print *, 'CopyInitialToD3STATE=', CopyInitialToD3STATE
      Print *, 'CopyD3STATEtoInitial=', CopyD3STATEtoInitial
      IF (CopyD3STATEtoInitial) THEN
        Print *, 'copying D3STATE to T'

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
        Print *, 'copying T to D3STATE'
        CALL COPY_T_to_D3STATE(LBi, UBi, LBj, UBj, UBk, UBt, ng, tile, 1, t)
      END IF
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
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
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
      integer iNode, i, j, k, idx, itrc, ibio
      integer iVar, iZ
      integer step
!
!  Assigning the STATE variables from the t array
!  ! We need to determine if the diagnostics need to be recomputed.
!
!      Print *, ' size(D3STATE,1) = ', size(D3STATE,1)
!      Print *, ' size(D3STATE,2) = ', size(D3STATE,2)
!      Print *, ' size(t,1)=', size(t,1)
!      Print *, ' size(t,2)=', size(t,2)
!      Print *, ' size(t,3)=', size(t,3)
!      Print *, ' size(t,4)=', size(t,4)
!      Print *, ' size(t,5)=', size(t,5)

!      Print *, 'Printing T average in biology_tile'
!      CALL Print_t_average(ng, tile, nstp)
!      CALL Print_t_average(ng, tile, nnew)
      PosMultiplier = PosMultiplier + 1
      IF (PosMultiplier == MULTIPLIER(ng)) THEN
        PosMultiplier = 0
        CALL SET_BOT_SURFINDICES(ng, tile)

!        Print *, 'Printing D3STATE before Source term integration'
!        CALL PRINT_AVERAGE_D3STATE(ng, tile)
        IF (SourceTermD3STATE) THEN
          CALL COPY_T_to_D3STATE(LBi, UBi, LBj, UBj, UBk, UBt, ng, tile, nstp, t)
!
!         Now the time stepping operations 
!         Right now we do Euler forward algorithm
!         It is probably needed to subdivide further the time interval
!         or to use a better integration method
!
          step = -1 ! not used
          call envforcing_bfm(ng, tile, step)
          call CalcVerticalExtinction( ) !     Compute extinction coefficient
          call EcologyDynamics           !     Compute reaction terms

          DO j=1,NO_D3_BOX_STATES
            IF (D3STATETYPE(j).ge.0) THEN
#ifndef EXPLICIT_SINK
              D3STATE(j,1:NO_BOXES) = D3STATE(j,1:NO_BOXES) +             &
     &          delt_bfm * D3SOURCE(j,1:NO_BOXES)
#else
              DO i=1,NO_BOXES
                 DO k=1,NO_D3_BOX_STATES
                    D3STATE(j,i) = D3STATE(j,i) +                         &
     &                   delt_bfm * (D3SOURCE(j,k,i) - D3SINK(j,k,i))
                 END DO
              END DO
#endif
            END IF
          END DO
!
!         Now copying back the field values
!
          CALL COPY_D3STATE_to_T(LBi, UBi, LBj, UBj, UBk, UBt, ng, tile, nnew, t)
!         Need to put code for the diagnostics. We do not put yet the dlux. Maybe never.
        END IF
      END IF
      END SUBROUTINE biology_tile
