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
      call set_bfm_fields_from_roms(ng, tile)
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
      DO idx=1,NO_BOXES_XY
         i = ListArrayWet(ng) % ListI(idx)
         j = ListArrayWet(ng) % ListJ(idx)
         lat = GRID(ng) % latp(i,j)
         ! The 
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
      SUBROUTINE INIT_BFM_SYSTEM_VARIABLE(ng, tile)
      USE mod_grid
      USE mod_param
      USE mod_biology
      USE api_bfm
      USE mem
      USE init_var_bfm_local, only : init_organic_constituents
      IMPLICIT NONE
      integer, intent(in) :: ng, tile
      logical, SAVE :: IsInitArray = .FALSE.
      logical, SAVE :: IsInitBFM = .FALSE.
      integer, SAVE :: NO_BOXES_XY_max = 0
      integer, SAVE :: NO_BOXES_Z_max = 0
      integer Nwetpoint, eProd
      integer idx, i, j
      integer idxBOT, idxSRF, eN, iNode
      integer istat
      integer bio_setup_loc
      integer, parameter :: namlst = 10
      integer, parameter :: file_id = 1453
      NAMELIST /PROC/ delt_bfm
# include "set_bounds.h"
      OPEN(file_id, FILE="bfm_input.nml")
      READ(file_id, NML = PROC)
      CLOSE(file_id)
      IF (.NOT. IsInitArray) THEN
        IsInitArray=.TRUE.
        allocate(AllocatedIJK(Ngrids))
        AllocatedIJK(:) = .FALSE.
        allocate(NO_BOXES_Z_arr(Ngrids), stat=istat)
        allocate(NO_BOXES_XY_arr(Ngrids), stat=istat)
        allocate(NO_BOXES_arr(Ngrids), stat=istat)
        allocate(ListArrayWet(Ngrids), stat=istat)
      END IF
      IF (.NOT. AllocatedIJK(ng)) THEN
        AllocatedIJK(ng) = .TRUE.
        NO_BOXES_Z_arr(ng) = N(ng)
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
        NO_BOXES_XY_arr(ng) = Nwetpoint
        IF (Nwetpoint .gt. NO_BOXES_XY_max) THEN
          NO_BOXES_XY_max = Nwetpoint
        END IF
        IF (N(ng) .gt. NO_BOXES_Z_max) THEN
          NO_BOXES_Z_max = N(ng)
        END IF
        NO_BOXES_XY = NO_BOXES_XY_max
        NO_BOXES_X  = NO_BOXES_XY_max
        NO_BOXES_Y  = 1
        NO_BOXES_Z  = NO_BOXES_Z_max
        !
        eProd = Nwetpoint * NO_BOXES_Z
        NO_BOXES_XY = Nwetpoint
        NO_BOXES_arr(ng) = eProd
        ListArrayWet(ng) % Nwetpoint = Nwetpoint
        allocate(ListArrayWet(ng) % ListI(Nwetpoint), stat=istat)
        allocate(ListArrayWet(ng) % ListJ(Nwetpoint), stat=istat)
! The arrays of bottom and surface
        allocate(BOTindices(Nwetpoint), stat=istat)
        allocate(SRFindices(Nwetpoint), stat=istat)
        eN = NO_BOXES_Z_arr(ng)
        DO iNode=1,Nwetpoint
          idxBOT = 1  + eN * (iNode-1)
          idxSRF = eN + eN * (iNode-1)
          BOTindices(iNode) = idxBOT
          SRFindices(iNode) = idxSRF
        END DO
! The WET arrays.
        idx=0
        DO j=Jstr-1,JendR
          DO i=Istr-1,IendR
#ifdef MASKING
            IF (GRID(ng) % rmask(i,j) .eq. 1) THEN
#endif
               idx = idx + 1
               ListArrayWet(ng) % ListI(idx) = i
               ListArrayWet(ng) % ListJ(idx) = j
#ifdef MASKING
            END IF
#endif
          END DO
        END DO
        
        IF ((NO_BOXES_XY_max .ne. 0).and.                               &
     &       (NO_BOXES_XY_max .gt. NO_BOXES_XY)) THEN
          Print *, 'NO_BOXES_XY_max=', NO_BOXES_XY_max
          Print *, 'NO_BOXES_XY=', NO_BOXES_XY
          Print *, 'We have a problem of allocation'
          STOP
        END IF
        NO_BOXES    = NO_BOXES_XY * NO_BOXES_Z
        NO_STATES   = NO_D3_BOX_STATES * NO_BOXES + NO_BOXES_XY
      END IF
      IF (.NOT. IsInitBFM) THEN
        ! Initialise the BFM with standalone settings
        call init_bfm(namlst)
        ! Initialise state variable names and diagnostics
        call set_var_info_bfm
        ! Allocate memory and give initial values
        ! to the pelagic system
        ! We need the variable well set for this to work.
        bio_setup_loc=1
        call init_var_bfm(bio_setup_loc)
        ! Initialize internal constitutents of functional groups
        call init_organic_constituents()
        ! Need to set up bfmtime adequately for the runs.
        ! Need also to set up the wind.
        call init_envforcing_bfm
        ! Read restart file (if flag)
        ! Overwrite previous initialization
        ! Initialise the diagnostic variables
        call CalcVerticalExtinction( )
        call CalcChlorophylla( )
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
!  Setting up the dimensions and initializing BFM if needed
!
      CALL INIT_BFM_SYSTEM_VARIABLE(ng, tile)
!
!  Assigning the STATE variables from the t array
!  ! We need to determine if the diagnostics need to be recomputed.
!
      Print *, ' size(D3STATE,1) = ', size(D3STATE,1)
      Print *, ' size(D3STATE,2) = ', size(D3STATE,2)
      Print *, ' size(t,1)=', size(t,1)
      Print *, ' size(t,2)=', size(t,2)
      Print *, ' size(t,3)=', size(t,3)
      Print *, ' size(t,4)=', size(t,4)
      Print *, ' size(t,5)=', size(t,5)
      DO iNode=1,NO_BOXES_XY
         i = ListArrayWet(ng) % ListI(iNode)
         j = ListArrayWet(ng) % ListJ(iNode)
         DO k=1,NO_BOXES_Z
            iZ = k
            idx = iZ + NO_BOXES_Z * (iNode-1)
!            Print *, '1: iZ=', iZ, ' iNode=', iNode, ' k=', k, ' idx=', idx
            DO iVar=stPelStateS, stPelStateE
               itrc = iVar - stPelStateS + 1
               ibio = idbio(itrc)
!               Print *, 'iVar=', iVar, ' itrc=', itrc, ' ibio=', ibio
!               Print *, 'i=', i, ' j=', j, ' k=', k, ' nstp=', nstp
               eVal = t(i, j, k, nstp, ibio)
!               Print *, ' eVal=', eVal
               D3STATE(itrc, idx) = eVal
            END DO
         END DO
      END DO
!
!     Now the time stepping operations 
!     Right now we do Euler forward algorithm
!     It is probably needed to subdivide further the time interval
!     or to use a better integration method
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
     &                 delt_bfm * (D3SOURCE(j,k,i) - D3SINK(j,k,i))
               END DO
            END DO
#endif
         END IF
      END DO
!
!     Now copying back the field values
!
      DO iNode=1,NO_BOXES_XY
         i = ListArrayWet(ng) % ListI(iNode)
         j = ListArrayWet(ng) % ListJ(iNode)
         DO k=1,NO_BOXES_Z
            iZ = k
            idx = iZ + NO_BOXES_Z * (iNode-1)
!            Print *, '2: iZ=', iZ, ' iNode=', iNode, ' k=', k, ' idx=', idx
            DO iVar=stPelStateS, stPelStateE
               itrc = iVar - stPelStateS + 1
               ibio = idbio(itrc)
!               Print *, 'iVar=', iVar, ' itrc=', itrc, ' ibio=', ibio
               t(i, j, k, nnew, ibio) = D3STATE(itrc, idx)
            END DO
         END DO
      END DO
! Need to put code for the diagnostics. We do not put yet the dlux. Maybe never.

      
      
      RETURN
      END SUBROUTINE biology_tile
