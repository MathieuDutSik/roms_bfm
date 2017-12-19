      SUBROUTINE read_BioPar (model, inp, out, Lwrite)
!
!svn $Id: hypoxia_srm_inp.h 851 2017-06-23 23:22:54Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This routine reads in Hypoxia Simple Respiration Model biological   !
!  model input parameters. They are specified in input script          !
!  "hypoxia_srm.in".                                                   !
!                                                                      !
!=======================================================================
!
      USE mod_param
      USE mod_parallel
      USE mod_biology
      USE mod_ncparam
      USE mod_scalars
      USE mod_grid
!
      implicit none
!
!  Imported variable declarations
!
      logical, intent(in) :: Lwrite
      integer, intent(in) :: model, inp, out
!
!  Local variable declarations.
!
      integer :: Npts, Nval
      integer :: iTrcStr, iTrcEnd
      integer :: i, ifield, igrid, itracer, itrc, ng, nline, status

      logical, dimension(NBT,Ngrids) :: Ltrc

      real(r8), dimension(NBT,Ngrids) :: Rbio

      real(r8), dimension(100) :: Rval

      character (len=40 ) :: KeyWord
      character (len=256) :: line
      character (len=256), dimension(200) :: Cval

      integer :: namlst=10
      integer eProd, idx, Nwetpoint
      integer NO_BOXES_Z_max, NO_BOXES_XY_max
# include "set_bounds.h"
!
!-----------------------------------------------------------------------
!  Initialize.
!-----------------------------------------------------------------------
!
      igrid=1                            ! nested grid counter
      itracer=0                          ! LBC tracer counter
      iTrcStr=1                          ! first LBC tracer to process
      iTrcEnd=NBT                        ! last  LBC tracer to process
      nline=0                            ! LBC multi-line counter
!
!-----------------------------------------------------------------------
!  Read in red tide (Stock et al., 2005; He et al., 2008) model
!  parameters.
!-----------------------------------------------------------------------
!
      
      allocate(NO_BOXES_Z_arr(Ngrids), stat=istat)
      allocate(NO_BOXES_XY_arr(Ngrids), stat=istat)
      allocate(NO_BOXES_arr(Ngrids), stat=istat)
      allocate(ListArrayWet(Ngrids), stat=istat)
      NO_BOXES_XY_max = 0
      NO_BOXES_Z_max = 0
      DO ng=1,Ngrids
         NO_BOXES_Z_arr(ng) = N(ng)
         Nwetpoint=0
         DO j=Jstr-1,JendR
            DO i=Istr-1,IendR
#ifdef MASKING
               IF (GRIDS(ng) % rmask(i,j) .eq. 1) THEN
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
         eProd = Nwetpoint * N(ng)
         NO_BOXES_arr(ng) = eProd
         ListArrayWet(ng) % Nwetpoint = Nwetpoint
         allocate(ListArrayWet(ng) % ListI(Nwetpoint), stat=istat)
         allocate(ListArrayWet(ng) % ListJ(Nwetpoint), stat=istat)
!
         idx=0
         DO j=Jstr-1,JendR
            DO i=Istr-1,IendR
#ifdef MASKING
               IF (GRIDS(ng) % rmask(i,j) .eq. 1) THEN
#endif
                  ListArrayWet(ng) % ListI(idx) = i
                  ListArrayWet(ng) % ListJ(idx) = j
#ifdef MASKING
               END IF
#endif
            END DO
         END DO
!
      END DO
!
      NO_BOXES_XY = NO_BOXES_XY_max
      NO_BOXES_X = NO_BOXES_XY_max
      NO_BOXES_Y = 1
      NO_BOXES_Z = NO_BOXES_Z_max
      NO_BOXES = NO_BOXES_XY * NO_BOXES_Z
      NO_STATES   = NO_D3_BOX_STATES * NO_BOXES + NO_BOXES_XY
      
! Initialise the BFM with standalone settings
      call init_bfm(namlst)
! Initialise state variable names and diagnostics
      call set_var_info_bfm
! Allocate memory and give initial values
! to the pelagic system
! We need the variable well set for this to work.
      call init_var_bfm(bio_setup)
! Initialize internal constitutents of functional groups
      call init_organic_constituents()

!     Need to set up bfmtime adequately for the runs.
!     Need also to set up the wind.
      call init_envforcing_bfm
!     Read restart file (if flag)
!     Overwrite previous initialization
!     Initialise the diagnostic variables
      call CalcVerticalExtinction( )
      call CalcChlorophylla( )

      

      END SUBROUTINE read_BioPar
