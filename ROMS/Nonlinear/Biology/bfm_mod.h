!
!svn $Id: hypoxia_srm_mod.h 851 2017-06-23 23:22:54Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Parameters for Hypoxia Simple Respiration Model:                    !
!                                                                      !
!  BioIter        Maximum number of iterations to achieve convergence  !
!                   of the nonlinear solution.                         !
!  ResRate        Total biological respiration rate (1/day).           !
!                                                                      !
!=======================================================================
!
      USE mod_param
!
      implicit none
!
!  Set biological tracer identification indices.
!
      integer, allocatable :: idbio(:)  ! Biological tracers

      integer :: iOxyg                  ! Dissolved oxygen concentration
      integer :: iPO4_                  ! Phosphate
      integer :: iNO3_                  ! Nitrate
      integer :: iNH4_                  ! Ammonium
      integer :: iSiOH                  ! Silicate
      integer :: iB1c_                  ! Aerobic or anaerobic bacteria
      integer :: iP1c_                  !
      integer :: iP1n_                  !
      integer :: iP1p_                  !
      integer :: iP1l_                  !
      integer :: iP1s_                  !
      integer :: iP2c_                  !
      integer :: iP2n_                  !
      integer :: iP2p_                  !
      integer :: iP2l_                  !
      integer :: iP3c_                  !
      integer :: iP3n_                  !
      integer :: iP3p_                  !
      integer :: iP3l_                  !
      integer :: iP4c_                  !
      integer :: iP4n_                  !
      integer :: iP4p_                  !
      integer :: iP4l_                  !
      integer :: iZ3c_                  !
      integer :: iZ4c_                  !
      integer :: iZ5c_                  !
      integer :: iZ6c_                  !
      integer :: iR1c_                  !
      integer :: iR1n_                  !
      integer :: iP1p_                  !
      integer :: iR2c_                  !
      integer :: iR6c_                  !
      integer :: iR6n_                  !
      integer :: iR6p_                  !
      integer :: iR6s_                  !
      integer :: iEIR_                  !
      integer :: iDIC_                  !
      integer :: iChlo                  !
      integer :: siP1_                  !
      integer :: siP2_                  !
      integer :: siP3_                  !
      integer :: siP4_                  !
      integer :: eiP1_                  !
      integer :: eiP2_                  !
      integer :: eiP3_                  !
      integer :: eiP4_                  !
      integer :: ruPTc                  !
      integer :: ruZTc                  !
!
!  Biological parameters.
!
      TYPE ARRAY_WET
      integer Nwetpoint
      integer, allocatable :: ListI(:)
      integer, allocatable :: ListJ(:)
      END TYPE ARRAY_WET
!
      integer, allocatable :: NO_BOXES_Z(:), NO_BOXES_XY(:), NO_BOXES(:)
      ARRAY_WET, allocatable :: ListArrayWet(:)
      
      CONTAINS

      SUBROUTINE initialize_biology
!
!=======================================================================
!                                                                      !
!  This routine sets several variables needed by the biology model.    !
!  It allocates and assigns biological tracers indices.                !
!                                                                      !
!=======================================================================
!
!  Local variable declarations
!
      integer :: i, ic
!
!-----------------------------------------------------------------------
!  Set number of biological tracers.
!-----------------------------------------------------------------------
!
      NBT=48

#if defined DIAGNOSTICS && defined DIAGNOSTICS_BIO
!
!-----------------------------------------------------------------------
!  Set sources and sinks biology diagnostic parameters.
!-----------------------------------------------------------------------
!
!  Set number of diagnostics terms.
!
      NDbio2d=0
      NDbio3d=0
!
!  Initialize biology diagnostic indices.
!
#endif
!
!-----------------------------------------------------------------------
!  Allocate various module variables.
!-----------------------------------------------------------------------
!
      IF (.not.allocated(BioIter)) THEN
        allocate ( BioIter(Ngrids) )
      END IF

      IF (.not.allocated(ResRate)) THEN
        allocate ( ResRate(Ngrids) )
      END IF
!
!  Allocate biological tracer vector.
!
      IF (.not.allocated(idbio)) THEN
        allocate ( idbio(NBT) )
      END IF

#if defined DIAGNOSTICS && defined DIAGNOSTICS_BIO
!
!  Allocate biological diagnostics vectors
!
      IF (.not.allocated(iDbio2)) THEN
        allocate ( iDbio2(NDbio2d) )
      END IF
#endif
!
!-----------------------------------------------------------------------
!  Initialize tracer identification indices.
!-----------------------------------------------------------------------
!
      ic=NAT+NPT+NCS+NNS
      DO i=1,NBT
        idbio(i)=ic+i
      END DO
      iOxyg = ic+1
      iPO4_ = ic+2
      iNO3_ = ic+3
      iNH4_ = ic+4
      iSiOH = ic+5
      iB1c_ = ic+6
      iP1c_ = ic+7
      iP1n_ = ic+8
      iP1p_ = ic+9
      iP1l_ = ic+10
      iP1s_ = ic+11
      iP2c_ = ic+12
      iP2n_ = ic+13
      iP2p_ = ic+14
      iP2l_ = ic+15
      iP3c_ = ic+16
      iP3n_ = ic+17
      iP3p_ = ic+18
      iP3l_ = ic+19
      iP4c_ = ic+20
      iP4n_ = ic+21
      iP4p_ = ic+22
      iP4l_ = ic+23
      iZ3c_ = ic+24
      iZ4c_ = ic+25
      iZ5c_ = ic+26
      iZ6c_ = ic+27
      iR1c_ = ic+28
      iR1n_ = ic+29
      iP1p_ = ic+30
      iR2c_ = ic+31
      iR6c_ = ic+32
      iR6n_ = ic+33
      iR6p_ = ic+34
      iR6s_ = ic+35
      iEIR_ = ic+36
      iDIC_ = ic+37
      iChlo = ic+38
      siP1_ = ic+39
      siP2_ = ic+40
      siP3_ = ic+41
      siP4_ = ic+42
      eiP1_ = ic+43
      eiP2_ = ic+44
      eiP3_ = ic+45
      eiP4_ = ic+46
      ruPTc = ic+47
      ruZTc = ic+48

      RETURN
      END SUBROUTINE initialize_biology
