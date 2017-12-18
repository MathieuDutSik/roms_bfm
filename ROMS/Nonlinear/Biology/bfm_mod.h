!
!svn $Id: hypoxia_srm_mod.h 851 2017-06-23 23:22:54Z arango $
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2017 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  Parameters for BFM model:                                           !
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

      integer :: iOxyg                  !  1: Dissolved oxygen concentration
      integer :: iPO4_                  !  2: Phosphate
      integer :: iNO3_                  !  3: Nitrate
      integer :: iNH4_                  !  4: Ammonium
      integer :: iO4n_                  !  5: 
      integer :: iSiOH                  !  6: Silicate
      integer :: iN6r_                  !  7:
      integer :: iB1c_                  !  8: Aerobic or anaerobic bacteria
      integer :: iB1n_                  !  9:
      integer :: iB1p_                  ! 10:
      integer :: iP1c_                  ! 11:
      integer :: iP1n_                  ! 12:
      integer :: iP1p_                  ! 13: 
      integer :: iP1l_                  ! 14:
      integer :: iP1s_                  ! 15: 
      integer :: iP2c_                  ! 16:
      integer :: iP2n_                  ! 17:
      integer :: iP2p_                  ! 18:
      integer :: iP2l_                  ! 19:
      integer :: iP3c_                  ! 20:
      integer :: iP3n_                  ! 21:
      integer :: iP3p_                  ! 22:
      integer :: iP3l_                  ! 23:
      integer :: iP4c_                  ! 24: 
      integer :: iP4n_                  ! 25:
      integer :: iP4p_                  ! 26:
      integer :: iP4l_                  ! 27: 
      integer :: iZ3c_                  ! 28:
      integer :: iZ3n_                  ! 29:
      integer :: iZ3p_                  ! 30:
      integer :: iZ4c_                  ! 31:
      integer :: iZ4n_                  ! 32:
      integer :: iZ4p_                  ! 33:
      integer :: iZ5c_                  ! 34:
      integer :: iZ5n_                  ! 35:
      integer :: iZ5p_                  ! 36:
      integer :: iZ6c_                  ! 37:
      integer :: iZ6n_                  ! 38:
      integer :: iZ6p_                  ! 39:
      integer :: iR1c_                  ! 40:
      integer :: iR1n_                  ! 41:
      integer :: iR1p_                  ! 42:
      integer :: iR2c_                  ! 43:
      integer :: iR3c_                  ! 44:
      integer :: iR6c_                  ! 45:
      integer :: iR6n_                  ! 46:
      integer :: iR6p_                  ! 47:
      integer :: iR6s_                  ! 48:
      integer :: iO3c_                  ! 49:
      integer :: iO3h_                  ! 50:
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
      iO4n_ = ic+5
      iSiOH = ic+6
      iN6r_ = ic+7
      iB1c_ = ic+8
      iB1n_ = ic+9
      iB1p_ = ic+10
      iP1c_ = ic+11
      iP1n_ = ic+12
      iP1p_ = ic+13
      iP1l_ = ic+14
      iP1s_ = ic+15
      iP2c_ = ic+16
      iP2n_ = ic+17
      iP2p_ = ic+18
      iP2l_ = ic+19
      iP3c_ = ic+20
      iP3n_ = ic+21
      iP3p_ = ic+22
      iP3l_ = ic+23
      iP4c_ = ic+24
      iP4n_ = ic+25
      iP4p_ = ic+26
      iP4l_ = ic+27
      iZ3c_ = ic+28
      iZ3n_ = ic+29
      iZ3p_ = ic+30
      iZ4c_ = ic+31
      iZ4n_ = ic+32
      iZ4p_ = ic+33
      iZ5c_ = ic+34
      iZ5n_ = ic+35
      iZ5p_ = ic+36
      iZ6c_ = ic+37
      iZ6n_ = ic+38
      iZ6p_ = ic+39
      iR1c_ = ic+40
      iR1n_ = ic+41
      iR1p_ = ic+42
      iR2c_ = ic+43
      iR3c_ = ic+44
      iR6c_ = ic+45
      iR6n_ = ic+46
      iR6p_ = ic+47
      iR6s_ = ic+48
      iO3c_ = ic+49
      iO3h_ = ic+50
      iEIR_ = ic+51
      iDIC_ = ic+52
      iChlo = ic+53
      siP1_ = ic+54
      siP2_ = ic+55
      siP3_ = ic+56
      siP4_ = ic+57
      eiP1_ = ic+58
      eiP2_ = ic+59
      eiP3_ = ic+60
      eiP4_ = ic+61
      ruPTc = ic+62
      ruZTc = ic+63

      RETURN
      END SUBROUTINE initialize_biology
