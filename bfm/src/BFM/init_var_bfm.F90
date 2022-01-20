#include"cppdefs_bfm.h"


MODULE init_var_bfm_local
   USE global_mem,ONLY: RLEN,ZERO
   IMPLICIT NONE

   real(RLEN),parameter :: nc_ratio_default = 0.0126_RLEN    ! Redfield
   real(RLEN),parameter :: pc_ratio_default = 0.7862e-3_RLEN ! Redfield
   real(RLEN),parameter :: lc_ratio_default = 0.03_RLEN      ! standard diatom value
   real(RLEN),parameter :: sc_ratio_default = 0.0145_RLEN    ! Redfield
   real(RLEN),parameter :: hc_ratio_default = ZERO


   character(len=*),parameter :: bfm_init_fname      = 'BFM_General.nml'
   character(len=*),parameter :: bfm_init_fname_ice  = 'BFM_General.nml'
   character(len=*),parameter :: bfm_init_fname_ben  = 'BFM_General.nml'

   CONTAINS

   subroutine init_constituents( c,n,p,l,s,h,nc,pc,lc,sc,hc)
     use global_mem, only: RLEN,ZERO
     IMPLICIT NONE
     real(RLEN),dimension(:),intent(in)             :: c
     real(RLEN),intent(in),optional                 :: nc,pc,lc,sc,hc
     real(RLEN),dimension(:),intent(inout),optional :: n,p,l,s,h
     real(RLEN)   :: nc_ratio,pc_ratio,lc_ratio,sc_ratio,hc_ratio

         nc_ratio = nc_ratio_default
         if (present(nc)) then
           if (nc>ZERO) nc_ratio = nc
         end if

         pc_ratio = pc_ratio_default
         if (present(pc)) then
           if (pc>ZERO) pc_ratio = pc
         end if

         lc_ratio = lc_ratio_default
         if (present(lc)) then
           if (lc>ZERO) lc_ratio = lc
         end if

         sc_ratio = sc_ratio_default
         if (present(sc)) then
           if (sc>ZERO) sc_ratio = sc
         end if

         hc_ratio = hc_ratio_default
         if (present(hc)) then
           if (hc>ZERO) hc_ratio = hc
         end if

         if (present(n)) then
           where (n==ZERO)
             n = nc_ratio*c
           end where
         end if
         if (present(p)) then
           where (p==ZERO)
             p = pc_ratio*c
           end where
         end if
         if (present(l)) then
           where (l==ZERO)
             l = lc_ratio*c
           end where
         end if
         if (present(s)) then
           where (s==ZERO)
             s = sc_ratio*c
           end where
         end if
         if (present(h)) then
           where (h==ZERO)
             h = hc_ratio*c
           end where
         end if
   end subroutine init_constituents


!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Initialise constituents of BFM variables
!
! !INTERFACE:
   subroutine init_organic_constituents(InitD3STATE)
!
! !DESCRIPTION:
!  Call the initialization of internal ratios for the living
! and non-living organic components
!
! !USES:
   use mem
   use mem_PelGlobal
   use mem_PelChem
   use mem_PelBac
   use mem_MesoZoo
   use mem_MicroZoo
   use mem_Phyto

   use constants, ONLY: HOURS_PER_DAY
   use mem_Param, ONLY: p_small,          &
                        CalcPelagicFlag,  &
                        CalcBenthicFlag,  &
                        CalcSeaiceFlag,   &
                        CalcPelChemistry, &
                        CalcTransportFlag
#ifdef INCLUDE_SEAICE
   use mem_SeaiceAlgae
   use mem_SeaiceBac
   use mem_SeaiceZoo
#endif
   IMPLICIT NONE
   logical, intent(in) :: InitD3STATE
!
! !INPUT PARAMETERS:
!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!
! !LOCAL VARIABLES:
   integer              :: i
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
   !---------------------------------------------
   ! Initialise other pelagic internal components
   ! with Redfield
   !---------------------------------------------
   IF (InitD3STATE) THEN
     do i = 1 , ( iiPelBacteria )
       call init_constituents( c=PelBacteria(i,iiC), &
         n=D3STATE(ppPelBacteria(i,iiN),:),  &
         p=D3STATE(ppPelBacteria(i,iiP),:),  &
         nc=p_qncPBA(i),  pc=p_qpcPBA(i) )
     end do
     do i = 1 , ( iiPhytoPlankton )
       call init_constituents( c=PhytoPlankton(i,iiC), &
         n=D3STATE(ppPhytoPlankton(i,iiN),:),  &
         p=D3STATE(ppPhytoPlankton(i,iiP),:),  &
         l=D3STATE(ppPhytoPlankton(i,iiL),:),  &
         s=D3STATE(ppPhytoPlankton(i,iiS),:),  &
         nc=p_qncPPY(i),  pc=p_qpcPPY(i),  lc=p_qlcPPY(i),  sc=p_qscPPY(i) )
     end do
     do i = 1 , ( iiMesoZooPlankton )
       call init_constituents( c=MesoZooPlankton(i,iiC), &
         n=D3STATE(ppMesoZooPlankton(i,iiN),:),  &
         p=D3STATE(ppMesoZooPlankton(i,iiP),:),  &
         nc=p_qncMEZ(i),  pc=p_qpcMEZ(i) )
     end do
     do i = 1 , ( iiMicroZooPlankton )
       call init_constituents( c=MicroZooPlankton(i,iiC), &
         n=D3STATE(ppMicroZooPlankton(i,iiN),:),  &
         p=D3STATE(ppMicroZooPlankton(i,iiP),:),  &
         nc=p_qncMIZ(i),  pc=p_qpcMIZ(i) )
     end do
     do i = 1 , ( iiPelDetritus )
       call init_constituents( c=PelDetritus(i,iiC), &
         n=D3STATE(ppPelDetritus(i,iiN),:),  &
         p=D3STATE(ppPelDetritus(i,iiP),:),  &
         s=D3STATE(ppPelDetritus(i,iiS),:) )
     end do
   END IF

#ifdef INCLUDE_BEN
   !---------------------------------------------
   ! Initialise other benthic internal components
   ! with Redfield
   !---------------------------------------------



#endif

#ifdef INCLUDE_SEAICE
   !---------------------------------------------
   ! Initialise other seaice internal components
   ! with Redfield
   !---------------------------------------------



#endif

end subroutine init_organic_constituents


END MODULE init_var_bfm_local


!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Initialise BFM variables
!
! !INTERFACE:
   subroutine init_var_bfm(setup, InitD3STATE)
!
! !DESCRIPTION:
!  Allocate BFM variables and give initial values of
!  parameters and state variables
!
! !USES:
#ifndef NOT_STANDALONE
   use api_bfm
   use global_mem
#endif
#ifdef BFM_GOTM
   use bio_var
   use bio_bfm
   use global_mem, ONLY: RLEN,ZERO,ONE
#endif
   use mem
   use mem_PelGlobal
   use mem_PelChem
   use mem_PelBac
   use mem_MesoZoo
   use mem_MicroZoo
   use mem_Phyto
   use init_var_bfm_local

   use constants, ONLY: HOURS_PER_DAY
   use mem_Param, ONLY: p_small,          &
                        CalcPelagicFlag,  &
                        CalcBenthicFlag,  &
                        CalcSeaiceFlag,   &
                        CalcPelChemistry, &
                        CalcTransportFlag

   use mem_Param, ONLY: AssignPelBenFluxesInBFMFlag
   use string_functions, ONLY: getseq_number,empty

#if defined INCLUDE_BEN && defined INCLUDE_BENPROFILES
   use mem_Param, ONLY: p_d_tot, p_sedlevels,p_sedsigma
#endif
#ifdef INCLUDE_SEAICE
   use mem_SeaiceAlgae
   use mem_SeaiceBac
   use mem_SeaiceZoo
#endif

   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in) :: setup
   logical, intent(in) :: InitD3STATE

!
! !REVISION HISTORY:
!  Original author(s): Marcello Vichi
!
! !LOCAL VARIABLES:
   integer,parameter    :: namlst=10
   integer              :: icontrol,i,j,n,Flun
   integer,parameter    :: NSAVE=120  ! Maximum no variables which can be saved
   character(len=64),dimension(NSAVE):: var_save
   character(len=64),dimension(NSAVE):: ave_save

!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! Definition of Initial Pelagic (D3) state variables
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

      real(RLEN) :: O2o0, N1p0, N3n0, N4n0, O4n0, N5s0, N6r0, B1c0,     &
     & P1c0, P1n0, P1p0, P1l0, P1s0, P2c0, P2n0, P2p0, P2l0, P3c0,      &
     & P3l0, P4c0, P4n0, P4p0, P4l0, Z3c0, Z3n0, Z3p0, Z4c0, Z4n0,      &
     & Z5n0, Z5p0, Z6c0, Z6n0, Z6p0, R1c0, R1n0, R1p0, R2c0, R3c0,      &
     & R6p0, R6s0, O3c0, O3h0, R6c0, R6n0, Z4p0, Z5c0, P3n0, P3p0,      &
     & B1n0, B1p0


       namelist /bfm_init_nml/ O2o0, N1p0, N3n0, N4n0, O4n0, N5s0,      &
     & B1n0, B1p0, P1c0, P1n0, P1p0, P1l0, P1s0, P2c0, P2n0, P2p0,      &
     & P3n0, P3p0, P3l0, P4c0, P4n0, P4p0, P4l0, Z3c0, Z3n0, Z3p0,      &
     & Z4p0, Z5c0, Z5n0, Z5p0, Z6c0, Z6n0, Z6p0, R1c0, R1n0, R1p0,      &
     & R6c0, R6n0, R6p0, R6s0, O3c0, O3h0, R2c0, R3c0, Z4c0, Z4n0,      &
     & P2l0, P3c0, N6r0, B1c0


   namelist /bfm_init_nml/ surface_flux_method,       &
                           bottom_flux_method,        &
                           n_surface_fluxes, InitVar

   namelist /bfm_save_nml/ var_save, ave_save

#ifdef INCLUDE_BEN



#endif

#ifdef INCLUDE_SEAICE



#endif


! COPYING
!
!   Copyright (C) 2015 BFM System Team (bfm_st@lists.cmcc.it)
!   Copyright (C) 2006 P. Ruardij and Marcello Vichi
!   (rua@nioz.nl, vichi@bo.ingv.it)
!
!   This program is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation;
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTEABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!EOP
!-----------------------------------------------------------------------
!BOC

   LEVEL2 'init_var_bfm'
   !---------------------------------------------
   ! Give zero initial values
   ! Overwritten by namelist parameters
   !---------------------------------------------
   surface_flux_method = -1
   bottom_flux_method = 0
   n_surface_fluxes = 1

   !---------------------------------------------
   ! Pelagic variables
   !---------------------------------------------

   O2o0 = _ZERO_
   N1p0 = _ZERO_
   N3n0 = _ZERO_
   N4n0 = _ZERO_
   O4n0 = _ZERO_
   N5s0 = _ZERO_
   N6r0 = _ZERO_
   B1c0 = _ZERO_
   B1n0 = _ZERO_
   B1p0 = _ZERO_
   P1c0 = _ZERO_
   P1n0 = _ZERO_
   P1p0 = _ZERO_
   P1l0 = _ZERO_
   P1s0 = _ZERO_
   P2c0 = _ZERO_
   P2n0 = _ZERO_
   P2p0 = _ZERO_
   P2l0 = _ZERO_
   P3c0 = _ZERO_
   P3n0 = _ZERO_
   P3p0 = _ZERO_
   P3l0 = _ZERO_
   P4c0 = _ZERO_
   P4n0 = _ZERO_
   P4p0 = _ZERO_
   P4l0 = _ZERO_
   Z3c0 = _ZERO_
   Z3n0 = _ZERO_
   Z3p0 = _ZERO_
   Z4c0 = _ZERO_
   Z4n0 = _ZERO_
   Z4p0 = _ZERO_
   Z5c0 = _ZERO_
   Z5n0 = _ZERO_
   Z5p0 = _ZERO_
   Z6c0 = _ZERO_
   Z6n0 = _ZERO_
   Z6p0 = _ZERO_
   R1c0 = _ZERO_
   R1n0 = _ZERO_
   R1p0 = _ZERO_
   R2c0 = _ZERO_
   R3c0 = _ZERO_
   R6c0 = _ZERO_
   R6n0 = _ZERO_
   R6p0 = _ZERO_
   R6s0 = _ZERO_
   O3c0 = _ZERO_
   O3h0 = _ZERO_


#ifdef INCLUDE_BEN
   !---------------------------------------------
   ! Benthic variables
   !---------------------------------------------



#endif

#ifdef INCLUDE_SEAICE
   !---------------------------------------------
   ! Seaice variables
   !---------------------------------------------



#endif

   !---------------------------------------------
   ! Initialize the structured array that
   ! defines if a variable is initialized with
   ! data. The namelist values override the
   ! assignment
   !---------------------------------------------
   InitVar = InputInfo(0, ZERO, "dummy.nc","dummy",ZERO,ZERO,ZERO,ZERO,.FALSE.,.FALSE.,.FALSE.)

   !---------------------------------------------
   ! Open and read the namelist
   !---------------------------------------------
   icontrol=0
   Print *, 'bfm_init_fname=', bfm_init_fname
   eFile = TRIM(BFM_Prefix_NML) // bfm_init_fname
   open(namlst,file=TRIM(eFile),action='read',status='old',err=100)
   var_save=""
   ave_save=""
   var_ave=.false.
   read(namlst,nml=bfm_save_nml,err=101)
   close(namlst)
   icontrol=1
100 if ( icontrol == 0 ) then
     LEVEL3 'I could not open ',trim(bfm_init_fname)
     LEVEL3 'The initial values of the BFM variables are set to ZERO'
     LEVEL3 'If thats not what you want you have to supply ',trim(bfm_init_fname)
     icontrol=1
  end if
101 if ( icontrol == 0 ) then
     FATAL 'I could not read bfm_save_nml'
     stop 'init_var_bfm'
   end if

   icontrol=0
   open(namlst,file=TRIM(eFile),action='read',status='old',err=102)
   read(namlst,nml=bfm_init_nml,err=103)
   close(namlst)
   icontrol=1
102 if ( icontrol == 0 ) then
     LEVEL3 'I could not open ',trim(bfm_init_fname)
     LEVEL3 'The initial values of the BFM variables are set to ZERO'
     LEVEL3 'If thats not what you want you have to supply ',trim(bfm_init_fname)
     icontrol=1
   end if
103 if ( icontrol == 0 ) then
     FATAL 'I could not read bfm_init_nml'
     stop 'init_var_bfm'
   end if


#ifdef INCLUDE_BEN
   icontrol=0
   open(namlst,file=bfm_init_fname_ben,action='read',status='old',err=104)
   read(namlst,nml=bfm_init_nml_ben,err=105)
   close(namlst)
   icontrol=1
104 if ( icontrol == 0 ) then
     LEVEL3 'I could not open ',trim(bfm_init_fname_ben)
     LEVEL3 'The initial values of the BFM variables are set to ZERO'
     LEVEL3 'If thats not what you want you have to supply ',trim(bfm_init_fname_ben)
     icontrol=1
   end if
105 if ( icontrol == 0 ) then
     FATAL 'I could not read bfm_init_nml_ben'
     stop 'init_var_bfm'
   end if

#endif

#ifdef INCLUDE_SEAICE
   icontrol=0
   open(namlst,file=bfm_init_fname_ice,action='read',status='old',err=106)
   read(namlst,nml=bfm_init_nml_ice,err=107)
   close(namlst)
   icontrol=1
106 if ( icontrol == 0 ) then
     LEVEL3 'I could not open ',trim(bfm_init_fname_ice)
     LEVEL3 'The initial values of the BFM variables are set to ZERO'
     LEVEL3 'If thats not what you want you have to supply ',trim(bfm_init_fname_ice)
     icontrol=1
   end if
107 if ( icontrol == 0 ) then
     FATAL 'I could not read bfm_init_nml_ice'
     stop 'init_var_bfm'
   end if

#endif

   !---------------------------------------------
   ! Check variable to be saved and
   ! set the corresponding flag value in var_ids
   !---------------------------------------------
   do i=1,NSAVE
      if (.NOT.empty(var_save(i))) then
            j=getseq_number(var_save(i),var_names,stEnd,.TRUE.)
            if ( j > 0 ) var_ids(j)=-1
      end if
      if ( .NOT.empty(var_save(i)) .AND. j==0 ) then
            STDERR 'Warning: variable ',trim(var_save(i)),' does not exist!'
      end if
   end do
   do i=1,NSAVE
      if (.NOT.empty(ave_save(i))) then
         j=getseq_number(ave_save(i),var_names,stEnd,.TRUE.)
         if ( .NOT.empty(ave_save(i)) .AND. j==0 ) then
            STDERR 'Warning: variable ',trim(ave_save(i)),' does not exist!'
         else if ( var_ids(j) <0 ) then
            STDERR 'Warning: Variable ',trim(ave_save(i)), &
               ' is already selected for output in var_save'
         else if ( j > 0 ) then
            var_ids(j) = -1
            var_ave(j) = .true.
            ave_ctl = .true.
         end if
      end if
   end do

#ifdef BFM_GOTM
   !---------------------------------------------
   ! Create pointers
   !---------------------------------------------
    call pointers_gotm_bfm()
#endif

   !---------------------------------------------
   ! Initialize BFM parameters
   !---------------------------------------------
   call Initialize

   !---------------------------------------------
   ! Initially set the number of sun hours
   ! equal to the number of hours in a day.
   !---------------------------------------------
   SUNQ = HOURS_PER_DAY

   !---------------------------------------------
   ! Initialise pelagic state variables
   ! also if using a benthic-only setup
   ! (for boundary conditions)
   !---------------------------------------------

     O2o = O2o0
     N1p = N1p0
     N3n = N3n0
     N4n = N4n0
     O4n = O4n0
     N5s = N5s0
     N6r = N6r0
     B1c = B1c0
     B1n = B1n0
     B1p = B1p0
     P1c = P1c0
     P1n = P1n0
     P1p = P1p0
     P1l = P1l0
     P1s = P1s0
     P2c = P2c0
     P2n = P2n0
     P2p = P2p0
     P2l = P2l0
     P3c = P3c0
     P3n = P3n0
     P3p = P3p0
     P3l = P3l0
     P4c = P4c0
     P4n = P4n0
     P4p = P4p0
     P4l = P4l0
     Z3c = Z3c0
     Z3n = Z3n0
     Z3p = Z3p0
     Z4c = Z4c0
     Z4n = Z4n0
     Z4p = Z4p0
     Z5c = Z5c0
     Z5n = Z5n0
     Z5p = Z5p0
     Z6c = Z6c0
     Z6n = Z6n0
     Z6p = Z6p0
     R1c = R1c0
     R1n = R1n0
     R1p = R1p0
     R2c = R2c0
     R3c = R3c0
     R6c = R6c0
     R6n = R6n0
     R6p = R6p0
     R6s = R6s0
     O3c = O3c0
     O3h = O3h0


#ifdef INCLUDE_BEN
   !---------------------------------------------
   ! Initialise benthic state variables
   ! also if using a ben-only setup
   ! (for boundary conditions)
   !---------------------------------------------



#endif

#ifdef INCLUDE_SEAICE
   !---------------------------------------------
   ! Initialise seaice state variables
   ! also if using a seaice-only setup
   ! (for boundary conditions)
   !---------------------------------------------



#endif

   !---------------------------------------------
   ! Check setup settings
   ! and finalize initialization
   !---------------------------------------------
   select case (setup)
      case (0)
         LEVEL2 "Fully coupled system, Pelagic, Benthic, Seaice"
      case (1) ! Pelagic only
         LEVEL2 "Pelagic-only setup (bio_setup=1), Switching off other systems"
         CalcBenthicFlag = 0
         CalcSeaiceFlag  = .FALSE.


      case (2) ! Benthic only
         LEVEL2 "Benthic-only setup (bio_setup=2), Switching off other systems"
         CalcPelagicFlag = .FALSE.
         CalcSeaiceFlag  = .FALSE.
         CalcPelBacteria = .FALSE.
         CalcPhytoPlankton = .FALSE.
         CalcMesoZooPlankton = .FALSE.
         CalcMicroZooPlankton = .FALSE.
         CalcPelDetritus = .FALSE.
         CalcInorganic = .FALSE.


      case (3) ! Pelagic-Benthic coupling
         LEVEL2 "Pelagic-Benthic setup (bio_setup=3), Switching off other systems"
         CalcSeaiceFlag  = .FALSE.

      case (4) ! SeaIce only
         LEVEL2 "Seaice-only setup (bio_setup=4), Switching off other systems"
         CalcPelagicFlag = .FALSE.
         CalcBenthicFlag = 0
         CalcPelBacteria = .FALSE.
         CalcPhytoPlankton = .FALSE.
         CalcMesoZooPlankton = .FALSE.
         CalcMicroZooPlankton = .FALSE.
         CalcPelDetritus = .FALSE.
         CalcInorganic = .FALSE.


      case (5) ! Pelagic-SeaIce coupling
         LEVEL2 "Pelagic-Seaice setup (bio_setup=5), Switching off other systems"
         CalcBenthicFlag = 0

   end select

#if defined INCLUDE_BEN
   !---------------------------------------------
   ! Check benthic model
   !---------------------------------------------
   select case (CalcBenthicFlag)
     case (0)
        LEVEL3 "Benthic model is: not used"
     case (1)
        LEVEL3 "Benthic model is: simple nutrient return"
     case (2)
        LEVEL3 "Benthic model is: benthos + intermediate nutrient return"
     case (3)
        LEVEL3 "Benthic model is: benthos + Ruardij & Van Raaphorst"
   end select
#endif

#if defined key_obcbfm
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  ! All variables are inizialized on obc default
  ! D3STATEOBC(:)=OBCSTATES
  ! Put D3STATEOBC(..)=NOOBCSTATES to exclude
  !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
    D3STATEOBC(:)=OBCSTATES
#endif

   !---------------------------------------------
   ! Check for transport flag
   !---------------------------------------------
#ifdef BFM_STANDALONE
   D3STATETYPE(:) = NOTRANSPORT
#else
   if (.NOT.CalcTransportFlag) D3STATETYPE(:) = NOTRANSPORT
#endif

   !----------------------------------------------------
   ! Zeroing of the switched off pelagic state variables
   !----------------------------------------------------

   do j = 1, iiPelBacteria
     if (.NOT.CalcPelBacteria(j)) then
       do i = 1,iiLastElement
         if ( ppPelBacteria(j,i) /= 0 ) then
           IF (InitD3STATE) THEN
             D3STATE(ppPelBacteria(j,i),:) = p_small
           END IF
           D3STATETYPE(ppPelBacteria(j,i)) = OFF
#if defined key_obcbfm
           D3STATEOBC(ppPelBacteria(j,i)) = NOOBCSTATES
#endif
         end if
       end do
     end if
   end do

   do j = 1, iiPhytoPlankton
     if (.NOT.CalcPhytoPlankton(j)) then
       do i = 1,iiLastElement
         if ( ppPhytoPlankton(j,i) /= 0 ) then
           IF (InitD3STATE) THEN
             D3STATE(ppPhytoPlankton(j,i),:) = p_small
           END IF
           D3STATETYPE(ppPhytoPlankton(j,i)) = OFF
#if defined key_obcbfm
           D3STATEOBC(ppPhytoPlankton(j,i)) = NOOBCSTATES
#endif
         end if
       end do
     end if
   end do

   do j = 1, iiMesoZooPlankton
     if (.NOT.CalcMesoZooPlankton(j)) then
       do i = 1,iiLastElement
         if ( ppMesoZooPlankton(j,i) /= 0 ) then
           IF (InitD3STATE) THEN
             D3STATE(ppMesoZooPlankton(j,i),:) = p_small
           END IF
           D3STATETYPE(ppMesoZooPlankton(j,i)) = OFF
#if defined key_obcbfm
           D3STATEOBC(ppMesoZooPlankton(j,i)) = NOOBCSTATES
#endif
         end if
       end do
     end if
   end do

   do j = 1, iiMicroZooPlankton
     if (.NOT.CalcMicroZooPlankton(j)) then
       do i = 1,iiLastElement
         if ( ppMicroZooPlankton(j,i) /= 0 ) then
           IF (InitD3STATE) THEN
             D3STATE(ppMicroZooPlankton(j,i),:) = p_small
           END IF
           D3STATETYPE(ppMicroZooPlankton(j,i)) = OFF
#if defined key_obcbfm
           D3STATEOBC(ppMicroZooPlankton(j,i)) = NOOBCSTATES
#endif
         end if
       end do
     end if
   end do

   do j = 1, iiPelDetritus
     if (.NOT.CalcPelDetritus(j)) then
       do i = 1,iiLastElement
         if ( ppPelDetritus(j,i) /= 0 ) then
           IF (InitD3STATE) THEN
             D3STATE(ppPelDetritus(j,i),:) = p_small
           END IF
           D3STATETYPE(ppPelDetritus(j,i)) = OFF
#if defined key_obcbfm
           D3STATEOBC(ppPelDetritus(j,i)) = NOOBCSTATES
#endif
         end if
       end do
     end if
   end do

   do j = 1, iiInorganic
     if (.NOT.CalcInorganic(j)) then
       do i = 1,iiLastElement
         if ( ppInorganic(j,i) /= 0 ) then
           IF (InitD3STATE) THEN
             D3STATE(ppInorganic(j,i),:) = p_small
           END IF
           D3STATETYPE(ppInorganic(j,i)) = OFF
#if defined key_obcbfm
           D3STATEOBC(ppInorganic(j,i)) = NOOBCSTATES
#endif
         end if
       end do
     end if
   end do



#ifdef INCLUDE_BEN
   !----------------------------------------------------
   ! Zeroing of the switched off benthic state variables
   !----------------------------------------------------



#endif


#ifdef INCLUDE_SEAICE
   !----------------------------------------------------
   ! Zeroing of the switched off seaice state variables
   !----------------------------------------------------



#endif


   !---------------------------------------------
   ! Write defined variables to stdout
   !---------------------------------------------
#ifdef BFM_PARALLEL
   Flun = LOGUNIT
#else
   Flun = stderr
#endif

   if (setup == 0 .OR. setup == 1 .OR. setup == 3 .OR. setup == 5 ) then
      LEVEL3 'Pelagic variables:'
      write(Flun,155) 'ID','Var','Unit','Long Name','Flag'
      do n=stPelStateS,stPelStateE
        write(Flun,156) n,trim(var_names(n)),trim(var_units(n)) &
          ,trim(var_long(n)),D3STATETYPE(n-stPelStateS+1)
      end do

   endif

#ifdef INCLUDE_BEN
   if (setup == 0 .OR. ( setup >= 2 .AND. setup <= 3 ) ) then
      LEVEL3 'Benthic variables:'
      write(Flun,155) 'ID','Var','Unit','Long Name','Flag'
      do n=stBenStateS,stBenStateE
        write(Flun,156) n,trim(var_names(n)),trim(var_units(n)) &
          ,trim(var_long(n)),D2STATETYPE_BEN(n-stBenStateS+1)
      end do

#ifdef INCLUDE_BENPROFILES
      !---------------------------------------------
      ! initialize the vertical grid for benthic
      ! nutrient profiles
      !---------------------------------------------
      LEVEL2 'Initialize the vertical grid for benthic profile diagnostics'
      LEVEL3 'Vertical sediment grid forced equal to model grid'
      p_sedlevels = NO_BOXES_Z_BEN
      call calc_sigma_depth(p_sedlevels,p_sedsigma,p_d_tot,seddepth)
      do n=1,p_sedlevels
         LEVEL3 n,seddepth(n)
      end do
#endif
   endif
#endif

#ifdef INCLUDE_SEAICE
   if (setup == 0 .OR. ( setup >= 4 .AND. setup <= 5 ) ) then
      LEVEL3 'Seaice variables:'
      write(Flun,155) 'ID','Var','Unit','Long Name','Flag'
      do n=stIceStateS,stIceStateE
        write(Flun,156) n,trim(var_names(n)),trim(var_units(n)) &
          ,trim(var_long(n)),D2STATETYPE_ICE(n-stIceStateS+1)
      end do

   endif
#endif

   return

155 FORMAT(10x,a4,1x,a5,1x,a12,1x,a40,1x,a10)
156 FORMAT(10x,i4,1x,a5,1x,a12,1x,a40,1x,i6)
   end subroutine init_var_bfm
!EOC

