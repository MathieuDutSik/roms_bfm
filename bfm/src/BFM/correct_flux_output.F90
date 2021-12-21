#include"cppdefs_bfm.h"
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE:routine to calculate sum of fluxes of D*SINK's and D*SOURCE's
!    into diagnostic variables defined as flux-variables in GlobalDefsBFM.model 
!
! !INTERFACE:
subroutine correct_flux_output(mode, nr0,zlev,nlev,out)
  !
  ! !DESCRIPTION:
  !     The information provided in the flux definition as given in GlobalDefsBFM.model
  !     is transferred into a number of array and scalars (flx_*). All
  !     elements of these arrays are defined in AllocateMem.F90. This
  !     information is used by this routines to calculates sums of fluxes
  !     from D3SINK, D2SINK, D3SOURCE, D2SOURCE, D3STATE and D2STATE ( See
  !     ModuleMem.F90 and AllocateMem.F90)
  !
  !     At input the value of zlev is 0 or 1. The value depend to which
  !      model BFM is coupled. The arrays in GOTM start at index 0. In most
  !     other coupled models at 1.
  !
  ! !USES:
  use constants, only: RLEN, ZERO, SEC_PER_DAY
  USE global_mem, ONLY: LOGUNIT

#ifdef NOPOINTERS
  use mem
#else

  use mem, only: NO_BOXES, PELBOTTOM,PELSURFACE,Depth, D3FLUX_FUNC, D3FLUX_MATRIX

#ifdef BFM_GOTM
  use bio_var, only: stPelStateS, stPelStateE
#else
  use api_bfm, only: stPelStateS, stPelStateE
#endif

#endif

#ifdef BFM_GOTM
  use bio_var, ONLY: SRFindices,BOTindices
#else
  use api_bfm, ONLY: SRFindices,BOTindices
#endif


  !
  ! !INPUT PARAMETERS:
  implicit none
  integer,intent(IN)                  ::mode
  integer,intent(IN)                  ::nr0
  integer,intent(IN)                  ::zlev
  integer,intent(IN)                  ::nlev
#ifdef BFM_GOTM
  real(RLEN),intent(OUT),dimension(0:nlev)  :: out
#else
  real(RLEN),intent(OUT),dimension(nlev)    :: out
#endif

  !
  ! !REVISION HISTORY:
  !  Original author(s): Piet Ruardij (NIOZ), Marcello Vichi (INGV)
  !

  !
  ! !LOCAL VARIABLES:
  integer      ::nr
  integer      ::i
  integer      ::k

  integer      ::idx_i, idx_j, origin, destination

#ifndef NOPOINTERS
  integer      ::idummy
#endif

#ifdef BFM_GOTM
  real(RLEN),dimension(0:NO_BOXES) :: hulp
#else
  real(RLEN),dimension(NO_BOXES)   :: hulp
#endif
  nr=nr0
  Print *, 'nr0=', nr0, ' RLEN=', RLEN, ' nlev=', nlev
  Print *, 'size(out)=', size(out)
  Print *, 'size(D3FLUX_FUNC)=', size(D3FLUX_FUNC,1), size(D3FLUX_FUNC,2)
  out(:) = D3FLUX_FUNC(nr0,:)
  Print *, '1 : max(out)=', maxval(out)
  do idx_i=stPelStateS, stPelStateE
       origin      = idx_i
       destination = idx_i
       Print *, 'allocated(D3FLUX_MATRIX(origin,destination)%p)=', allocated(D3FLUX_MATRIX(origin,destination)%p)
       if( allocated( D3FLUX_MATRIX(origin,destination)%p ) ) then
        do idx_j=1, SIZE(D3FLUX_MATRIX(origin,destination)%p)
    if( ABS(D3FLUX_MATRIX(origin,destination)%p(idx_j)) .eq. nr0 ) then
       if( D3FLUX_MATRIX(origin,destination)%dir(idx_j) == 0  ) then ! "A->B" => (out flow) => flux < ZERO => D3SINK
            out(BOTindices) = out(BOTindices) -                         &
     &    SIGN( 1, D3FLUX_MATRIX(origin,destination)%p(idx_j) ) *       &
     &     max(ZERO, PELBOTTOM(origin,:)) / Depth(BOTindices)
           out(SRFindices) = out(SRFindices) -                          &
     &    SIGN( 1, D3FLUX_MATRIX(origin,destination)%p(idx_j) ) *       &
     &     max(ZERO, PELSURFACE(origin,:)) / Depth(SRFindices)
       else ! "A<-B" => (in flow) => flux > ZERO => D3SOURCE
            out(BOTindices) = out(BOTindices) +                         &
     &    SIGN( 1, D3FLUX_MATRIX(origin,destination)%p(idx_j) ) *       &
     &     min(ZERO, PELBOTTOM(origin,:)) / Depth(BOTindices)
            out(SRFindices) = out(SRFindices) +                         &
     &    SIGN( 1, D3FLUX_MATRIX(origin,destination)%p(idx_j) ) *       &
     &     min(ZERO, PELSURFACE(origin,:)) / Depth(SRFindices)
              endif
           end if
        end do
     end if
  end do
  Print *, '2 : max(out)=', maxval(out)
end subroutine correct_flux_output
