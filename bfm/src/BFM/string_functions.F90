!$Id: util.F90,v 1.1 2005-06-27 10:54:33 kbk Exp $
#include"cppdefs_bfm.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: string_functions --- 
!
! !INTERFACE:
   MODULE string_functions
!
! !DESCRIPTION:
!   String functions to manupulate with arrays of strings.
! !USES:
   IMPLICIT NONE

   public empty, index_trim, getseq_number, replace_char
!
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: to test on empty string
!
! !INTERFACE:
     logical function empty(string)
!
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
     character(len=*)         ::string
!
! !REVISION HISTORY:
!  Original author(s): Piet Ruardij
!
!EOP
!
! !LOCAL VARIABLES:
!-------------------------------------------------------------------------
!BOC
     empty=lle(string,' ')
     return
     end function empty
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  index_trim
!
! !INTERFACE:
     integer function index_trim(string,text)
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
     character(len=*)             :: string*(*),text*(*)
!-------------------------------------------------------------------------
!BOC
     if (.not.empty(text)) then
        index_trim=index(string,text(1:len_trim(text)))
     else
        index_trim=0
     endif
     return
     end function index_trim
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  getseq_number
!
! !INTERFACE:
     integer function getseq_number(string,name,n,exact,start)
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
     implicit none
     integer,intent(IN)                       :: n
     character(len=*),dimension(n),intent(IN) :: name
     character(len=*),intent(IN)              :: string
     logical,intent(IN)                       :: exact
     integer,intent(IN),optional              :: start
!
!EOP
!
! !LOCAL VARIABLES:
     integer                 :: i,j,k
!-------------------------------------------------------------------------
!BOC
     k=1
     if ( present(start) ) k=start
     do i=k,n
       j=index_trim(name(i),string)
       if (j == 1) then
         if (exact ) THEN
           if (name(i) /= string) j=0
         endif
         if (j == 1) then
           getseq_number=i
           return
         endif
       endif
     enddo
     getseq_number=0
     return
     end function getseq_number
!EOC

!-----------------------------------------------------------------------
!EOC
!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE:  replace_char
!
! !INTERFACE:
     subroutine replace_char(str,tar,rep)
! !DESCRIPTION:
!
! !INPUT PARAMETERS:
     implicit none
    character(LEN=*), intent(INOUT) :: str
    character(LEN=*), intent(IN)    :: tar, rep
!
!EOP
!
! !LOCAL VARIABLES:
     integer                 :: times
!-------------------------------------------------------------------------
!BOC
     times = scan(str, tar)
     do while ( times .ne. 0 )
        str(times:times) = rep
        times = scan(str, tar)
     end do
   end subroutine replace_char
!EOC

!-----------------------------------------------------------------------

end module string_functions

!-----------------------------------------------------------------------
