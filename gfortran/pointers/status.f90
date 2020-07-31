!------------------------------------------------------------------------------
!
!  Program name:
!
!    pointer_status
!
!  Purpose:
!
!    This program shows the status of pointers at different stages
!
!  Note:
!
!    The disassociation of p did not affect q even though they were 
!    both pointing at the same object. After being nullified, p can 
!    be associated again either with the same or a different object later. 
!    The last line just illustrates that a NULLIFY statement can 
!    have more than one pointer argument. 
!
!------------------------------------------------------------------------------

PROGRAM array_pointers

  IMPLICIT NONE 
  REAL, POINTER :: p, q 	! undefined association status
  REAL, TARGET :: t1 = 3.4, t2 = 4.5 

  p => t1			! p points to t1
  PRINT *, "After p => t1, ASSOCIATED(p) = ", ASSOCIATED(p)

  q => t1			! q also points to t1
  PRINT *, "After q => t1, ASSOCIATED(q) = ", ASSOCIATED(q)

  PRINT *, "ASSOCIATED(p, q) = ", ASSOCIATED(p, q)         
  PRINT *, "but ASSOCIATED(p, t2) = ", ASSOCIATED(p, t2)

  NULLIFY(p)			
  PRINT *, "After NULLIFY(p), ASSOCIATED(p) = ", ASSOCIATED(p)

  PRINT *, "ASSOCIATED(p, q) = ", ASSOCIATED(p, q)            
  ! The disassociation of p did not affect q even though 
  ! they were both pointing at the same object.
  PRINT *, "ASSOCIATED(q, t1) = ", ASSOCIATED(q, t1)          

  ! After being nullified, p can be associated again either 
  ! with the same or different object later
  p => t2			! p points to t2
  PRINT *, "After p => t2, ASSOCIATED(p) = ", ASSOCIATED(p) 	

  ! a NULLIFY statement can have more than one pointer argument
  NULLIFY(p, q)

END PROGRAM array_pointers
