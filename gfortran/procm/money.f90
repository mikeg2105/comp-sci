!-------------------------------------------------------------------------
!
!  File name:
!
!    money.f90
!
!  Purpose:
!
!    This file contains the code fragements, from which a complete module 
!    money_module can be constructed that allows you to run the program 
!    money_main.
!
!-------------------------------------------------------------------------

! Module:

MODULE moneytype
IMPLICIT NONE

  TYPE money
    INTEGER :: pounds, pence
  END TYPE money

  INTERFACE OPERATOR (+)
    MODULE PROCEDURE addmoney
  END INTERFACE

  INTERFACE OPERATOR (-)
    MODULE PROCEDURE negatemoney, subtractmoney
  END INTERFACE

  CONTAINS

  FUNCTION addmoney(a,b)
    TYPE (money) :: addmoney
    TYPE (money), INTENT(IN) :: a,b
    INTEGER :: carry, temppence
    temppence = a%pence + b%pence
    carry = 0
      IF (temppence>100) THEN
        temppence = temppence - 100
        carry = 1
      END IF
    addmoney%pounds = a%pounds + b%pounds + carry
    addmoney%pence = temppence
  END FUNCTION addmoney

  FUNCTION negatemoney(a)
    TYPE (money) :: negatemoney
    TYPE (money), INTENT(IN) :: a
    negatemoney%pounds = -a%pounds
    negatemoney%pence = -a%pence
  END FUNCTION negatemoney

  FUNCTION subtractmoney(a,b)
    TYPE (money) :: subtractmoney
    TYPE (money), INTENT(IN) :: a,b
    INTEGER :: temppound, temppence, carry
    temppence = a%pence - b%pence
    temppound = a%pounds - b%pounds
      IF ((temppence<0).AND.(temppound>0)) THEN
        temppence = 100 + temppence
        temppound = temppound - 1
      ELSE IF ((temppence>0).AND.(temppound<0)) THEN
        temppence = temppence - 100
        temppound = temppound + 1
      END IF
    subtractmoney%pence = temppence
    subtractmoney%pounds = temppound
  END FUNCTION subtractmoney

END MODULE moneytype

PROGRAM test
USE moneytype
IMPLICIT NONE
  
  TYPE (money) :: pay, spend
  TYPE (money) :: add, neg, sub
  pay=money(2,36)
  spend=money(1,35)


  add=pay + spend
  neg= - pay
  sub=pay-spend

  WRITE(*,*) add, neg, sub

END PROGRAM test


    
   




  