!
! An example solution to exercise 2:
!
! Write out a program which asks the user to enter;
!       (1) His first name
!       (2) Year of his birth
!       (3) His surname
!     The program then goes on to say;
!
!     Hello <firstname> <surname> YOU ARE <nn> years old.  
!
!
       program myage
       character*12 :: name , surname 
       integer :: birth , age 
       
         print*, 'Please Enter your firstname:'
         read(*,'(A)' ) name
         
         print*,'Please enter your year of birth'
         read*,birth
         age = 2005-birth
         
         print*, 'Please Enter your surname:'
         read(*,'(A)' ) surname
       
         print *, ' Hello  ' , name , surname , ' you are ', age , ' years old'
         end
         
