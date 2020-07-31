 
 (1) (a)Write a program to calculate the sum, the product and
        the quotient of two numbers and print it out.
        Experiment with what happens when a divide by zero is attempted.
     (b) If you have not already done so, write the program such that
         it asks the user to enter two numbers reads the numbers in
         and then prints out the results with some explanation
         such as 'the product is' ....
  (2) Write a program which asks the user to enter;
       (1) His first name
       (2) Year of his birth
       (3) His surname
       The program should then go on to say;

           Hello <firstname> <surname> YOU ARE <nn> years old.  



  (3) The area of a triangle is given by the formula:
    
       AREA  = square root of  ( S (S-A)(S-B)(S-C) )

        where;  S = 0.5 * ( A + B + C )

      Write a program which reads the values for the three
      sides from the keyboard and prints out the numeric value of the area.
      USE : READ and PRINT statements.

  (4) Write a program which prints out a conversion table from
      Centigrade to Fahrenheit,   from 0C to 21C inclusive.
       The formula for the conversion is:

            F = 32 + 1.8C
       USE : DO LOOPS

       (b) Modify this program so that it asks the user to enter
           the range of the table ( start temp. and stop temp. )
           and the increment C degrees. Such as in steps of 5C
           before printing out the table.
       
    

  (5) Taylors series expansion for sin(x) is given by the formula

          sinx = x - x^3/3!  + x^5/5! - x7/7! + ......
       where x is in radians.
         While programming you could simply use the expression
          SIN(X) which returns the correct value, by using the
          above formula without you being aware of it.
       As an exercise write a program which calculates sin(x).
       The problem is that there are infinite numbers of terms 
       in this series making it impossible to calculate sin(x) 
       without putting a limit on the number of terms. 
       One option may be to limit the series to say 6 terms.
       USE: DO LOOPS
       Another option is to note that the absolute value of 
       each consecutive term get smaller and smaller.
       Therefore we could decide to terminate the calculations
       once the absolute value of the next term in the series 
       becomes smaller than a pre defined limit (EPSILON).
       USE : DO WHILE LOOPS
       

  (6) Solution of the roots of a quadratic equation ;

        ax^2 + bx + c = 0 is given by the formula;

                -b +- sqrt( b^2 - 4ac )    
        x1,2 =  -----------------------
                          2a 

         Write a program to;
          - read in the three parameters a,b and c
          - calculate the roots if they exist and print them out
            otherwise print out that there are no real roots. 
        Run the program with different values of input a,b and c.             
        USE: IF STATEMENTS.

  
		 

  (7) The formula for the mean and the standard deviation of a sample of
      (n) data points (x) are given as ;

                  x1+x2+ .... xn
         xmean = -----------------
                       n

        
                           (x1-xmean)^2 + (x2-xmean)^2 + ...  (xn-xmean)^2 
         xstd =square root( ---------------------------------------------- )
                                          (n-1) 

      Read a fixed number of data values from the keyboard
      and calculate their mean and standard deviation and 
      print them out on the screen




  (8) Modify the program written for exercise (7) so that; 
     (a) The user is prompted for the data filename,
     (b) All the data values in that file ( one number per line )
         are read until the end of the file is reached,
     (c) A new file is created and the normalised data values
         are written into that file.
                                            (x-xmean)
         Use the formula; normalised(x) =  ----------
                                               xstd

  (9a) Example of subroutines and functions:
      Write the program (7) above in a modular way so that it has 
      the following parts: 
      - Main Program
      - An external procedure to read data from the terminal
      - An external procedure to calculate the meanvalue 
      - An external procedure to calculate the standard deviation
      The advantage of coding the program in this manner is that the
      your procedures can be called as many times as needed from many
      places as necessary with different parameters. This allows you to 
      create a library of routines for you to use in the same program
      or from other programs you may write in the future.
 
  (9b) Write the same program by using internal procedures.  
       Internal procedures are easier to write due to the scoping
       rules that makes all variables in the containing procedure
       (or main program) automatically accessible from within the 
       internal procedure. I recommend you use internal procedures
       to implement localised tasks that are not likely to be repeated
       elsewhere. 
       


 (10) Have a look at the example solutions directory to check 
     solutions listed there against your programs. 
     Your programs need not look like the example solutions. 
     Rarely two people come up with the same programming
     solution to the same problem. However, you will have 
     to judge your program against following criteria; 
     - Is your program easy to read.
     - Will it be easy to modify later for small alterations
       such as to increase the data size.
     - Is it efficient 'avoiding unnecessary or repeated calculations' !
     - Is it wasteful with data 'too many variables or large arrays 
       used where not necessary'.

 (11) We sometimes do not know the size of the data we are expected to
      handle until we start running the program. 
      Exercise 8 is a good example of such a problem, where we start 
      reading data from a file without any prior knowledge of how large 
      or small the data file is going to be.
      This is where allocatable arrays come handy.
      Ex8b.f95 is a modified solution to exercise (8) that uses an 
      allocatable array to hold the data values. 




