! Exercise 4. Conversion table- from Celcius to Fahrenheid .
!
        program conv 
        integer ::  i 
        real    :: fahrenheit , celcius
                            
        print *, ' CELCIUS  :  FAHRENHEID' 
        
           do i = 0 , 21         
            celcius = i 
            fahrenheit = 32 + 1.8*celcius
            print * , celcius , fahrenheit
           enddo  
        
        end
        
 