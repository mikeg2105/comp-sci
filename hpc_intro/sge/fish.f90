        program fisheries
        real :: fstock , price , catch , rate , fstocklow
        real :: stocknew , pricenew , catchnew , ratenew , totcatch
        integer :: boats , time , boatsnew
        
! 
!  initial conditions.
!
           
        totcatch =0.0
        fstock = 10000.0
        fstocklow =5000.0  
        boats = 100        
!        rate = 0.01  
        rate =0.015 
        rate0 = rate
        catch = boats*2.0
        price = 1.0/catch 
        time = 0      
        write (*,101) 
101     format( ' Fisheries Simulation Program' ///)
        write (*,102) fstock
102     format(' Initial Stock (tons) : ', G12.2 )  
        write(*,104) 
104     format ( 'Enter enforcement introduction stock limit (tons)' )
        read(*,*) fstocklow  
        
10      time = time + 1
        
        stocknew = fstock*(1.0 + rate ) - catch
        pricenew = 1.0/catch
        if ( fstock .GE. fstocklow ) then
          boatsnew =  fstock/100 - price*10 - catch/100 
        else if ( boats.gt. 20 ) then
           boatsnew = boats/2.0 
        endif
        
        if (fstock .GE. fstocklow ) then
           catchnew = boats*2.0 
        else
           catchnew = boats*2.0*(1.0- (fstocklow-fstock)/fstocklow )
        endif
        
        if (fstock .GE.fstocklow ) then  
           ratenew =rate0 
        else 
           ratenew = rate0*(1.0 - (sqrt((fstocklow-fstock)**2) )/fstocklow )
        endif 
  
          
        
! 
!  update variables 
!        
        fstock =stocknew
        rate = ratenew
        boats = boatsnew
        catch = catchnew
        price =pricenew   
        totcatch = totcatch + catch
        write(*,*) time , fstock , boats ,catch  
               
        if (fstock.le.0 .or. boats.le.0 .OR.time .ge. 1000 .or. &
        & fstock .gt. 1.0E12 )then
           write (*,112) totcatch/time , time 
 112     format( ' Average catch :', F16.2 , 'tons  over ', I4 ,' months' ) 
         
           stop
        endif   
        go to 10
        
        end
        
        
        
        
        
        
        
        
        
        
