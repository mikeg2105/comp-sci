program optimpath
!
! (c) Written by D.Savas on 12th of April 2012 
!  Email: D.Savas@sheffield.ac.uk 
!
! Calculates the optimal path thro' a connected graph with symmetric distances.
! Used as an example of solving the problem of going from LOCATION-X to LOCATION-Y
! using a predefined road network.
! The distance graph is symmetric in the sense that the distance from A to B 
! is the same as the distance from B to A. 
! The connection graph is expected to be an ASCII file with the following format.
!  FIRST THREE LINES : Captions ETC. ( Ignored by the program )
!  Next lines come in pairs. 
! If there are N-nodes in the graph there will be N pairs of lines.
!  LINE 1 :  Number of the Node  ,  Name of the node in 'quotes'  , Number of Connections to it.
!  LINE 2 : ( Integer - Real ) number pairs for each connection representing the connected node and distance.
!  LINE 3 : BLANK_LINE to make data more readable. 
! For example :
!  Line 1 :    1 , 'SMALLVILLE'  ,  3
!  Line 2 :    4 , 6.7   ,  3 , 9.8   ,  11 ,  2.6 
!  Line 3 : BLANK_LINE 
!  
!
! 
!
integer , parameter :: NC = 20   ! number of cities/nodes
real DISTMAT( NC,NC ) 
character*16 CNAME(NC) 
character*64 hop(NC,NC) 
real , parameter :: ALARGE= 1.0E30 
LOGICAL STATE_CHANGE
INTEGER IFAIL 

! Distance matrix shows the distance between each city.  Originally we assume no roads
! so the distance is very very large.

    DISTMAT= ALARGE
!!!	HOP = '*' 
! distance from city i to itself is 0 of course.
! also no hops to get to itself 
    DO i = 1 , NC
      DISTMAT(I,I) = 0.0 
      HOP(I,I) = ' '
    ENDDO 

    call readata( CNAME , DISTMAT , HOP , NC , IFAIL )
    IF ( IFAIL.NE.0 ) STOP
!
! CALCULATE THE DISTANCE MATRIX AND CONNECTIVITY ( IN HOP ARRAY ) 

10  STATE_CHANGE =.FALSE.
    do i = 1 , NC-1
      do j = i , nc
!!!        if(DISTMAT(I,J).EQ. ALARGE ) THEN 
! no direct connection between i and j . So search for indirect connection .
! however commented out because it is possible to have a one hop route that is shorter than direct route  
            DO K = 1 , NC 
              IF( DISTMAT(J,K) .NE. ALARGE .AND. DISTMAT(I,K) .NE. ALARGE ) THEN
                 DISTANCE = DISTMAT(J,K) + DISTMAT(I,K) 
                 IF( DISTANCE .LT. DISTMAT(I,J) ) THEN 
                   STATE_CHANGE = .TRUE.
                   DISTMAT(I,J) = DISTANCE   
                   DISTMAT(J,I) = DISTMAT(I,J) 
				    CALL ADDHOP( I, J , K, NC , HOP  ) 
                 ENDIF
              ENDIF 
            ENDDO 
!!!        endif 
          
       enddo
    enddo    
   
! CHECK TO SEE IF THERE HAS BEEN ANY IMPROVEMENT THE LAST TIME.
   if( STATE_CHANGE ) GO TO 10
 12  write(*,*) 'To finish enter  0  for start and stop cities.' 
    write ( *,*) ' Enter journey start city number followed by journey end city number' 
    read( *,* ) ISTART , IEND 
    If ( istart.eq. 0 .or. iend .eq. 0 ) stop

    DIST = DISTMAT( ISTART , IEND ) 
    WRITE(*,101 ) DIST 
 101 FORMAT( ' A ROUTE IS FOUND WITH DISTANCE: ', F14.3 ) 
    CALL PRINTR ( ISTART , IEND , NC , HOP ) 
          
    write(*,*) ' '
    go to 12
!   write(*,*) DISTMAT 
  
   end
   


   subroutine readata(  CNAME , DISTMAT , HOP , NC, IFAIL )
   integer      NC , IFAIL
   character*16 CNAME(NC) , CITYNAME
   character*64 hop(NC,NC) 
   real DISTMAT( NC,NC ), DIST(NC)
   integer numroads , neighb(NC)  
   IFAIL = 1 

    open(unit=12,file='road-data.txt', STATUS='OLD', IOSTAT = IFAIL )
    IF ( IFAIL .NE. 0 ) THEN 
      WRITE(*,*) 'ERROR: Can not open the road data file: road-data.txt.' 
      STOP

    ENDIF 
    read(12,199) cityname
    
 199 format( //A ) 
 ! Read all road distance information and update the distance matrix 
    do i=1, 20, 1
        read(12, *, ERR= 30 ) II  ,CITYNAME , NUMROADS 
        read(12, *, ERR=30 ) ( neighb(j) , DIST(j) , j=1,NUMROADS )  
        read(12,* )
        do j = 1 , NUMROADS  
          DISTMAT( II , NEIGHB(J))  = DIST(J)
		  HOP( II,NEIGHB(J) ) = ' ' 
		  
        enddo

        CNAME(II) = CITYNAME
    end do
    close(11)
    IFAIL= 0 
    RETURN
 30 IFAIL = 1
    RETURN
   end
   
   subroutine addhop( I,J, K , NC , HOP  )
   integer I,J,K , NC 
   character*64 HOPP , HOP(NC,NC)
! increments the hop array by adding an extra hop to current hops.
! the current hops are all BLANK for direct connections and '*' for unconnected ones.   
! For node I to node j connection add K into the hop.
!
   
    ICC = ICHAR('A') -1 
    LL= LEN_TRIM( HOP(I,K) )  
    LL2 = LEN_TRIM( HOP(J,K) )
 
    HOPP = HOP(I,K)(1:LL)//CHAR(ICC+K)// HOP(J,K)(1:LL2) 
	 
    HOP(I,J) = HOPP 
    HOP(J,I) = HOPP
	RETURN
    END


    SUBROUTINE  PRINTR ( ISTART , IEND , NC, HOP ) 
    INTEGER ISTART , IEND , NC , IS , IE , LL
    CHARACTER*64 HOP(NC,NC) , HOPP
        IF ( ISTART.GT. IEND ) THEN
         IS = IEND 
         IE = ISTART
    ELSE
         IS= ISTART
         IE = IEND 
    ENDIF 
    WRITE(*, '(A)' , ADVANCE='NO' ) 'ROUTE IS: ' 
    WRITE(*, '(I3)' , ADVANCE='NO' ) ISTART 
    LL= LEN_TRIM( HOP(IS,IE) )   
    HOPP = HOP(IS,IE)(1:LL) 
    DO K = 1 , LL
     KK = K 
     IF( ISTART.GT.IEND ) KK = LL+1-K   
     ICC = ICHAR( HOPP(KK:KK) ) - ICHAR('A' ) + 1
     WRITE(*, '(I3)' , ADVANCE='NO' ) ICC 
   ENDDO 
      WRITE(*, '(I3)' , ADVANCE='NO' ) IEND
    RETURN
    END