       SUBROUTINE OutputFields_Cell(IO,NI,NJ,X,Y,U,V,P,R)
         IMPLICIT NONE

         INTEGER NI,NJ,IO
         REAL(8),DIMENSION(NI,NJ):: X,Y
         REAL(8),DIMENSION(0:NI,0:NJ)::U,V,P,R
       
         Write(IO,*) 'VARIABLES = "X", "Y", "U", "V", "P", "R"' 
         Write(IO,*) 'ZONE I=',NI,', J=',NJ,', DATAPACKING=BLOCK, VARLOCATION=([3-20]=CELLCENTERED)'
         Write(IO,'(100E25.16)') X(1:NI,1:NJ) 
         Write(IO,'(100E25.16)') Y(1:NI,1:NJ)
         Write(IO,'(100E25.16)') U(1:NI-1,1:NJ-1)
         Write(IO,'(100E25.16)') V(1:NI-1,1:NJ-1)
         Write(IO,'(100E25.16)') P(1:NI-1,1:NJ-1)
		 Write(IO,'(100E25.16)') R(1:NI-1,1:NJ-1)

       END SUBROUTINE 
