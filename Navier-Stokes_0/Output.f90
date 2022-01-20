       SUBROUTINE OutputFields_Cell(IO,NI,NJ,X,Y,U,V,P)
         IMPLICIT NONE

         INTEGER NI,NJ,IO
         REAL(8),DIMENSION(NI,NJ):: X,Y
         REAL(8),DIMENSION(0:NI,0:NJ)::U,V,P
       
         Write(IO,*) 'VARIABLES = "X", "Y", "U", "V", "P"' 
         Write(IO,*) 'ZONE I=',NI,', J=',NJ,', DATAPACKING=BLOCK, VARLOCATION=([3-20]=CELLCENTERED)'
         Write(IO,'(100E25.16)') X(1:NI,1:NJ) 
         Write(IO,'(100E25.16)') Y(1:NI,1:NJ)
         Write(IO,'(100E25.16)') U(1:NI-1,1:NJ-1)
         Write(IO,'(100E25.16)') V(1:NI-1,1:NJ-1)
         Write(IO,'(100E25.16)') P(1:NI-1,1:NJ-1)

       END SUBROUTINE
	   
	   SUBROUTINE OutputFields_Node(IO,NI,NJ,X,Y,U,V,P)
         IMPLICIT NONE
		 
         INTEGER NI,NJ,IO
         REAL(8),DIMENSION(1:NI,1:NJ):: X,Y
         REAL(8),DIMENSION(1:NI,1:NJ):: U,V,P
		 
		 !print*, U

         Write(IO,*) 'VARIABLES = "X", "Y", "U", "V", "P"' 
         Write(IO,*) 'ZONE I=',NI,', J=',NJ, ', DATAPACKING=BLOCK'
         Write(IO,'(100E25.16)') X 
         Write(IO,'(100E25.16)') Y
         Write(IO,'(100E25.16)') U
         Write(IO,'(100E25.16)') V
         Write(IO,'(100E25.16)') P

       END SUBROUTINE

