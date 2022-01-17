!************************************************************************************************
       SUBROUTINE Output_Fields(IO,NI,NJ,X,Y,U,V,P)
         IMPLICIT NONE
 
         INTEGER NI,NJ,IO
         REAL,DIMENSION(NI,NJ):: X,Y
         REAL,DIMENSION(NI,NJ):: U,V,P

         Write(IO,*) 'VARIABLES = "X", "Y", "U", "V", "P"' 
         Write(IO,*) 'ZONE I=',NI,', J=',NJ, ', DATAPACKING=BLOCK'
         Write(IO,'(100E25.16)') X(1:NI,1:NJ) 
         Write(IO,'(100E25.16)') Y(1:NI,1:NJ)
         Write(IO,'(100E25.16)') U(1:NI,1:NJ)
         Write(IO,'(100E25.16)') V(1:NI,1:NJ)
         Write(IO,'(100E25.16)') P(1:NI,1:NJ)

       END  SUBROUTINE 
	   
	   SUBROUTINE Output_Cf(NI, NJ, X, Cf, Cf_B)
		integer NI, NJ, i, j
		REAL,DIMENSION(NI,NJ):: X
		
		REAL, DIMENSION(NI) :: Cf, Cf_B
		
		open (1, file = "Cf.plt")
		
		Write(1,*) 'VARIABLES = "Re_x", "Cf", "Cf_B"' 
		
		DO i=1,NI
          write(1, *), X(i, 1), Cf(i), Cf_B(i)
        ENDDO
		close(1)
	   end subroutine

!************************************************************************************************
