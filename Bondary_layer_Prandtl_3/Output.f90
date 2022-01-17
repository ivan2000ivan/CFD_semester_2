!************************************************************************************************
       SUBROUTINE Output_Fields(IO,NI,NJ,X,Y,U,V,P,R)
         IMPLICIT NONE
 
         INTEGER NI,NJ,IO
         real(8),DIMENSION(NI,NJ):: X,Y
         real(8),DIMENSION(NI,NJ):: U,V,P,R

         Write(IO,*) 'VARIABLES = "X", "Y", "U", "V", "P", "R"' 
         Write(IO,*) 'ZONE I=',NI,', J=',NJ, ', DATAPACKING=BLOCK'
         Write(IO,'(100E25.16)') X(1:NI,1:NJ) 
         Write(IO,'(100E25.16)') Y(1:NI,1:NJ)
         Write(IO,'(100E25.16)') U(1:NI,1:NJ)
         Write(IO,'(100E25.16)') V(1:NI,1:NJ)
         Write(IO,'(100E25.16)') P(1:NI,1:NJ)
		 Write(IO,'(100E25.16)') R(1:NI,1:NJ)

       END  SUBROUTINE 
	   
	   SUBROUTINE Output_Gnuplot(NI, NJ, X, Cf, Cf_B)
		integer NI, NJ, i, j
		real(8),DIMENSION(NI,NJ):: X
		
		real(8), DIMENSION(NI) :: Cf, Cf_B
		
		open (1, file = "Cf.txt")
		DO i=1,NI
          write(1, *), X(i, 1), Cf(i)
        ENDDO
		close(1)
		
		open (1, file = "Cf_B.txt")
		DO i=1,NI
          write(1, *), X(i, 1), Cf_B(i)
        ENDDO
		close(1)
	   end subroutine

!************************************************************************************************
