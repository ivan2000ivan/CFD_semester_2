        Program Pr
        Implicit none
 
        INTEGER, PARAMETER:: IO = 12 ! input-output unit
        INTEGER I,J,NI, NJ, NITER
        REAL(8) L,H,U0,MU,NU,R0,P0
        REAL(8) dx,dy,CFL,Uref,EPS, Umax
        REAL(8),ALLOCATABLE :: X_Node(:,:),Y_Node(:,:)
        REAL(8),ALLOCATABLE :: X_Cell(:,:),Y_Cell(:,:)
        REAL(8),ALLOCATABLE :: P_c(:,:),U_c(:,:),V_c(:,:), U(:,:), V(:,:), P(:,:), Cf(:), Cf_B(:), U_n(:,:), V_n(:,:), P_n(:,:)
 
        write(*,*) 'Read input file'
        open(IO,FILE='Input.txt')
        read(IO,*) L
        read(IO,*) H
        read(IO,*) NI
        read(IO,*) NJ
        read(IO,*) NITER
        read(IO,*) CFL
        read(IO,*) Uref
        read(IO,*) EPS
 
        read(IO,*) U0
        read(IO,*) MU
        read(IO,*) R0
        read(IO,*) P0
        CLOSE(IO)

        allocate(X_Node(NI,NJ)) ! mesh nodes X-coordinates
        allocate(Y_Node(NI,NJ)) ! mesh nodes Y-coordinates
        allocate(X_Cell(0:NI,0:NJ)) ! cell centers X-coordinates
        allocate(Y_Cell(0:NI,0:NJ)) ! cell centers Y-coordinates

!------ Cell-centered variables
        allocate(U_c(0:NI,0:NJ))   ! Velocity U
        allocate(V_c(0:NI,0:NJ))   ! Velocity V
        allocate(P_c(0:NI,0:NJ))   ! Pressure
		
		allocate(U(NI,NJ))  
        allocate(V(NI,NJ))   
        allocate(P(NI,NJ))
		
		allocate(U_n(1:NI,1:NJ))  
        allocate(V_n(1:NI,1:NJ))   
        allocate(P_n(1:NI,1:NJ))
		
		allocate(Cf(NI), Cf_B(NI))

        dx=L/(NI-1)
        dy=H/(NJ-1)

!------ Coordinate of nodes
        DO I=1,NI
          DO J=1,NJ
            X_Node(I,J)=(I-1)*dx
            Y_Node(I,J)=(J-1)*dy
          END DO
        END DO

!------ Coordinate of cell centers
        X_Cell(0,1:NJ)=-dx/2
        Y_Cell(0,1:NJ)=Y_Node(1,1:NJ)+dy/2
        X_Cell(1:NI,0)=X_Node(1:NI,1)+dx/2
        Y_Cell(1:NI,0)=-dy/2
        DO I=1,NI
          DO J=1,NJ
            X_Cell(I,J)=X_Node(I,J)+dx/2
            Y_Cell(I,J)=Y_Node(I,J)+dy/2
          END DO
        END DO

!----------------- Parameters ------------------------

        NU=MU/R0

        write(*,*)'L= ', L, 'NI= ', NI, 'dx= ', dx
        write(*,*)'H= ', H, 'NJ= ', NJ, 'dy= ', dy
        write(*,*)'ReL= ', U0*L/NU
        pause    

!----------------- Initial fields -----------------------------

        DO I=0,NI
          DO J=0,NJ
            U_c(I,J)=U0
            V_c(I,J)=1.0e-5
            P_c(I,J)=0.0
          ENDDO
        ENDDO
		

!---------------- Solve Prandtl equations ---------------------
 
        write(*,*) 'Solve Navier-Stokes equations' 
        Call NavierStokes(X_Node, Y_Node, U_c, V_c, P_c, dx, dy, U0, NI, NJ, NU, EPS, CFL, NITER)
		
		DO I=1,NI
          DO J=1,NJ
            U(I, J) = (U_c(I, J) + U_c(I, J-1)) / 2.0
			V(I, J) = (V_c(I, J) + V_c(I, J-1)) / 2.0
			P(I, J) = (P_c(I, J) + P_c(I, J-1)) / 2.0
          ENDDO
        ENDDO
		
		DO I=1,NI
          DO J=1,NJ
			U_n(I, J) = (U_c(I-1, J-1) + U_c(I-1, J) + U_c(I, J-1) + U_c(I, J)) / 4.0
			V_n(I, J) = (V_c(I-1, J-1) + V_c(I-1, J) + V_c(I, J-1) + V_c(I, J)) / 4.0
			P_n(I, J) = (P_c(I-1, J-1) + P_c(I-1, J) + P_c(I, J-1) + P_c(I, J)) / 4.0
			!print*, U_n(I,J)
          ENDDO
		  !print*, ""
        ENDDO
		
		do i = 1, NI
			Cf(i) = -MU * (3 * U_n(i, 1) - 4 * U_n(i, 2) + U_n(i, 3)) / (2 * dy)
			Umax = ((U_n(i, 0))**2 + (V_n(i, 0))**2)**0.5
			do j = 1, NJ
				Umax = max(Umax, ((U_n(i, j))**2 + (V_n(i, j))**2)**0.5)
			enddo
			Cf(i) = 2 * Cf(i) / (R0 * (Umax)**2)
			Cf_B(I) = 0.664 / sqrt(U0 * X_Node(I, 1) / NU)
			!print*, Umax
		enddo
                       
 !----------------- Output data ------------------------------
 
        write(*,*) 'Output data cell (Navier-Stokes)' 
        Open(IO,FILE='data.plt')
        !Call OutputFields_Cell(IO,NI,NJ,X_Node,Y_Node,U_c,V_c,P_c)
		Call OutputFields_Node(IO,NI,NJ,X_Node,Y_Node,U_n,V_n,P_n)
        Close(IO)
		
		open (1, file = "Cf.plt")
		
		Write(1,*) 'VARIABLES = "Re_x", "Cf", "Cf_B"' 
		
		DO i=1,NI
          write(1, *), X_Node(i, 1), Cf(i), Cf_B(i)
        ENDDO
		close(1)
     
        END PROGRAM
   
