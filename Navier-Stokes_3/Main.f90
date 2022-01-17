        Program Pr
        Implicit none
 
        INTEGER, PARAMETER:: IO = 12 ! input-output unit
        INTEGER I,J,NI, NJ, NITER
        REAL(8) L,H,U0,MU,NU,R0,P0, h_V, U1, U2
        REAL(8) dx,dy,CFL,Uref,EPS
        REAL(8),ALLOCATABLE :: X_Node(:,:),Y_Node(:,:)
        REAL(8),ALLOCATABLE :: X_Cell(:,:),Y_Cell(:,:)
        REAL(8),ALLOCATABLE :: P_c(:,:),U_c(:,:),V_c(:,:), R_c(:,:)
 
        write(*,*) 'Read input file'
        open(IO,FILE='Input.txt')
        read(IO,*) L
        read(IO,*) H
		read(IO,*) h_V
        read(IO,*) NI
        read(IO,*) NJ
        read(IO,*) NITER
        read(IO,*) CFL
        read(IO,*) Uref
        read(IO,*) EPS
 
        read(IO,*) U1
		read(IO, *) U2
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
		allocate(R_c(0:NI,0:NJ))   ! Density

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
		U0 = (U1 + U2) / 2.0

        write(*,*)'L= ', L, 'NI= ', NI, 'dx= ', dx
        write(*,*)'H= ', H, 'NJ= ', NJ, 'dy= ', dy
        write(*,*)'ReL= ', U0*L/NU
		print*, "P0 = ", P0
		print*, "R0 = ", R0
        pause    

!----------------- Initial fields -----------------------------

        DO I=0,NI
          DO J=0,NJ
			if (Y_Cell(i,j) < h_V) then
				U_c(i, j) = U1
			else
				U_c(i, j) = U2
			endif
            V_c(I,J)=0.0
            P_c(I,J)=P0
			R_c(I,J)=R0
          ENDDO
        ENDDO
		

!---------------- Solve Prandtl equations ---------------------
 
        write(*,*) 'Solve Navier-Stokes equations' 
        Call NavierStokes(X_Node, Y_Node, U_c, V_c, P_c, R_c, dx, dy, U0, NI, NJ, MU, EPS, CFL, NITER, P0, R0, h_V, U1, U2, Y_Cell)
                       
 !----------------- Output data ------------------------------
 
        write(*,*) 'Output data cell (Navier-Stokes)' 
        Open(IO,FILE='data.plt')
        Call OutputFields_Cell(IO,NI,NJ,X_Node,Y_Node,U_c,V_c,P_c,R_c)
        Close(IO)
		
		deallocate (P_c, U_c, V_c, R_c, X_Node, Y_Node, X_Cell, Y_Cell)
        END PROGRAM
   
