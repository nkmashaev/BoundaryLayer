PROGRAM NAVIERSTOKES
	IMPLICIT NONE
	CHARACTER(*), PARAMETER 				::	INPUT_FILE = "input.txt", &
												OUTPUT_FILE = "output.plt", &
												RESIDUALS_FILE = "residuals.dat"
	INTEGER(4), PARAMETER					::	COMMON_IO = 101
	REAL(8), DIMENSION(:), ALLOCATABLE		::	NODE_ARR_X, &
												NODE_ARR_Y
	REAL(8), DIMENSION(:, :), ALLOCATABLE	::	P_MATR, &
												U_MATR, &
												V_MATR
	REAL(8)									::	START_X, &
												END_X, &
												START_Y, &
												END_Y, &
												DX, &
												DY, &
												MU, &
												U0, &
												V0, &
												P0, &
												CFL, &
												CF_ANALYT, &
												CF_THEOR, &
												RE_X
	INTEGER(4)								::	NUMB_X, &
												NUMB_Y, &
												NUMB_I, &
												NUMB_J, &
												I, &
												J, &
												MAX_ITER

		OPEN(COMMON_IO, FILE=INPUT_FILE)
		CALL READ_INPUT_FILE(COMMON_IO, START_X, END_X, NUMB_X, &
						START_Y, END_Y, NUMB_Y, U0, V0, &
						P0, MU, CFL, MAX_ITER)
		CLOSE(COMMON_IO)
		NUMB_I = NUMB_Y + 1
		NUMB_J = NUMB_X + 1
		
		ALLOCATE(NODE_ARR_X(NUMB_X))
		ALLOCATE(NODE_ARR_Y(NUMB_Y))
		ALLOCATE(U_MATR(NUMB_I, NUMB_J))
		ALLOCATE(V_MATR(NUMB_I, NUMB_J))
		ALLOCATE(P_MATR(NUMB_I, NUMB_J))
		DX = CREATE_DIM(NODE_ARR_X, START_X, END_X, NUMB_X)
		DY = CREATE_DIM(NODE_ARR_Y, START_Y, END_Y, NUMB_Y)
		DO J = 1, NUMB_J, 1
			DO I = 1, NUMB_I, 1
				P_MATR(I, J) = P0
				V_MATR(I, J) = V0
				U_MATR(I, J) = U0
			ENDDO
		ENDDO
		
		OPEN(COMMON_IO, FILE=RESIDUALS_FILE)
		CALL NS_SOLVER(	COMMON_IO, NUMB_I, NUMB_J, MAX_ITER, &
						DX, DY, MU, U0, V0, P0, CFL, &
						P_MATR, U_MATR, V_MATR)			
		CLOSE(COMMON_IO)
		
		OPEN(COMMON_IO, FILE=OUTPUT_FILE)
		CALL WRITE_TO_TECPLOT(COMMON_IO, NUMB_X, NUMB_Y, &
						NODE_ARR_X, NODE_ARR_Y, P_MATR, &
						U_MATR, V_MATR)
		CLOSE(COMMON_IO)
		
		OPEN(COMMON_IO, FILE='cf.dat')
		DO J=2, NUMB_X, 1
			RE_X = U0 * (J - 0.5D0) * DX / MU
			CF_THEOR = 0.664D0 / SQRT(RE_X)
			CF_ANALYT = 4.0D0 * MU * U_MATR(2, J) / (DY * U0**2)
			WRITE(COMMON_IO, *) RE_X , CF_THEOR, CF_ANALYT
		ENDDO
		CLOSE(COMMON_IO)
		
		DEALLOCATE(P_MATR)
		DEALLOCATE(V_MATR)
		DEALLOCATE(U_MATR)
		DEALLOCATE(NODE_ARR_Y)
		DEALLOCATE(NODE_ARR_X)
					
	CONTAINS
		SUBROUTINE READ_INPUT_FILE(	IO, &
									START_X, &
									END_X, &
									NUMB_X, &
									START_Y, &
									END_Y, &
									NUMB_Y, &
									U0, &
									V0, &
									P0, &
									MU, &
									CFL, &
									MAX_ITER)
									
		!READ THE INPUT FILE WITH INITIALIZATION PARAMETERS
		!IO - THE BUFFER ID
		!START_X, END_X, NUMB_X - A DEFINITION OF X-DIMENSION
		!START_Y, END_Y, NUMB_Y - A DEFINITION OF Y-DIMENSION
		!U0, V0, P0	- INPUT PARAMETERS
		!DENSITY - CONSTANT DENSITY OF FLUID
		!MU - CONSTANT DYNAMIC VISCOSITY COEFFICIENT
		!CFL - COURANT NUMBER
		!MAX_ITER - MAX NUMBER OF ITERATION 
		
			IMPLICIT NONE
			INTEGER(4), INTENT(IN)	::	IO
			INTEGER(4), INTENT(OUT)	::	NUMB_X, &
										NUMB_Y, &
										MAX_ITER
			REAL(8), INTENT(OUT)	::	START_X, &
										END_X, &
										START_Y, &
										END_Y, &
										U0, &
										V0, &
										P0, &
										MU, &
										CFL
			READ(IO, *)	START_X, END_X, NUMB_X
			READ(IO, *) START_Y, END_Y, NUMB_Y
			READ(IO, *) U0, V0
			READ(IO, *) P0
			READ(IO, *) MU
			READ(IO, *) CFL
			READ(IO, *) MAX_ITER
		
		END SUBROUTINE READ_INPUT_FILE
		
		
		
		REAL(8)  FUNCTION CREATE_DIM(	NODE_ARR_XI, &
										START_XI, &
										END_XI, &
										NUMB_XI)
		! NODE_ARR_XI - CREATE XI-COORDINATES
		! START_XI, END_XI - SECTION DEFINITION
		! NUMB_XI - NUMBER OF NODES
			IMPLICIT NONE
			REAL(8), DIMENSION(:), INTENT(OUT) 	:: 	NODE_ARR_XI
			INTEGER(4), INTENT(IN) 				:: 	NUMB_XI
			REAL(8), INTENT(IN) 				:: 	START_XI, &
													END_XI
			REAL(8)								::	DXI
			INTEGER(4)							::	I
			
			DXI = (END_XI - START_XI) / (NUMB_XI - 1)
			DO I = 1, NUMB_XI, 1
				NODE_ARR_XI(I) = START_XI + (I - 1) * DXI
			ENDDO
			CREATE_DIM = DXI
			
		END FUNCTION CREATE_DIM
		
		
		
		SUBROUTINE NS_SOLVER(	IO, &
									NUMB_I, &
									NUMB_J, &
									MAX_ITER, &
									DX, &
									DY, &
									MU, &
									U0, &
									V0, &
									P0, &
									CFL, &
									P_MATR, &
									U_MATR, &
									V_MATR)		
		! SOLVE A SYSTEM OF NAVIER STOKES EQUATIONS
		! IO - THE BUFFER ID
		! NUMB_I, NUMB_J - DIMENSION OF SYSTEM
		! MAX_ITER - ITERATION LIMITER
		! DX, DY - COORDIANTES STEPS
		! MU, - DYNAMIC VISCOSITY COEFFICIENT
		! U0, V0 - INPUT COMPONENTS OF VELOCITY
		! P0 - INPUT PRESSURE
		! CFL - COURANT NUMBER
		! P_MATR - PRESSURE
		! U_MATR - X-VELOCITY
		! V_MATR - Y_VELOCITY
		
			IMPLICIT NONE
			INTEGER(4), INTENT(IN)						::	IO, &
															NUMB_I, &
															NUMB_J, &
															MAX_ITER
			REAL(8), DIMENSION(:, :), INTENT(IN OUT)	::	U_MATR, &
															V_MATR, &
															P_MATR
			REAL(8), INTENT(IN)							::	DX, &
															DY, &
															MU, &
															U0, &
															V0, &
															P0, &
															CFL
			INTEGER(4)									::	I, &
															J, &
															N
			REAL(8)										::	UIP, &
															VIP, &
															UIM, &
															VIM, &
															UJP, &
															VJP, &
															UJM, &
															VJM, &
															UIPF, &
															VIPF, &
															UIMF, &
															VIMF, &
															UJPF, &
															VJPF, &
															UJMF, &
															VJMF, &
															PIPF, &
															PIMF, &
															PJPF, &
															PJMF, &
															CONV1, &
															CONV2, &
															PRS, &
															DIFF1, &
															DIFF2, &
															RES_C_MAX, &
															RES_X_MAX, &
															RES_Y_MAX, &
															MU_X, &
															MU_Y, &
															DT, &
															V_NORM, & 
															RES_C, &
															RES_X, & 
															RES_Y, &
															A
			
			V_NORM = SQRT(U0 * U0 + V0 * V0)
			A = V_NORM * V_NORM
			DT = CFL * (DX * DY) / (2.0D0 * (DX + DY) * V_NORM)
			MU_X = MU / (DX * DX)
			MU_Y = MU / (DY * DY)
			
			DO N = 1, MAX_ITER, 1
			
				RES_C_MAX = 0.0D0
				RES_X_MAX = 0.0D0
				RES_Y_MAX = 0.0D0
				DO J = 2, NUMB_J - 1, 1
					DO I = 2, NUMB_I - 1, 1
						UIP = 0.5D0 * (U_MATR(I, J) + U_MATR(I + 1, J))
						VIP = 0.5D0 * (V_MATR(I, J) + V_MATR(I + 1, J))
						UIM = 0.5D0 * (U_MATR(I, J) + U_MATR(I - 1, J))
						VIM = 0.5D0 * (V_MATR(I, J) + V_MATR(I - 1, J))
						UJP = 0.5D0 * (U_MATR(I, J) + U_MATR(I, J + 1))
						VJP = 0.5D0 * (V_MATR(I, J) + V_MATR(I, J + 1))
						UJM = 0.5D0 * (U_MATR(I, J) + U_MATR(I, J - 1))
						VJM = 0.5D0 * (V_MATR(I, J) + V_MATR(I, J - 1))
						
						IF (UJP >= 0.0D0) THEN
							UJPF = U_MATR(I, J)
							VJPF = V_MATR(I, J)
							PJPF = P_MATR(I, J + 1)
						ELSE
							UJPF = U_MATR(I, J + 1)
							VJPF = V_MATR(I, J + 1)
							PJPF = P_MATR(I, J)
						ENDIF
						
						IF (UJM >= 0.0D0) THEN
							UJMF = U_MATR(I, J - 1)
							VJMF = V_MATR(I, J - 1)
							PJMF = P_MATR(I, J)
						ELSE
							UJMF = U_MATR(I, J)
							VJMF = V_MATR(I, J)
							PJMF = P_MATR(I, J - 1)
						ENDIF
						
						IF (VIP >= 0.0D0) THEN
							UIPF = U_MATR(I, J)
							VIPF = V_MATR(I, J)
							PIPF = P_MATR(I + 1, J)
						ELSE
							UIPF = U_MATR(I + 1, J)
							VIPF = V_MATR(I + 1, J)
							PIPF = P_MATR(I, J)
						ENDIF
						
						IF (VIM >= 0.0D0) THEN
							UIMF = U_MATR(I - 1, J)
							VIMF = V_MATR(I - 1, J)
							PIMF = P_MATR(I, J)
						ELSE
							UIMF = U_MATR(I, J)
							VIMF = V_MATR(I, J)
							PIMF = P_MATR(I - 1, J)
						ENDIF
						
						CONV1 = (UJPF - UJMF) / DX
						IF (I == 2) THEN
							CONV2 = VIPF / DY
						ELSE IF (I == NUMB_I - 1) THEN
							CONV2 = -VIMF / DY
						ELSE
							CONV2 = (VIPF - VIMF) / DY
						ENDIF
						RES_C = CONV1 + CONV2
						RES_C_MAX = MAX(RES_C_MAX, ABS(RES_C))
						P_MATR(I, J) = P_MATR(I, J) - DT * A * RES_C
						
						CONV1 = (UJPF * UJP - UJMF * UJM) / DX
						CONV2 = (UIPF * VIP - UIMF * VIM) / DY
						PRS = (PJPF - PJMF) / DX
						DIFF1 = MU_X * (U_MATR(I, J + 1) - 2.0D0 * U_MATR(I, J) + U_MATR(I, J - 1))
						DIFF2 = MU_Y * (U_MATR(I + 1, J) - 2.0D0 * U_MATR(I, J) + U_MATR(I - 1, J))
						RES_X = CONV1 + CONV2 + PRS - DIFF1 - DIFF2
						RES_X_MAX = MAX(RES_X_MAX, ABS(RES_X))
						U_MATR(I, J) = U_MATR(I, J) - DT * RES_X					
						
						CONV1 = (VJPF * UJP - VJMF * UJM) / DX
						CONV2 = (VIPF * VIP - VIMF * VIM) / DY
						PRS = (PIPF - PIMF) / DY
						DIFF1 = MU_X * (V_MATR(I, J + 1) - 2.0D0 * V_MATR(I, J) + V_MATR(I, J - 1))
						DIFF2 = MU_Y * (V_MATR(I + 1, J) - 2.0D0 * V_MATR(I, J) + V_MATR(I - 1, J))
						RES_Y = CONV1 + CONV2 + PRS - DIFF1 - DIFF2
						RES_Y_MAX = MAX(RES_Y_MAX, ABS(RES_Y))
						V_MATR(I, J) = V_MATR(I, J) - DT * RES_Y
					ENDDO
				ENDDO

				DO I = 2, NUMB_I - 1, 1
					U_MATR(I, 1) = 2.0D0 * U0 - U_MATR(I, 2)
					V_MATR(I, 1) = 2.0D0 * V0 - V_MATR(I, 2)
					P_MATR(I, 1) = P_MATR(I, 2)
					
					U_MATR(I, NUMB_J) = U_MATR(I, NUMB_J - 1)
					V_MATR(I, NUMB_J) = V_MATR(I, NUMB_J - 1)
					P_MATR(I, NUMB_J) = 2.0D0 * P0 - P_MATR(I, NUMB_J - 1)
				ENDDO
				
				DO J = 1, NUMB_J, 1
					U_MATR(1, J) = -U_MATR(2, J)
					V_MATR(1, J) = -V_MATR(2, J)
					P_MATR(1, J) = P_MATR(2, J)
					
					!SYMMETRY
					U_MATR(NUMB_I, J) = U_MATR(NUMB_I - 1, J)
					V_MATR(NUMB_I, J) = V_MATR(NUMB_I - 1, J)
					P_MATR(NUMB_I, J) = P_MATR(NUMB_I - 1, J)
					!OUTER FLOW CONDITION
					!U_MATR(NUMB_I, J) = 2.0D0 * U0 - U_MATR(NUMB_I - 1, J)
					!V_MATR(NUMB_I, J) = 2.0D0 * V0 - V_MATR(NUMB_I - 1, J)
					!P_MATR(NUMB_I, J) = P_MATR(NUMB_I - 1, J)
				ENDDO
				
				WRITE(*, *) N, RES_C_MAX, RES_X_MAX, RES_Y_MAX
				WRITE(IO, *) N, RES_C_MAX, RES_X_MAX, RES_Y_MAX
			ENDDO
			
		END SUBROUTINE NS_SOLVER
		
		
		
		SUBROUTINE WRITE_TO_TECPLOT(IO, 		&
									NUMB_X,		&
									NUMB_Y,		&
									NODE_ARR_X,	&
									NODE_ARR_Y,	&
									P_MATR,		&
									U_MATR,		&
									V_MATR)
		! OUTPUT DATA TO TECPLOT IN A FINITE VOLUME FORMAT
		! IO - THE BUFFER ID
		! NUMB_X, NUMB_Y - NUMBER OF NODES (NUMB_Y x NUMB_X)
		! NODE_ARR_X - X-COORDINATES OF NODES
		! NODE_ARR_Y - Y-COORDINATES OF NODES
		! P_MATR - PRESSURE
		! U_MATR - X-COMPONENT OF THE VELOCITY VECTOR
		! V_MATR - Y-COMPONENT OF THE VELOCITY VECTOR
			INTEGER(4), INTENT(IN)				::	IO
			INTEGER(4), INTENT(IN)				::	NUMB_X, &
													NUMB_Y
			REAL(8), DIMENSION(:), INTENT(IN)	::	NODE_ARR_X, &
													NODE_ARR_Y
			REAL(8), DIMENSION(:, :), INTENT(IN)::	U_MATR, &
													V_MATR, &
													P_MATR
			INTEGER(4)							::	I, J
													
			WRITE(IO, *)	'VARIABLES = "X", "Y", "P", "U", "V"'
			WRITE(IO, *)	'ZONE I=', NUMB_Y, ', J=', NUMB_X, ', DATAPACKING=BLOCK, VARLOCATION=([3-5]=CELLCENTERED)'
			WRITE(IO, '(100F14.7)')	((NODE_ARR_X(J), I=1, NUMB_Y), J=1, NUMB_X)
			WRITE(IO, '(100F14.7)')	((NODE_ARR_Y(I), I=1, NUMB_Y), J=1, NUMB_X)
			WRITE(IO, '(100F14.7)') P_MATR(2:NUMB_Y, 2:NUMB_X)
			WRITE(IO, '(100F14.7)') U_MATR(2:NUMB_Y, 2:NUMB_X)
			WRITE(IO, '(100F14.7)') V_MATR(2:NUMB_Y, 2:NUMB_X)
		END SUBROUTINE WRITE_TO_TECPLOT
END PROGRAM NAVIERSTOKES