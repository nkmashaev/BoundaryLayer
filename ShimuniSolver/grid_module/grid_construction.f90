SUBMODULE (GRID_MODULE) GRID_CONSTRUCTION
	
	CONTAINS
	
		REAL(8) MODULE FUNCTION CREATE_DIM1(	NODE_ARR_XI, &
												CELL_ARR_XI, &
												START_XI, &
												END_XI, &
												NUMB_XI)
			IMPLICIT NONE
			REAL(8), DIMENSION(:), INTENT(OUT) 	:: 	NODE_ARR_XI, &
													CELL_ARR_XI
			INTEGER(4), INTENT(IN) 				:: 	NUMB_XI
			REAL(8), INTENT(IN) 				:: 	START_XI, &
													END_XI
			REAL(8)								::	DXI
			INTEGER(4)							::	I
			
			DXI = (END_XI - START_XI) / (NUMB_XI - 1)
			DO I = 1, NUMB_XI, 1
				NODE_ARR_XI(I) = START_XI + (I - 1) * DXI
				CELL_ARR_XI(I) = START_XI + (I - 1.5D0) * DXI
			ENDDO
			CELL_ARR_XI(NUMB_XI + 1) = END_XI + 0.5D0 * DXI
			
			CREATE_DIM1 = DXI
			
		END FUNCTION CREATE_DIM1

		REAL(8) MODULE FUNCTION CREATE_DIM2(	NODE_ARR_XI, &
												START_XI, &
												END_XI, &
												NUMB_XI)
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
			CREATE_DIM2 = DXI
			
		END FUNCTION CREATE_DIM2
		
END SUBMODULE GRID_CONSTRUCTION