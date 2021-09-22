SUBMODULE (TECPLOT_MODULE) TECPLOT_OUTPUT
	CONTAINS
	
		MODULE SUBROUTINE FV_OUTPUT1(IO, 		&
									NUMB_X,		&
									NUMB_Y,		&
									NODE_ARR_X,	&
									NODE_ARR_Y,	&
									P_ARR,		&
									U_MATR,		&
									V_MATR)
			INTEGER(4), INTENT(IN)				::	IO
			INTEGER(4), INTENT(IN)				::	NUMB_X, &
													NUMB_Y
			REAL(8), DIMENSION(:), INTENT(IN)	::	NODE_ARR_X, &
													NODE_ARR_Y, &
													P_ARR
			REAL(8), DIMENSION(:, :), INTENT(IN)::	U_MATR, &
													V_MATR
			INTEGER(4)							::	I, J
													
			WRITE(IO, *)	'VARIABLES = "X", "Y", "P", "U", "V"'
			WRITE(IO, *)	'ZONE I=', NUMB_Y, ', J=', NUMB_X, ', DATAPACKING=BLOCK, VARLOCATION=([3-6]=CELLCENTERED)'
			WRITE(IO, '(100F14.7)')	((NODE_ARR_X(J), I=1, NUMB_Y), J=1, NUMB_X)
			WRITE(IO, '(100F14.7)')	((NODE_ARR_Y(I), I=1, NUMB_Y), J=1, NUMB_X)
			WRITE(IO, '(100F14.7)') ((P_ARR(J), I=2, NUMB_Y), J=2, NUMB_X)
			WRITE(IO, '(100F14.7)') U_MATR(2:NUMB_Y, 2:NUMB_X)
			WRITE(IO, '(100F14.7)') V_MATR(2:NUMB_Y, 2:NUMB_X)
		END SUBROUTINE FV_OUTPUT1
		
		MODULE SUBROUTINE FV_OUTPUT2(IO, 	&
									NUMB_X,		&
									NUMB_Y,		&
									NODE_ARR_X,	&
									NODE_ARR_Y,	&
									P_MATR,		&
									U_MATR,		&
									V_MATR)
			INTEGER(4), INTENT(IN)				::	IO
			INTEGER(4), INTENT(IN)				::	NUMB_X, &
													NUMB_Y
			REAL(8), DIMENSION(:), INTENT(IN)	::	NODE_ARR_X, &
													NODE_ARR_Y
			REAL(8), DIMENSION(:, :), INTENT(IN)::	U_MATR, &
													V_MATR,	&
													P_MATR
			WRITE(IO, *)	'VARIABLES = "X", "Y", "P", "U", "V"'
			WRITE(IO, *)	'ZONE I=', NUMB_Y, ', J=', NUMB_X, ', DATAPACKING=BLOCK, VARLOCATION=([3-6]=CELLCENTERED)'
			WRITE(IO, '(100F14.7)')	((NODE_ARR_X(J), I=1, NUMB_Y), J=1, NUMB_X)
			WRITE(IO, '(100F14.7)')	((NODE_ARR_Y(I), I=1, NUMB_Y), J=1, NUMB_X)
			WRITE(IO, '(100F14.7)') P_MATR(2:NUMB_Y, 2:NUMB_X)
			WRITE(IO, '(100F14.7)') U_MATR(2:NUMB_Y, 2:NUMB_X)
			WRITE(IO, '(100F14.7)') V_MATR(2:NUMB_Y, 2:NUMB_X)
		END SUBROUTINE FV_OUTPUT2
		
END SUBMODULE TECPLOT_OUTPUT