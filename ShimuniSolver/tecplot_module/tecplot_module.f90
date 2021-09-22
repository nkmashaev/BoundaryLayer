MODULE TECPLOT_MODULE
	IMPLICIT NONE
	
	INTERFACE
		
		MODULE SUBROUTINE FV_OUTPUT1(IO, 	&
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
		END SUBROUTINE FV_OUTPUT2
		
	END INTERFACE
	
	INTERFACE FV_OUTPUT
		MODULE PROCEDURE FV_OUTPUT1, FV_OUTPUT2
	END INTERFACE FV_OUTPUT

END MODULE TECPLOT_MODULE