MODULE load_data
!
USE constants
USE observations
!
IMPLICIT NONE
!
CONTAINS
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE input_model ( sel_X )
	!
    ! Reads the model's specification from text file
    !
	IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(OUT) :: sel_X(num_tot_X,num_states)
	! 
	! Declaring local variables
	!
	INTEGER :: open_err								! Open error code
	INTEGER :: n, i_x, i_tot_x, i					! Loop index
	REAL(8) :: tot_x(num_N,num_tot_X)				! Matrix of available explanatory variables
	REAL(8) :: tot_xn(num_tot_X-1)
	!
	! Beginning execution
	!
	OPEN(UNIT = unit_model,FILE = file_model,STATUS = 'old',ACTION = 'read',IOSTAT=open_err)
	IF (open_err .NE. 0) THEN
		WRITE(*,*) 'PROBLEM: Unable to open file MODEL.TXT'
		STOP
	END IF
	!
	! First, skip the first row (containing names of states variables)
	!
	READ(unit_model,1001)
	1001 FORMAT()
	!
	! Second, read the following num_tot_K rows, num_states columns
	!
	DO n = 1, num_tot_X
        READ(unit_model,*) sel_X(n,:)
	    PRINT*, 'n = ', n
    END DO
	!
	CLOSE(UNIT = unit_model)
    !
    num_X_s = SUM(sel_X(:,1))
    num_X_z = SUM(sel_X(:,2))
    num_X_y = SUM(sel_X(:,3))
    !
    ! Ending execution and returning control
    !
	END SUBROUTINE input_model
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
 !   SUBROUTINE input_data ( theta_inf, theta_sup )
	!!
	!IMPLICIT NONE
	!! 
	!! Declaring dummy variables
	!!
	!REAL(8), INTENT(OUT) :: theta_inf(num_theta)	! Lower bounds for generating completely 
	!												! random starting points
	!REAL(8), INTENT(OUT) :: theta_sup(num_theta)	! Upper bounds for generating completely
	!												! random starting points
	!! 
	!! Declaring local variables
	!!
	!INTEGER :: open_err								! Open error code
	!INTEGER :: n, i_x, i_tot_x, i					! Loop index
	!REAL(8) :: tot_x(num_N,num_tot_X)				! Matrix of available explanatory variables
	!REAL(8) :: tot_xn(num_tot_X-1)
	!!
	!! Beginning execution
	!!
	!! Opening data file
	!!
	!OPEN ( UNIT=unit_data, FILE=file_data, STATUS='old', ACTION='read', IOSTAT=open_err )
	!IF (open_err .NE. 0) THEN
	!	WRITE(*,*) 'PROBLEM: Unable to open file	DATI.TXT'
	!	STOP
	!END IF
	!!
	!! Reading data
	!!
	!! First, skip the first row (containing names)
	!!
	!READ (unit_data,1001)
	!1001 FORMAT ()
	!!
	!! Second, read the following num_N rows
	!!
 !   tot_x(:,1) = 1.d0
	!DO n = 1, num_N
 !       READ (unit_data,*) d(n), ab(n), as(n), ar(n), c(n), tot_xn, peso(n)
	!    PRINT*, 'n = ', n, ' ; peso = ', peso(n)
 !       tot_x(n,2:) = tot_xn
 !   END DO
	!!
	!CLOSE ( UNIT=unit_data )
 !   normpeso = peso/REAL(SUM(peso))
	!!
 !   ! Select explanatory variables in pib
 !   !
 !   i_x = 0
 !   DO i_tot_x = 1, num_tot_X
 !       IF (sel_X(i_tot_x,1) .EQ. 1) THEN
 !           i_x = i_x+1
 !           x_pib(:,i_x) = tot_x(:,i_tot_x)
 !       END IF
 !   END DO
	!!
 !   ! Select explanatory variables in pis
 !   !
 !   i_x = 0
 !   DO i_tot_x = 1, num_tot_X
 !       IF (sel_X(i_tot_x,2) .EQ. 1) THEN
 !           i_x = i_x+1
 !           x_pis(:,i_x) = tot_x(:,i_tot_x)
 !       END IF
 !   END DO
	!!
 !   ! Select explanatory variables in xib
 !   !
 !   i_x = 0
 !   DO i_tot_x = 1, num_tot_X
 !       IF (sel_X(i_tot_x,3) .EQ. 1) THEN
 !           i_x = i_x+1
 !           x_xib(:,i_x) = tot_x(:,i_tot_x)
 !       END IF
 !   END DO
	!!
 !   ! Select explanatory variables in xis
 !   !
 !   i_x = 0
 !   DO i_tot_x = 1, num_tot_X
 !       IF (sel_X(i_tot_x,4) .EQ. 1) THEN
 !           i_x = i_x+1
 !           x_xis(:,i_x) = tot_x(:,i_tot_x)
 !       END IF
 !   END DO
	!!
 !   ! Select explanatory variables in rho
 !   !
 !   i_x = 0
 !   DO i_tot_x = 1, num_tot_X
 !       IF (sel_X(i_tot_x,5) .EQ. 1) THEN
 !           i_x = i_x+1
 !           x_rho(:,i_x) = tot_x(:,i_tot_x)
 !       END IF
 !   END DO
	!!
 !   ! Select explanatory variables in kappa
 !   !
 !   i_x = 0
 !   DO i_tot_x = 1, num_tot_X
 !       IF (sel_X(i_tot_x,6) .EQ. 1) THEN
 !           i_x = i_x+1
 !           x_kappa(:,i_x) = tot_x(:,i_tot_x)
 !       END IF
 !   END DO
	!!
 !   ! Select explanatory variables in omegab
 !   !
 !   i_x = 0
 !   DO i_tot_x = 1, num_tot_X
 !       IF (sel_X(i_tot_x,7) .EQ. 1) THEN
 !           i_x = i_x+1
 !           x_omegab(:,i_x) = tot_x(:,i_tot_x)
 !       END IF
 !   END DO
	!!
 !   ! Select explanatory variables in omegas
 !   !
 !   i_x = 0
 !   DO i_tot_x = 1, num_tot_X
 !       IF (sel_X(i_tot_x,8) .EQ. 1) THEN
 !           i_x = i_x+1
 !           x_omegas(:,i_x) = tot_x(:,i_tot_x)
 !       END IF
 !   END DO
 !   !
	!! Computing bounds for random starting values of the parameters
	!!
 !   theta_inf = -1.d0
 !   theta_sup = 1.d0
	!!
 !   ! Defining identity matrix(num_theta) to be used to compute derivatives
 !   !
 !   DO i = 1, num_psi
 !       idmat(i,:) = 0.d0
 !       idmat(i,i) = 1.d0
 !   END DO    
 !   !
 !   ! Ending execution and returning control
 !   !
	!END SUBROUTINE input_data
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
	!SUBROUTINE input_names ( )
	!!
	!IMPLICIT NONE
	!! 
	!! Declaring local variables
	!!
	!INTEGER :: open_err								! Open error code
	!INTEGER :: i_x, i_tot_x							! Loop index
	!CHARACTER(len=10) :: tot_names(num_tot_X)		! Available explanatory variables' names
	!!
	!! Beginning execution
	!!
	!! Opening names file
	!!
	!OPEN ( UNIT=unit_names, FILE=file_names, STATUS='old', ACTION='read', IOSTAT=open_err )
	!IF (open_err .NE. 0) THEN
	!	WRITE(*,*) 'PROBLEM: Unable to open file	NOMI.TXT'
	!	STOP
	!END IF
	!!
	!! Reading names
	!!
	!DO i_tot_x = 1, num_tot_X
	!	READ (unit_names,*) tot_names(i_tot_x)
	!END DO
	!!
	!CLOSE ( UNIT=unit_names )
	!!
 !   ! Select names of explanatory variables in pib
 !   !
 !   i_x = 0
 !   DO i_tot_x = 1, num_tot_X
 !       IF (sel_X(i_tot_x,1) .EQ. 1) THEN
 !           i_x = i_x+1
 !           names_pib(i_x) = tot_names(i_tot_x)
 !       END IF
 !   END DO
	!!
 !   ! Select names of explanatory variables in pis
 !   !
 !   i_x = 0
 !   DO i_tot_x = 1, num_tot_X
 !       IF (sel_X(i_tot_x,2) .EQ. 1) THEN
 !           i_x = i_x+1
 !           names_pis(i_x) = tot_names(i_tot_x)
 !       END IF
 !   END DO
	!!
 !   ! Select names of explanatory variables in xib
 !   !
 !   i_x = 0
 !   DO i_tot_x = 1, num_tot_X
 !       IF (sel_X(i_tot_x,3) .EQ. 1) THEN
 !           i_x = i_x+1
 !           names_xib(i_x) = tot_names(i_tot_x)
 !       END IF
 !   END DO
	!!
 !   ! Select names of explanatory variables in xis
 !   !
 !   i_x = 0
 !   DO i_tot_x = 1, num_tot_X
 !       IF (sel_X(i_tot_x,4) .EQ. 1) THEN
 !           i_x = i_x+1
 !           names_xis(i_x) = tot_names(i_tot_x)
 !       END IF
 !   END DO
	!!
 !   ! Select names of explanatory variables in rho
 !   !
 !   i_x = 0
 !   DO i_tot_x = 1, num_tot_X
 !       IF (sel_X(i_tot_x,5) .EQ. 1) THEN
 !           i_x = i_x+1
 !           names_rho(i_x) = tot_names(i_tot_x)
 !       END IF
 !   END DO
	!!
 !   ! Select names of explanatory variables in kappa
 !   !
 !   i_x = 0
 !   DO i_tot_x = 1, num_tot_X
 !       IF (sel_X(i_tot_x,6) .EQ. 1) THEN
 !           i_x = i_x+1
 !           names_kappa(i_x) = tot_names(i_tot_x)
 !       END IF
 !   END DO
	!!
 !   ! Select names of explanatory variables in omegab
 !   !
 !   i_x = 0
 !   DO i_tot_x = 1, num_tot_X
 !       IF (sel_X(i_tot_x,7) .EQ. 1) THEN
 !           i_x = i_x+1
 !           names_omegab(i_x) = tot_names(i_tot_x)
 !       END IF
 !   END DO
	!!
 !   ! Select names of explanatory variables in omegas
 !   !
 !   i_x = 0
 !   DO i_tot_x = 1, num_tot_X
 !       IF (sel_X(i_tot_x,8) .EQ. 1) THEN
 !           i_x = i_x+1
 !           names_omegas(i_x) = tot_names(i_tot_x)
 !       END IF
 !   END DO
	!!
	!! Ending execution and returning control
	!!
	!END SUBROUTINE input_names
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
END MODULE load_data