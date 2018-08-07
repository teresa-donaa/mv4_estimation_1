MODULE maxlik
!
USE constants
USE observations
USE printing_routines
USE utilities
!
IMPLICIT NONE
!
CONTAINS
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE estimate ( i_stime, theta, obj, grad, task )
    !
    IMPLICIT NONE
    ! 
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: i_stime
    REAL(8), INTENT(INOUT) :: theta(num_theta)
    REAL(8), INTENT(OUT) :: obj
    REAL(8), INTENT(OUT) :: grad(num_theta)
    CHARACTER(len=60), INTENT(OUT) :: task
    !
    ! Declaring local variables
    !
    INTEGER :: isave(44), iter, n
    INTEGER, PARAMETER :: m = 20, iprint = -1
    INTEGER, ALLOCATABLE :: nbd(:), iwa(:)
    REAL(8) :: dsave(29)
    REAL(8), ALLOCATABLE :: wa(:), l(:), u(:)
    CHARACTER(len=60) :: csave
    LOGICAL :: lsave(4)
    !
    ! Beginning execution
    !
    n = num_theta
    ALLOCATE ( nbd(n), iwa(3*n) )
    ALLOCATE ( wa(2*m*n + 5*n + 11*m*m + 8*m), l(n), u(n) )
    nbd = 0
    task = 'START'
    iter = 0
    DO WHILE ((task(1:2) .EQ. 'FG') .OR. (task .EQ. 'NEW_X') .OR. (task .EQ. 'START')) 
        !
        CALL setulb ( n, m, theta, l, u, nbd, obj, grad, factr, pgtol, &
            wa, iwa, task, iprint, csave, lsave, isave, dsave )
        ! 
        IF (task(1:5) .EQ. 'NEW_X') THEN
            WRITE(*,19) i_stime, iter, obj
19          FORMAT('Stima ', I4, ', Iterazione: ', I5, ' ; -LL/N: ', ES15.8)
            iter = iter+1
        END IF
        IF (task(1:2) .eq. 'FG') CALL loglik_fct ( n, theta, obj, grad )
        IF (error_FLAG) RETURN
        !
    END DO
    DEALLOCATE ( nbd, iwa )
    DEALLOCATE ( wa, l, u )
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE estimate
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE loglik_fct ( nn, theta, fc, gc ) 
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: nn
    REAL(8), INTENT(IN) :: theta(nn)        
    REAL(8), INTENT(OUT) :: fc
    REAL(8), INTENT(OUT) :: gc(nn)
    !
    ! Declaring local variables
    !
    INTEGER :: n
    REAL(8) :: psi(num_psi), pdpsi(num_psi), p, pd(num_theta)
    !
    ! Beginning execution
    !
    ! Initializing the loglikelihood and its gradient
    !
    fc = 0.d0
    gc = 0.d0
    !
    ! Beginning loop over observations
    !
    IF (to0 .EQ. 1) CALL open_write_file(unit_loglik,file_loglik)
    DO n = 1, num_N
        !
        ! Computing individual regime probability and associated gradient
        !
        CALL compute_psi_parameters(x_m_s(n,:),x_m_q(n,:),x_m_y(n,:),x_sigma_s(n,:), &
            theta,psi)
!@SP        CALL individual_contribution(psi,d(n),as(n),riskav(n),horizon(n),p)
        CALL INDIVIDUAL_CONTRIBUTION_DV(psi,idmat,d(n),as(n),riskav(n),horizon(n), &
            p,pdpsi,num_psi)
        !
        ! Adding individual contribution to total loglikelihood
        !
        IF ((p .LE. 0.d0) .OR. (ISNAN(p))) THEN
            !
!@SP       
OPEN(UNIT = unit_theta,FILE = file_theta)
WRITE(unit_theta,10) theta
10 FORMAT(<num_theta>(ES30.23, 1X))
CLOSE(UNIT = unit_theta)
CALL INDIVIDUAL_CONTRIBUTION(psi,d(n),as(n),riskav(n),horizon(n),p)
!@SP
            error_flag = .TRUE.
            fc = 1.d10
            RETURN
            !
        ELSE
            !
            fc = fc-normpeso(n)*LOG(p)
            !
            pd = MATMUL(matrix_dtheta_dpsiprime(x_m_s(n,:),x_m_q(n,:),x_m_y(n,:),x_sigma_s(n,:), &
                            theta),pdpsi)
            gc = gc-normpeso(n)*pd/p
            !
        END IF
        !
        IF (to0 .EQ. 1) THEN
            !
            WRITE(unit_loglik,1) n, d(n), as(n), riskav(n), horizon(n), &
                p, fc
1           FORMAT(I6, 1X, I3, 1X, ES15.8, 1X, <2>(I3, 1X), <2>(ES15.8, 1X))
            !
        END IF
        !
    END DO
    CLOSE(UNIT=unit_loglik)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE loglik_fct
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    FUNCTION loglik_politope ( theta ) 
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: theta(num_theta)        
    !
    ! Declaring local variables
    !
    INTEGER :: n
    REAL(8) :: psi(num_psi), p
    !
    ! Declaring function's type
    !
    REAL(8) :: loglik_politope
    !
    ! Beginning execution
    !
    ! Initializing the loglikelihood and its gradient
    !
    loglik_politope = 0.d0
    !
    ! Beginning loop over observations
    !
    IF (to0 .EQ. 1) CALL open_write_file(unit_loglik,file_loglik)
    DO n = 1, num_N
        !
        ! Computing individual regime probability and associated gradient
        !
        CALL compute_psi_parameters(x_m_s(n,:),x_m_q(n,:),x_m_y(n,:),x_sigma_s(n,:), &
            theta,psi)
        CALL individual_contribution(psi,d(n),as(n),riskav(n),horizon(n),p)
        !
        ! Adding individual contribution to total loglikelihood
        !
        IF ((p .LE. 0.d0) .OR. (ISNAN(p))) THEN
            !
!@SP       
OPEN(UNIT = unit_theta,FILE = file_theta)
WRITE(unit_theta,10) theta
10 FORMAT(<num_theta>(ES30.23, 1X))
CLOSE(UNIT = unit_theta)
CALL INDIVIDUAL_CONTRIBUTION(psi,d(n),as(n),riskav(n),horizon(n),p)
!@SP
            error_flag = .TRUE.
            loglik_politope = 1.d10
            RETURN
            !
        ELSE
            !
            loglik_politope = loglik_politope-normpeso(n)*LOG(p)
            !
        END IF
        !
        IF (to0 .EQ. 1) THEN
            !
            WRITE(unit_loglik,1) n, d(n), as(n), riskav(n), horizon(n), &
                p, loglik_politope
1           FORMAT(I6, 1X, I3, 1X, ES15.8, 1X, <2>(I3, 1X), <2>(ES15.8, 1X))
            !
        END IF
        !
    END DO
    CLOSE(UNIT=unit_loglik)
    !
    ! Ending execution and returning control
    !
    END FUNCTION loglik_politope
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE maxlik
