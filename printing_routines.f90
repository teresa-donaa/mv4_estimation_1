MODULE printing_routines
!
USE constants
USE observations
USE utilities
!
IMPLICIT NONE
!
CONTAINS
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE open_read_file ( unit_number, file_name )
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: unit_number
    CHARACTER(len=20), INTENT(IN) :: file_name
    !
    ! Declaring local variables
    !
    INTEGER :: open_err                 ! Open file error code
    !
    ! Beginning execution
    !
    OPEN ( UNIT=unit_number, FILE=file_name, STATUS='old', IOSTAT=open_err )
    IF (open_err .NE. 0) THEN
        WRITE(*,1234) file_name
        STOP
    END IF
    1234 FORMAT ('Unable to open input file ', A20)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE open_read_file 
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE open_write_file ( unit_number, file_name )
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: unit_number
    CHARACTER(len=20), INTENT(IN) :: file_name
    !
    ! Declaring local variables
    !
    INTEGER :: open_err                 ! Open file error code
    !
    ! Beginning execution
    !
    OPEN ( UNIT=unit_number, FILE=file_name, STATUS='replace', IOSTAT=open_err )
    IF (open_err .NE. 0) THEN
        WRITE(*,1234) file_name
        STOP
    END IF
    1234 FORMAT ('Unable to open input file ', A20)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE open_write_file 
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE print_res ( i_stime, objf, theta, task, grad )
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    INTEGER, INTENT(IN) :: i_stime      ! Estimation trial number
    REAL(8), INTENT(IN) :: objf         ! Latent criterion function at the optimum
    REAL(8), INTENT(IN) :: theta(num_theta)
    CHARACTER(len=60), INTENT(IN) :: task
    REAL(8), INTENT(IN) :: grad(num_theta)
    !
    ! Declaring local variables
    !
    REAL(8) :: psi(num_psi)
    REAL(8) :: m_s, m_q, m_z, m_y, sigma_s, sigma_z, sigma_y, rho_sz, rho_sy, rho_zy
    REAL(8) :: gamma, kappa, q, delta_z(num_L-1), delta_y(num_H-1)
    !
    ! Beginning execution
    !
    CALL compute_psi_parameters(x_m_s(1,:),x_m_q(1,:),x_m_y(1,:),x_sigma_s(1,:), &
        theta,psi)
    CALL compute_b_parameters(psi,m_s,m_q,m_z,m_y, &
            sigma_s,sigma_z,sigma_y,rho_sz,rho_sy,rho_zy, &
            gamma,kappa,q,delta_z,delta_y)
    !
    ! Print res file
    !
    WRITE (unit_res,35) i_stime, objf, theta, &
        sigma_z, sigma_y, rho_sz, rho_sy, rho_zy, &
        delta_z, delta_y, &
        error_flag, task, grad
35  FORMAT ( I3, " # ", ES25.18, " #", <num_theta>(1X, ES25.18), &
        " #", <5>(1X, ES15.8), &
        " #", <num_L-1>(1X, ES15.8), &
        " #", <num_H-1>(1X, ES15.8), &
        " # ", L5, " # ", A60, " # ", <num_theta>(ES15.8,1X))
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE print_res 
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!    SUBROUTINE print_final_results ( theta, stderr, objf, grad, vmat )
!    !
!    ! Computes final results and writes them on file
!    !
!    IMPLICIT NONE
!    !
!    ! Declaring dummy variables
!    !
!    REAL(8), INTENT(IN) :: theta(num_theta)             ! Estimates of parameters
!    REAL(8), INTENT(IN) :: stderr(num_theta)            ! Estimated asymptotic standard errors
!    REAL(8), INTENT(IN) :: objf                         ! Optimized latent criterion 
!    REAL(8), INTENT(IN) :: grad(num_theta)              ! Gradient at theta_star
!    REAL(8), INTENT(IN) :: vmat(num_theta,num_theta)    ! Estimated variance matrix of the parameters
!    !
!    ! Declaring local variables
!    !
!    INTEGER :: i_par, i_L, ind                          ! Loop indexes
!    !
!    ! Beginning execution
!    !
!    ! Computing statistics for the parameter estimates
!    !
!    ! m_pib
!    !
!    WRITE (unit_fin_res,3525) 
!    3525 FORMAT ( /, &
!        '-----------------------------------------------------------------', /, &
!        '    Parameter     Estimate       AsySE1  AsyT-Ratio1     Gradient', /, &
!        '-----------------------------------------------------------------', /, &
!        '                           M_PI B                                ', /, &
!        '-----------------------------------------------------------------' )
!    !
!    DO i_par = 1, num_X_pib
!        WRITE (unit_fin_res, 9999) names_pib(i_par), theta(i_par), stderr(i_par), &
!            theta(i_par)/stderr(i_par), grad(i_par)
!        9999 FORMAT ( 1X, A12, 4(1X, ES12.5) )
!    END DO
!    !
!    ! m_pis
!    !
!    WRITE (unit_fin_res,35251) 
!    35251 FORMAT ( &
!        '-----------------------------------------------------------------', /, &
!        '                           M_PI S                                ', /, &
!        '-----------------------------------------------------------------' )
!    !
!    ind = num_X_pib
!    DO i_par = ind+1, ind+num_X_pis
!        WRITE (unit_fin_res, 9999) names_pis(i_par-ind), theta(i_par), stderr(i_par), &
!            theta(i_par)/stderr(i_par), grad(i_par)
!    END DO
!    !
!    ! xib
!    !
!    WRITE (unit_fin_res,3526) 
!    3526 FORMAT ( &
!        '-----------------------------------------------------------------', /, &
!        '                            XI B                                 ', /, &
!        '-----------------------------------------------------------------' )
!    !
!    ind = num_X_pib+num_X_pis
!    DO i_par = ind+1, ind+num_X_xib
!        WRITE (unit_fin_res, 9999) names_xib(i_par-ind), theta(i_par), stderr(i_par), &
!            theta(i_par)/stderr(i_par), grad(i_par)
!    END DO
!    !
!    ! xis
!    !
!    WRITE (unit_fin_res,35261) 
!    35261 FORMAT ( &
!        '-----------------------------------------------------------------', /, &
!        '                            XI S                                 ', /, &
!        '-----------------------------------------------------------------' )
!    !
!    ind = num_X_pib+num_X_pis+num_X_xib
!    DO i_par = ind+1, ind+num_X_xis
!        WRITE (unit_fin_res, 9999) names_xis(i_par-ind), theta(i_par), stderr(i_par), &
!            theta(i_par)/stderr(i_par), grad(i_par)
!    END DO
!    !
!    ! rho
!    !
!    WRITE (unit_fin_res,35262) 
!    35262 FORMAT ( &
!        '-----------------------------------------------------------------', /, &
!        '                             RHO                                 ', /, &
!        '-----------------------------------------------------------------' )
!    !
!    ind = num_X_pib+num_X_pis+num_X_xib+num_X_xis
!    DO i_par = ind+1, ind+num_X_rho
!        WRITE (unit_fin_res, 9999) names_rho(i_par-ind), theta(i_par), stderr(i_par), &
!            theta(i_par)/stderr(i_par), grad(i_par)
!    END DO
!    !
!    ! kappa
!    !
!    WRITE (unit_fin_res,35263) 
!    35263 FORMAT ( &
!        '-----------------------------------------------------------------', /, &
!        '                              KAPPA                              ', /, &
!        '-----------------------------------------------------------------' )
!    !
!    ind = num_X_pib+num_X_pis+num_X_xib+num_X_xis+num_X_rho
!    DO i_par = ind+1, ind+num_X_kappa
!        WRITE (unit_fin_res, 9999) names_kappa(i_par-ind), theta(i_par), &
!            stderr(i_par), theta(i_par)/stderr(i_par), grad(i_par)
!    END DO
!    !
!    ! omegab
!    !
!    WRITE (unit_fin_res,89911) 
!    89911 FORMAT ( &
!        '-----------------------------------------------------------------', /, &
!        '                            OMEGA B                              ', /, &
!        '-----------------------------------------------------------------' )
!    ind = num_theta_beta 
!    DO i_par = ind+1, ind+num_X_omegab
!        WRITE (unit_fin_res, 9999) names_omegab(i_par-ind), theta(i_par), stderr(i_par), &
!            theta(i_par)/stderr(i_par), grad(i_par)
!    END DO
!    !
!    ! omegas
!    !
!    WRITE (unit_fin_res,89912) 
!    89912 FORMAT ( &
!        '-----------------------------------------------------------------', /, &
!        '                            OMEGA S                              ', /, &
!        '-----------------------------------------------------------------' )
!    ind = num_theta_beta+num_X_omegab
!    DO i_par = ind+1, ind+num_X_omegas
!        WRITE (unit_fin_res, 9999) names_omegas(i_par-ind), theta(i_par), stderr(i_par), &
!            theta(i_par)/stderr(i_par), grad(i_par)
!    END DO
!    !
!    ! sigma_pis
!    !
!    WRITE (unit_fin_res,89913) 
!    89913 FORMAT ( &
!        '-----------------------------------------------------------------', /, &
!        '                          SIGMA PIS                              ', /, &
!        '-----------------------------------------------------------------' )
!    ind = num_theta_beta+num_theta_omega+1
!    WRITE (unit_fin_res, 1991) EXP(theta(ind)), stderr(ind), &
!        EXP(theta(ind))/stderr(ind), grad(ind)
!1991 FORMAT ( 1X, 'sigma pis   ', 4(1X, ES12.5) )
!    !
!    ! sigma_pis
!    !
!    IF (switch_Sigma1 .EQ. 1) THEN
!        WRITE (unit_fin_res,89914) 
!        89914 FORMAT ( &
!            '-----------------------------------------------------------------', /, &
!            '                          R PIB PIS                              ', /, &
!            '-----------------------------------------------------------------' )
!        ind = num_theta_beta+num_theta_omega+num_theta_sigma
!        WRITE (unit_fin_res, 1992) max_rpibpis*(2.d0/greek_pi*DATAN(theta(ind))), stderr(ind), &
!            max_rpibpis*(2.d0/greek_pi*DATAN(theta(ind)))/stderr(ind), grad(ind)
!1992    FORMAT ( 1X, 'r pib pis   ', 4(1X, ES12.5) )
!    END IF
!    !
!    ! delta1
!    !
!    WRITE (unit_fin_res,8887) 
!    8887 FORMAT ( &
!        '-----------------------------------------------------------------', /, &
!        '                              DELTA                              ', /, &
!        '-----------------------------------------------------------------' )
!    ind = num_theta_beta+num_theta_omega+num_theta_sigma 
!    WRITE (unit_fin_res, 93991) num_theta_delta, EXP(theta(ind+1)), stderr(ind+1), &
!        EXP(theta(ind+1))/stderr(ind+1), grad(ind+1)
!    93991 FORMAT ( 1X, '    delta(', I1, ')', 4(1X, ES12.5) )
!    !
!    WRITE (unit_fin_res,3524) 
!    !
!    WRITE (unit_fin_res, 9997) -num_N*objf
!    9997 FORMAT ( 1X, '  Total objf', 1X, ES12.5 )
!    !
!    WRITE (unit_fin_res, 9992) -objf
!    9992 FORMAT ( 1X, '   Avg. objf', 1X, ES12.5 )
!    !
!    WRITE (unit_fin_res,3524) 
!    3524 FORMAT ( &
!        '---------------------------------------------------------------------------------------------------------------------' )
!    !
!    WRITE (unit_fin_res,14121) (vmat(ind,:), ind = 1, num_theta)
!    14121 FORMAT (1X, &
!        1X, 'Estimated asymptotic variance matrix of theta:', /, &
!        <num_theta>(1X,ES12.5))
!    !
!    WRITE (unit_fin_res,3524) 
!    !
!    ! Ending execution and returning control
!    !
!    END SUBROUTINE print_final_results
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
END MODULE printing_routines
