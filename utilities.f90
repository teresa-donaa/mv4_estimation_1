MODULE utilities
!
USE constants
USE observations
!
IMPLICIT NONE
!
CONTAINS
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE compute_psi_parameters ( xn_m_s, xn_m_z, xn_m_y, xn_sigma_s, &
        theta, psi )
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: xn_m_s(num_X_m_s)    
    REAL(8), INTENT(IN) :: xn_m_z(num_X_m_z)    
    REAL(8), INTENT(IN) :: xn_m_y(num_X_m_y)
    REAL(8), INTENT(IN) :: xn_sigma_s(num_X_sigma_s)
    REAL(8), INTENT(IN) :: theta(num_theta)     ! Deep parameters
    REAL(8), INTENT(OUT) :: psi(num_psi)  
    !
    ! Declaring local parameters
    !
    INTEGER :: il, iu
    REAL(8) :: arg(num_tot_X)
    !
    ! Beginning execution
    !
    ! Metaparameters
    !
    ! m_s
    !
    arg = 0.d0
    il = 1
    iu = num_X_m_s
    arg(:num_X_m_s) = theta(il:iu)*xn_m_s
    psi(1) = SUM(arg)
    !
    ! m_z
    !
    arg = 0.d0
    il = iu+1
    iu = iu+num_X_m_z
    arg(:num_X_m_z) = theta(il:iu)*xn_m_z
    psi(2) = SUM(arg)
    !
    ! m_y
    !
    arg = 0.d0
    il = iu+1
    iu = iu+num_X_m_y
    arg(:num_X_m_y) = theta(il:iu)*xn_m_y
    psi(num_psi_m) = SUM(arg)
    !
    ! sigma_s
    !
    arg = 0.d0
    il = iu+1
    iu = iu+num_X_sigma_s
    arg(:num_X_sigma_s) = theta(il:iu)*xn_sigma_s
    psi(num_psi_m+num_psi_sigma) = SUM(arg)
    ! 
    ! delta_z & delta_y
    !
    psi(num_psi_m+num_psi_sigma+1:) = theta(num_theta_beta+num_theta_sigma_s+1:)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE compute_psi_parameters
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    !SUBROUTINE compute_b_parameters ( psi, mpis, mpis, xib, xis, rho, kappa, &
    !    omegab, omegas, delta, sigmapib, sigmapis, rpibpis, sigmaz )
    !
    ! Declaring dummy variables
    !
    !REAL(8), INTENT(IN) :: psi(num_psi)
    !REAL(8), INTENT(OUT) :: mpib
    !REAL(8), INTENT(OUT) :: mpis
    !REAL(8), INTENT(OUT) :: xib
    !REAL(8), INTENT(OUT) :: xis
    !REAL(8), INTENT(OUT) :: rho
    !REAL(8), INTENT(OUT) :: kappa                   
    !REAL(8), INTENT(OUT) :: omegab
    !REAL(8), INTENT(OUT) :: omegas
    !REAL(8), INTENT(OUT) :: delta(num_theta_delta)
    !REAL(8), INTENT(OUT) :: sigmapib
    !REAL(8), INTENT(OUT) :: sigmapis
    !REAL(8), INTENT(OUT) :: rpibpis
    !REAL(8), INTENT(OUT) :: sigmaz
    !
    ! Beginning execution
    !
    ! Unconditional mean parameters
    !
    !mpib = psi(1)
    !mpis = psi(2)
    !xib = psi(3)
    !xis = psi(4)
    !rho = 2.d0/pi*ATAN(psi(5))
    !kappa = EXP(psi(num_psi_m))
    !
    !omegab = EXP(psi(num_psi_m+1))
    !omegas = EXP(psi(num_psi_m+num_psi_omega))
    !
    !sigmapib = 1.d0
    !sigmapis = EXP(psi(num_psi_m+num_psi_omega+1))
    !IF (switch_Sigma1 .EQ. 1) THEN
    !    rpibpis = max_rpibpis*(2.d0/pi*ATAN(psi(num_psi_m+num_psi_omega+num_psi_sigma)))
    !ELSE IF (switch_Sigma2 .EQ. 1) THEN
    !    rpibpis = 0.d0
    !END IF
    !
    !delta(1) = EXP(psi(num_psi))
    !
    !sigmaz = 1.d0
    !
    ! Ending execution and returning control
    !
    !END SUBROUTINE compute_b_parameters
!! 
!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!! 
!    FUNCTION matrix_dtheta_dpsiprime ( xn_pib, xn_pis, xn_xib, xn_xis, &
!        xn_rho, xn_kappa, xn_omegab, xn_omegas, theta )
!    ! 
!    ! Computes the matrix d theta / d psi',
!    ! dove theta è il vettore dei parametri rispetto ai quali si minimizza - log(L)
!    ! e psi è il vettore dei parametri rispetto ai quali si deriva
!    ! 
!    ! Declaring dummy variables
!    ! 
!    REAL(8), INTENT(IN) :: xn_pib(num_X_pib)    
!    REAL(8), INTENT(IN) :: xn_pis(num_X_pis)    
!    REAL(8), INTENT(IN) :: xn_xib(num_X_xib)
!    REAL(8), INTENT(IN) :: xn_xis(num_X_xis)
!    REAL(8), INTENT(IN) :: xn_rho(num_X_rho)
!    REAL(8), INTENT(IN) :: xn_kappa(num_X_kappa)
!    REAL(8), INTENT(IN) :: xn_omegab(num_X_omegab)
!    REAL(8), INTENT(IN) :: xn_omegas(num_X_omegas)
!    REAL(8), INTENT(IN) :: theta(num_theta)
!    !
!    ! Declaring local variables
!    !
!    INTEGER :: il, iu
!    ! 
!    ! Declaring function type
!    ! 
!    REAL(8) :: matrix_dtheta_dpsiprime(num_theta,num_psi)
!    ! 
!    ! Beginning execution
!    ! 
!    ! Initialising matrix_dtheta_dpsiprime
!    ! 
!    matrix_dtheta_dpsiprime = 0.d0
!    !
!    ! betapib
!    !
!    il = 1
!    iu = num_X_pib
!    matrix_dtheta_dpsiprime(il:iu,1) = xn_pib
!    !
!    ! betapis
!    !
!    il = iu+1
!    iu = iu+num_X_pis
!    matrix_dtheta_dpsiprime(il:iu,2) = xn_pis
!    !
!    ! betaxib
!    !
!    il = iu+1
!    iu = iu+num_X_xib
!    matrix_dtheta_dpsiprime(il:iu,3) = xn_xib
!    !
!    ! betaxis
!    !
!    il = iu+1
!    iu = iu+num_X_xis
!    matrix_dtheta_dpsiprime(il:iu,4) = xn_xis
!    !
!    ! betarho
!    !
!    il = iu+1
!    iu = iu+num_X_rho
!    matrix_dtheta_dpsiprime(il:iu,5) = xn_rho
!    !
!    ! betakappa
!    !
!    il = iu+1
!    iu = num_theta_beta
!    matrix_dtheta_dpsiprime(il:iu,num_psi_m) = xn_kappa
!    !
!    ! omegab
!    !
!    il = iu+1
!    iu = iu+num_X_omegab
!    matrix_dtheta_dpsiprime(il:iu,num_psi_m+1) = xn_omegab
!    !
!    ! omegas
!    !
!    il = iu+1
!    iu = num_theta_beta+num_theta_omega
!    matrix_dtheta_dpsiprime(il:iu,num_psi_m+num_psi_omega) = xn_omegas
!    !
!    ! sigmapis
!    !
!    iu = iu+1
!    matrix_dtheta_dpsiprime(iu,num_psi_m+num_psi_omega+1) = 1.d0
!    !
!    ! rpibpis
!    !
!    IF (switch_Sigma1 .EQ. 1) THEN
!        iu = num_theta_beta+num_theta_omega+num_theta_sigma
!        matrix_dtheta_dpsiprime(iu,num_psi_m+num_psi_omega+num_psi_sigma) = 1.d0
!    END IF
!    !
!    ! delta
!    !
!    iu = num_theta
!    matrix_dtheta_dpsiprime(iu,num_psi) = 1.d0
!    ! 
!    ! Ending execution and returning control
!    ! 
!    END FUNCTION matrix_dtheta_dpsiprime
!! 
!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!! 
!    FUNCTION matrix_dthetaF_dthetaprime ( theta )
!    ! 
!    ! Computes the matrix d theta_F / d theta',
!    ! dove theta_F è il vettore dei parametri "finali" in output (il vero delta_1)
!    ! e theta è il vettore dei parametri rispetto ai quali si minimizza - log(L)
!    ! 
!    ! Declaring dummy variables
!    ! 
!    REAL(8), INTENT(IN) :: theta(num_theta)
!    ! 
!    ! Declaring function type
!    ! 
!    REAL(8) :: matrix_dthetaF_dthetaprime(num_theta,num_theta)
!    !
!    ! Declaring local variables
!    !
!    INTEGER :: i, iu
!    ! 
!    ! Beginning execution
!    ! 
!    ! Initialising matrix_dtheta_dbprime
!    ! 
!    matrix_dthetaF_dthetaprime = 0.d0
!    !
!    ! Diagonale uguale a 1, a parte gli ultimi tre elementi in basso a destra
!    !
!    DO i = 1, num_theta-3
!        matrix_dthetaF_dthetaprime(i,i) = 1.d0
!    END DO
!    !
!    ! sigmapis
!    !
!    iu = num_theta_beta+num_theta_omega+1
!    matrix_dthetaF_dthetaprime(iu,iu) = EXP(theta(iu))
!    !
!    ! rpibpis
!    !
!    IF (switch_Sigma1 .EQ. 1) THEN
!        iu = num_theta_beta+num_theta_omega+num_theta_sigma
!        matrix_dthetaF_dthetaprime(iu,iu) = 2.d0*max_rpibpis/(pi*(1.d0+theta(iu)**2))
!    END IF
!    !
!    ! delta1
!    !
!    iu = num_theta
!    matrix_dthetaF_dthetaprime(iu,iu) = EXP(theta(iu))
!    ! 
!    ! Ending execution and returning control
!    ! 
!    END FUNCTION matrix_dthetaF_dthetaprime
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
END MODULE utilities
