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
    SUBROUTINE compute_psi_parameters ( xn_m_s, xn_m_q, xn_m_y, xn_sigma_s, &
        theta, psi )
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: xn_m_s(num_X_m_s)    
    REAL(8), INTENT(IN) :: xn_m_q(num_X_m_q)    
    REAL(8), INTENT(IN) :: xn_m_y(num_X_m_y)
    REAL(8), INTENT(IN) :: xn_sigma_s(num_X_sigma_s)
    REAL(8), INTENT(IN) :: theta(num_theta)     ! Deep parameters
    REAL(8), INTENT(OUT) :: psi(num_psi)  
    !
    ! Declaring local parameters
    !
    INTEGER :: il, iu, ipsi
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
    ! m_q
    !
    arg = 0.d0
    il = iu+1
    iu = iu+num_X_m_q
    arg(:num_X_m_q) = theta(il:iu)*xn_m_q
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
    psi(num_psi_m+1) = SUM(arg)
    ipsi = num_psi_m+1
    !
    ! sigma_z
    !
    IF (switch_sigma_z .EQ. 1) THEN
        !
        il = iu+1
        iu = il
        ipsi = ipsi+1
        psi(ipsi) = theta(il)
        !
    END IF
    !
    ! sigma_y
    !
    IF (switch_sigma_y .EQ. 1) THEN
        !
        il = iu+1
        iu = il
        ipsi = ipsi+1
        psi(ipsi) = theta(il)
        !
    END IF
    !
    ! rho_sz
    !
    IF (switch_rho_sz .EQ. 1) THEN
        !
        il = iu+1
        iu = il
        ipsi = ipsi+1
        psi(ipsi) = theta(il)
        !
    END IF
    !
    ! rho_sy
    !
    IF (switch_rho_sy .EQ. 1) THEN
        !
        il = iu+1
        iu = il
        ipsi = ipsi+1
        psi(ipsi) = theta(il)
        !
    END IF
    !
    ! rho_zy
    !
    IF (switch_rho_zy .EQ. 1) THEN
        !
        il = iu+1
        iu = il
        ipsi = ipsi+1
        psi(ipsi) = theta(il)
        !
    END IF
    ! 
    ! delta_z & delta_y
    !
    psi(num_psi_m+num_psi_sigma+num_psi_rho+1:) = &
        theta(num_theta_beta+num_theta_sigma_s+num_theta_sigma_z_y+num_theta_rho+1:)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE compute_psi_parameters
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE compute_b_parameters ( psi, m_s, m_q, m_z, m_y, &
        sigma_s, sigma_z, sigma_y, rho_sz, rho_sy, rho_zy, &
        gamma, kappa, q, delta_z, delta_y )
    !
    ! Compute unconditional parameters
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: psi(num_psi)
    REAL(8), INTENT(OUT) :: m_s
    REAL(8), INTENT(OUT) :: m_q
    REAL(8), INTENT(OUT) :: m_z
    REAL(8), INTENT(OUT) :: m_y
    REAL(8), INTENT(OUT) :: sigma_s
    REAL(8), INTENT(OUT) :: sigma_z
    REAL(8), INTENT(OUT) :: sigma_y
    REAL(8), INTENT(OUT) :: rho_sz
    REAL(8), INTENT(OUT) :: rho_sy
    REAL(8), INTENT(OUT) :: rho_zy
    REAL(8), INTENT(OUT) :: gamma
    REAL(8), INTENT(OUT) :: kappa
    REAL(8), INTENT(OUT) :: q
    REAL(8), INTENT(OUT) :: delta_z(num_delta_z)
    REAL(8), INTENT(OUT) :: delta_y(num_delta_y)
    !
    ! Declaring local variables
    !
    INTEGER :: ipsi, il, iu, i
    !
    ! Beginning 
    !
    ! Beginning execution
    ! 
    ! means
    !
    m_s = psi(1)
    m_q = psi(2)
    m_y = psi(3)
    !
    ! risk aversion and planning horizon
    !
    q = EXP(psi(2))
    kappa = EXP(psi(3))
    gamma = (q+kappa/apriori_T)/(1.d0+kappa/apriori_T)
    m_z = LOG(gamma)
    !
    ! sigma_s
    !
    sigma_s = minimum_sigma_s+EXP(psi(num_psi_m+1))
    ipsi = num_psi_m+1
    !
    ! sigma_z
    !
    IF (switch_sigma_z .EQ. 1) THEN
        !
        ipsi = ipsi+1
        sigma_z = EXP(psi(ipsi))
        !
    ELSE IF (switch_sigma_z .EQ. 0) THEN
        !
        sigma_z = 1.d0
        !
    END IF
    !
    ! sigma_y
    !
    IF (switch_sigma_y .EQ. 1) THEN
        !
        ipsi = ipsi+1
        sigma_y = EXP(psi(ipsi))
        !
    ELSE IF (switch_sigma_y .EQ. 0) THEN
        !
        sigma_y = 1.d0
        !
    END IF
    !
    ! rho_sz
    !
    IF (switch_rho_sz .EQ. 1) THEN
        !
        ipsi = ipsi+1
        rho_sz = twooverpi*ATAN(psi(ipsi))
        !
    ELSE IF (switch_rho_sz .EQ. 0) THEN
        !
        rho_sz = 0.d0
        !
    END IF
    !
    ! rho_sy
    !
    IF (switch_rho_sy .EQ. 1) THEN
        !
        ipsi = ipsi+1
        rho_sy = twooverpi*ATAN(psi(ipsi))
        !
    ELSE IF (switch_rho_sy .EQ. 0) THEN
        !
        rho_sy = 0.d0
        !
    END IF
    !
    ! rho_zy
    !
    IF (switch_rho_zy .EQ. 1) THEN
        !
        ipsi = ipsi+1
        rho_zy = rho_sz*rho_sy+SQRT((1.d0-rho_sz**2)*(1.d0-rho_sy**2))*twooverpi*ATAN(psi(ipsi))
        !
    END IF
    !
    ! delta_z
    !
    il = ipsi
    delta_z(1) = 0.d0
    delta_z(2:) = EXP(psi(il+1:il+num_theta_delta_z))
    DO i = 3, num_delta_z
        !
        delta_z(i) = delta_z(i)+delta_z(i-1)
        !
    END DO
    !
    ! delta_y
    !
    il = il+num_theta_delta_z
    delta_y(1) = 0.d0
    delta_y(2:) = EXP(psi(il+1:il+num_theta_delta_y))
    DO i = 3, num_delta_y
        !
        delta_y(i) = delta_y(i)+delta_y(i-1)
        !
    END DO
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE compute_b_parameters
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE compute_thetaF_parameters ( theta, thetaF )
    !
    ! Computes the theta_F parameters from the theta parameters, 
    ! where theta_F is the "final" parameters vector on output 
    ! and theta is the vector of parameters w.r.t. which we minimize -log(L)
    !
    IMPLICIT NONE
    !
    ! Declaring dummy variables
    !
    REAL(8), INTENT(IN) :: theta(num_theta)
    REAL(8), INTENT(OUT) :: thetaF(num_theta)
    !
    ! Declaring local variables
    !
    INTEGER :: i, iu
    REAL(8) :: rho_sz, rho_sy
    ! 
    ! Beginning execution
    ! 
    ! thetaF = theta for the parameters in the means and in sigma_s
    !
    DO i = 1, num_theta_beta+num_theta_sigma_s
        !
        thetaF(i) = theta(i)
        !
    END DO
    iu = num_theta_beta+num_theta_sigma_s
    !
    ! sigma_z
    !
    IF (switch_sigma_z .EQ. 1) THEN
        !
        iu = iu+1
        thetaF(iu) = EXP(theta(iu))
        !
    END IF
    !
    ! sigma_y
    !
    IF (switch_sigma_y .EQ. 1) THEN
        !
        iu = iu+1
        thetaF(iu) = EXP(theta(iu))
        !
    END IF
    !
    ! rho_sz
    !
    IF (switch_rho_sz .EQ. 1) THEN
        !
        iu = iu+1
        rho_sz = twooverpi*ATAN(theta(iu))
        thetaF(iu) = rho_sz
        !
    ELSE IF (switch_rho_sz .EQ. 0) THEN
        !
        rho_sz = 0.d0
        !
    END IF
    !
    ! rho_sy
    !
    IF (switch_rho_sy .EQ. 1) THEN
        !
        iu = iu+1
        rho_sy = twooverpi*ATAN(theta(iu))
        thetaF(iu) = rho_sy
        !
    ELSE IF (switch_rho_sy .EQ. 0) THEN
        !
        rho_sy = 0.d0
        !
    END IF
    !
    ! rho_zy
    !
    IF (switch_rho_zy .EQ. 1) THEN
        !
        iu = iu+1
        thetaF(iu) = rho_sz*rho_sy+SQRT((1.d0-rho_sz**2)*(1.d0-rho_sy**2))*twooverpi*ATAN(theta(iu))
        !
    END IF
    !
    ! delta_z
    !
    thetaF(iu+1:iu+num_theta_delta_z) = EXP(theta(iu+1:iu+num_theta_delta_z))
    DO i = iu+2, iu+num_theta_delta_z
        !
        thetaF(i) = thetaF(i)+thetaF(i-1)
        !
    END DO
    iu = iu+num_theta_delta_z
    !
    ! delta_y
    !
    thetaF(iu+1:iu+num_theta_delta_y) = EXP(theta(iu+1:iu+num_theta_delta_y))
    DO i = iu+2, iu+num_theta_delta_y
        !
        thetaF(i) = thetaF(i)+thetaF(i-1)
        !
    END DO
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE compute_thetaF_parameters
! 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 
    FUNCTION matrix_dtheta_dpsiprime ( xn_m_s, xn_m_q, xn_m_y, xn_sigma_s, theta )
    ! 
    ! Computes the matrix d theta / d psi',
    ! where theta are the parameters w.r.t. which we minimize -log(L),
    ! and psi are the parameters w.r.t. which we differentiate
     !
	IMPLICIT NONE
	!
    ! Declaring dummy variables
    ! 
    REAL(8), INTENT(IN) :: xn_m_s(num_X_m_s)    
    REAL(8), INTENT(IN) :: xn_m_q(num_X_m_q)    
    REAL(8), INTENT(IN) :: xn_m_y(num_X_m_y)
    REAL(8), INTENT(IN) :: xn_sigma_s(num_X_sigma_s)
    REAL(8), INTENT(IN) :: theta(num_theta)     ! Deep parameters
    !
    ! Declaring local variables
    !
    INTEGER :: il, iu, ipsi, i
    ! 
    ! Declaring function type
    ! 
    REAL(8) :: matrix_dtheta_dpsiprime(num_theta,num_psi)
    ! 
    ! Beginning execution
    ! 
    ! Initialising matrix_dtheta_dpsiprime
    ! 
    matrix_dtheta_dpsiprime = 0.d0
    !
    ! beta_m_s
    !
    il = 1
    iu = num_X_m_s
    matrix_dtheta_dpsiprime(il:iu,1) = xn_m_s
    !
    ! beta_m_q
    !
    il = iu+1
    iu = iu+num_X_m_q
    matrix_dtheta_dpsiprime(il:iu,2) = xn_m_q
    !
    ! beta_m_y
    !
    il = iu+1
    iu = iu+num_X_m_y
    matrix_dtheta_dpsiprime(il:iu,3) = xn_m_y
    !
    ! beta_sigma_s
    !
    il = iu+1
    iu = iu+num_X_sigma_s
    matrix_dtheta_dpsiprime(il:iu,4) = xn_sigma_s
    ipsi = 4
    !
    ! sigma_z
    !
    IF (switch_sigma_z .EQ. 1) THEN
        !
        il = iu+1
        iu = il
        ipsi = ipsi+1
        matrix_dtheta_dpsiprime(iu,ipsi) = 1.d0
        !
    END IF
    !
    ! sigma_y
    !
    IF (switch_sigma_y .EQ. 1) THEN
        !
        il = iu+1
        iu = il
        ipsi = ipsi+1
        matrix_dtheta_dpsiprime(iu,ipsi) = 1.d0
        !
    END IF
    !
    ! rho_sz
    !
    IF (switch_rho_sz .EQ. 1) THEN
        !
        il = iu+1
        iu = il
        ipsi = ipsi+1
        matrix_dtheta_dpsiprime(iu,ipsi) = 1.d0
        !
    END IF
    !
    ! rho_sy
    !
    IF (switch_rho_sy .EQ. 1) THEN
        !
        il = iu+1
        iu = il
        ipsi = ipsi+1
        matrix_dtheta_dpsiprime(iu,ipsi) = 1.d0
        !
    END IF
    !
    ! rho_zy
    !
    IF (switch_rho_zy .EQ. 1) THEN
        !
        il = iu+1
        iu = il
        ipsi = ipsi+1
        matrix_dtheta_dpsiprime(iu,ipsi) = 1.d0
        !
    END IF
    !
    ! delta_z
    !
    DO i = 1, num_theta_delta_z
        !
        il = iu+1
        iu = il
        ipsi = ipsi+1
        matrix_dtheta_dpsiprime(iu,ipsi) = 1.d0
        !
    END DO
    !
    ! delta_y
    !
    DO i = 1, num_theta_delta_y
        !
        il = iu+1
        iu = il
        ipsi = ipsi+1
        matrix_dtheta_dpsiprime(iu,ipsi) = 1.d0
        !
    END DO
    ! 
    ! Ending execution and returning control
    ! 
    END FUNCTION matrix_dtheta_dpsiprime
! 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 
    FUNCTION matrix_dthetaF_dthetaprime ( theta )
    ! 
    ! Computes the matrix d theta_F / d theta',
    ! where theta_F is the "final" parameters vector on output 
    ! and theta is the vector of parameters w.r.t. which we minimize -log(L)
    !
    IMPLICIT NONE
    ! 
    ! Declaring dummy variables
    ! 
    REAL(8), INTENT(IN) :: theta(num_theta)
    ! 
    ! Declaring function type
    ! 
    REAL(8) :: matrix_dthetaF_dthetaprime(num_theta,num_theta)
    !
    ! Declaring local variables
    !
    INTEGER :: i, iu
    REAL(8) :: rho_sz, rho_sy, rho_zy
    ! 
    ! Beginning execution
    ! 
    ! Initialising matrix_dtheta_dbprime
    ! 
    matrix_dthetaF_dthetaprime = 0.d0
    !
    ! Unit diagonal block for the parameters in the means and in sigma_s
    !
    DO i = 1, num_theta_beta+num_theta_sigma_s
        !
        matrix_dthetaF_dthetaprime(i,i) = 1.d0
        !
    END DO
    iu = num_theta_beta+num_theta_sigma_s
    !
    ! sigma_z
    !
    IF (switch_sigma_z .EQ. 1) THEN
        !
        iu = iu+1
        matrix_dthetaF_dthetaprime(iu,iu) = EXP(theta(iu))
        !
    END IF
    !
    ! sigma_y
    !
    IF (switch_sigma_y .EQ. 1) THEN
        !
        iu = iu+1
        matrix_dthetaF_dthetaprime(iu,iu) = EXP(theta(iu))
        !
    END IF
    !
    ! Correlation parameters
    !
    IF ((switch_rho_sz .EQ. 1) .AND. (switch_rho_sy .EQ. 1)) THEN
        !
        rho_sz = twooverpi*ATAN(theta(iu+1))
        rho_sy = twooverpi*ATAN(theta(iu+2))
        rho_zy = rho_sz*rho_sy+SQRT((1.d0-rho_sz**2)*(1.d0-rho_sy**2))*twooverpi*ATAN(theta(iu+3))
        matrix_dthetaF_dthetaprime(iu+1,iu+1) = twooverpi/(1.d0+theta(iu+1)**2)
        matrix_dthetaF_dthetaprime(iu+3,iu+1) = twooverpi/(1.d0+theta(iu+1)**2)* &
            (rho_sy-rho_sz*(rho_zy-rho_sz*rho_sy)/(1.d0-rho_sz**2))
        matrix_dthetaF_dthetaprime(iu+2,iu+2) = twooverpi/(1.d0+theta(iu+2)**2)
        matrix_dthetaF_dthetaprime(iu+3,iu+2) = twooverpi/(1.d0+theta(iu+2)**2)* &
            (rho_sz-rho_sy*(rho_zy-rho_sz*rho_sy)/(1.d0-rho_sy**2))
        matrix_dthetaF_dthetaprime(iu+3,iu+3) = twooverpi/(1.d0+theta(iu+3)**2)* &
            SQRT((1.d0-rho_sz**2)*(1.d0-rho_sy**2))
        iu = iu+3
        !
    ELSE IF ((switch_rho_sz .EQ. 1) .AND. (switch_rho_sy .EQ. 0)) THEN
        !
        rho_sz = twooverpi*ATAN(theta(iu+1))
        rho_zy = SQRT(1.d0-rho_sz**2)*twooverpi*ATAN(theta(iu+2))
        matrix_dthetaF_dthetaprime(iu+1,iu+1) = twooverpi/(1.d0+theta(iu+1)**2)
        matrix_dthetaF_dthetaprime(iu+2,iu+1) = -twooverpi/(1.d0+theta(iu+1)**2)* &
            rho_sz*rho_zy/(1.d0-rho_sz**2)
        matrix_dthetaF_dthetaprime(iu+2,iu+2) = twooverpi/(1.d0+theta(iu+2)**2)* &
            SQRT(1.d0-rho_sz**2)
        iu = iu+2
        !
    ELSE IF ((switch_rho_sz .EQ. 0) .AND. (switch_rho_sy .EQ. 1)) THEN
        !
        rho_sy = twooverpi*ATAN(theta(iu+1))
        rho_zy = SQRT(1.d0-rho_sy**2)*twooverpi*ATAN(theta(iu+2))
        matrix_dthetaF_dthetaprime(iu+1,iu+1) = twooverpi/(1.d0+theta(iu+1)**2)
        matrix_dthetaF_dthetaprime(iu+2,iu+1) = -twooverpi/(1.d0+theta(iu+1)**2)* &
            rho_sy*rho_zy/(1.d0-rho_sy**2)
        matrix_dthetaF_dthetaprime(iu+2,iu+2) = twooverpi/(1.d0+theta(iu+2)**2)* &
            SQRT(1.d0-rho_sy**2)
        iu = iu+2
        !
    ELSE IF ((switch_rho_sz .EQ. 0) .AND. (switch_rho_sy .EQ. 0)) THEN
        !
        rho_zy = twooverpi*ATAN(theta(iu+1))
        matrix_dthetaF_dthetaprime(iu+1,iu+1) = twooverpi/(1.d0+theta(iu+1)**2)
        !
    END IF
    !
    ! delta_z
    !
    DO i = iu+1, iu+num_theta_delta_z
        !
        matrix_dthetaF_dthetaprime(i:iu+num_theta_delta_z,i) = EXP(theta(i))
        !
    END DO
    iu = iu+num_theta_delta_z
    !
    ! delta_y
    !
    DO i = iu+1, iu+num_theta_delta_y
        !
        matrix_dthetaF_dthetaprime(i:iu+num_theta_delta_y,i) = EXP(theta(i))
        !
    END DO
    ! 
    ! Ending execution and returning control
    ! 
    END FUNCTION matrix_dthetaF_dthetaprime
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
END MODULE utilities
