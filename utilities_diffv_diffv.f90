!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.13 (r6666M) -  1 Mar 2018 15:30
!
!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.13 (r6666M) -  1 Mar 2018 15:30
!
MODULE UTILITIES_DIFFV_DIFFV
!
  USE CONSTANTS
  USE OBSERVATIONS
  USE DIFFSIZES
  USE DIFFSIZES
!  Hint: nbdirsmax should be the maximum number of differentiation directions
  IMPLICIT NONE

CONTAINS
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
  SUBROUTINE COMPUTE_PSI_PARAMETERS(xn_m_s, xn_m_q, xn_m_y, xn_sigma_s, &
&   theta, psi)
    IMPLICIT NONE
!
! Ending execution and returning control
!
!
! Declaring dummy variables
!
    REAL*8, INTENT(IN) :: xn_m_s(num_x_m_s)
    REAL*8, INTENT(IN) :: xn_m_q(num_x_m_q)
    REAL*8, INTENT(IN) :: xn_m_y(num_x_m_y)
    REAL*8, INTENT(IN) :: xn_sigma_s(num_x_sigma_s)
! Deep parameters
    REAL*8, INTENT(IN) :: theta(num_theta)
    REAL*8, INTENT(OUT) :: psi(num_psi)
!
! Declaring local parameters
!
    INTEGER :: il, iu, ipsi
    REAL*8 :: arg(num_tot_x)
    INTRINSIC SUM
!
! Beginning execution
!
! Metaparameters
!
! m_s
!
    arg = 0.d0
    il = 1
    iu = num_x_m_s
    arg(:num_x_m_s) = theta(il:iu)*xn_m_s
    psi(1) = SUM(arg)
!
! m_q
!
    arg = 0.d0
    il = iu + 1
    iu = iu + num_x_m_q
    arg(:num_x_m_q) = theta(il:iu)*xn_m_q
    psi(2) = SUM(arg)
!
! m_y
!
    arg = 0.d0
    il = iu + 1
    iu = iu + num_x_m_y
    arg(:num_x_m_y) = theta(il:iu)*xn_m_y
    psi(num_psi_m) = SUM(arg)
!
! sigma_s
!
    arg = 0.d0
    il = iu + 1
    iu = iu + num_x_sigma_s
    arg(:num_x_sigma_s) = theta(il:iu)*xn_sigma_s
    psi(num_psi_m+1) = SUM(arg)
    ipsi = num_psi_m + 1
!
! sigma_z
!
    IF (switch_sigma_z .EQ. 1) THEN
!
      il = iu + 1
      iu = il
      ipsi = ipsi + 1
      psi(ipsi) = theta(il)
!
    END IF
!
! sigma_y
!
    IF (switch_sigma_y .EQ. 1) THEN
!
      il = iu + 1
      iu = il
      ipsi = ipsi + 1
      psi(ipsi) = theta(il)
!
    END IF
!
! rho_sz
!
    IF (switch_rho_sz .EQ. 1) THEN
!
      il = iu + 1
      iu = il
      ipsi = ipsi + 1
      psi(ipsi) = theta(il)
!
    END IF
!
! rho_sy
!
    IF (switch_rho_sy .EQ. 1) THEN
!
      il = iu + 1
      iu = il
      ipsi = ipsi + 1
      psi(ipsi) = theta(il)
!
    END IF
!
! rho_zy
!
    IF (switch_rho_zy .EQ. 1) THEN
!
      il = iu + 1
      iu = il
      ipsi = ipsi + 1
      psi(ipsi) = theta(il)
!
    END IF
! 
! delta_z & delta_y
!
    psi(num_psi_m+num_psi_sigma+num_psi_rho+1:) = theta(num_theta_beta+&
&     num_theta_sigma_s+num_theta_sigma_z_y+num_theta_rho+1:)
  END SUBROUTINE COMPUTE_PSI_PARAMETERS
!  Differentiation of compute_b_parameters_dv in forward (tangent) mode (with options multiDirectional i4 dr8 r4):
!   variations   of useful results: rho_zy qd rho_zyd rho_szd m_s
!                q m_y m_z sigma_yd sigma_s sigma_y sigma_z sigma_sd
!                delta_y delta_z delta_yd sigma_zd rho_syd delta_zd
!                rho_sy rho_sz m_zd
!   with respect to varying inputs: psi
!  Differentiation of compute_b_parameters in forward (tangent) mode (with options multiDirectional i4 dr8 r4):
!   variations   of useful results: rho_zy m_s q m_y m_z sigma_s
!                sigma_y sigma_z delta_y delta_z rho_sy rho_sz
!   with respect to varying inputs: psi
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
  SUBROUTINE COMPUTE_B_PARAMETERS_DV_DV(psi, psid0, psid, m_s, m_sd0, &
&   m_sd, m_q, m_z, m_zd0, m_zd, m_zdd, m_y, m_yd0, m_yd, sigma_s, &
&   sigma_sd0, sigma_sd, sigma_sdd, sigma_z, sigma_zd0, sigma_zd, &
&   sigma_zdd, sigma_y, sigma_yd0, sigma_yd, sigma_ydd, rho_sz, rho_szd0&
&   , rho_szd, rho_szdd, rho_sy, rho_syd0, rho_syd, rho_sydd, rho_zy, &
&   rho_zyd0, rho_zyd, rho_zydd, gamma, kappa, q, qd0, qd, qdd, delta_z&
&   , delta_zd0, delta_zd, delta_zdd, delta_y, delta_yd0, delta_yd, &
&   delta_ydd, nbdirs, nbdirs0)
    USE DIFFSIZES
    IMPLICIT NONE
!
!
! Ending execution and returning control
!
!
! Declaring dummy variables
!
    REAL*8, INTENT(IN) :: psi(num_psi)
    REAL*8, INTENT(IN) :: psid0(nbdirsmax0, num_psi)
    REAL*8, INTENT(IN) :: psid(nbdirsmax, num_psi)
    REAL*8, INTENT(OUT) :: m_s
    REAL*8, DIMENSION(nbdirsmax0), INTENT(OUT) :: m_sd0
    REAL*8, DIMENSION(nbdirsmax), INTENT(OUT) :: m_sd
    REAL*8, INTENT(OUT) :: m_q
    REAL*8, INTENT(OUT) :: m_z
    REAL*8, DIMENSION(nbdirsmax0), INTENT(OUT) :: m_zd0
    REAL*8, DIMENSION(nbdirsmax), INTENT(OUT) :: m_zd
    REAL*8, DIMENSION(nbdirsmax0, nbdirsmax), INTENT(OUT) :: m_zdd
    REAL*8, INTENT(OUT) :: m_y
    REAL*8, DIMENSION(nbdirsmax0), INTENT(OUT) :: m_yd0
    REAL*8, DIMENSION(nbdirsmax), INTENT(OUT) :: m_yd
    REAL*8, INTENT(OUT) :: sigma_s
    REAL*8, DIMENSION(nbdirsmax0), INTENT(OUT) :: sigma_sd0
    REAL*8, DIMENSION(nbdirsmax), INTENT(OUT) :: sigma_sd
    REAL*8, DIMENSION(nbdirsmax0, nbdirsmax), INTENT(OUT) :: sigma_sdd
    REAL*8, INTENT(OUT) :: sigma_z
    REAL*8, DIMENSION(nbdirsmax0), INTENT(OUT) :: sigma_zd0
    REAL*8, DIMENSION(nbdirsmax), INTENT(OUT) :: sigma_zd
    REAL*8, DIMENSION(nbdirsmax0, nbdirsmax), INTENT(OUT) :: sigma_zdd
    REAL*8, INTENT(OUT) :: sigma_y
    REAL*8, DIMENSION(nbdirsmax0), INTENT(OUT) :: sigma_yd0
    REAL*8, DIMENSION(nbdirsmax), INTENT(OUT) :: sigma_yd
    REAL*8, DIMENSION(nbdirsmax0, nbdirsmax), INTENT(OUT) :: sigma_ydd
    REAL*8, INTENT(OUT) :: rho_sz
    REAL*8, DIMENSION(nbdirsmax0), INTENT(OUT) :: rho_szd0
    REAL*8, DIMENSION(nbdirsmax), INTENT(OUT) :: rho_szd
    REAL*8, DIMENSION(nbdirsmax0, nbdirsmax), INTENT(OUT) :: rho_szdd
    REAL*8, INTENT(OUT) :: rho_sy
    REAL*8, DIMENSION(nbdirsmax0), INTENT(OUT) :: rho_syd0
    REAL*8, DIMENSION(nbdirsmax), INTENT(OUT) :: rho_syd
    REAL*8, DIMENSION(nbdirsmax0, nbdirsmax), INTENT(OUT) :: rho_sydd
    REAL*8, INTENT(OUT) :: rho_zy
    REAL*8, DIMENSION(nbdirsmax0), INTENT(OUT) :: rho_zyd0
    REAL*8, DIMENSION(nbdirsmax), INTENT(OUT) :: rho_zyd
    REAL*8, DIMENSION(nbdirsmax0, nbdirsmax), INTENT(OUT) :: rho_zydd
    REAL*8, INTENT(OUT) :: gamma
    REAL*8, DIMENSION(nbdirsmax0) :: gammad0
    REAL*8, DIMENSION(nbdirsmax) :: gammad
    REAL*8, DIMENSION(nbdirsmax0, nbdirsmax) :: gammadd
    REAL*8, INTENT(OUT) :: kappa
    REAL*8, DIMENSION(nbdirsmax0) :: kappad0
    REAL*8, DIMENSION(nbdirsmax) :: kappad
    REAL*8, DIMENSION(nbdirsmax0, nbdirsmax) :: kappadd
    REAL*8, INTENT(OUT) :: q
    REAL*8, DIMENSION(nbdirsmax0), INTENT(OUT) :: qd0
    REAL*8, DIMENSION(nbdirsmax), INTENT(OUT) :: qd
    REAL*8, DIMENSION(nbdirsmax0, nbdirsmax), INTENT(OUT) :: qdd
    REAL*8, INTENT(OUT) :: delta_z(num_delta_z)
    REAL*8, INTENT(OUT) :: delta_zd0(nbdirsmax0, num_delta_z)
    REAL*8, INTENT(OUT) :: delta_zd(nbdirsmax, num_delta_z)
    REAL*8, INTENT(OUT) :: delta_zdd(nbdirsmax0, nbdirsmax, num_delta_z)
    REAL*8, INTENT(OUT) :: delta_y(num_delta_y)
    REAL*8, INTENT(OUT) :: delta_yd0(nbdirsmax0, num_delta_y)
    REAL*8, INTENT(OUT) :: delta_yd(nbdirsmax, num_delta_y)
    REAL*8, INTENT(OUT) :: delta_ydd(nbdirsmax0, nbdirsmax, num_delta_y)
!
! Declaring local variables
!
    INTEGER :: ipsi, il, iu, i
    INTRINSIC EXP
    INTRINSIC LOG
    INTRINSIC ATAN
    INTRINSIC SQRT
    REAL*8 :: arg1
    REAL*8, DIMENSION(nbdirsmax0) :: arg1d0
    REAL*8, DIMENSION(nbdirsmax) :: arg1d
    REAL*8, DIMENSION(nbdirsmax0, nbdirsmax) :: arg1dd
    REAL*8 :: result1
    REAL*8, DIMENSION(nbdirsmax0) :: result1d0
    REAL*8, DIMENSION(nbdirsmax) :: result1d
    REAL*8, DIMENSION(nbdirsmax0, nbdirsmax) :: result1dd
    INTEGER :: nd
    INTEGER :: nbdirs
    REAL*8 :: result10
    REAL*8, DIMENSION(nbdirsmax0) :: result10d
    INTEGER :: nd0
    INTEGER :: nbdirs0
    q = EXP(psi(2))
    kappa = EXP(psi(3))
    DO nd0=1,nbdirs0
      qd0(nd0) = psid0(nd0, 2)*EXP(psi(2))
      kappad0(nd0) = psid0(nd0, 3)*EXP(psi(3))
      gammad0(nd0) = ((qd0(nd0)+kappad0(nd0)/apriori_t)*(1.d0+kappa/&
&       apriori_t)-(q+kappa/apriori_t)*kappad0(nd0)/apriori_t)/(1.d0+&
&       kappa/apriori_t)**2
    END DO
    gamma = (q+kappa/apriori_t)/(1.d0+kappa/apriori_t)
    DO nd0=1,nbdirs0
      qdd(nd0, :) = 0.0_8
      sigma_sdd(nd0, :) = 0.0_8
      m_zdd(nd0, :) = 0.0_8
      kappadd(nd0, :) = 0.0_8
      gammadd(nd0, :) = 0.0_8
    END DO
    DO nd=1,nbdirs
!
! Beginning 
!
! Beginning execution
! 
! means
!
      m_sd(nd) = psid(nd, 1)
      m_yd(nd) = psid(nd, 3)
      qd(nd) = psid(nd, 2)*EXP(psi(2))
      kappad(nd) = psid(nd, 3)*EXP(psi(3))
      gammad(nd) = ((qd(nd)+kappad(nd)/apriori_t)*(1.d0+kappa/apriori_t)&
&       -(q+kappa/apriori_t)*kappad(nd)/apriori_t)/(1.d0+kappa/apriori_t&
&       )**2
      DO nd0=1,nbdirs0
!
! risk aversion and planning horizon
!
        qdd(nd0, nd) = psid(nd, 2)*psid0(nd0, 2)*EXP(psi(2))
        kappadd(nd0, nd) = psid(nd, 3)*psid0(nd0, 3)*EXP(psi(3))
        gammadd(nd0, nd) = (((qdd(nd0, nd)+kappadd(nd0, nd)/apriori_t)*(&
&         1.d0+kappa/apriori_t)+(qd(nd)+kappad(nd)/apriori_t)*kappad0(&
&         nd0)/apriori_t-((qd0(nd0)+kappad0(nd0)/apriori_t)*kappad(nd)+(&
&         q+kappa/apriori_t)*kappadd(nd0, nd))/apriori_t)*(1.d0+kappa/&
&         apriori_t)**2-((qd(nd)+kappad(nd)/apriori_t)*(1.d0+kappa/&
&         apriori_t)-(q+kappa/apriori_t)*kappad(nd)/apriori_t)*2*(1.d0+&
&         kappa/apriori_t)*kappad0(nd0)/apriori_t)/((1.d0+kappa/&
&         apriori_t)**2)**2
        m_zdd(nd0, nd) = (gammadd(nd0, nd)*gamma-gammad(nd)*gammad0(nd0)&
&         )/gamma**2
!
! sigma_s
!
        sigma_sdd(nd0, nd) = psid(nd, num_psi_m+1)*psid0(nd0, num_psi_m+&
&         1)*EXP(psi(num_psi_m+1))
      END DO
      m_zd(nd) = gammad(nd)/gamma
      sigma_sd(nd) = psid(nd, num_psi_m+1)*EXP(psi(num_psi_m+1))
    END DO
    DO nd0=1,nbdirs0
      m_sd0(nd0) = psid0(nd0, 1)
      m_yd0(nd0) = psid0(nd0, 3)
      m_zd0(nd0) = gammad0(nd0)/gamma
      sigma_sd0(nd0) = psid0(nd0, num_psi_m+1)*EXP(psi(num_psi_m+1))
    END DO
    m_s = psi(1)
    m_q = psi(2)
    m_y = psi(3)
    m_z = LOG(gamma)
    sigma_s = minimum_sigma_s+EXP(psi(num_psi_m+1))
    ipsi = num_psi_m + 1
!
! sigma_z
!
    IF (switch_sigma_z .EQ. 1) THEN
!
      ipsi = ipsi + 1
      DO nd0=1,nbdirs0
        sigma_zdd(nd0, :) = 0.0_8
      END DO
      DO nd=1,nbdirs
        DO nd0=1,nbdirs0
          sigma_zdd(nd0, nd) = psid(nd, ipsi)*psid0(nd0, ipsi)*EXP(psi(&
&           ipsi))
        END DO
        sigma_zd(nd) = psid(nd, ipsi)*EXP(psi(ipsi))
      END DO
      DO nd0=1,nbdirs0
        sigma_zd0(nd0) = psid0(nd0, ipsi)*EXP(psi(ipsi))
      END DO
      sigma_z = EXP(psi(ipsi))
!
    ELSE IF (switch_sigma_z .EQ. 0) THEN
!
      sigma_z = 1.d0
!
      DO nd=1,nbdirs
        sigma_zd(nd) = 0.0_8
      END DO
      DO nd0=1,nbdirs0
        sigma_zd0(nd0) = 0.0_8
        sigma_zdd(nd0, :) = 0.0_8
      END DO
    ELSE
      DO nd=1,nbdirs
        sigma_zd(nd) = 0.0_8
      END DO
      DO nd0=1,nbdirs0
        sigma_zd0(nd0) = 0.0_8
        sigma_zdd(nd0, :) = 0.0_8
      END DO
    END IF
!
! sigma_y
!
    IF (switch_sigma_y .EQ. 1) THEN
!
      ipsi = ipsi + 1
      DO nd0=1,nbdirs0
        sigma_ydd(nd0, :) = 0.0_8
      END DO
      DO nd=1,nbdirs
        DO nd0=1,nbdirs0
          sigma_ydd(nd0, nd) = psid(nd, ipsi)*psid0(nd0, ipsi)*EXP(psi(&
&           ipsi))
        END DO
        sigma_yd(nd) = psid(nd, ipsi)*EXP(psi(ipsi))
      END DO
      DO nd0=1,nbdirs0
        sigma_yd0(nd0) = psid0(nd0, ipsi)*EXP(psi(ipsi))
      END DO
      sigma_y = EXP(psi(ipsi))
!
    ELSE IF (switch_sigma_y .EQ. 0) THEN
!
      sigma_y = 1.d0
!
      DO nd=1,nbdirs
        sigma_yd(nd) = 0.0_8
      END DO
      DO nd0=1,nbdirs0
        sigma_ydd(nd0, :) = 0.0_8
        sigma_yd0(nd0) = 0.0_8
      END DO
    ELSE
      DO nd=1,nbdirs
        sigma_yd(nd) = 0.0_8
      END DO
      DO nd0=1,nbdirs0
        sigma_ydd(nd0, :) = 0.0_8
        sigma_yd0(nd0) = 0.0_8
      END DO
    END IF
!
! rho_sz
!
    IF (switch_rho_sz .EQ. 1) THEN
!
      ipsi = ipsi + 1
      DO nd0=1,nbdirs0
        rho_szdd(nd0, :) = 0.0_8
      END DO
      DO nd=1,nbdirs
        DO nd0=1,nbdirs0
          rho_szdd(nd0, nd) = -(twooverpi*psid(nd, ipsi)*2*psi(ipsi)*&
&           psid0(nd0, ipsi)/(1.0+psi(ipsi)**2)**2)
        END DO
        rho_szd(nd) = twooverpi*psid(nd, ipsi)/(1.0+psi(ipsi)**2)
      END DO
      DO nd0=1,nbdirs0
        rho_szd0(nd0) = twooverpi*psid0(nd0, ipsi)/(1.0+psi(ipsi)**2)
      END DO
      rho_sz = twooverpi*ATAN(psi(ipsi))
!
    ELSE IF (switch_rho_sz .EQ. 0) THEN
!
      rho_sz = 0.d0
!
      DO nd=1,nbdirs
        rho_szd(nd) = 0.0_8
      END DO
      DO nd0=1,nbdirs0
        rho_szdd(nd0, :) = 0.0_8
        rho_szd0(nd0) = 0.0_8
      END DO
    ELSE
      DO nd=1,nbdirs
        rho_szd(nd) = 0.0_8
      END DO
      DO nd0=1,nbdirs0
        rho_szdd(nd0, :) = 0.0_8
        rho_szd0(nd0) = 0.0_8
      END DO
    END IF
!
! rho_sy
!
    IF (switch_rho_sy .EQ. 1) THEN
!
      ipsi = ipsi + 1
      DO nd0=1,nbdirs0
        rho_sydd(nd0, :) = 0.0_8
      END DO
      DO nd=1,nbdirs
        DO nd0=1,nbdirs0
          rho_sydd(nd0, nd) = -(twooverpi*psid(nd, ipsi)*2*psi(ipsi)*&
&           psid0(nd0, ipsi)/(1.0+psi(ipsi)**2)**2)
        END DO
        rho_syd(nd) = twooverpi*psid(nd, ipsi)/(1.0+psi(ipsi)**2)
      END DO
      DO nd0=1,nbdirs0
        rho_syd0(nd0) = twooverpi*psid0(nd0, ipsi)/(1.0+psi(ipsi)**2)
      END DO
      rho_sy = twooverpi*ATAN(psi(ipsi))
!
    ELSE IF (switch_rho_sy .EQ. 0) THEN
!
      rho_sy = 0.d0
!
      DO nd=1,nbdirs
        rho_syd(nd) = 0.0_8
      END DO
      DO nd0=1,nbdirs0
        rho_sydd(nd0, :) = 0.0_8
        rho_syd0(nd0) = 0.0_8
      END DO
    ELSE
      DO nd=1,nbdirs
        rho_syd(nd) = 0.0_8
      END DO
      DO nd0=1,nbdirs0
        rho_sydd(nd0, :) = 0.0_8
        rho_syd0(nd0) = 0.0_8
      END DO
    END IF
!
! rho_zy
!
    IF (switch_rho_zy .EQ. 1) THEN
!
      ipsi = ipsi + 1
      arg1 = (1.d0-rho_sz**2)*(1.d0-rho_sy**2)
      DO nd0=1,nbdirs0
        arg1d0(nd0) = -(2*rho_sz*rho_szd0(nd0)*(1.d0-rho_sy**2)) - (1.d0&
&         -rho_sz**2)*2*rho_sy*rho_syd0(nd0)
        IF (arg1 .EQ. 0.0) THEN
          result1d0(nd0) = 0.0_8
        ELSE
          result1d0(nd0) = arg1d0(nd0)/(2.0*SQRT(arg1))
        END IF
      END DO
      result1 = SQRT(arg1)
      DO nd0=1,nbdirs0
        rho_zydd(nd0, :) = 0.0_8
        arg1dd(nd0, :) = 0.0_8
        result1dd(nd0, :) = 0.0_8
      END DO
      DO nd=1,nbdirs
        DO nd0=1,nbdirs0
          arg1dd(nd0, nd) = 2*(rho_sz*rho_szd(nd)*2*rho_sy*rho_syd0(nd0)&
&           ) - 2*((rho_szd0(nd0)*rho_szd(nd)+rho_sz*rho_szdd(nd0, nd))*&
&           (1.d0-rho_sy**2)) - 2*((1.d0-rho_sz**2)*(rho_syd0(nd0)*&
&           rho_syd(nd)+rho_sy*rho_sydd(nd0, nd))) + 2*(2*rho_sz*&
&           rho_szd0(nd0)*rho_sy*rho_syd(nd))
        END DO
        arg1d(nd) = -(2*rho_sz*rho_szd(nd)*(1.d0-rho_sy**2)) - (1.d0-&
&         rho_sz**2)*2*rho_sy*rho_syd(nd)
        IF (arg1 .EQ. 0.0) THEN
          DO nd0=1,nbdirs0
            result1dd(nd0, nd) = 0.0_8
          END DO
          result1d(nd) = 0.0_8
        ELSE
          result10 = SQRT(arg1)
          DO nd0=1,nbdirs0
            IF (arg1 .EQ. 0.0) THEN
              result10d(nd0) = 0.0_8
            ELSE
              result10d(nd0) = arg1d0(nd0)/(2.0*SQRT(arg1))
            END IF
            result1dd(nd0, nd) = (arg1dd(nd0, nd)*2.0*result10-arg1d(nd)&
&             *2.0*result10d(nd0))/(2.0*result10)**2
          END DO
          result1d(nd) = arg1d(nd)/(2.0*result10)
        END IF
        DO nd0=1,nbdirs0
          rho_zydd(nd0, nd) = rho_szdd(nd0, nd)*rho_sy + rho_szd(nd)*&
&           rho_syd0(nd0) + rho_szd0(nd0)*rho_syd(nd) + rho_sz*rho_sydd(&
&           nd0, nd) + twooverpi*(result1dd(nd0, nd)*ATAN(psi(ipsi))+&
&           result1d(nd)*psid0(nd0, ipsi)/(1.0+psi(ipsi)**2)+(psid(nd, &
&           ipsi)*result1d0(nd0)*(1.0+psi(ipsi)**2)-result1*psid(nd, &
&           ipsi)*2*psi(ipsi)*psid0(nd0, ipsi))/(1.0+psi(ipsi)**2)**2)
        END DO
        rho_zyd(nd) = rho_szd(nd)*rho_sy + rho_sz*rho_syd(nd) + &
&         twooverpi*(result1d(nd)*ATAN(psi(ipsi))+result1*psid(nd, ipsi)&
&         /(1.0+psi(ipsi)**2))
      END DO
      DO nd0=1,nbdirs0
        rho_zyd0(nd0) = rho_szd0(nd0)*rho_sy + rho_sz*rho_syd0(nd0) + &
&         twooverpi*(result1d0(nd0)*ATAN(psi(ipsi))+result1*psid0(nd0, &
&         ipsi)/(1.0+psi(ipsi)**2))
      END DO
      rho_zy = rho_sz*rho_sy + result1*twooverpi*ATAN(psi(ipsi))
!
    ELSE
      DO nd=1,nbdirs
        rho_zyd(nd) = 0.0_8
      END DO
      DO nd0=1,nbdirs0
        rho_zyd0(nd0) = 0.0_8
        rho_zydd(nd0, :) = 0.0_8
      END DO
    END IF
!
! delta_z
!
    il = ipsi
    delta_z(1) = 0.d0
    DO nd0=1,nbdirs0
      delta_zdd(nd0, :, :) = 0.0_8
    END DO
    DO nd=1,nbdirs
      DO nd0=1,nbdirs0
        delta_zdd(nd0, nd, :) = 0.0_8
        delta_zdd(nd0, nd, 2:) = psid(nd, il+1:il+num_theta_delta_z)*&
&         psid0(nd0, il+1:il+num_theta_delta_z)*EXP(psi(il+1:il+&
&         num_theta_delta_z))
      END DO
      delta_zd(nd, :) = 0.0_8
      delta_zd(nd, 2:) = psid(nd, il+1:il+num_theta_delta_z)*EXP(psi(il+&
&       1:il+num_theta_delta_z))
    END DO
    DO nd0=1,nbdirs0
      delta_zd0(nd0, :) = 0.0_8
      delta_zd0(nd0, 2:) = psid0(nd0, il+1:il+num_theta_delta_z)*EXP(psi&
&       (il+1:il+num_theta_delta_z))
    END DO
    delta_z(2:) = EXP(psi(il+1:il+num_theta_delta_z))
    DO i=3,num_delta_z
      DO nd=1,nbdirs
        DO nd0=1,nbdirs0
!
          delta_zdd(nd0, nd, i) = delta_zdd(nd0, nd, i) + delta_zdd(nd0&
&           , nd, i-1)
        END DO
        delta_zd(nd, i) = delta_zd(nd, i) + delta_zd(nd, i-1)
      END DO
      DO nd0=1,nbdirs0
        delta_zd0(nd0, i) = delta_zd0(nd0, i) + delta_zd0(nd0, i-1)
      END DO
      delta_z(i) = delta_z(i) + delta_z(i-1)
    END DO
!
!
! delta_y
!
    il = il + num_theta_delta_z
    delta_y(1) = 0.d0
    DO nd0=1,nbdirs0
      delta_ydd(nd0, :, :) = 0.0_8
    END DO
    DO nd=1,nbdirs
      DO nd0=1,nbdirs0
        delta_ydd(nd0, nd, :) = 0.0_8
        delta_ydd(nd0, nd, 2:) = psid(nd, il+1:il+num_theta_delta_y)*&
&         psid0(nd0, il+1:il+num_theta_delta_y)*EXP(psi(il+1:il+&
&         num_theta_delta_y))
      END DO
      delta_yd(nd, :) = 0.0_8
      delta_yd(nd, 2:) = psid(nd, il+1:il+num_theta_delta_y)*EXP(psi(il+&
&       1:il+num_theta_delta_y))
    END DO
    DO nd0=1,nbdirs0
      delta_yd0(nd0, :) = 0.0_8
      delta_yd0(nd0, 2:) = psid0(nd0, il+1:il+num_theta_delta_y)*EXP(psi&
&       (il+1:il+num_theta_delta_y))
    END DO
    delta_y(2:) = EXP(psi(il+1:il+num_theta_delta_y))
    DO i=3,num_delta_y
      DO nd=1,nbdirs
        DO nd0=1,nbdirs0
!
          delta_ydd(nd0, nd, i) = delta_ydd(nd0, nd, i) + delta_ydd(nd0&
&           , nd, i-1)
        END DO
        delta_yd(nd, i) = delta_yd(nd, i) + delta_yd(nd, i-1)
      END DO
      DO nd0=1,nbdirs0
        delta_yd0(nd0, i) = delta_yd0(nd0, i) + delta_yd0(nd0, i-1)
      END DO
      delta_y(i) = delta_y(i) + delta_y(i-1)
    END DO
  END SUBROUTINE COMPUTE_B_PARAMETERS_DV_DV
!  Differentiation of compute_b_parameters in forward (tangent) mode (with options multiDirectional i4 dr8 r4):
!   variations   of useful results: rho_zy m_s q m_y m_z sigma_s
!                sigma_y sigma_z delta_y delta_z rho_sy rho_sz
!   with respect to varying inputs: psi
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
  SUBROUTINE COMPUTE_B_PARAMETERS_DV(psi, psid, m_s, m_sd, m_q, m_z, &
&   m_zd, m_y, m_yd, sigma_s, sigma_sd, sigma_z, sigma_zd, sigma_y, &
&   sigma_yd, rho_sz, rho_szd, rho_sy, rho_syd, rho_zy, rho_zyd, gamma, &
&   kappa, q, qd, delta_z, delta_zd, delta_y, delta_yd, nbdirs)
    USE DIFFSIZES
    IMPLICIT NONE
!
!
! Ending execution and returning control
!
!
! Declaring dummy variables
!
    REAL*8, INTENT(IN) :: psi(num_psi)
    REAL*8, INTENT(IN) :: psid(nbdirsmax, num_psi)
    REAL*8, INTENT(OUT) :: m_s
    REAL*8, DIMENSION(nbdirsmax), INTENT(OUT) :: m_sd
    REAL*8, INTENT(OUT) :: m_q
    REAL*8, INTENT(OUT) :: m_z
    REAL*8, DIMENSION(nbdirsmax), INTENT(OUT) :: m_zd
    REAL*8, INTENT(OUT) :: m_y
    REAL*8, DIMENSION(nbdirsmax), INTENT(OUT) :: m_yd
    REAL*8, INTENT(OUT) :: sigma_s
    REAL*8, DIMENSION(nbdirsmax), INTENT(OUT) :: sigma_sd
    REAL*8, INTENT(OUT) :: sigma_z
    REAL*8, DIMENSION(nbdirsmax), INTENT(OUT) :: sigma_zd
    REAL*8, INTENT(OUT) :: sigma_y
    REAL*8, DIMENSION(nbdirsmax), INTENT(OUT) :: sigma_yd
    REAL*8, INTENT(OUT) :: rho_sz
    REAL*8, DIMENSION(nbdirsmax), INTENT(OUT) :: rho_szd
    REAL*8, INTENT(OUT) :: rho_sy
    REAL*8, DIMENSION(nbdirsmax), INTENT(OUT) :: rho_syd
    REAL*8, INTENT(OUT) :: rho_zy
    REAL*8, DIMENSION(nbdirsmax), INTENT(OUT) :: rho_zyd
    REAL*8, INTENT(OUT) :: gamma
    REAL*8, DIMENSION(nbdirsmax) :: gammad
    REAL*8, INTENT(OUT) :: kappa
    REAL*8, DIMENSION(nbdirsmax) :: kappad
    REAL*8, INTENT(OUT) :: q
    REAL*8, DIMENSION(nbdirsmax), INTENT(OUT) :: qd
    REAL*8, INTENT(OUT) :: delta_z(num_delta_z)
    REAL*8, INTENT(OUT) :: delta_zd(nbdirsmax, num_delta_z)
    REAL*8, INTENT(OUT) :: delta_y(num_delta_y)
    REAL*8, INTENT(OUT) :: delta_yd(nbdirsmax, num_delta_y)
!
! Declaring local variables
!
    INTEGER :: ipsi, il, iu, i
    INTRINSIC EXP
    INTRINSIC LOG
    INTRINSIC ATAN
    INTRINSIC SQRT
    REAL*8 :: arg1
    REAL*8, DIMENSION(nbdirsmax) :: arg1d
    REAL*8 :: result1
    REAL*8, DIMENSION(nbdirsmax) :: result1d
    INTEGER :: nd
    INTEGER :: nbdirs
    REAL*8 :: result10
    q = EXP(psi(2))
    kappa = EXP(psi(3))
    gamma = (q+kappa/apriori_t)/(1.d0+kappa/apriori_t)
    DO nd=1,nbdirs
!
! Beginning 
!
! Beginning execution
! 
! means
!
      m_sd(nd) = psid(nd, 1)
      m_yd(nd) = psid(nd, 3)
!
! risk aversion and planning horizon
!
      qd(nd) = psid(nd, 2)*EXP(psi(2))
      kappad(nd) = psid(nd, 3)*EXP(psi(3))
      gammad(nd) = ((qd(nd)+kappad(nd)/apriori_t)*(1.d0+kappa/apriori_t)&
&       -(q+kappa/apriori_t)*kappad(nd)/apriori_t)/(1.d0+kappa/apriori_t&
&       )**2
      m_zd(nd) = gammad(nd)/gamma
!
! sigma_s
!
      sigma_sd(nd) = psid(nd, num_psi_m+1)*EXP(psi(num_psi_m+1))
    END DO
    m_s = psi(1)
    m_q = psi(2)
    m_y = psi(3)
    m_z = LOG(gamma)
    sigma_s = minimum_sigma_s+EXP(psi(num_psi_m+1))
    ipsi = num_psi_m + 1
!
! sigma_z
!
    IF (switch_sigma_z .EQ. 1) THEN
!
      ipsi = ipsi + 1
      DO nd=1,nbdirs
        sigma_zd(nd) = psid(nd, ipsi)*EXP(psi(ipsi))
      END DO
      sigma_z = EXP(psi(ipsi))
!
    ELSE IF (switch_sigma_z .EQ. 0) THEN
!
      sigma_z = 1.d0
!
      DO nd=1,nbdirs
        sigma_zd(nd) = 0.0_8
      END DO
    ELSE
      DO nd=1,nbdirs
        sigma_zd(nd) = 0.0_8
      END DO
    END IF
!
! sigma_y
!
    IF (switch_sigma_y .EQ. 1) THEN
!
      ipsi = ipsi + 1
      DO nd=1,nbdirs
        sigma_yd(nd) = psid(nd, ipsi)*EXP(psi(ipsi))
      END DO
      sigma_y = EXP(psi(ipsi))
!
    ELSE IF (switch_sigma_y .EQ. 0) THEN
!
      sigma_y = 1.d0
!
      DO nd=1,nbdirs
        sigma_yd(nd) = 0.0_8
      END DO
    ELSE
      DO nd=1,nbdirs
        sigma_yd(nd) = 0.0_8
      END DO
    END IF
!
! rho_sz
!
    IF (switch_rho_sz .EQ. 1) THEN
!
      ipsi = ipsi + 1
      DO nd=1,nbdirs
        rho_szd(nd) = twooverpi*psid(nd, ipsi)/(1.0+psi(ipsi)**2)
      END DO
      rho_sz = twooverpi*ATAN(psi(ipsi))
!
    ELSE IF (switch_rho_sz .EQ. 0) THEN
!
      rho_sz = 0.d0
!
      DO nd=1,nbdirs
        rho_szd(nd) = 0.0_8
      END DO
    ELSE
      DO nd=1,nbdirs
        rho_szd(nd) = 0.0_8
      END DO
    END IF
!
! rho_sy
!
    IF (switch_rho_sy .EQ. 1) THEN
!
      ipsi = ipsi + 1
      DO nd=1,nbdirs
        rho_syd(nd) = twooverpi*psid(nd, ipsi)/(1.0+psi(ipsi)**2)
      END DO
      rho_sy = twooverpi*ATAN(psi(ipsi))
!
    ELSE IF (switch_rho_sy .EQ. 0) THEN
!
      rho_sy = 0.d0
!
      DO nd=1,nbdirs
        rho_syd(nd) = 0.0_8
      END DO
    ELSE
      DO nd=1,nbdirs
        rho_syd(nd) = 0.0_8
      END DO
    END IF
!
! rho_zy
!
    IF (switch_rho_zy .EQ. 1) THEN
!
      ipsi = ipsi + 1
      arg1 = (1.d0-rho_sz**2)*(1.d0-rho_sy**2)
      result1 = SQRT(arg1)
      DO nd=1,nbdirs
        arg1d(nd) = -(2*rho_sz*rho_szd(nd)*(1.d0-rho_sy**2)) - (1.d0-&
&         rho_sz**2)*2*rho_sy*rho_syd(nd)
        IF (arg1 .EQ. 0.0) THEN
          result1d(nd) = 0.0_8
        ELSE
          result10 = SQRT(arg1)
          result1d(nd) = arg1d(nd)/(2.0*result10)
        END IF
        rho_zyd(nd) = rho_szd(nd)*rho_sy + rho_sz*rho_syd(nd) + &
&         twooverpi*(result1d(nd)*ATAN(psi(ipsi))+result1*psid(nd, ipsi)&
&         /(1.0+psi(ipsi)**2))
      END DO
      rho_zy = rho_sz*rho_sy + result1*twooverpi*ATAN(psi(ipsi))
!
    ELSE
      DO nd=1,nbdirs
        rho_zyd(nd) = 0.0_8
      END DO
    END IF
!
! delta_z
!
    il = ipsi
    delta_z(1) = 0.d0
    DO nd=1,nbdirs
      delta_zd(nd, :) = 0.0_8
      delta_zd(nd, 2:) = psid(nd, il+1:il+num_theta_delta_z)*EXP(psi(il+&
&       1:il+num_theta_delta_z))
    END DO
    delta_z(2:) = EXP(psi(il+1:il+num_theta_delta_z))
    DO i=3,num_delta_z
      DO nd=1,nbdirs
!
        delta_zd(nd, i) = delta_zd(nd, i) + delta_zd(nd, i-1)
      END DO
      delta_z(i) = delta_z(i) + delta_z(i-1)
    END DO
!
!
! delta_y
!
    il = il + num_theta_delta_z
    delta_y(1) = 0.d0
    DO nd=1,nbdirs
      delta_yd(nd, :) = 0.0_8
      delta_yd(nd, 2:) = psid(nd, il+1:il+num_theta_delta_y)*EXP(psi(il+&
&       1:il+num_theta_delta_y))
    END DO
    delta_y(2:) = EXP(psi(il+1:il+num_theta_delta_y))
    DO i=3,num_delta_y
      DO nd=1,nbdirs
!
        delta_yd(nd, i) = delta_yd(nd, i) + delta_yd(nd, i-1)
      END DO
      delta_y(i) = delta_y(i) + delta_y(i-1)
    END DO
  END SUBROUTINE COMPUTE_B_PARAMETERS_DV
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
  SUBROUTINE COMPUTE_B_PARAMETERS(psi, m_s, m_q, m_z, m_y, sigma_s, &
&   sigma_z, sigma_y, rho_sz, rho_sy, rho_zy, gamma, kappa, q, delta_z, &
&   delta_y)
    IMPLICIT NONE
!
!
! Ending execution and returning control
!
!
! Declaring dummy variables
!
    REAL*8, INTENT(IN) :: psi(num_psi)
    REAL*8, INTENT(OUT) :: m_s
    REAL*8, INTENT(OUT) :: m_q
    REAL*8, INTENT(OUT) :: m_z
    REAL*8, INTENT(OUT) :: m_y
    REAL*8, INTENT(OUT) :: sigma_s
    REAL*8, INTENT(OUT) :: sigma_z
    REAL*8, INTENT(OUT) :: sigma_y
    REAL*8, INTENT(OUT) :: rho_sz
    REAL*8, INTENT(OUT) :: rho_sy
    REAL*8, INTENT(OUT) :: rho_zy
    REAL*8, INTENT(OUT) :: gamma
    REAL*8, INTENT(OUT) :: kappa
    REAL*8, INTENT(OUT) :: q
    REAL*8, INTENT(OUT) :: delta_z(num_delta_z)
    REAL*8, INTENT(OUT) :: delta_y(num_delta_y)
!
! Declaring local variables
!
    INTEGER :: ipsi, il, iu, i
    INTRINSIC EXP
    INTRINSIC LOG
    INTRINSIC ATAN
    INTRINSIC SQRT
    REAL*8 :: arg1
    REAL*8 :: result1
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
    gamma = (q+kappa/apriori_t)/(1.d0+kappa/apriori_t)
    m_z = LOG(gamma)
!
! sigma_s
!
    sigma_s = minimum_sigma_s+EXP(psi(num_psi_m+1))
    ipsi = num_psi_m + 1
!
! sigma_z
!
    IF (switch_sigma_z .EQ. 1) THEN
!
      ipsi = ipsi + 1
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
      ipsi = ipsi + 1
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
      ipsi = ipsi + 1
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
      ipsi = ipsi + 1
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
      ipsi = ipsi + 1
      arg1 = (1.d0-rho_sz**2)*(1.d0-rho_sy**2)
      result1 = SQRT(arg1)
      rho_zy = rho_sz*rho_sy + result1*twooverpi*ATAN(psi(ipsi))
!
    END IF
!
! delta_z
!
    il = ipsi
    delta_z(1) = 0.d0
    delta_z(2:) = EXP(psi(il+1:il+num_theta_delta_z))
    DO i=3,num_delta_z
!
      delta_z(i) = delta_z(i) + delta_z(i-1)
    END DO
!
!
! delta_y
!
    il = il + num_theta_delta_z
    delta_y(1) = 0.d0
    delta_y(2:) = EXP(psi(il+1:il+num_theta_delta_y))
    DO i=3,num_delta_y
!
      delta_y(i) = delta_y(i) + delta_y(i-1)
    END DO
  END SUBROUTINE COMPUTE_B_PARAMETERS
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
  SUBROUTINE COMPUTE_THETAF_PARAMETERS(theta, thetaf)
    IMPLICIT NONE
!
!
! Ending execution and returning control
!
!
! Declaring dummy variables
!
    REAL*8, INTENT(IN) :: theta(num_theta)
    REAL*8, INTENT(OUT) :: thetaf(num_theta)
!
! Declaring local variables
!
    INTEGER :: i, iu
    REAL*8 :: rho_sz, rho_sy
    INTRINSIC EXP
    INTRINSIC ATAN
    INTRINSIC SQRT
    REAL*8 :: arg1
    REAL*8 :: result1
! 
! Beginning execution
! 
! thetaF = theta for the parameters in the means and in sigma_s
!
    DO i=1,num_theta_beta+num_theta_sigma_s
!
      thetaf(i) = theta(i)
    END DO
!
    iu = num_theta_beta + num_theta_sigma_s
!
! sigma_z
!
    IF (switch_sigma_z .EQ. 1) THEN
!
      iu = iu + 1
      thetaf(iu) = EXP(theta(iu))
!
    END IF
!
! sigma_y
!
    IF (switch_sigma_y .EQ. 1) THEN
!
      iu = iu + 1
      thetaf(iu) = EXP(theta(iu))
!
    END IF
!
! rho_sz
!
    IF (switch_rho_sz .EQ. 1) THEN
!
      iu = iu + 1
      rho_sz = twooverpi*ATAN(theta(iu))
      thetaf(iu) = rho_sz
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
      iu = iu + 1
      rho_sy = twooverpi*ATAN(theta(iu))
      thetaf(iu) = rho_sy
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
      iu = iu + 1
      arg1 = (1.d0-rho_sz**2)*(1.d0-rho_sy**2)
      result1 = SQRT(arg1)
      thetaf(iu) = rho_sz*rho_sy + result1*twooverpi*ATAN(theta(iu))
!
    END IF
!
! delta_z
!
    thetaf(iu+1:iu+num_theta_delta_z) = EXP(theta(iu+1:iu+&
&     num_theta_delta_z))
    DO i=iu+2,iu+num_theta_delta_z
!
      thetaf(i) = thetaf(i) + thetaf(i-1)
    END DO
!
    iu = iu + num_theta_delta_z
!
! delta_y
!
    thetaf(iu+1:iu+num_theta_delta_y) = EXP(theta(iu+1:iu+&
&     num_theta_delta_y))
    DO i=iu+2,iu+num_theta_delta_y
!
      thetaf(i) = thetaf(i) + thetaf(i-1)
    END DO
  END SUBROUTINE COMPUTE_THETAF_PARAMETERS
! 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 
  FUNCTION MATRIX_DTHETA_DPSIPRIME(xn_m_s, xn_m_q, xn_m_y, xn_sigma_s, &
&   theta)
    IMPLICIT NONE
!
! 
! Ending execution and returning control
! 
!
! Declaring dummy variables
! 
    REAL*8, INTENT(IN) :: xn_m_s(num_x_m_s)
    REAL*8, INTENT(IN) :: xn_m_q(num_x_m_q)
    REAL*8, INTENT(IN) :: xn_m_y(num_x_m_y)
    REAL*8, INTENT(IN) :: xn_sigma_s(num_x_sigma_s)
! Deep parameters
    REAL*8, INTENT(IN) :: theta(num_theta)
!
! Declaring local variables
!
    INTEGER :: il, iu, ipsi, i
! 
! Declaring function type
! 
    REAL*8 :: matrix_dtheta_dpsiprime(num_theta, num_psi)
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
    iu = num_x_m_s
    matrix_dtheta_dpsiprime(il:iu, 1) = xn_m_s
!
! beta_m_q
!
    il = iu + 1
    iu = iu + num_x_m_q
    matrix_dtheta_dpsiprime(il:iu, 2) = xn_m_q
!
! beta_m_y
!
    il = iu + 1
    iu = iu + num_x_m_y
    matrix_dtheta_dpsiprime(il:iu, 3) = xn_m_y
!
! beta_sigma_s
!
    il = iu + 1
    iu = iu + num_x_sigma_s
    matrix_dtheta_dpsiprime(il:iu, 4) = xn_sigma_s
    ipsi = 4
!
! sigma_z
!
    IF (switch_sigma_z .EQ. 1) THEN
!
      il = iu + 1
      iu = il
      ipsi = ipsi + 1
      matrix_dtheta_dpsiprime(iu, ipsi) = 1.d0
!
    END IF
!
! sigma_y
!
    IF (switch_sigma_y .EQ. 1) THEN
!
      il = iu + 1
      iu = il
      ipsi = ipsi + 1
      matrix_dtheta_dpsiprime(iu, ipsi) = 1.d0
!
    END IF
!
! rho_sz
!
    IF (switch_rho_sz .EQ. 1) THEN
!
      il = iu + 1
      iu = il
      ipsi = ipsi + 1
      matrix_dtheta_dpsiprime(iu, ipsi) = 1.d0
!
    END IF
!
! rho_sy
!
    IF (switch_rho_sy .EQ. 1) THEN
!
      il = iu + 1
      iu = il
      ipsi = ipsi + 1
      matrix_dtheta_dpsiprime(iu, ipsi) = 1.d0
!
    END IF
!
! rho_zy
!
    IF (switch_rho_zy .EQ. 1) THEN
!
      il = iu + 1
      iu = il
      ipsi = ipsi + 1
      matrix_dtheta_dpsiprime(iu, ipsi) = 1.d0
!
    END IF
!
! delta_z
!
    DO i=1,num_theta_delta_z
!
      il = iu + 1
      iu = il
      ipsi = ipsi + 1
      matrix_dtheta_dpsiprime(iu, ipsi) = 1.d0
    END DO
!
!
! delta_y
!
    DO i=1,num_theta_delta_y
!
      il = iu + 1
      iu = il
      ipsi = ipsi + 1
      matrix_dtheta_dpsiprime(iu, ipsi) = 1.d0
    END DO
  END FUNCTION MATRIX_DTHETA_DPSIPRIME
! 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 
  FUNCTION MATRIX_DTHETAF_DTHETAPRIME(theta)
    IMPLICIT NONE
!
! 
! Ending execution and returning control
! 
! 
! Declaring dummy variables
! 
    REAL*8, INTENT(IN) :: theta(num_theta)
! 
! Declaring function type
! 
    REAL*8 :: matrix_dthetaf_dthetaprime(num_theta, num_theta)
!
! Declaring local variables
!
    INTEGER :: i, iu
    REAL*8 :: rho_sz, rho_sy, rho_zy
    INTRINSIC EXP
    INTRINSIC ATAN
    INTRINSIC SQRT
    REAL*8 :: arg1
    REAL*8 :: result1
! 
! Beginning execution
! 
! Initialising matrix_dtheta_dbprime
! 
    matrix_dthetaf_dthetaprime = 0.d0
!
! Unit diagonal block for the parameters in the means and in sigma_s
!
    DO i=1,num_theta_beta+num_theta_sigma_s
!
      matrix_dthetaf_dthetaprime(i, i) = 1.d0
    END DO
!
    iu = num_theta_beta + num_theta_sigma_s
!
! sigma_z
!
    IF (switch_sigma_z .EQ. 1) THEN
!
      iu = iu + 1
      matrix_dthetaf_dthetaprime(iu, iu) = EXP(theta(iu))
!
    END IF
!
! sigma_y
!
    IF (switch_sigma_y .EQ. 1) THEN
!
      iu = iu + 1
      matrix_dthetaf_dthetaprime(iu, iu) = EXP(theta(iu))
!
    END IF
!
! Correlation parameters
!
    IF (switch_rho_sz .EQ. 1 .AND. switch_rho_sy .EQ. 1) THEN
!
      rho_sz = twooverpi*ATAN(theta(iu+1))
      rho_sy = twooverpi*ATAN(theta(iu+2))
      arg1 = (1.d0-rho_sz**2)*(1.d0-rho_sy**2)
      result1 = SQRT(arg1)
      rho_zy = rho_sz*rho_sy + result1*twooverpi*ATAN(theta(iu+3))
      matrix_dthetaf_dthetaprime(iu+1, iu+1) = twooverpi/(1.d0+theta(iu+&
&       1)**2)
      matrix_dthetaf_dthetaprime(iu+3, iu+1) = twooverpi/(1.d0+theta(iu+&
&       1)**2)*(rho_sy-rho_sz*(rho_zy-rho_sz*rho_sy)/(1.d0-rho_sz**2))
      matrix_dthetaf_dthetaprime(iu+2, iu+2) = twooverpi/(1.d0+theta(iu+&
&       2)**2)
      matrix_dthetaf_dthetaprime(iu+3, iu+2) = twooverpi/(1.d0+theta(iu+&
&       2)**2)*(rho_sz-rho_sy*(rho_zy-rho_sz*rho_sy)/(1.d0-rho_sy**2))
      arg1 = (1.d0-rho_sz**2)*(1.d0-rho_sy**2)
      result1 = SQRT(arg1)
      matrix_dthetaf_dthetaprime(iu+3, iu+3) = twooverpi/(1.d0+theta(iu+&
&       3)**2)*result1
      iu = iu + 3
!
    ELSE IF (switch_rho_sz .EQ. 1 .AND. switch_rho_sy .EQ. 0) THEN
!
      rho_sz = twooverpi*ATAN(theta(iu+1))
      arg1 = 1.d0 - rho_sz**2
      result1 = SQRT(arg1)
      rho_zy = result1*twooverpi*ATAN(theta(iu+2))
      matrix_dthetaf_dthetaprime(iu+1, iu+1) = twooverpi/(1.d0+theta(iu+&
&       1)**2)
      matrix_dthetaf_dthetaprime(iu+2, iu+1) = -(twooverpi/(1.d0+theta(&
&       iu+1)**2)*rho_sz*rho_zy/(1.d0-rho_sz**2))
      arg1 = 1.d0 - rho_sz**2
      result1 = SQRT(arg1)
      matrix_dthetaf_dthetaprime(iu+2, iu+2) = twooverpi/(1.d0+theta(iu+&
&       2)**2)*result1
      iu = iu + 2
!
    ELSE IF (switch_rho_sz .EQ. 0 .AND. switch_rho_sy .EQ. 1) THEN
!
      rho_sy = twooverpi*ATAN(theta(iu+1))
      arg1 = 1.d0 - rho_sy**2
      result1 = SQRT(arg1)
      rho_zy = result1*twooverpi*ATAN(theta(iu+2))
      matrix_dthetaf_dthetaprime(iu+1, iu+1) = twooverpi/(1.d0+theta(iu+&
&       1)**2)
      matrix_dthetaf_dthetaprime(iu+2, iu+1) = -(twooverpi/(1.d0+theta(&
&       iu+1)**2)*rho_sy*rho_zy/(1.d0-rho_sy**2))
      arg1 = 1.d0 - rho_sy**2
      result1 = SQRT(arg1)
      matrix_dthetaf_dthetaprime(iu+2, iu+2) = twooverpi/(1.d0+theta(iu+&
&       2)**2)*result1
      iu = iu + 2
!
    ELSE IF (switch_rho_sz .EQ. 0 .AND. switch_rho_sy .EQ. 0) THEN
!
      rho_zy = twooverpi*ATAN(theta(iu+1))
      matrix_dthetaf_dthetaprime(iu+1, iu+1) = twooverpi/(1.d0+theta(iu+&
&       1)**2)
!
    END IF
!
! delta_z
!
    DO i=iu+1,iu+num_theta_delta_z
!
      matrix_dthetaf_dthetaprime(i:iu+num_theta_delta_z, i) = EXP(theta(&
&       i))
    END DO
!
    iu = iu + num_theta_delta_z
!
! delta_y
!
    DO i=iu+1,iu+num_theta_delta_y
!
      matrix_dthetaf_dthetaprime(i:iu+num_theta_delta_y, i) = EXP(theta(&
&       i))
    END DO
  END FUNCTION MATRIX_DTHETAF_DTHETAPRIME
END MODULE UTILITIES_DIFFV_DIFFV
