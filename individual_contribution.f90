SUBROUTINE individual_contribution ( psi, dn, asn, riskavn, horizonn, p )
!
USE constants
USE observations
USE utilities
USE sp
!
IMPLICIT NONE
! 
! Declaring dummy variables
!   
REAL(8), INTENT(IN) :: psi(num_psi)         ! Vector of metaparameters
INTEGER, INTENT(IN) :: dn                   ! Observed regime dummy
REAL(8), INTENT(IN) :: asn
INTEGER, INTENT(IN) :: riskavn              ! Observed self-reported risk attitude
INTEGER, INTENT(IN) :: horizonn             ! Observed self-reported planning horizon
REAL(8), INTENT(OUT) :: p                   ! Computed regime probability
!
! Declaring external routines
!
REAL(8), EXTERNAL :: TVTL                   ! Trivariate normal cdf
REAL(8), EXTERNAL :: BVND                   ! Bivariate normal survival sdf
REAL(8), EXTERNAL :: PHID                   ! Univariate normal cdf
! 
! Declaring local variables
! 
REAL(8) :: m_s, m_q, m_z, m_y 
REAL(8) :: sigma_s, sigma_z, sigma_y
REAL(8) :: rho_sz, rho_sy, rho_zy
REAL(8) :: gamma, kappa, q
REAL(8) :: delta_z(num_delta_z), delta_y(num_delta_y)
REAL(8) :: m_zs, m_ys, psi_z, psi_y, rho_psi_zy
REAL(8) :: p_as, p_zy, eps
REAL(8) :: limit3(3,4), corr3(3), pcomp(4,2)
! 
! Beginning execution
! 
! Computing the b parameters
! 
CALL compute_b_parameters(psi,m_s,m_q,m_z,m_y, &
        sigma_s,sigma_z,sigma_y,rho_sz,rho_sy,rho_zy, &
        gamma,kappa,q,delta_z,delta_y)
pcomp = 0.d0
limit3 = 0.d0
! 
! Computing individual contribution to the loglikelihood
! 
SELECT CASE (dn)
    !
    CASE (1)   ! 0 < as < 1          
        ! 
        ! density of a_s
        ! 
        p_as = pdf_N01((q*asn-m_s)/sigma_s)/sigma_s
        !
        ! probability of self-reported risk aversion and planning horizon
        !
        m_zs = m_z+rho_sz*sigma_z/sigma_s*(m_s-q*asn)
        m_ys = m_y+rho_sy*sigma_y/sigma_s*(m_s-q*asn)
        psi_z = sigma_z*SQRT(1.d0-rho_sz**2)
        psi_y = sigma_y*SQRT(1.d0-rho_sy**2)
        rho_psi_zy = (rho_zy-rho_sz*rho_sy)/(SQRT(1.d0-rho_sz**2)*SQRT(1.d0-rho_sy**2))
        !
        IF ((riskavn .EQ. 1) .AND. (horizonn .EQ. 1)) THEN
            !
            limit3(2,1) = (delta_z(1)-m_zs)/psi_z
            limit3(3,1) = (delta_y(1)-m_ys)/psi_y
            pcomp(1,1) = BVND(-limit3(2,1),-limit3(3,1),rho_psi_zy)
            !
        ELSE IF ((riskavn .EQ. 1) .AND. (horizonn .EQ. 2)) THEN
            !
            limit3(2,1) = (delta_z(1)-m_zs)/psi_z
            limit3(3,1) = (delta_y(2)-m_ys)/psi_y
            limit3(2,3) = (delta_z(1)-m_zs)/psi_z
            limit3(3,3) = (delta_y(1)-m_ys)/psi_y
            pcomp(1,1) = BVND(-limit3(2,1),-limit3(3,1),rho_psi_zy) 
            pcomp(3,1) = BVND(-limit3(2,3),-limit3(3,3),rho_psi_zy)
            !
        ELSE IF ((riskavn .EQ. 1) .AND. (horizonn .EQ. 3)) THEN
            !
            limit3(2,1) = (delta_z(1)-m_zs)/psi_z
            limit3(3,1) = (delta_y(3)-m_ys)/psi_y
            limit3(2,3) = (delta_z(1)-m_zs)/psi_z
            limit3(3,3) = (delta_y(2)-m_ys)/psi_y
            pcomp(1,1) = BVND(-limit3(2,1),-limit3(3,1),rho_psi_zy)
            pcomp(3,1) = BVND(-limit3(2,3),-limit3(3,3),rho_psi_zy)
            !
        ELSE IF ((riskavn .EQ. 1) .AND. (horizonn .EQ. 4)) THEN
            !
            limit3(2,1) = (delta_z(1)-m_zs)/psi_z
            limit3(2,3) = (delta_z(1)-m_zs)/psi_z
            limit3(3,3) = (delta_y(3)-m_ys)/psi_y
            pcomp(1,1) = PHID(limit3(2,1)) 
            pcomp(3,1) = BVND(-limit3(2,3),-limit3(3,3),rho_psi_zy)
            !
        ELSE IF ((riskavn .EQ. 2) .AND. (horizonn .EQ. 1)) THEN
            !
            limit3(2,1) = (delta_z(2)-m_zs)/psi_z
            limit3(3,1) = (delta_y(1)-m_ys)/psi_y
            limit3(2,2) = (delta_z(1)-m_zs)/psi_z
            limit3(3,2) = (delta_y(1)-m_ys)/psi_y
            pcomp(1,1) = BVND(-limit3(2,1),-limit3(3,1),rho_psi_zy) 
            pcomp(2,1) = BVND(-limit3(2,2),-limit3(3,2),rho_psi_zy)
            !
        ELSE IF ((riskavn .EQ. 2) .AND. (horizonn .EQ. 2)) THEN
            !
            limit3(2,1) = (delta_z(2)-m_zs)/psi_z
            limit3(3,1) = (delta_y(2)-m_ys)/psi_y
            limit3(2,2) = (delta_z(1)-m_zs)/psi_z
            limit3(3,2) = (delta_y(2)-m_ys)/psi_y
            limit3(2,3) = (delta_z(2)-m_zs)/psi_z
            limit3(3,3) = (delta_y(1)-m_ys)/psi_y
            limit3(2,4) = (delta_z(1)-m_zs)/psi_z
            limit3(3,4) = (delta_y(1)-m_ys)/psi_y
            pcomp(1,1) = BVND(-limit3(2,1),-limit3(3,1),rho_psi_zy) 
            pcomp(2,1) = BVND(-limit3(2,2),-limit3(3,2),rho_psi_zy) 
            pcomp(3,1) = BVND(-limit3(2,3),-limit3(3,3),rho_psi_zy) 
            pcomp(4,1) = BVND(-limit3(2,4),-limit3(3,4),rho_psi_zy) 
            !
        ELSE IF ((riskavn .EQ. 2) .AND. (horizonn .EQ. 3)) THEN
            !
            limit3(2,1) = (delta_z(2)-m_zs)/psi_z
            limit3(3,1) = (delta_y(3)-m_ys)/psi_y
            limit3(2,2) = (delta_z(1)-m_zs)/psi_z
            limit3(3,2) = (delta_y(3)-m_ys)/psi_y
            limit3(2,3) = (delta_z(2)-m_zs)/psi_z
            limit3(3,3) = (delta_y(2)-m_ys)/psi_y
            limit3(2,4) = (delta_z(1)-m_zs)/psi_z
            limit3(3,4) = (delta_y(2)-m_ys)/psi_y
            pcomp(1,1) = BVND(-limit3(2,1),-limit3(3,1),rho_psi_zy) 
            pcomp(2,1) = BVND(-limit3(2,2),-limit3(3,2),rho_psi_zy)
            pcomp(3,1) = BVND(-limit3(2,3),-limit3(3,3),rho_psi_zy) 
            pcomp(4,1) = BVND(-limit3(2,4),-limit3(3,4),rho_psi_zy) 
            !
        ELSE IF ((riskavn .EQ. 2) .AND. (horizonn .EQ. 4)) THEN
            !
            limit3(2,1) = (delta_z(2)-m_zs)/psi_z
            limit3(2,2) = (delta_z(1)-m_zs)/psi_z
            limit3(2,3) = (delta_z(2)-m_zs)/psi_z
            limit3(3,3) = (delta_y(3)-m_ys)/psi_y
            limit3(2,4) = (delta_z(1)-m_zs)/psi_z
            limit3(3,4) = (delta_y(3)-m_ys)/psi_y
            pcomp(1,1) = PHID(limit3(2,1)) 
            pcomp(2,1) = PHID(limit3(2,2)) 
            pcomp(3,1) = BVND(-limit3(2,3),-limit3(3,3),rho_psi_zy) 
            pcomp(4,1) = BVND(-limit3(2,4),-limit3(3,4),rho_psi_zy)
            !
        ELSE IF ((riskavn .EQ. 3) .AND. (horizonn .EQ. 1)) THEN
            !
            limit3(3,1) = (delta_y(1)-m_ys)/psi_y
            limit3(2,2) = (delta_z(2)-m_zs)/psi_z
            limit3(3,2) = (delta_y(1)-m_ys)/psi_y
            pcomp(1,1) = PHID(limit3(3,1)) 
            pcomp(2,1) = BVND(-limit3(2,2),-limit3(3,2),rho_psi_zy) 
            !
        ELSE IF ((riskavn .EQ. 3) .AND. (horizonn .EQ. 2)) THEN
            !
            limit3(3,1) = (delta_y(2)-m_ys)/psi_y
            limit3(2,2) = (delta_z(2)-m_zs)/psi_z
            limit3(3,2) = (delta_y(2)-m_ys)/psi_y
            limit3(3,3) = (delta_y(1)-m_ys)/psi_y
            limit3(2,4) = (delta_z(2)-m_zs)/psi_z
            limit3(3,4) = (delta_y(1)-m_ys)/psi_y
            pcomp(1,1) = PHID(limit3(3,1)) 
            pcomp(2,1) = BVND(-limit3(2,2),-limit3(3,2),rho_psi_zy) 
            pcomp(3,1) = PHID(limit3(3,3)) 
            pcomp(4,1) = BVND(-limit3(2,4),-limit3(3,4),rho_psi_zy)
            !
        ELSE IF ((riskavn .EQ. 3) .AND. (horizonn .EQ. 3)) THEN
            !
            limit3(3,1) = (delta_y(3)-m_y)/sigma_y
            limit3(2,2) = (delta_z(2)-m_z)/sigma_z
            limit3(3,2) = (delta_y(3)-m_y)/sigma_y
            limit3(3,3) = (delta_y(2)-m_y)/sigma_y
            limit3(2,4) = (delta_z(2)-m_z)/sigma_z
            limit3(3,4) = (delta_y(2)-m_y)/sigma_y
            pcomp(1,1) = PHID(limit3(3,1)) 
            pcomp(2,1) = BVND(-limit3(2,2),-limit3(3,2),rho_psi_zy) 
            pcomp(3,1) = PHID(limit3(3,3)) 
            pcomp(4,1) = BVND(-limit3(2,4),-limit3(3,4),rho_psi_zy)
            !
        ELSE IF ((riskavn .EQ. 3) .AND. (horizonn .EQ. 4)) THEN
            !
            limit3(2,2) = (delta_z(2)-m_z)/sigma_z
            limit3(3,3) = (delta_y(3)-m_y)/sigma_y
            limit3(2,4) = (delta_z(2)-m_z)/sigma_z
            limit3(3,4) = (delta_y(3)-m_y)/sigma_y
            pcomp(1,1) = 1.d0 
            pcomp(2,1) = PHID(limit3(2,2))
            pcomp(3,1) = PHID(limit3(3,3)) 
            pcomp(4,1) = BVND(-limit3(2,4),-limit3(3,4),rho_psi_zy)
            !
        END IF
        !
        p_zy = pcomp(1,1)-pcomp(2,1)-pcomp(3,1)+pcomp(4,1)
        ! 
    CASE (2)   ! as = 0
        ! 
        ! joint probability
        ! 
        p_zy = 1.d0
        corr3 = (/ rho_sz, rho_sy, rho_zy /)
        limit3(1,:) = -m_s/sigma_s
        eps = 1.D-13
        !
        IF ((riskavn .EQ. 1) .AND. (horizonn .EQ. 1)) THEN
            !
            limit3(2,1) = (delta_z(1)-m_z)/sigma_z
            limit3(3,1) = (delta_y(1)-m_y)/sigma_y
            pcomp(1,1) = TVTL(0,limit3(:,1),corr3,eps)
            !
        ELSE IF ((riskavn .EQ. 1) .AND. (horizonn .EQ. 2)) THEN
            !
            limit3(2,1) = (delta_z(1)-m_z)/sigma_z
            limit3(3,1) = (delta_y(2)-m_y)/sigma_y
            limit3(2,3) = (delta_z(1)-m_z)/sigma_z
            limit3(3,3) = (delta_y(1)-m_y)/sigma_y
            pcomp(1,1) = TVTL(0,limit3(:,1),corr3,eps)
            pcomp(3,1) = TVTL(0,limit3(:,3),corr3,eps)
            !
        ELSE IF ((riskavn .EQ. 1) .AND. (horizonn .EQ. 3)) THEN
            !
            limit3(2,1) = (delta_z(1)-m_z)/sigma_z
            limit3(3,1) = (delta_y(3)-m_y)/sigma_y
            limit3(2,3) = (delta_z(1)-m_z)/sigma_z
            limit3(3,3) = (delta_y(2)-m_y)/sigma_y
            pcomp(1,1) = TVTL(0,limit3(:,1),corr3,eps)
            pcomp(3,1) = TVTL(0,limit3(:,3),corr3,eps)
            !
        ELSE IF ((riskavn .EQ. 1) .AND. (horizonn .EQ. 4)) THEN
            !
            limit3(2,1) = (delta_z(1)-m_z)/sigma_z
            limit3(2,3) = (delta_z(1)-m_z)/sigma_z
            limit3(3,3) = (delta_y(3)-m_y)/sigma_y
            pcomp(1,1) = BVND(-limit3(1,1),-limit3(2,1),rho_sz) 
            pcomp(3,1) = TVTL(0,limit3(:,3),corr3,eps)
            !
        ELSE IF ((riskavn .EQ. 2) .AND. (horizonn .EQ. 1)) THEN
            !
            limit3(2,1) = (delta_z(2)-m_z)/sigma_z
            limit3(3,1) = (delta_y(1)-m_y)/sigma_y
            limit3(2,2) = (delta_z(1)-m_z)/sigma_z
            limit3(3,2) = (delta_y(1)-m_y)/sigma_y
            pcomp(1,1) = TVTL(0,limit3(:,1),corr3,eps)
            pcomp(2,1) = TVTL(0,limit3(:,2),corr3,eps)
            !
        ELSE IF ((riskavn .EQ. 2) .AND. (horizonn .EQ. 2)) THEN
            !
            limit3(2,1) = (delta_z(2)-m_z)/sigma_z
            limit3(3,1) = (delta_y(2)-m_y)/sigma_y
            limit3(2,2) = (delta_z(1)-m_z)/sigma_z
            limit3(3,2) = (delta_y(2)-m_y)/sigma_y
            limit3(2,3) = (delta_z(2)-m_z)/sigma_z
            limit3(3,3) = (delta_y(1)-m_y)/sigma_y
            limit3(2,4) = (delta_z(1)-m_z)/sigma_z
            limit3(3,4) = (delta_y(1)-m_y)/sigma_y
            pcomp(1,1) = TVTL(0,limit3(:,1),corr3,eps)
            pcomp(2,1) = TVTL(0,limit3(:,2),corr3,eps) 
            pcomp(3,1) = TVTL(0,limit3(:,3),corr3,eps)
            pcomp(4,1) = TVTL(0,limit3(:,4),corr3,eps)
            !
        ELSE IF ((riskavn .EQ. 2) .AND. (horizonn .EQ. 3)) THEN
            !
            limit3(2,1) = (delta_z(2)-m_z)/sigma_z
            limit3(3,1) = (delta_y(3)-m_y)/sigma_y
            limit3(2,2) = (delta_z(1)-m_z)/sigma_z
            limit3(3,2) = (delta_y(3)-m_y)/sigma_y
            limit3(2,3) = (delta_z(2)-m_z)/sigma_z
            limit3(3,3) = (delta_y(2)-m_y)/sigma_y
            limit3(2,4) = (delta_z(1)-m_z)/sigma_z
            limit3(3,4) = (delta_y(2)-m_y)/sigma_y
            pcomp(1,1) = TVTL(0,limit3(:,1),corr3,eps)
            pcomp(2,1) = TVTL(0,limit3(:,2),corr3,eps) 
            pcomp(3,1) = TVTL(0,limit3(:,3),corr3,eps)
            pcomp(4,1) = TVTL(0,limit3(:,4),corr3,eps)
            !
        ELSE IF ((riskavn .EQ. 2) .AND. (horizonn .EQ. 4)) THEN
            !
            limit3(2,1) = (delta_z(2)-m_z)/sigma_z
            limit3(2,2) = (delta_z(1)-m_z)/sigma_z
            limit3(2,3) = (delta_z(2)-m_z)/sigma_z
            limit3(3,3) = (delta_y(3)-m_y)/sigma_y
            limit3(2,4) = (delta_z(1)-m_z)/sigma_z
            limit3(3,4) = (delta_y(3)-m_y)/sigma_y
            pcomp(1,1) = BVND(-limit3(1,1),-limit3(2,1),rho_sz) 
            pcomp(2,1) = BVND(-limit3(1,2),-limit3(2,2),rho_sz) 
            pcomp(3,1) = TVTL(0,limit3(:,3),corr3,eps)
            pcomp(4,1) = TVTL(0,limit3(:,4),corr3,eps)
            !
        ELSE IF ((riskavn .EQ. 3) .AND. (horizonn .EQ. 1)) THEN
            !
            limit3(3,1) = (delta_y(1)-m_y)/sigma_y
            limit3(2,2) = (delta_z(2)-m_z)/sigma_z
            limit3(3,2) = (delta_y(1)-m_y)/sigma_y
            pcomp(1,1) = BVND(-limit3(1,1),-limit3(3,1),rho_sy) 
            pcomp(2,1) = TVTL(0,limit3(:,2),corr3,eps)
            !
        ELSE IF ((riskavn .EQ. 3) .AND. (horizonn .EQ. 2)) THEN
            !
            limit3(3,1) = (delta_y(2)-m_y)/sigma_y
            limit3(2,2) = (delta_z(2)-m_z)/sigma_z
            limit3(3,2) = (delta_y(2)-m_y)/sigma_y
            limit3(3,3) = (delta_y(1)-m_y)/sigma_y
            limit3(2,4) = (delta_z(2)-m_z)/sigma_z
            limit3(3,4) = (delta_y(1)-m_y)/sigma_y
            pcomp(1,1) = BVND(-limit3(1,1),-limit3(3,1),rho_sy) 
            pcomp(2,1) = TVTL(0,limit3(:,2),corr3,eps)
            pcomp(3,1) = BVND(-limit3(1,3),-limit3(3,3),rho_sy) 
            pcomp(4,1) = TVTL(0,limit3(:,4),corr3,eps)
            !
        ELSE IF ((riskavn .EQ. 3) .AND. (horizonn .EQ. 3)) THEN
            !
            limit3(3,1) = (delta_y(3)-m_y)/sigma_y
            limit3(2,2) = (delta_z(2)-m_z)/sigma_z
            limit3(3,2) = (delta_y(3)-m_y)/sigma_y
            limit3(3,3) = (delta_y(2)-m_y)/sigma_y
            limit3(2,4) = (delta_z(2)-m_z)/sigma_z
            limit3(3,4) = (delta_y(2)-m_y)/sigma_y
            pcomp(1,1) = BVND(-limit3(1,1),-limit3(3,1),rho_sy) 
            pcomp(2,1) = TVTL(0,limit3(:,2),corr3,eps)
            pcomp(3,1) = BVND(-limit3(1,3),-limit3(3,3),rho_sy) 
            pcomp(4,1) = TVTL(0,limit3(:,4),corr3,eps)
            !
        ELSE IF ((riskavn .EQ. 3) .AND. (horizonn .EQ. 4)) THEN
            !
            limit3(2,2) = (delta_z(2)-m_z)/sigma_z
            limit3(3,3) = (delta_y(3)-m_y)/sigma_y
            limit3(2,4) = (delta_z(2)-m_z)/sigma_z
            limit3(3,4) = (delta_y(3)-m_y)/sigma_y
            pcomp(1,1) = 1.d0 
            pcomp(2,1) = BVND(-limit3(1,2),-limit3(2,2),rho_sz) 
            pcomp(3,1) = BVND(-limit3(1,3),-limit3(3,3),rho_sy) 
            pcomp(4,1) = TVTL(0,limit3(:,4),corr3,eps)
            !
        END IF
        !
        p_as = pcomp(1,1)-pcomp(2,1)-pcomp(3,1)+pcomp(4,1)
        ! 
    CASE (3)   ! as = 1
        ! 
        ! joint probability
        ! 
        p_zy = 1.d0
        corr3 = (/ rho_sz, rho_sy, rho_zy /)
        limit3(1,:) = (q-m_s)/sigma_s
        eps = 1.D-13
        !
        IF ((riskavn .EQ. 1) .AND. (horizonn .EQ. 1)) THEN
            !
            limit3(2,1) = (delta_z(1)-m_z)/sigma_z
            limit3(3,1) = (delta_y(1)-m_y)/sigma_y
            pcomp(1,1) = BVND(-limit3(2,1),-limit3(3,1),rho_zy) 
            pcomp(1,2) = TVTL(0,limit3(:,1),corr3,eps)
            !
        ELSE IF ((riskavn .EQ. 1) .AND. (horizonn .EQ. 2)) THEN
            !
            limit3(2,1) = (delta_z(1)-m_z)/sigma_z
            limit3(3,1) = (delta_y(2)-m_y)/sigma_y
            limit3(2,3) = (delta_z(1)-m_z)/sigma_z
            limit3(3,3) = (delta_y(1)-m_y)/sigma_y
            pcomp(1,1) = BVND(-limit3(2,1),-limit3(3,1),rho_zy) 
            pcomp(3,1) = BVND(-limit3(2,3),-limit3(3,3),rho_zy) 
            pcomp(1,2) = TVTL(0,limit3(:,1),corr3,eps)
            pcomp(3,2) = TVTL(0,limit3(:,3),corr3,eps)
            !
        ELSE IF ((riskavn .EQ. 1) .AND. (horizonn .EQ. 3)) THEN
            !
            limit3(2,1) = (delta_z(1)-m_z)/sigma_z
            limit3(3,1) = (delta_y(3)-m_y)/sigma_y
            limit3(2,3) = (delta_z(1)-m_z)/sigma_z
            limit3(3,3) = (delta_y(2)-m_y)/sigma_y
            pcomp(1,1) = BVND(-limit3(2,1),-limit3(3,1),rho_zy) 
            pcomp(3,1) = BVND(-limit3(2,3),-limit3(3,3),rho_zy) 
            pcomp(1,2) = TVTL(0,limit3(:,1),corr3,eps)
            pcomp(3,2) = TVTL(0,limit3(:,3),corr3,eps)
            !
        ELSE IF ((riskavn .EQ. 1) .AND. (horizonn .EQ. 4)) THEN
            !
            limit3(2,1) = (delta_z(1)-m_z)/sigma_z
            limit3(2,3) = (delta_z(1)-m_z)/sigma_z
            limit3(3,3) = (delta_y(3)-m_y)/sigma_y
            pcomp(1,1) = PHID(limit3(2,1)) 
            pcomp(3,1) = BVND(-limit3(2,3),-limit3(3,3),rho_zy) 
            pcomp(1,2) = BVND(-limit3(1,1),-limit3(2,1),rho_sz)
            pcomp(3,2) = TVTL(0,limit3(:,3),corr3,eps)
            !
        ELSE IF ((riskavn .EQ. 2) .AND. (horizonn .EQ. 1)) THEN
            !
            limit3(2,1) = (delta_z(2)-m_z)/sigma_z
            limit3(3,1) = (delta_y(1)-m_y)/sigma_y
            limit3(2,2) = (delta_z(1)-m_z)/sigma_z
            limit3(3,2) = (delta_y(1)-m_y)/sigma_y
            pcomp(1,1) = BVND(-limit3(2,1),-limit3(3,1),rho_zy) 
            pcomp(2,1) = BVND(-limit3(2,2),-limit3(3,2),rho_zy) 
            pcomp(1,2) = TVTL(0,limit3(:,1),corr3,eps)
            pcomp(2,2) = TVTL(0,limit3(:,2),corr3,eps)
            !
        ELSE IF ((riskavn .EQ. 2) .AND. (horizonn .EQ. 2)) THEN
            !
            limit3(2,1) = (delta_z(2)-m_z)/sigma_z
            limit3(3,1) = (delta_y(2)-m_y)/sigma_y
            limit3(2,2) = (delta_z(1)-m_z)/sigma_z
            limit3(3,2) = (delta_y(2)-m_y)/sigma_y
            limit3(2,3) = (delta_z(2)-m_z)/sigma_z
            limit3(3,3) = (delta_y(1)-m_y)/sigma_y
            limit3(2,4) = (delta_z(1)-m_z)/sigma_z
            limit3(3,4) = (delta_y(1)-m_y)/sigma_y
            pcomp(1,1) = BVND(-limit3(2,1),-limit3(3,1),rho_zy) 
            pcomp(2,1) = BVND(-limit3(2,2),-limit3(3,2),rho_zy) 
            pcomp(3,1) = BVND(-limit3(2,3),-limit3(3,3),rho_zy) 
            pcomp(4,1) = BVND(-limit3(2,4),-limit3(3,4),rho_zy) 
            pcomp(1,2) = TVTL(0,limit3(:,1),corr3,eps)
            pcomp(2,2) = TVTL(0,limit3(:,2),corr3,eps) 
            pcomp(3,2) = TVTL(0,limit3(:,3),corr3,eps)
            pcomp(4,2) = TVTL(0,limit3(:,4),corr3,eps)
            !
        ELSE IF ((riskavn .EQ. 2) .AND. (horizonn .EQ. 3)) THEN
            !
            limit3(2,1) = (delta_z(2)-m_z)/sigma_z
            limit3(3,1) = (delta_y(3)-m_y)/sigma_y
            limit3(2,2) = (delta_z(1)-m_z)/sigma_z
            limit3(3,2) = (delta_y(3)-m_y)/sigma_y
            limit3(2,3) = (delta_z(2)-m_z)/sigma_z
            limit3(3,3) = (delta_y(2)-m_y)/sigma_y
            limit3(2,4) = (delta_z(1)-m_z)/sigma_z
            limit3(3,4) = (delta_y(2)-m_y)/sigma_y
            pcomp(1,1) = BVND(-limit3(2,1),-limit3(3,1),rho_zy) 
            pcomp(2,1) = BVND(-limit3(2,2),-limit3(3,2),rho_zy) 
            pcomp(3,1) = BVND(-limit3(2,3),-limit3(3,3),rho_zy) 
            pcomp(4,1) = BVND(-limit3(2,4),-limit3(3,4),rho_zy) 
            pcomp(1,2) = TVTL(0,limit3(:,1),corr3,eps)
            pcomp(2,2) = TVTL(0,limit3(:,2),corr3,eps) 
            pcomp(3,2) = TVTL(0,limit3(:,3),corr3,eps)
            pcomp(4,2) = TVTL(0,limit3(:,4),corr3,eps)
            !
        ELSE IF ((riskavn .EQ. 2) .AND. (horizonn .EQ. 4)) THEN
            !
            limit3(2,1) = (delta_z(2)-m_z)/sigma_z
            limit3(2,2) = (delta_z(1)-m_z)/sigma_z
            limit3(2,3) = (delta_z(2)-m_z)/sigma_z
            limit3(3,3) = (delta_y(3)-m_y)/sigma_y
            limit3(2,4) = (delta_z(1)-m_z)/sigma_z
            limit3(3,4) = (delta_y(3)-m_y)/sigma_y
            pcomp(1,1) = PHID(limit3(2,1)) 
            pcomp(2,1) = PHID(limit3(2,2)) 
            pcomp(3,1) = BVND(-limit3(2,3),-limit3(3,3),rho_zy) 
            pcomp(4,1) = BVND(-limit3(2,4),-limit3(3,4),rho_zy) 
            pcomp(1,2) = BVND(-limit3(1,1),-limit3(2,1),rho_sz) 
            pcomp(2,2) = BVND(-limit3(1,2),-limit3(2,2),rho_sz) 
            pcomp(3,2) = TVTL(0,limit3(:,3),corr3,eps)
            pcomp(4,2) = TVTL(0,limit3(:,4),corr3,eps)
            !
        ELSE IF ((riskavn .EQ. 3) .AND. (horizonn .EQ. 1)) THEN
            !
            limit3(3,1) = (delta_y(1)-m_y)/sigma_y
            limit3(2,2) = (delta_z(2)-m_z)/sigma_z
            limit3(3,2) = (delta_y(1)-m_y)/sigma_y
            pcomp(1,1) = PHID(limit3(3,1)) 
            pcomp(2,1) = BVND(-limit3(2,2),-limit3(3,2),rho_zy) 
            pcomp(1,2) = BVND(-limit3(1,1),-limit3(3,1),rho_sy) 
            pcomp(2,2) = TVTL(0,limit3(:,2),corr3,eps)
            !
        ELSE IF ((riskavn .EQ. 3) .AND. (horizonn .EQ. 2)) THEN
            !
            limit3(3,1) = (delta_y(2)-m_y)/sigma_y
            limit3(2,2) = (delta_z(2)-m_z)/sigma_z
            limit3(3,2) = (delta_y(2)-m_y)/sigma_y
            limit3(3,3) = (delta_y(1)-m_y)/sigma_y
            limit3(2,4) = (delta_z(2)-m_z)/sigma_z
            limit3(3,4) = (delta_y(1)-m_y)/sigma_y
            pcomp(1,1) = PHID(limit3(3,1)) 
            pcomp(2,1) = BVND(-limit3(2,2),-limit3(3,2),rho_zy) 
            pcomp(3,1) = PHID(limit3(3,3)) 
            pcomp(4,1) = BVND(-limit3(2,4),-limit3(3,4),rho_zy)
            pcomp(1,2) = BVND(-limit3(1,1),-limit3(3,1),rho_sy) 
            pcomp(2,2) = TVTL(0,limit3(:,2),corr3,eps)
            pcomp(3,2) = BVND(-limit3(1,3),-limit3(3,3),rho_sy) 
            pcomp(4,2) = TVTL(0,limit3(:,4),corr3,eps)
            !
        ELSE IF ((riskavn .EQ. 3) .AND. (horizonn .EQ. 3)) THEN
            !
            limit3(3,1) = (delta_y(3)-m_y)/sigma_y
            limit3(2,2) = (delta_z(2)-m_z)/sigma_z
            limit3(3,2) = (delta_y(3)-m_y)/sigma_y
            limit3(3,3) = (delta_y(2)-m_y)/sigma_y
            limit3(2,4) = (delta_z(2)-m_z)/sigma_z
            limit3(3,4) = (delta_y(2)-m_y)/sigma_y
            pcomp(1,1) = PHID(limit3(3,1)) 
            pcomp(2,1) = BVND(-limit3(2,2),-limit3(3,2),rho_zy) 
            pcomp(3,1) = PHID(limit3(3,3)) 
            pcomp(4,1) = BVND(-limit3(2,4),-limit3(3,4),rho_zy) 
            pcomp(1,2) = BVND(-limit3(1,1),-limit3(3,1),rho_sy)
            pcomp(2,2) = TVTL(0,limit3(:,2),corr3,eps)
            pcomp(3,2) = BVND(-limit3(1,3),-limit3(3,3),rho_sy) 
            pcomp(4,2) = TVTL(0,limit3(:,4),corr3,eps)
            !
        ELSE IF ((riskavn .EQ. 3) .AND. (horizonn .EQ. 4)) THEN
            !
            limit3(2,2) = (delta_z(2)-m_z)/sigma_z
            limit3(3,3) = (delta_y(3)-m_y)/sigma_y
            limit3(2,4) = (delta_z(2)-m_z)/sigma_z
            limit3(3,4) = (delta_y(3)-m_y)/sigma_y
            pcomp(1,1) = 1.d0 
            pcomp(2,1) = PHID(limit3(2,2)) 
            pcomp(3,1) = PHID(limit3(3,3)) 
            pcomp(4,1) = BVND(-limit3(2,4),-limit3(3,4),rho_zy) 
            pcomp(1,2) = PHID(limit3(1,1))
            pcomp(2,2) = BVND(-limit3(1,2),-limit3(2,2),rho_sz) 
            pcomp(3,2) = BVND(-limit3(1,3),-limit3(3,3),rho_sy) 
            pcomp(4,2) = TVTL(0,limit3(:,4),corr3,eps)
            !
        END IF
        !
        p_as = pcomp(1,1)-pcomp(2,1)-pcomp(3,1)+pcomp(4,1)-(pcomp(1,2)-pcomp(2,2)-pcomp(3,2)+pcomp(4,2))
        ! 
END SELECT
!
! Evaluating the integrand
!
p = p_as*p_zy+minimum_p
! 
! Ending subroutine and returning execution
! 
END SUBROUTINE individual_contribution
