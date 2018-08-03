MODULE observations
!
! Collects the data
!
USE constants
!
IMPLICIT NONE
!
! Explanatory variables 
!
INTEGER :: sel_X(num_tot_X,num_states)
INTEGER :: num_X_m_s, num_X_m_q, num_X_m_y, num_X_sigma_s
!
! Parameters
!
INTEGER :: num_theta_beta, num_theta_sigma_s, num_theta
REAL(8), ALLOCATABLE :: theta_inf(:), theta_sup(:), theta(:), grad(:)
!
! Data
!
INTEGER, DIMENSION(num_N), SAVE :: d, riskav, horizon
REAL(8), DIMENSION(num_N), SAVE :: ab, as, peso, normpeso
REAL(8), ALLOCATABLE :: x_m_s(:,:)
REAL(8), ALLOCATABLE :: x_m_q(:,:)
REAL(8), ALLOCATABLE :: x_m_y(:,:)
REAL(8), ALLOCATABLE :: x_sigma_s(:,:)
REAL(8) :: idmat(num_psi,num_psi)                       ! Used for analytic differentiation
!
! Names
!
CHARACTER(len = 15), ALLOCATABLE :: names_m_s(:)	
CHARACTER(len = 15), ALLOCATABLE :: names_m_q(:)	
CHARACTER(len = 15), ALLOCATABLE :: names_m_y(:)
CHARACTER(len = 15), ALLOCATABLE :: names_sigma_s(:)
!
! Other
!
LOGICAL :: error_flag
!
    CONTAINS
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE CreateVariables ( )
    !
    ! Allocates arrays whose size is known at runtime only
    !
    IMPLICIT NONE
    !
    ! Beginning execution
    !
    ALLOCATE(theta_inf(num_theta),theta_sup(num_theta),theta(num_theta),grad(num_theta), &
        x_m_s(num_N,num_X_m_s),x_m_q(num_N,num_X_m_q),x_m_y(num_N,num_X_m_y),x_sigma_s(num_N,num_X_sigma_s), &
        names_m_s(num_X_m_s),names_m_q(num_X_m_q),names_m_y(num_X_m_y),names_sigma_s(num_X_sigma_s))
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE CreateVariables
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    SUBROUTINE DestroyVariables ( )
    !
    ! Deallocates arrays whose size is known at runtime only
    !
    IMPLICIT NONE
    !
    ! Beginning execution
    !
    DEALLOCATE(theta_inf,theta_sup,theta,grad, &
        x_m_s,x_m_q,x_m_y,x_sigma_s, &
        names_m_s,names_m_q,names_m_y,names_sigma_s)
    !
    ! Ending execution and returning control
    !
    END SUBROUTINE DestroyVariables
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Ending module
!
END MODULE observations
