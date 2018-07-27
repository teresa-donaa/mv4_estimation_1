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
INTEGER :: num_X_s
INTEGER :: num_X_z
INTEGER :: num_X_y
!
!INTEGER, SAVE :: d(num_N)						
!INTEGER, SAVE :: c(num_N)						
!REAL(8), SAVE :: ab(num_N)						
!REAL(8), SAVE :: as(num_N)						
!REAL(8), SAVE :: ar(num_N)
!REAL(8), SAVE :: x_pib(num_N,num_X_pib)
!REAL(8), SAVE :: x_pis(num_N,num_X_pis)
!REAL(8), SAVE :: x_xib(num_N,num_X_xib)
!REAL(8), SAVE :: x_xis(num_N,num_X_xis)
!REAL(8), SAVE :: x_rho(num_N,num_X_rho)
!REAL(8), SAVE :: x_kappa(num_N,num_X_kappa)     
!REAL(8), SAVE :: x_omegab(num_N,num_X_omegab)	
!REAL(8), SAVE :: x_omegas(num_N,num_X_omegas)
!!
!REAL(8), SAVE :: peso(num_N)                    ! Sampling weights
!REAL(8), SAVE :: normpeso(num_N)                ! Normalized sampling weights
!!
!CHARACTER(len=10), SAVE :: names_pib(num_X_pib)	
!CHARACTER(len=10), SAVE :: names_pis(num_X_pis)	
!CHARACTER(len=10), SAVE :: names_xib(num_X_xib)
!CHARACTER(len=10), SAVE :: names_xis(num_X_xis)
!CHARACTER(len=10), SAVE :: names_rho(num_X_rho)
!CHARACTER(len=10), SAVE :: names_kappa(num_X_kappa)	
!CHARACTER(len=10), SAVE :: names_omegab(num_X_omegab)	
!CHARACTER(len=10), SAVE :: names_omegas(num_X_omegas)	
!!
!REAL*8, SAVE :: idmat(num_psi,num_psi)          ! Used for analytic differentiation
!LOGICAL, SAVE :: error_flag
!
! Ending module
!
END MODULE observations
