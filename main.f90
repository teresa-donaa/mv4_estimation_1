PROGRAM main
!
USE constants
USE observations
USE printing_routines
USE load_data
USE starting_points
USE maxlik
!USE UTILITIES_DV_DV
!USE asymptotic_variance
!
IMPLICIT NONE
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Declaring variables
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
CHARACTER (len = 12) :: date_char(3)        ! Not used
INTEGER :: itime(8)                         ! Starting date
INTEGER :: seed(2)                          ! Seed for r.n. generation
INTEGER :: i_stime                          ! Estimation trial loop index
REAL(8) :: llcheck
REAL(8) :: objf                             ! Loglikelihood function at the optimum
!REAL(8) :: grad(num_theta)                  ! Gradient at optimal theta
!CHARACTER(len=60) :: task
!REAL(8) :: v1mat(num_theta,num_theta)       ! Estimated variance matrix of the parameters
!REAL(8) :: stderr1(num_theta)               ! Estimated asymptotic standard errors
!!
!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!! Beginning execution
!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!
! Loading the model's specification
!
CALL input_model(sel_X)
CALL CreateVariables()
!
IF (compute_var_as .EQ. 0) THEN
    !
    ! Opening output files
    !
    CALL open_write_file(unit_res,file_res)
    !
    ! Loading data
    !
    CALL input_data(theta_inf,theta_sup)
    !
    ! Initializing random number generator
    !
    CALL date_and_time(date_char(1),date_char(2),date_char(3),itime)
    seed(1) = itime(7)*itime(8)
    seed(2) = itime(5)*itime(6)
    !
    ! Starting loop 
    !
    DO i_stime = 1, num_stime
        !
        ! Creating random starting values of the parameters 
        !
        CALL admissible_starting_point(seed,theta,theta_inf,theta_sup,llcheck)
        !
        ! Estimation 
        !
        error_flag = .FALSE.
!@SP
CALL loglik_fct(num_theta,theta,objf,grad)
!@SP
!        IF (to0 .EQ. 1) CALL loglik_fct(num_theta,theta,objf,grad)
!        IF (to0 .NE. 1) CALL estimate(i_stime,theta,objf,grad,task)
!        !
!        ! Printing intermediate trace output 
!        !
!        CALL print_res(i_stime,objf,theta,task,grad)
!        ! 
    END DO 
    !
    ! Closes optimization stage
    !
    CLOSE(UNIT=unit_res)
    !
END IF
!!
!! Computing final statistics
!!
!IF (compute_var_as .EQ. 1) THEN
!    !
!    ! Reading parameters estimates
!    !
!    CALL admissible_starting_point(seed,theta,theta_inf,theta_sup,llcheck)
!    !
!    ! Opening and reading other input files
!    !
!    CALL input_names()
!    CALL input_data(theta_inf,theta_sup)
!    CALL open_write_file(unit_fin_res,file_fin_res)
!    CALL open_write_file(unit_vartheta,file_vartheta)
!    !
!    ! Computing and writing asymptotic standard errors
!    !
!    CALL compute_asymptotic_variance(theta,llcheck,objf,grad,v1mat,stderr1)
!    CALL print_final_results(theta,stderr1,objf,grad,v1mat)
!    CLOSE(UNIT=unit_fin_res)
!    CLOSE(UNIT=unit_vartheta)
!    !
!END IF
!
CALL DestroyVariables()
!
! End of execution
!
END PROGRAM main
