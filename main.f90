PROGRAM main
!
USE constants
USE observations
USE printing_routines
USE load_data
USE starting_points
USE maxlik
USE simplex
!USE UTILITIES_DV_DV
USE asymptotic_variance
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
CHARACTER(len=60) :: task
REAL(8) :: min_objf
INTEGER :: iter
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Beginning execution
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
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
!@SP
!    CALL date_and_time(date_char(1),date_char(2),date_char(3),itime)
!    seed(1) = itime(7)*itime(8)
!    seed(2) = itime(5)*itime(6)
    seed = 1
!@SP
    !
    ! Starting loop 
    !
    min_objf = 1.d3
    DO i_stime = 1, num_stime
        !
        ! Creating random starting values of the parameters 
        !
        CALL admissible_starting_point(seed,theta,theta_inf,theta_sup,llcheck)
        !
        ! Estimation 
        !
        error_flag = .FALSE.
        IF (to0 .EQ. 1) CALL loglik_fct(num_theta,theta,objf,grad)
        IF (to0 .NE. 1) THEN
            !
            ! POLITOPE Estimation 
            !
            IF (switch_politope .EQ. 1) THEN
                !
                CALL open_write_file(unit_politope,file_politope)
                CALL politope(i_stime,loglik_politope,theta,pert_theta,1,100,objf,iter)
                CLOSE(UNIT=unit_politope)
                min_objf = MIN(objf,min_objf)
                !
            END IF
            !
            ! L-BFGS Estimation 
            !
            IF ((switch_lbfgs .EQ. 1) .AND. (objf .LE. min_objf+2.d0)) THEN
                !
                CALL estimate(i_stime,theta,objf,grad,task)
                min_objf = MIN(objf,min_objf)
                !
            END IF
            !
        END IF
        !
        ! Printing intermediate trace output 
        !
        CALL print_res(i_stime,objf,theta,task,grad)
        ! 
    END DO 
    !
    ! Closes optimization stage
    !
    CLOSE(UNIT=unit_res)
    !
END IF
!
! Computing final statistics
!
IF (compute_var_as .EQ. 1) THEN
    !
    ! Reading parameters estimates
    !
    CALL admissible_starting_point(seed,theta,theta_inf,theta_sup,llcheck)
    !
    ! Opening and reading other input files
    !
    CALL input_names()
    CALL input_data(theta_inf,theta_sup)
    CALL open_write_file(unit_fin_res,file_fin_res)
    CALL open_write_file(unit_varthetaQML,file_varthetaQML)
    !
    ! Computing and writing asymptotic standard errors
    !
    CALL compute_asymptotic_variance(theta,llcheck,thetaF,objf,grad,varQMLmat,stderrQML)
    CALL print_final_results(thetaF,stderrQML,objf,grad,varQMLmat)
    CLOSE(UNIT = unit_fin_res)
    CLOSE(UNIT = unit_varthetaQML)
    !
END IF
!
CALL DestroyVariables()
!
! End of execution
!
END PROGRAM main
