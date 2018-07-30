MODULE constants
!
IMPLICIT NONE
! 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Selecting optimization method and other switches
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
INTEGER, PARAMETER :: to0 = 0           ! Modalità di test
INTEGER, PARAMETER :: to1 = 1           ! Prime stime
INTEGER, PARAMETER :: to2 = 0           ! Seconde stime
INTEGER, PARAMETER :: to3 = 0           ! Varianza
!
!REAL(8), PARAMETER :: max_rpibpis = 1.d0
!INTEGER, PARAMETER :: switch_Sigma1 = 0 ! = 1: nonzero rpibpis, 
!                                        ! constrained between (-max_rpibpis , +max_rpibpis)
!INTEGER, PARAMETER :: switch_Sigma2 = 1 ! = 1: rpibpis = 0
!!
!INTEGER, PARAMETER :: num_stime = to0*1+to1*100+to2*1+to3*0
!                                        ! Total number of completed estimation trials
INTEGER, PARAMETER :: compute_var_as = to0*0+to1*0+to2*0+to3*1      
                                        ! Switch to compute the asymptotic variance matrix
! 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Declaring sample characteristics and selecting explanatory variables
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! SCF 1995 - 2016 data
!
INTEGER, PARAMETER :: num_N = 24811     ! mv4 data 95-16
INTEGER, PARAMETER :: num_tot_X = 42    ! Total number of available explanatory variables
INTEGER, PARAMETER :: num_L = 3         ! Degrees of self-reported risk aversion
INTEGER, PARAMETER :: num_H = 4         ! Degrees of self-reported planning horizon
INTEGER, PARAMETER :: num_states = 4    ! Number of state variables: (m_s,m_z,m_y,sigma_s)
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Defining the number of parameters
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! theta
!
INTEGER, PARAMETER :: num_theta_delta_z = num_L-2
INTEGER, PARAMETER :: num_theta_delta_y = num_H-2
!
! psi
!
INTEGER, PARAMETER :: num_psi_m = 3       
INTEGER, PARAMETER :: num_psi_sigma = 1   
INTEGER, PARAMETER :: num_psi = num_psi_m+num_psi_sigma+num_theta_delta_z+num_theta_delta_y
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Declaring optimization switches
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 
REAL(8), PARAMETER :: factr = (to0+to1)*1.d+7+(to2+to3)*1.d+1
REAL(8), PARAMETER :: pgtol = (to0+to1)*1.0d-4+(to2+to3)*1.d-5
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Declaring parameters about the output files
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! 
CHARACTER(len=30), PARAMETER :: file_data = 'data_1.txt'
CHARACTER(len=30), PARAMETER :: file_names = 'names_data_1.txt'
CHARACTER(len=30), PARAMETER :: file_model = 'model.txt'
CHARACTER(len=30), PARAMETER :: file_theta = 'theta.txt'
!
CHARACTER(len=20), PARAMETER :: file_res = 'res.txt'
!CHARACTER(len=20), PARAMETER :: file_res_to1 = 'res_to1.txt'
!CHARACTER(len=20), PARAMETER :: file_res_to2 = 'res_to2.txt'
!CHARACTER(len=20), PARAMETER :: file_fin_res = 'fin_res.txt'
!CHARACTER(len=20), PARAMETER :: file_vartheta = 'vartheta.txt'
!CHARACTER(len=20), PARAMETER :: file_imat = 'imat.txt'
!CHARACTER(len=20), PARAMETER :: file_jmat = 'jmat.txt'
!CHARACTER(len=20), PARAMETER :: file_invjmat = 'invjmat.txt'
!CHARACTER(len=20), PARAMETER :: file_dl = 'dl.txt'
!CHARACTER(len=20), PARAMETER :: file_loglik = 'loglik.txt'
!
INTEGER, PARAMETER :: unit_data = 1             
INTEGER, PARAMETER :: unit_names = 2            
INTEGER, PARAMETER :: unit_model = 3
INTEGER, PARAMETER :: unit_theta = 4
!
INTEGER, PARAMETER :: unit_res = 21             
!INTEGER, PARAMETER :: unit_res_to1 = 22         
!INTEGER, PARAMETER :: unit_res_to2 = 23         
!INTEGER, PARAMETER :: unit_fin_res = 24         
!INTEGER, PARAMETER :: unit_vartheta = 25        
!INTEGER, PARAMETER :: unit_imat = 26            
!INTEGER, PARAMETER :: unit_jmat = 27            
!INTEGER, PARAMETER :: unit_invjmat = 28     
!INTEGER, PARAMETER :: unit_dl = 29      
!INTEGER, PARAMETER :: unit_loglik = 30
! 
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Declaring constants
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
REAL(8), PARAMETER :: minimum_p = 1.d-10                                ! Minimum probability value
REAL(8), PARAMETER :: pi = 3.14159265358979323846264338328d0            ! Pi
REAL(8), PARAMETER :: invsqrtpi = 0.564189583547756286948079451561d0    ! 1/Sqrt[Pi]
REAL(8), PARAMETER :: sqrt2 = 1.41421356237309504880168872421d0         ! Sqrt[2]
REAL(8), PARAMETER :: invsqrt2pi = 0.398942280401432677939946059934d0   ! 1/Sqrt[2*Pi]
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
END MODULE constants
