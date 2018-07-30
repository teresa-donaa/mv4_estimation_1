SUBROUTINE individual_contribution ( psi, dn, asn, riskavn, horizonn, p )
!
USE constants
USE observations
!USE UTILITIES_DV_DV
!USE SP_DV_DV
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
! Declaring local variables
! 
!REAL(8) :: mpib, mpis, xib, xis, rho, kappa, omegab, omegas, sigmapib, sigmapis, rpibpis
!REAL(8) :: sigmaz, delta(num_theta_delta)
!REAL(8) :: pd, pc
!REAL(8) :: mtaub, mtaus, sigma2taub, sigma2taus, sigmataubtaus, logkappa, rho2
!REAL(8) :: den, den2, mab, mas, sigma2ab, sigma2as, sigmaabas, sigmaab, sigmaas, rabas
!REAL(8) :: kos, kob, qsrb, qbrs, arg, argb, args, mb, ms, sigma2b, sigma2s, sigmabs, rbs
!REAL(8) :: taub, taus, sigmataubas, sigmaabtaus, sigmaabs
! 
! Beginning execution
! 
! Computing the b parameters
! 
!CALL compute_b_parameters(psi,mpib,mpis,xib,xis,rho,kappa, &
!    omegab,omegas,delta,sigmapib,sigmapis,rpibpis,sigmaz)
!mtaub = mpib-kappa*arn*xib
!mtaus = mpis-kappa*arn*xis
!sigma2taub = sigmapib**2
!sigma2taus = sigmapis**2
!sigmataubtaus = rpibpis*sigmapib*sigmapis
!logkappa = LOG(kappa)
!rho2 = rho**2
!! 
!! Computing individual contribution to the loglikelihood - PORTFOLIO CHOICE
!! 
!SELECT CASE (dn)
!    CASE (1)
!        ! 
!        ! Computing the density 
!        ! 
!        den = kappa*(1.d0-rho2)
!        den2 = den**2
!        mab = (mtaub-rho*mtaus)/(den*omegab)
!        mas = (mtaus-rho*mtaub)/(den*omegas)
!        sigma2ab = (sigma2taub-2.d0*rho*sigmataubtaus+rho2*sigma2taus)/(den2*omegab**2)
!        sigma2as = (sigma2taus-2.d0*rho*sigmataubtaus+rho2*sigma2taub)/(den2*omegas**2)
!        sigmaabas = ((1.d0+rho2)*sigmataubtaus-rho*(sigma2taub+sigma2taus))/(den2*omegab*omegas)
!        sigmaab = SQRT(sigma2ab)
!        sigmaas = SQRT(sigma2as)
!        rabas = sigmaabas/(sigmaab*sigmaas)
!        pd = pdf_biv_N01((abn-mab)/sigmaab,(asn-mas)/sigmaas,rabas)/(sigmaab*sigmaas)
!        ! 
!    CASE (2)
!        ! 
!        ! Computing the product of density and probability 
!        ! 
!        kos = kappa*omegas
!        mas = mtaus/kos+rho*omegab*arn/omegas
!        sigma2as = sigma2taus/kos**2
!        sigmataubas = sigmataubtaus/kos
!        taus = kappa*(omegas*asn-omegab*rho*arn)
!        arg = (rho*taus-kappa*omegab*(1.d0-rho2)*arn-(mtaub+sigmataubas/sigma2as*(asn-mas)))/ &
!            SQRT(sigma2taub-sigmataubas**2/sigma2as)
!        pd = 1.d0/SQRT(sigma2as)*pdf_N01((asn-mas)/SQRT(sigma2as))*cdf_N01(arg,.FALSE.)
!        ! 
!    CASE (3)
!        ! 
!        ! Computing the product of density and probability 
!        ! 
!        kob = kappa*omegab
!        mab = mtaub/kob
!        sigma2ab = sigma2taub/kob**2
!        sigmaabtaus = sigmataubtaus/kob
!        taub = kob*abn
!        arg = (rho*taub-(mtaus+sigmaabtaus/sigma2ab*(abn-mab)))/SQRT(sigma2taus-sigmaabtaus**2/sigma2ab)
!        pd = 1.d0/SQRT(sigma2ab)*pdf_N01((abn-mab)/SQRT(sigma2ab))*cdf_N01(arg,.FALSE.)
!        ! 
!    CASE (4)
!        ! 
!        ! Computing the product of density and probability 
!        ! 
!        den = kappa*(omegab**2-2.d0*rho*omegab*omegas+omegas**2)
!        qsrb = omegas-rho*omegab
!        qbrs = omegab-rho*omegas
!        mab = (omegab*mtaub-omegas*mtaus+kappa*omegas*qsrb*(1-arn))/den
!        ms = qsrb*mtaub+qbrs*mtaus
!        sigma2ab = (omegab**2*sigma2taub+omegas**2*sigma2taus-2.d0*omegab*omegas*sigmataubtaus)/den**2
!        sigma2s = qsrb**2*sigma2taub+qbrs**2*sigma2taus+2.d0*qsrb*qbrs*sigmataubtaus
!        sigmaabs = (omegab*qsrb*sigma2taub-omegas*qbrs*sigma2taus+(omegab**2-omegas**2)*sigmataubtaus)/den
!        arg = (kappa*omegab*omegas*(1.d0-rho2)*(1.d0-arn)-(ms+sigmaabs/sigma2ab*(abn-mab)))/ &
!            SQRT(sigma2s-sigmaabs**2/sigma2ab)
!        pd = 1.d0/SQRT(sigma2ab)*pdf_N01((abn-mab)/SQRT(sigma2ab))*cdf_N01(arg,.TRUE.)
!        ! 
!    CASE (5)
!        ! 
!        ! Computing the probability 
!        ! 
!        argb = (-kappa*omegab*arn-mtaub)/SQRT(sigma2taub)
!        args = (-kappa*omegab*rho*arn-mtaus)/SQRT(sigma2taus)
!        rbs = sigmataubtaus/SQRT(sigma2taub*sigma2taus)
!        pd = cdf_biv_N01(argb,args,rbs)
!        ! 
!    CASE (6)
!        ! 
!        ! Computing the probability 
!        ! 
!        mb = omegab*mtaub-omegas*mtaus
!        ms = -mtaus
!        sigma2b = omegab**2*sigma2taub+omegas**2*sigma2taus-2.d0*omegab*omegas*sigmataubtaus
!        sigma2s = sigma2taus
!        sigmabs = -omegab*sigmataubtaus+omegas*sigma2taus
!        argb = (-kappa*omegas*(omegas-rho*omegab)-kappa*omegab*(omegab-rho*omegas)*arn-mb)/SQRT(sigma2b)
!        args = (-kappa*omegas+kappa*omegab*rho*arn-ms)/SQRT(sigma2s)
!        rbs = sigmabs/SQRT(sigma2b*sigma2s)
!        pd = cdf_biv_N01(argb,args,rbs)
!        ! 
!    CASE (7)
!        ! 
!        ! Computing the probability 
!        ! 
!        mb = -mtaub
!        ms = -omegab*mtaub+omegas*mtaus
!        sigma2b = sigma2taub
!        sigma2s = omegab**2*sigma2taub+omegas**2*sigma2taus-2.d0*omegab*omegas*sigmataubtaus
!        sigmabs = omegab*sigma2taub-omegas*sigmataubtaus
!        argb = (-kappa*omegab*(1.d0-arn)-mb)/SQRT(sigma2b)
!        args = (-kappa*omegab*(omegab-rho*omegas)*(1.d0-arn)-ms)/SQRT(sigma2s)
!        rbs = sigmabs/SQRT(sigma2b*sigma2s)
!        pd = cdf_biv_N01(argb,args,rbs)
!        ! 
!END SELECT
!! 
!! Computing individual contribution to the loglikelihood - SELF-REPORTED RISK AVERSION
!! 
!SELECT CASE (cn)
!    ! 
!    ! 1 "Take substantial financial risks expecting to earn substantial returns" or
!    !   "Take above average financial risks expecting to earn above average returns"
!    ! 
!    CASE (1)
!        pc = cdf_N01(-logkappa/sigmaz,.FALSE.)
!    ! 
!    ! 2 "Take average financial risks expecting to earn average returns"
!    ! 
!    CASE (2)
!        pc = cdf_N01((delta(1)-logkappa)/sigmaz,.FALSE.)-cdf_N01(-logkappa/sigmaz,.FALSE.)
!    ! 
!    ! 3 "Not willing to take any financial risks"
!    ! 
!    CASE (3)
!        pc = cdf_N01((delta(1)-logkappa)/sigmaz,.TRUE.)
!    ! 
!END SELECT
!!
!! Evaluating the integrand
!!
!p = pd*pc+minimum_p
! 
! Ending subroutine and returning execution
! 
END SUBROUTINE individual_contribution
