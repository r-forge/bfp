      function laplace(nObs, nCoefs, linPred, L)
      implicit none
!
!	  Computes the eighth-order Laplace approximation in the logistic
!	  regression case.
!
!	  Input:
!
!	  nObs		number of observations
!	  nCoefs	number of coefficients
!	  linPred	linear predictor vector of length nObs
!	  L		nCoefs x nObs matrix containing products of the
!	  		the covariance Cholesky factor with the design columns
!
	  integer nObs, nCoefs
	  double precision, dimension(0:(nObs - 1)) :: linPred
	  double precision, dimension(0:(nCoefs - 1),0:(nObs - 1)) :: L
!
	  double precision laplace
          double precision ddot
          external ddot
!
	  integer i, j
	  double precision eT4, eT3T3, eT6,
     $         eT3T5, eT4T4, eT8,
     $         mu, m2, oneMinus12m2, m3m3, m4m4, m3m5, m6, m8,
     $         lii2, lii3, lii4, lij, liiljjlij
!
	  double precision, dimension(0:(nObs - 1)) :: m3
      double precision, dimension(0:(nObs - 1)) :: m4
	  double precision, dimension(0:(nObs - 1)) :: m5
	  double precision, dimension(0:(nObs - 1)) :: lProds
!
	  eT4 = 0.0d0
	  eT3T3 = 0.0d0
	  eT6 = 0.0d0
	  eT3T5 = 0.0d0
	  eT4T4 = 0.0d0
	  eT8 = 0.0d0
!     
          do i = 0,(nObs - 1)
             mu = 1.0d0 / (1.0d0 + exp(- linPred(i)))
             m2 = mu * (1.0d0 - mu)
             oneMinus12m2 = 1.0d0 - 12.0d0 * m2
             m3(i) = m2 * (1.0d0 - 2.0d0 * mu)
             m3m3 = m3(i) * m3(i)
             m4(i) = m2 * (1.0d0 - 6.0d0 * m2)
             m4m4 = m4(i) * m4(i)
             m5(i) = m3(i) * oneMinus12m2
             m3m5 = m3(i) * m5(i)
             m6 = m4(i) * oneMinus12m2 - 12.0d0 * m3m3
             m8 = m6 * oneMinus12m2 - 48.0d0 * m3m5 - 36.0d0 * m4m4
!     
             lProds(i) = ddot(nCoefs,
     $                        L(:, i), 1,
     $                        L(:, i), 1)
             lii2 = lProds(i) * lProds(i)
             lii3 = lii2 * lProds(i)
             lii4 = lii2 * lii2
!
             eT4 = eT4 + m4(i) * lii2
             eT3T3 = eT3T3 + m3m3 * lii3
             eT6 = eT6 + m6 * lii3
             eT3T5 = eT3T5 + m3m5 * lii4
             eT4T4 = eT4T4 + m4m4 * lii4
             eT8 = eT8 + m8 * lii4
!     
             do j = 0, (i - 1)
             	lij = ddot(nCoefs,
     $                     L(:, i), 1,
     $                     L(:, j), 1)
                liiljjlij = lProds(i) * lProds(j) * lij
                eT3T3 = eT3T3 + 2.0d0 * m3(i) * m3(j) * liiljjlij
                eT3T5 = eT3T5 + (m3(i) * m5(j) * lProds(j) +
     $               m3(j) * m5(i) * lProds(i)) * liiljjlij
                eT4T4 = eT4T4 + 2.0d0 * m4(i) * m4(j) * liiljjlij * lij
             end do
!     
	  end do
!
          laplace =  - 1.0d0 / 8.0d0 * eT4 +
     $               5.0d0 / 24.0d0 * eT3T3 -
     $               1.0d0 / 48.0d0 * eT6 +
     $               7.0d0 / 96.0d0 * eT3T5 +
     $               35.0d0 / 384.0d0 * eT4T4 -
     $               1.0d0 / 384.0d0 * eT8
!     
          end function laplace


