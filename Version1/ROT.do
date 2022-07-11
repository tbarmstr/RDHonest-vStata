* add this after IKBW_fit ===================

  class RDData scalar ROTBW_fit(class RDData scalar df, real matrix kernC,| string kernel, real scalar boundary) {
		/* Calculate bandwidth for sharp RD based on local linear regression
		using method by in Fan and Gijbels (1996, Chapter 4.2) */
		/* note: no option for order != 1 */

		if(args() < 3)  kernel = "triangular"

		real scalar order, Nm, Np, N, kernInd
		real scalar h1, f0, m3, m3sigma, B, V
		real vector col1, col2, col3, col4, col5, m3coeff, tempY
		real matrix s, X, m3mat

		X = df.Xm\df.Xp
		order = 1
		Nm = length(df.Xm)
		Np = length(df.Xp)
		N = Nm + Np
		
		/* specify boundary: do we need this? */
		if(args() < 4)  {
			if (min(X) >= 0 | max(X) <= 0) {
				boundary = 1
			}
			else {
				boundary = 0
			}	
		}

		/* kernels must be specified as strings */

		/* translate kernel strings into kernel indices */
		kernInd = ((strpos(kernel,"uni")>0) ? 0 :
			(strpos(kernel,"tri")>0 ? 1 : 2))

		/* access kernel constants */
		/* select kernC row with matching kernel index, order == 1 and
			boundary == TRUE */
		s = select(kernC,kernC[.,1]:==kernInd)
		s = select(s,s[.,2]:==1)
		s = select(s,s[.,3]:==boundary)
		s = s'
		/* s is now a colvector, so we can avoid complications with variance() */
		
		
		/* estimate f_X(0) */
		h1 = 1.843 *min( (sqrt(variance(X)), (qtile(X, 0.75)-qtile(X, 0.25))/1.349) ) / (N^(1/5))
		f0 = sum(abs(X) :<= h1)/(2 * N * h1)
		
		/* run a linear regression of y on x^p, p=0,1,2,3,4 (0 to 1+3), then extract
		coefficients of x^(1+2) */
		col1 = X :^ 0; col2 = X; col3 = X :^ 2; col4 = X :^ 3; col5 = X :^ 4;
		m3mat = col1, col2, col3, col4, col5
		m3coeff = invsym(m3mat'm3mat)*m3mat'(df.Ym\df.Yp)
		/* maybe use cross() for speed
		m3coeff = invsym( cross(m3mat,m3mat) )*( cross(m3mat,(df.Ym\df.Yp)) )
		*/
		m3 = m3coeff[4]
		m3sigma = variance( (df.Ym\df.Yp) - m3mat*m3 )
		
		/* calculate results */
		B = m3 * s[6]
		V = m3sigma * s[9] / f0
		
		return (V/(B^2 * 2 * 2 * N))^(1/5)
	}
	
* add this to accessory functions ============

  // calculate the percentile of a vector
  real scalar qtile(real vector xs, real scalar qt) {
    
	real vector tmp
    real scalar idx
    
	// qt takes value from 0 to 1
    if ((qt <= 0) | (qt >= 1)) {
        _error(3300)
    }
    
    tmp = sort(xs, 1)
    idx = rows(xs)*qt
    
    if (mod(idx, 1)) {
        return(tmp[ceil(idx)])
    }
    else {
        return(((tmp[idx] + tmp[idx + 1])/2))
    }
  }

  *==============================================
  * Rule of thumb choosing M
  *==============================================
  
  // Rule of thumb for choosing M (perhaps should be put into a class here)
  class RDData scalar MROT_fit(class RDData scalar df) {
		/* Use global quadratic regression to estimate a bound on the second derivative */
		/* note: this is only for RD, nor FRD yet */

		real scalar M1, M2, m3, m3sigma
		real vector col1, col2, col3, col4, col5, m3coeff
		real matrix s, m3mat
		
		/* Step 1: Estimate global polynomial regression */
		// below cutoff
		col1 = df.Xm :^ 0; col2 = df.Xm; col3 = df.Xm :^ 2; col4 = df.Xm :^ 3; col5 = df.Xm :^ 4;
		m3mat = col1, col2, col3, col4, col5
		m3coeff = invsym(m3mat'm3mat)*m3mat'(df.Ym)
		M1 = MROT_aux_endpoint(m3coeff, df.Xm)

		// above cutoff
		col1 = df.Xp :^ 0; col2 = df.Xp; col3 = df.Xp :^ 2; col4 = df.Xp :^ 3; col5 = df.Xp :^ 4;
		m3mat = col1, col2, col3, col4, col5
		m3coeff = invsym(m3mat'm3mat)*m3mat'(df.Yp)
		M2 = MROT_aux_endpoint(m3coeff, df.Xp)

		return max((M1, M2))
	}
	

  /* Auxilary function for MROT_fit: maximum occurs either at endpoints, or else 
  at the extremum, -r1[4]/(4*r1[5]), if the extremum is in the support */
  real scalar MROT_aux_endpoint(real vector coef, real vector X){
		
		real scalar Fmin, Fmax, Fe
		
		if (abs(coef[5])<=1e-10) {
			Fe = .
		}
		else {
			Fe = -coef[4]/(4*coef[5])
		}
		
		Fmin = abs( 2*coef[3] + 6*min(X)*coef[4] + 12*min(X)*coef[5])
		Fmax = abs( 2*coef[3] + 6*max(X)*coef[4] + 12*max(X)*coef[5])
		
		if ( min(X)<Fe & max(X)>Fe ) {
			M = max(abs( 2*coef[3] + 6*Fe*coef[4] + 12*Fe*coef[5]), max((Fmin,Fmax)))
		}
		else {
			M = max((Fmin,Fmax))
		}
  }