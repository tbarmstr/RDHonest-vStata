*! rdhonest v1.0.0 KHLee MAY-27-2019

program define rdhonest, eclass
	version 14.2
	syntax varlist(min=1) [if] [in] ///
	 [,m(real) c(real 0) hp(real 0) kernel(string) hm(real 0) opt_criterion(string) ///
	 bw_equal(integer 1) alpha(real 0.05) beta(real 0.8) se_method(string) ///
         j(integer 3) sclass(string) order(integer 1) se_initial(string)]
  * SZ: set m() as optional as well, and MROT as default 

  marksample touse
  preserve
  qui keep if `touse'
  local y `1'
  // disp "`y'"
  local x = subinstr("`2'",",","",.) // strip running variable of commas
  // disp "`x'"

  // set more off for now
  set more off

  /* set defaults */
  if (`hm' == 0 & `bw_equal' == 1) local hm = `hp' // default to equate hm and hp
  if ("`kernel'" == "") local kernel "tri"
  if ("`opt_criterion'" == "") local opt_criterion "MSE"
  if ("`se_method'" == "") local se_method "nn"
  if ("`sclass'" == "") local sclass "H"
  if ("`se_initial'" == "") local se_initial "IKEHW"

  // "hacky" method to force Stata to recognize
  // RDData, RDOptions and RDResults
  qui do rdhonest_build.do

  set matastrict on
  mata {
    Y = st_data(.,("`y'"),0);	X = st_data(.,("`x'"),0)

    /* initialize data frame and map Y and X in */
    df = RDDataPrep(X,Y,`c')

    /* printf("`se_method', `j', `sclass', `order', `se_initial'") */

    /* initialize options */
    m = `m'
    bw_equal = `bw_equal'
    kernel = "`kernel'"
    hp = `hp';  hm = `hm'
    opt_criterion = "`opt_criterion'"
    alpha = `alpha'
    beta = `beta'
    se_method = "`se_method'"
    j = `j'
    sclass = "`sclass'"
    order =  `order'
    se_initial = "`se_initial'"

    opt = RDOptionsPrep(m,bw_equal,kernel,hp,hm,opt_criterion,
        alpha,beta,se_method,j,sclass,order,se_initial)

    ret = RDResults()
  }

  // load kernel constants into memory. Note: uniform == 0, triangular == 1,
  // epanechnikov == 2.
  use kernC, replace
  local kernArgs kernel order boundary mu0 mu1 mu2 mu3 mu4 nu0 nu1 nu2 nu3 ///
    nu4 pi0 pi1 pi2 pi3 pi4 pmse

  mata {
    kernC = st_data(.,"`kernArgs'")

    // run RDHonest.fit procedure

    if("`kernel'"=="optimal") _error("RDTOpt_fit not implemented yet.")
    /* ret = RDTOpt_fit(df,opt,kernC) */
    else ret = RDHonest_fit(df,opt,kernC,1)
    // output passed from Mata to Stata in the body of RDHonest_fit
  }
  // pass output from Stata to the user
  restore
  di ""
  di "Call: `y' is dependent variable; `x' is running variable."
  di "Inference by se_method: " se_method
  di "Estimate: " estimate
  di "Maximum Bias: " bias
  di "Std. Error: " stdErr
  di ""
  di "Confidence intervals:"
  di "("lEqLimit ", " uEqLimit "), " "(" lLimit ", Inf), " "(-Inf, " uLimit ")"
  di ""
  di "Bandwidth below cutoff: " hm
  di "Bandwidth above cutoff: " hp
  di "Number of effective observations: " effObs

  // grep: complete this bunch of output, emulating RDHonest implemented in R

  ereturn clear
/*
  ereturn scalar ATE = ATE
  ereturn scalar EffObs = EffObs
  ereturn scalar bias = bias
  ereturn scalar sd = sd

  mata mata clear
*/

end

// mata code begins here

version 14.2
mata:
  class RDData {
    /* this is the RD data frame object, containing the cutoff,
    the X and Y matrices before and after the cutoff, and eventually
    the supplied_vars.

    We work in normalized terms, i.e. X - cutoff. */

    /* variables */
    real matrix Xm,Xp,Ym,Yp,sigma2m,sigma2p
    real scalar cutoff

    /* functions */
    void setup()
  }

  void RDData::setup(real matrix X, real matrix Y, real scalar c,| real vector sigma2) {

    real matrix s

    /* assign matrices Xm, Xp, Ym, Yp */
    this.cutoff = c
    this.Xm = select(X,X:<c);	this.Xp = select(X,X:>=c)
    this.Xm = this.Xm :- c;	this.Xp = this.Xp :- c /* normalize X matrices */
    this.Ym = select(Y,X:<c);	this.Yp = select(Y,X:>=c)

    /* initialize supplied variances */
    if(args()==4) {
      if(length(sigma2)==length(X)) {
        this.sigma2m = select(sigma2,X:<c)
        this.sigma2p = select(sigma2,X:>=c)
      }
      else _error("length of supplied variances does not match data")
    }
    else {
      /* explicitly list variances as being missing */
      this.sigma2m = J(length(this.Xm),1,.)
      this.sigma2p = J(length(this.Xp),1,.)
    }

    /* sort data */
    s = sort((this.Xm,this.Ym,this.sigma2m),1)
    this.Xm = s[.,1];   this.Ym = s[.,2];  this.sigma2m = s[.,3]
    s = sort((this.Xp,this.Yp,this.sigma2p),1)
    this.Xp = s[.,1];   this.Yp = s[.,2];  this.sigma2p = s[.,3]

  }
  class RDOptions {
		/* this class stores all the options input by the user */

		/* variables */
		real scalar m, bw_equal, hp, hm, alpha, beta, j, order
		string kernel, opt_criterion, se_method, sclass, se_initial

		/* functions */
		void setup()
	}

  void RDOptions::setup(real scalar u_m ,| real scalar u_bw_equal,
	 	string u_kernel, real scalar u_hp, real scalar u_hm,
		string u_opt_criterion, real scalar u_alpha, real scalar u_beta,
		string u_se_method, real scalar u_j, string u_sclass,
		real scalar u_order, string u_se_initial) {

		/* assign values and defaults */
		m = u_m
		if(args() >= 2) bw_equal = (u_bw_equal > 0 ? 1 : 0)
                else bw_equal = 1

		if(args() >= 3 & u_kernel != "") kernel = u_kernel
		else kernel = "tri"

		if(args() >= 4 & u_hp != .)	hp = u_hp
		else hp = . /* need to explicitly specify bandwidth is missing */

		hm = ((bw_equal == 1 | (args() >= 5 & u_hm == .)) ? hp : u_hm)

		if(args() >= 6 & u_opt_criterion != .)	opt_criterion = u_opt_criterion
		else opt_criterion = ""
		/* need to explicitly specify opt_criterion is missing */

		if(args() >= 7)         alpha = u_alpha
		else alpha = 0.05
		if(args() >= 8)         beta = u_beta
		else beta = 0.8
		if(args() >= 9)         se_method = u_se_method
		else se_method = "nn"
		if(args() >= 10)        j = u_j
		else j = 3
		if(args() >= 11)	sclass = u_sclass
		else sclass = "H"
		if(args() >= 12)	order = u_order
		else order = 1
		if(args() >= 13)	se_initial = u_se_initial
		else se_initial = "IKEHW"

	}

  class LPRegOutput {
		/* contains output from LPReg */
		real scalar theta /* treatment effect */
		real vector sigma2 /* squared residuals */
		real vector var /* variance */
		real vector w /* observation weights */
		real scalar eo /* effective observations */
		real matrix Gamma /* to recompute weight matrix */

		real vector wf() /* weight function, akin to Kolesar's w() */
	}

  real vector LPRegOutput::wf(real matrix x, real scalar h,
		string kernel, real scalar order, real matrix u_Gamma){

		real matrix R,kernWeights, w
		R = J(1,order+1,x) :^ (0..order)
		kernWeights = EqKernWeight(x/h,kernel,0,0)
		w = ((invsym(u_Gamma) * (kernWeights :* R)')[1,])'
		return(w)
	}

  class RDLPregOutput extends LPRegOutput {

		public:
			real vector wm /* weights for negative observations */
			real vector wp /* weights for positive observations */
			real vector sigma2m, sigma2p /* squared residuals */
			real vector se /* standard error vector */
			real matrix GammaM, GammaP /* variance weight matrices */

		private:
			real vector sigma2
			real vector var
			real vector w
			real matrix Gamma
	}

  class RDData scalar RDDataPrep(real matrix X, real matrix Y, real scalar c ,| real vector sigma2) {
    /* constructs df */
    class RDData scalar df
    if(args()==4) df.setup(X,Y,c,sigma2)
    else df.setup(X,Y,c)
    return(df)
  }

  class RDOptions scalar RDOptionsPrep(real scalar u_m, real scalar u_bw_equal,
    string u_kernel, real scalar u_hp, real scalar u_hm,
    string u_opt_criterion, real scalar u_alpha, real scalar u_beta,
    string u_se_method, real scalar u_j, string u_sclass,
    real scalar u_order, string u_se_initial) {
    /* constructs opt */
    class RDOptions scalar opt
    opt.setup(u_m, u_bw_equal, u_kernel, u_hp, u_hm, u_opt_criterion, u_alpha,
      u_beta, u_se_method, u_j, u_sclass, u_order, u_se_initial)
    return(opt)
  }

  class RDData scalar IKBW_fit(class RDData scalar df, real matrix kernC,| string kernel,
		real scalar verbose) {
		/* Calculate bandwidth for sharp RD based on local linear regression
		using method by Imbens and Kalyanaraman (ReStud 2012) */
		/* note: no option for order != 1 */

		if(args() < 3)  kernel = "triangular"
		if(args() < 4)  verbose = 0

		real scalar order, Nm, Np, N, kernInd, cons, varm, varp, h2m, h2p, m2m, m2p
		real scalar h1, f0, m3, rm, rp
		real vector col1, col2, col3, col4, col5, m3coeff, tempY
		real matrix s,X, m3mat, tempX, m2mat

		X = df.Xm\df.Xp
		order = 1
		Nm = length(df.Xm)
		Np = length(df.Xp)
		N = Nm + Np

		/* kernels must be specified as strings */

		/* translate kernel strings into kernel indices */
		kernInd = ((strpos(kernel,"uni")>0) ? 0 :
			(strpos(kernel,"tri")>0 ? 1 : 2))

		/* access kernel constants */
		/* select kernC row with matching kernel index, order == 1 and
			boundary == TRUE */
		s = select(kernC,kernC[.,1]:==kernInd)
		s = select(s,s[.,2]:==1)
		s = select(s,s[.,3]:==1)
		s = s'
		/* s is now a colvector, so we can avoid complications with variance() */

		cons = (s[9]/(s[6]^2))^(1/5)
		df = RDPrelimVar(df, kernC, "Silverman")
		h1 = 1.84 * sqrt(variance(X))/(N^(1/5))
		f0 = sum(abs(X) :<= h1)/(2 * N * h1)
		varm = df.sigma2m[1]; varp = df.sigma2p[1]
		/* run a linear regression of y on 1{x >= 0}, x, x^2 and x^3, then extract
		coefficients of x^3 */
		col1 = X :^ 0; col2 = (X :>= 0); col3 = X; col4 = X :^ 2;  col5 = X :^ 3
		m3mat = col1, col2, col3, col4, col5
		m3coeff = invsym(m3mat'm3mat)*m3mat'(df.Ym\df.Yp)
		m3 = 6 * m3coeff[5]

		h2m = 7200^(1/7) * (varm/(f0 * m3^2))^(1/7) * Nm^(-1/7)
		h2p = 7200^(1/7) * (varp/(f0 * m3^2))^(1/7) * Np^(-1/7)

		/* run a linear regression of y_m on x_m and x_m^2, amongst obs that lie
		within the bandwidth */

		tempY = select(df.Ym,df.Xm :>= -h2m)
		tempX = select(df.Xm,df.Xm :>= -h2m)
		col1 = tempX :^ 0; col2 = tempX; col3 = tempX :^ 2
		m2mat = col1, col2, col3
		m2m = 2 * (invsym(m2mat'm2mat) * m2mat'tempY)[3]

		/* same procedure: y_p on x_p and x_p^2, within bandwidth */

		tempY = select(df.Yp,df.Xp :<= h2p)
		tempX = select(df.Xp,df.Xp :<= h2p)
		col1 = tempX :^ 0; col2 = tempX; col3 = tempX :^ 2
		m2mat = col1, col2, col3
		m2p = 2 * (invsym(m2mat'm2mat) * m2mat'tempY)[3]

		rm = 2160 * varm/(sum(df.Xm :>= -h2m) * h2m^4)
		rp = 2160 * varp/(sum(df.Xp :<= h2p) * h2p^4)

		/* output parameters if verbose */
		if (verbose) {
			printf("\n h1: %9.0g \n N_{-}: %9.0g \n  N_{+}: %9.0g", h1, Nm, Np)
			printf("\n f0: %9.0g \n sigma^2_{+}(0): %9.0g^2 \n sigma^2_{-}(0): %9.0g^2",
				f0, sqrt(varp), sqrt(varm))
			printf("\n m3: %9.0g \n h_{2,+}: %9.0g \n h_{2,-}: %9.0g",m3,h2p,h2m)
			printf("\n m^{(2)}_{+}: %9.0g \n m^{(2)}_{-}: %9.0g",m2p,m2m)
			printf("\n r_{+}: %9.0g \n r_{-}: %9.0g \n\n",rp,rm)
		}
		return(cons * ((varp + varm)/(f0 * N * ((m2p - m2m)^2 + rm + rp)))^(1/5))
	}

  class RDData scalar RDPrelimVar(class RDData scalar df, real matrix kernC,| string se_initial) {
		/* Computes preliminary estimates of variances, to be used in
		optimal bandwidth calculations. */

		real matrix X, matm, matp
		class RDLPregOutput scalar r1
		real scalar h1

		/* set defaults */
		if(args()==2) se_initial = "IKEHW"

		/* concatenate Xm and Xp and compute "rule of thumb" bandwidth */
		X = df.Xm\df.Xp
		h1 = max((1.84*sqrt(variance(X))/(sum(length(X)))^(1/5),
			uniqrows(df.Xp)[2],uniqrows(abs(df.Xm))[2]))

		if(strpos(se_initial,"Silverman") > 0 & strpos(se_initial,"SilvermanNN") == 0) {
			matm = select(df.Ym,df.Xm :>= -h1)
			df.sigma2m = J(length(df.Xm),1,variance(matm))
			matp = select(df.Yp,df.Xp :<= h1)
			df.sigma2p = J(length(df.Xp),1,variance(matp))
		}
		else if(strpos(se_initial,"SilvermanNN") > 0) {
			r1 = RDLPreg(df,h1,"uniform",1,h1,"nn")
			df.sigma2m = J(length(df.Xm),1,mean(r1.sigma2m))
			df.sigma2p = J(length(df.Xp),1,mean(r1.sigma2p))
		}
		else if(strpos(se_initial,"IKdemeaned") > 0) {
			h1 = IKBW_fit(df, kernC)
			r1 = RDLPreg(df,h1,"tri",1,h1,"demeaned")
			df.sigma2m = J(length(df.Xm),1,mean(r1.sigma2m))
			df.sigma2p = J(length(df.Xp),1,mean(r1.sigma2p))
		}
		else if(strpos(se_initial,"IKEHW") > 0) {
			h1 = IKBW_fit(df, kernC)
			r1 = RDLPreg(df,h1,"tri",1,h1,"EHW")
			df.sigma2m = J(length(df.Xm),1,mean(r1.sigma2m))
			df.sigma2p = J(length(df.Xp),1,mean(r1.sigma2p))
		}
		else if(strpos(se_initial,"NN") > 0 & strpos(se_initial,"SilvermanNN") == 0) {
			df.sigma2p = sigmaNN(df.Xp,df.Yp,3)
			df.sigma2m = sigmaNN(df.Xm,df.Ym,3)
		}
		else {
			_error("Unknown method for estimating initial variance.")
		}

		return(df)
	}

	/* v2 stata impl. of LPReg */

  class LPRegOutput scalar LPReg (real matrix X, real matrix Y, real scalar h,
		string kernel, real scalar order, string se_method, real scalar j,
		| real matrix u_sigma2) {

    /* grep: need to debug for sclass H order 2 */

		real scalar ord, effObs, v
		real matrix R, W, Gamma, beta, hsigma2, dsigma2, nsigma2
		class LPRegOutput scalar output

		ord = order

                /* printf("The bandwidth for LPRegOutput is %9.0g \n", h) */

		R = J(1,ord+1,X) :^ (0..ord)
		W = EqKernWeight(X/h,kernel,0,0)             
                /* st_matrix("W",W) */
                /* stata("matlist W") */
		Gamma = quadcross(R,W :* R)

		if(sum(select(W,W:>0)) <= ord | h == 0| det(Gamma)==0) {
			printf("{red}Error: bandwidth too small or singular matrix.\n")
			exit(error(3352))
		}

		beta = invsym(Gamma)* quadcross(R,W :* Y)
		hsigma2 = (Y - R*beta) :^2

		/* Robust variance-based formulae */

    /* printf("LPReg: Now checking for variance types \n") */

		if(args()==8 & strpos(se_method,"supplied_var") > 0) {
      /* printf("Entering LPReg supplied variance conditional \n") */
                        v = EHW(Gamma,W,R,u_sigma2)
                }

		if(strpos(se_method,"EHW") > 0 | strpos(se_method,"ehw") > 0) {
			v = EHW(Gamma,W,R,hsigma2)
		}

		if(strpos(se_method,"demeaned") > 0) {
			dsigma2 = (Y:-beta[1]) :^2
			v = EHW(Gamma,W,R,dsigma2)
		}

		if(strpos(se_method,"nn") > 0 | strpos(se_method,"NN") > 0 ) {
			nsigma2 = sigmaNN(X,Y,j)
			v = EHW(Gamma,W,R,nsigma2)
		}

    /* catch error when no variance is provided */
    if(v[1]==.) _error("No variance was computed by LPReg")

		effObs = 1/sum(linWeight(X,h,kernel,ord,Gamma):^2)

		/* return output */
		output.theta = beta[1]
		output.sigma2 = hsigma2
		output.var = v
		output.w = W
		output.eo = effObs
		output.Gamma = Gamma
		return(output)
	}

	/* Implementation of RDLPreg in mata */
  class RDLPregOutput scalar RDLPreg(class RDData scalar df, real scalar hp,
	|string kernel, real scalar order, real scalar hm, string se_method,
		real scalar no_warning, real scalar j) {

    /* grep: need to debug for sclass H order 2 */

		/* check for errors */
		if (hp <= 0) _error("Non-positive bandwidth hp")

	 	/* set defaults */
		if (args() <= 2 | kernel == "") kernel = "tri"
		if (args() <= 3 | order == .) order = 1
	  if (args() <= 4 | hm == .)	hm = hp  /* default to equate hm and hp */
	  if (args() <= 5 | se_method == "") se_method = "nn"
		if (args() <= 6 | no_warning == .) no_warning = 0
		if (args() <= 7 | j == .) j = 3

		/* variable declarations */
		real scalar plugin /* to be implemented later */
		real vector Wm, Wp, Ym, Yp, sigma2m, sigma2p
		real matrix Xm, Xp
		class LPRegOutput scalar rm, rp
		class RDLPregOutput scalar output

		/* set kernel weights */
		Wm = (hm <= 0 ? 0 :* df.Xm : EqKernWeight(df.Xm/hm,kernel,0,0))
		Wp = (hp <= 0 ? 0 :* df.Xp : EqKernWeight(df.Xp/hp,kernel,0,0))

		/* variance calculations are faster if we only keep data with positive
	 	kernel weights */

		Xp = select(df.Xp,Wp:>0);	Xm = select(df.Xm,Wm:>0)
		Yp = select(df.Yp,Wp:>0);	Ym = select(df.Ym,Wm:>0)
    sigma2p = select(df.sigma2p,Wp:>0); sigma2m = select(df.sigma2m,Wm:>0)

		if ((length(Xm) < 3*order | length(Xp) < 3*order) & !no_warning) {
			printf("{red}Warning: Too few observations to compute RD estimates.\n")
			printf("{red}Only %g control and %g treated units with positive weights.\n",
			length(Xm),length(Xp))
		}

		if ((length(uniqrows(Xm)) <= order | length(uniqrows(Xp)) <= order) & !no_warning) {
			printf("{red}Warning: Too few distinct values to compute RD estimates.\n")
			printf("{red}Only %g unique control and %g unique treated values, \n",
			length(uniqrows(Xm)),length(uniqrows(Xp)))
			printf("{red}with positive weights, for the running variable. \n")
		}

    if(strpos(se_method,"supplied_var") > 0) {
      /* diagnostics */
      /* printf("{red}entering RDLPreg supplied variance conditional \n")
      printf("truncated sigma2m has first entry %9.0g \n", sigma2m[1]) */
      rm = LPReg(Xm,Ym,hm,kernel,order,se_method,j,sigma2m)
  		rp = LPReg(Xp,Yp,hp,kernel,order,se_method,j,sigma2p)
      /* printf("rm.var has value %9.0g \n", rm.var[1]) */
    }
    else {
      rm = LPReg(Xm,Ym,hm,kernel,order,se_method,j)
  		rp = LPReg(Xp,Yp,hp,kernel,order,se_method,j)
    }

		output.theta = rp.theta - rm.theta
		/* implement plugin asymptotic variance later */
		plugin = .
		output.se = (rp.var+rm.var,plugin) :^ 0.5
		output.eo = rm.eo + rp.eo
		output.wm = rm.w;			output.wp = rp.w
		output.sigma2m = rm.sigma2;		output.sigma2p = rp.sigma2
		output.GammaM = rm.Gamma;		output.GammaP = rp.Gamma

		return(output)
	}

	/* define struct for results from RDHonest_fit */

  struct RDResults {
		real scalar estimate
		real scalar lff /* not sure what this is */
		real scalar bias, sd, lower, upper, hl, eo, hp, hm, naive
	}

	/* implement RDHonest_fit */

  struct RDResults scalar RDHonest_fit (class RDData scalar df,
		class RDOptions scalar opt, real matrix kernC, real scalar returnFlag) {
    /* returnFlag is either 0 or 1. 1 means we send to Stata the results of
    the computation, while 0 means we do not. */

		/* declare variables */
		class RDLPregOutput scalar r1
		real scalar z, hp, hm
		real vector wp, wm, hVec
		struct RDResults scalar results

		hp = opt.hp
		hm = opt.hm

		/* initial se estimate */

    if ((mean(df.sigma2p :== .) == 1 | mean(df.sigma2m :== .) == 1)
      & (strpos(opt.se_method,"supplied_var") > 0 | hp == .)) {
      printf("Generating initial variance estimates via RDPrelimVar \n")
      df = RDPrelimVar(df,kernC,opt.se_initial)
    }

		/* deal with missing bandwidths */
    if(opt.hp <= 0) {
        // hacky fix: previously opt.hp == .
      printf("Bandwidth (hp) missing or invalid. Running RDOptBW_fit \n")
      hVec = RDOptBW_fit(df,opt,kernC)
      /* will need: m, kernel, opt_criterion, bw.equal, alpha,
        beta, sclass, order */
      hp = hVec[1]
      hm = hVec[2]
		}

		/* Suppress warnings about too few observations */
		r1 = RDLPreg(df,hp, opt.kernel, opt.order, hm, opt.se_method, 1, opt.j)
		wp = r1.wf(df.Xp,hp,opt.kernel,opt.order,r1.GammaP)
		wm = r1.wf(df.Xm,hm,opt.kernel,opt.order,r1.GammaM)

		/* If bandwidths are too small, set bias to machine maximum */
		if (sum(wp :>= 0) == 0 | sum(wm :>= 0) == 0) {
			results.bias = results.sd = results.upper = results.hl = sqrt(1.0X+4e/10)
			results.lower = -results.upper
		}
		else {
			results.sd = r1.se[1]
			if (opt.sclass == "T") results.bias = opt.m/2 *
				(sum(abs(wp :* (df.Xp :^2))) + sum(abs(wm :* (df.Xm :^2))))
			else if(opt.sclass == "H" & opt.order == 1) results.bias = -opt.m/2 *
				(sum(wp :* df.Xp :^2) + sum(wm :* df.Xm :^2))
			else if(opt.sclass == "H" & opt.order == 2) {
				/* Need to find numerically. Unlike in the R code, we will use
                                brute-force grid search */
        // noisy debug
        // printf("Numerically minimizing for sclass H, order 2 \n")
        real scalar objP, objM
        objP = gridSearchOptim(&bp(),0,opt.hp,wp,df)
        // printf("Found objP, which is %9.0g \n", objP)
        objM = gridSearchOptim(&bm(),0,opt.hm,wm,df)
        // printf("Found objM, which is %9.0g \n", objM)
        /* need to debug */
        results.bias = -(opt.m/2) * (objP + objM)
        // printf("Found bias for sclass H, order 2 \n")
        // printf("results.bias is %9.0g \n", results.bias)
        // grep: gives bias of zero, which is incongruent with result from R
			}
			else {
				_error("Don't know how to compute bias for specified sclass and order.")
			}

                        results.lower = r1.theta - results.bias - invnormal(1-opt.alpha) * results.sd
			results.upper = r1.theta + results.bias + invnormal(1-opt.alpha) * results.sd
			/* find hl directly using CVb */
			results.hl = (CVb(results.bias/results.sd, opt.alpha)) * results.sd
      // note to self: hl indicates "distance" from estimate to lower bound
      // and from estimate to upper bound of symmetric 95% honest CI
		}

		/* Finally, calculate coverage of naive CIs */
		z = invnormal(1-opt.alpha/2)
		results.naive = normal(z-results.bias/results.sd) - normal(-z-results.bias/results.sd)

		/* input known values into results */
		results.estimate = r1.theta
		results.hp = hp
		results.hm = hm
		results.eo = r1.eo
    // printf("Number of effective observations is %9.0g \n", results.eo)

    /* if this is the final iteration of RDHonest_fit,
      send results back to Stata as numerical scalars */
    if (returnFlag == 1)  {
      // *** numerical scalars ***
      st_numscalar("estimate",results.estimate)
      st_numscalar("bias",results.bias)
      st_numscalar("stdErr",results.sd)
      st_numscalar("hp",results.hp)
      st_numscalar("hm",results.hm)
      st_numscalar("effObs",results.eo)
      // one-sided confidence interval limits
      st_numscalar("lLimit",results.lower)
      st_numscalar("uLimit",results.upper)
      // symmetric two-sided confidence interval limits
      st_numscalar("lEqLimit",results.estimate - results.hl)
      st_numscalar("uEqLimit",results.estimate + results.hl)
      // *** strings from specified options ***
      st_strscalar("se_method",opt.se_method)
      // grep: add string scalar outputs as required
    }

		return(results)
	}

  // accessory functions for sclass H and order 2 bias calculations
  real scalar bp(real scalar b, class RDData scalar df, real vector wp) {

    real scalar i, val
    real vector pointMaxSq

    pointMaxSq = J(length(df.Xp),1,0)
    for(i = 1; i <= length(pointMaxSq); i++) pointMaxSq[i] = max((df.Xp[i]-b,0))^2

    val = -2 * sum(wp :* pointMaxSq)
    return(val)
    // to debug
  }

  real scalar bm(real scalar b, class RDData scalar df, real vector wm) {

    real scalar i, val
    real vector pointMaxSq

    pointMaxSq = J(length(df.Xm),1,0)
    for(i = 1; i <= length(pointMaxSq); i++) pointMaxSq[i] = max((abs(df.Xm[i])-b,0))^2

    val = -2 * sum(wm :* pointMaxSq)
    return(val)
  }

  /* code for RDOptBW and accessory functions */

  /* begin with objective function for optimizing bandwidth when bw_equal == 0 */

  void RDOptBWEvalNotEqual (real scalar todo, real vector b, class RDData scalar df,
    class RDOptions scalar opt, real matrix kernC, val, grad, hess) {

    struct RDResults scalar r
    class RDOptions scalar optLocal

    /* copy local options from opt but set hp, hm and se_method differently */
    /* printf("current values of hp and hm are %9.0g and %9.0g \n", b[1], b[2]) */

    optLocal.setup(opt.m, opt.bw_equal, opt.kernel, abs(b[1]), abs(b[2]), opt.opt_criterion,
      opt.alpha, opt.beta, "supplied_var", opt.j, opt.sclass, opt.order, opt.se_initial)

    /* the vector h is supposed to have two dimensions: h[1] is hp and h[2] is hm */
    r = RDHonest_fit(df,optLocal,kernC,0)
    
    if (strpos(optLocal.opt_criterion,"OCI") > 0) val = 2*r.bias + r.sd*(
      invnormal(1-optLocal.alpha) + invnormal(optLocal.beta))
    else if (strpos(optLocal.opt_criterion,"MSE") > 0) val = (r.bias)^2 + (r.sd)^2
    else if (strpos(optLocal.opt_criterion,"FLCI") > 0) val = r.hl
    else {
      printf("optimality criterion %s not implemented yet \n", optLocal.opt_criterion)
      _error("exiting RDOptBWEvalNotEqual()")
    }
  }

  /* now write objective function for optimizing bandwidth when bw_equal == 1 */
  void RDOptBWEval (real scalar todo, real scalar b, class RDData scalar df,
    class RDOptions scalar opt, real matrix kernC, real scalar hmin, real scalar hmax,
    val, grad, hess) {

    struct RDResults scalar r
    class RDOptions scalar optLocal

    /* copy local options from opt but set hp, hm and se_method differently */
    /* printf("current value of hp is %9.0g \n", b) */

    /* check to make sure bandwidth is feasible, else kicks the function back into range */
    if (b < hmin) b = hmin + 0.2*(hmax-hmin)/2
    else if (b > hmax) b = hmax - 0.2*(hmax-hmin)/2

    optLocal.setup(opt.m, opt.bw_equal, opt.kernel, abs(b), abs(b), opt.opt_criterion,
      opt.alpha, opt.beta, "supplied_var", opt.j, opt.sclass, opt.order, opt.se_initial)

    r = RDHonest_fit(df,optLocal,kernC,0)
    if (strpos(optLocal.opt_criterion,"OCI") > 0) val = 2*r.bias + r.sd*(
      invnormal(1-optLocal.alpha) + invnormal(optLocal.beta))
    else if (strpos(optLocal.opt_criterion,"MSE") > 0) val = (r.bias)^2 + (r.sd)^2
    else if (strpos(optLocal.opt_criterion,"FLCI") > 0) val = r.hl
    else {
      printf("optimality criterion %s not implemented yet \n", optLocal.opt_criterion)
      _error("exiting RDOptBWEval()")
         }
  }

  /* objective function for golden section search */

  real scalar RDOptBWGss (real scalar h, class RDData scalar df,
    class RDOptions scalar opt, real matrix kernC) {

    struct RDResults scalar r
    class RDOptions scalar optLocal
    real scalar val

    /* To help debug: bw_equal = TRUE trivially */
    optLocal.setup(opt.m, opt.bw_equal, opt.kernel, h, h, opt.opt_criterion,
      opt.alpha, opt.beta, "supplied_var", opt.j, opt.sclass, opt.order, opt.se_initial)

    r = RDHonest_fit(df,optLocal,kernC,0)

    if (strpos(optLocal.opt_criterion,"OCI") > 0) val = 2*r.bias + r.sd*(
      invnormal(1-optLocal.alpha) + invnormal(optLocal.beta))
    else if (strpos(optLocal.opt_criterion,"MSE") > 0) val = (r.bias)^2 + (r.sd)^2
    else if (strpos(optLocal.opt_criterion,"FLCI") > 0) val = r.hl
    else {
      printf("optimality criterion %s not implemented yet \n", optLocal.opt_criterion)
      _error("exiting RDOptBWGss()")
    }
    return(val)
  }

  // now implement function to fit optimal bandwidth

  real vector RDOptBW_fit (class RDData scalar df, class RDOptions scalar opt,
    real matrix kernC) {

    transmorphic S /* to optimize function */
    real scalar h0, hmin, hmax, hp, hm // ROT guess, minimal and maximal h
    real vector h_opt, supp
    real matrix X

    /* First check if sigma2 is supplied */
    if (mean(df.sigma2p :== .) == 1
      | mean(df.sigma2m :== .) == 1) df = RDPrelimVar(df,kernC,opt.se_initial)

    /* compute Silverman rule of thumb bandwidth */
    X = df.Xm\df.Xp
    h0 = 1.06 * sqrt(variance(X)) * length(X)^(-1/5)

    /* initialize bandwidth optimization */
    S = optimize_init()
    optimize_init_argument(S,1,df)
    optimize_init_argument(S,2,opt)
    optimize_init_argument(S,3,kernC)
    optimize_init_which(S, "min")

    if (opt.bw_equal == 0) {
      optimize_init_evaluator(S,&RDOptBWEvalNotEqual())
      optimize_init_technique(S,"bfgs")
      optimize_init_params(S, (max(df.Xp)/2,max(abs(df.Xm))/2))
      h_opt = optimize(S)
      hp = h_opt[1]
      hm = h_opt[2]
    }
    else {
      hmin = max((uniqrows(df.Xp)[opt.order+1],uniqrows(abs(df.Xm))[opt.order+1]))
      hmax = max(abs(df.Xp\df.Xm))
      /* Optimize piecewise-constant function using modified golden section search.
      Careful: criterion may not be unimodal (but appears to be so, even for
       triangular kernel) */
      if (strpos(opt.kernel,"uni") > 0) {
         /* printf("Entering golden section search \n") */
         supp = sort(uniqrows(df.Xp\abs(df.Xm)),1)
         /* changed select(supp,supp :> hmin) */
         hp = hm = gss(&RDOptBWGss(), select(supp,supp :>= hmin), df, opt, kernC)
         /* printf("Golden section bandwidth obtained is %9.0g \n",hp) */
      }
      else {
        // set up parameters of optimization problem
        optimize_init_argument(S,4,hmin)
        optimize_init_argument(S,5,hmax)
        optimize_init_evaluator(S,&RDOptBWEval())
        optimize_init_params(S,h0)
        optimize_init_technique(S,"bfgs")
        h_opt = optimize(S)
        hp = hm = h_opt
      }
    }
    /* A REMPLIR */
    return((hp,hm))
  }

	/*  *** ACCESSORY FUNCTIONS *** */
  // 

  // modified golden section optimizer to use with mata objective functions
  real scalar gss(pointer scalar f, real vector xs, class RDData scalar df,
    class RDOptions scalar opt, real matrix kernC) {

    real scalar gr, a, b, c, d, i, h_opt, i_opt, val
    real vector supp

    /* st_matrix("xs",xs)
    stata("matlist xs") */
    
    gr = (sqrt(5)+1)/2
    a = 1
    b = length(xs)
    c = round(b - (b - a) / gr)
    d = round(a + (b - a) / gr)

    while (b-a > 100) {
      if ((*f)(xs[c],df,opt,kernC) < (*f)(xs[d],df,opt,kernC)) b = d
      else a = c
      // recompute both c and d
      c = round(b - (b - a) / gr)
      d = round(a + (b - a) / gr)
    }
    supp = xs[a..b]
    
    /* st_matrix("supp",supp) */
    /* stata("matlist supp") */
    /* printf("The support we get is from %9.0g to %9.0g \n",a,b) */

    // iteratively apply f to support points to find optimal bandwidth

    h_opt = 1.0X+4e
    for (i=1; i<=length(supp); i++) {
      val = (*f)(supp[i],df,opt,kernC)
      /* printf("The new value of val is %9.0g \n",val) */
      if(val <= h_opt) {
        h_opt = val
        i_opt = i
      }
      /* printf("The optimal index we get is %9.0g \n",i_opt) */
    }
    return(supp[i_opt])
  }

// implement CVb with noncentral chi-squared distribution

real scalar CVb(real scalar B,| real scalar alpha) {
    if (args()==1) alpha = 0.05
    
    // declare variables
    real scalar root
    
    // check for null bias/alpha out of range
    if(B == .) return(.)
    if(B < 0 | alpha <= 0 | alpha >= 1) _error("B negative or alpha out of range")

    root = sqrt(invnchi2(1,B^2,1-alpha))
    return(root)
    }

// cv is depreciated. objective function and original implementation of cv follow
  void cvObj(todo, c, real scalar B, real scalar alpha, val, grad, hess) {
    val = (normal(c-B) - normal(-c-B) - (1-alpha))^2
  }

  real scalar cv(real scalar B,| real scalar alpha) {
    // noisy debug
    // printf("running function cv \n")
    // set default value of alpha
    if (args()==1) alpha = 0.05

    // declare variables
    real scalar root
    transmorphic S

    // check for null bias/alpha out of range
    if(B == .) return(.)
    if(B < 0 | alpha <= 0 | alpha >= 1) _error("B negative or alpha out of range")

    // initialize optimization
    S = optimize_init()
    optimize_init_evaluator(S,&cvObj())
    optimize_init_argument(S,1,B)
    optimize_init_argument(S,2,alpha)
    optimize_init_params(S,B)
    optimize_init_technique(S, "bfgs")
    optimize_init_which(S, "min")
    root = optimize(S)

		return(root)
	}

	/* EHW function */

  real scalar EHW(real matrix Gamma, real matrix W, real matrix R,
		real matrix sigma2) {
		return(quadcross((sigma2 :* W :* R)*invsym(Gamma),
			(W :* R) * invsym(Gamma))[1,1])
	}

	/* to implement weight function, thinking of estimator as linear */

  real vector linWeight(x,h,kernel,order,Gamma) {
		real matrix R,kernWeights, w
		R = J(1,order+1,x) :^ (0..order)
		kernWeights = EqKernWeight(x/h,kernel,0,0)
		w = ((invsym(Gamma) * (kernWeights :* R)')[1,])'
		return(w)
	}

	/* enhance rdrobust_kweight() to EqKern() standards by
		including order and boundary parameters */

  real vector EqKernWeight(real matrix Xh, string kernel,
		 real scalar order, real scalar boundary)
	{
		/* input Xh as X/h */
		real matrix u,w
		u = Xh
		if (boundary == 0 & order == 0) {
			if (kernel=="epanechnikov"|kernel=="epa") {
				w = (0.75:*(1:-u:^2):*(abs(u):<=1))
			}
			else if(kernel=="uniform"|kernel=="uni") {
				w = (0.5:*(abs(u):<=1))
			}
			else if(kernel=="triangular"|kernel=="tri"){
				w = (1:-abs(u)):*(abs(u):<=1)
			}
			else {
				printf("{red}Error: Unsupported kernel type\n")
				exit(error(198))
			}
		}
		else if (boundary == 0 & order == 1) { /* check if this is correct */
			if (kernel=="epanechnikov"|kernel=="epa") {
				w = (0.75:*(1:-u:^2):*(abs(u):<=1))
			}
			else if(kernel=="uniform"|kernel=="uni") {
				w = (0.5:*(abs(u):<=1))
			}
			else if(kernel=="triangular"|kernel=="tri"){
				w = (1:-abs(u)):*(abs(u):<=1)
			}
			else {
				printf("{red}Error: Unsupported kernel type\n")
				exit(error(198))
			}
		}
		else if (boundary == 0 & order == 2) { /* check if this is correct */
			if (kernel=="epanechnikov"|kernel=="epa") {
				w = (15/32) :* (3:-7:*u:^2) :* (1:- u:^2):*(abs(u):<=1)
			}
			else if(kernel=="uniform"|kernel=="uni") {
				w = ((9 :- 15 :* u:^2) / 8) :* (abs(u):<=1)
			}
			else if(kernel=="triangular"|kernel=="tri"){
				w = 6/7 :* (2:-5:*u:^2) :* (1:-abs(u)) :* (abs(u):<=1)
			}
			else {
				printf("{red}Error: Unsupported kernel type\n")
				exit(error(198))
			}
		}
		else if (boundary == 1 & order == 0) {
			if (kernel=="epanechnikov"|kernel=="epa") {
				w = 1.5 :* (1 :- u:^2):*(abs(u):<=1)
			}
			else if(kernel=="uniform"|kernel=="uni") {
				w = (abs(u):<=1)
			}
			else if(kernel=="triangular"|kernel=="tri"){
				w = 2 * (1 :- u) :* (abs(u):<=1)
			}
			else {
				printf("{red}Error: Unsupported kernel type \n")
				exit(error(198))
			}
		}
		else if (boundary == 1 & order == 1) {
			if (kernel=="epanechnikov"|kernel=="epa") {
				w = 6/19 * (16 :- 30*u) :* (1 :- u:^2) :* (abs(u):<=1)
			}
			else if(kernel=="uniform"|kernel=="uni") {
				w = (4 :- 6 :* u) :* (abs(u):<=1)
			}
			else if(kernel=="triangular"|kernel=="tri"){
				w = 6 * (1 :- 2 :* u) :* (1 :- u) :* (abs(u):<=1)
			}
			else {
				printf("{red}Error: Unsupported kernel type\n")
				exit(error(198))
			}
		}
		else if (boundary == 1 & order == 2) {
			if (kernel=="epanechnikov"|kernel=="epa") {
				w = 1/8 * (85 :- 400*u :+ 385*u:^2) :* (1 :- u:^2) :* (abs(u):<=1)
			}
			else if(kernel=="uniform"|kernel=="uni") {
				w = (9 :- 36*u :+ 30*u:^2) :* (abs(u):<=1)
			}
			else if(kernel=="triangular"|kernel=="tri"){
				w = 12 * (1 :- 5*u :+ 5*u:^2) :* (1 :- u) :* (abs(u):<=1)
			}
			else {
				printf("{red}Error: Unsupported kernel type\n")
				exit(error(198))
			}
		}
		else {
			printf("{red}Error: Invalid order or boundary specified.\n")
			exit(error(198))
		}
		return(w)
	}


	/* generate KNN standard errors */

  real vector sigmaNN(real matrix X, real matrix Y, real scalar j) {
    /* assume X is sorted */
    real scalar n,k,d,Jk
    real matrix sigma2,s,ind,YTemp

    n = length(X)
    sigma2 = J(n,1,.)

    for(k = 1; k <= n; k++) {
      s = max((k-j,1))..k

      if (length(s)==1) {
        /* avoid error when trying to join a nonexistent vector
        with one that exists */
        d = (sort(abs(X[(k+1)..min((k+j,n))]:-X[k]),1))[j]
      }
      else if (k == n) {
        /* avoid error when trying to join a nonexistent vector
        with one that exists */
        d = (sort(abs((X[s[1..length(s)-1]]):-X[k]),1))[j]
      }
      else {
        d = (sort(abs((X[s[1..length(s)-1]]\X[(k+1)..min((k+j,n))]):-X[k]),1))[j]
      }
      ind = (abs(X :- X[k]) :<= d)
      ind[k] = 0
      Jk = sum(ind)
      /* compute mean of Y with positive indices */
      YTemp = (sort((ind,Y),-1))[1..Jk,2]
      sigma2[k] = Jk/(Jk+1) * (Y[k] - mean(YTemp))^2
    }
    return(sigma2)
  }

  // brute force grid search over the desired interval
  real scalar gridSearchOptim(pointer scalar f, real scalar lower, real scalar upper,
    real vector w, class RDData scalar df) {

    real scalar delta, k, j
    real vector obj, arg, optMat
    real matrix s

    // k is the number of grid points
    k = 1000
    // attempt to replicate the seq() function in R
    delta = (upper-lower)/(k-1)
    s = ((0::(k-1)) :* delta) :+ lower

    // initialize pointwise objective functions and arg-optimizers
    arg = 1::k
    obj = J(k,1,0)
    for(j = 1; j <= k; j++) {
      // optimize pointwise
      obj[j] = (*f)(s[j],df,w)
    }
    // sort entries by size and find minimizing index
    optMat = sort((arg,obj),2)

    return(optMat[1,2])
  }
end
