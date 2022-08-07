*! version 1  
*! CYG




program define rdhonestV1, eclass

version 14.2

    	syntax [anything(name=0)] [if] [in]  ///
		   [,m(numlist >0 min=1)  ///
		   c(real 0) /// 
		   h(real 0) ///
		   kernel(string) ///
		   opt_criterion(string) ///
		   alpha(real 0.05) beta(real 0.8) ///
		   se_method(string) ///
           j(integer 3) ///
		   cluster(varname) sigma2(varlist) ///
		   t0(real 0) ///
		   weight(varname)] ///
		   //order(integer 1)///]
	
               /*  find y x treatment and covariates  */
	
	//       DepVar runingVar TreatmentVar CoVar 
	// RD:    y        runv                covar 
    // FRD:   y        runv      trtmnt    covar 
	
*********************************************
//          sub-programs from ivregress ados 
//          extracting endog exog iv and y
//          _iv_parse 
*********************************************
	
    _iv_parse `0'
	
   	local y `s(lhs)'
	local endog `s(endog)'
	local exog `s(exog)'
	local exog1 : word 1 of `exog'
	cap local exog2 : subinstr local exog "`exog1'" ""
	local inst `s(inst)'
	local 0 `s(zero)'
	
	// which case
	// sharp rdd---Srd = 1
	
	local Srd  `=("`inst'"=="")'
		
	// 
	local runv `=cond(`Srd',"`exog1'","`inst'")'
	local covar `=cond(`Srd',"`exog2'","`exog'")'
	local trtmnt `endog'
	
	// No covariates---NoCovar =1
	local NoCovar `=("`covar'"=="")'
	
	// Use-suuplied PreVar --- SelfPreVar=1
	// sigma2 can be a matrix
	local SelfPreVar `="`sigma2'" != ""'
	
	// cluster specified --- Cluster=1
	local Cluster `="`cluster'" != ""'
	
	// weights of obs provided --- Weight=1
	local Weight `="`weight'" != ""'	
	
		   
		              /* mark the data sample*/
						   
	      	 
  marksample touse
  markout `touse' `y' `runv' `covar' `trtmnt'  `cluster' `sigma2' `weight', strok
  preserve
  qui keep if `touse'
  
  
  // set more off for now
  set more off

  
                        /* set defaults for string options*/
 
  if ("`kernel'" == "") local kernel "tri"
  if ("`opt_criterion'" == "") local opt_criterion "MSE"
  if ("`se_method'" == "") local se_method "nn"
  if ("`m'" == "") local m = -1
   
                         
  
   
	                       /* Mata */
 set matastrict on
  
  
  mata{
  	
    Y = st_data(.,("`y'"),0);
	X = st_data(.,("`runv'"),0)	
	if (!`Srd') Y = st_data(.,("`y'","`trtmnt'"),0);
	if (!`NoCovar') Z = st_data(.,("`covar'"),0);	
	
	Class = (`Srd'>0 ? "srd":"frd")
	
	if (`SelfPreVar') sigma2 = st_data(.,("`sigma2'"),0); 
	if (!`SelfPreVar') sigma2 = J(rows(Y),cols(Y)^2,.);
	if (`Cluster') cluster = st_data(.,("`cluster'"),0);
	if (!`Cluster') cluster = J(rows(Y),1,1);
	if (`Weight') weight = st_data(.,("`weights'")); 
	if (!`Weight') weight = J(rows(Y),1,1);
  
    /* initialize data frame and map Y and X in */
    
	df = RDDataPrep(X,Y,`c',sigma2,weight,cluster,Class)

    /* printf("`se_method', `j', `sclass', `order', `se_initial'") */

    /* initialize options */
	m=(`=subinstr("`m'"," ",",",.)')
	kernel = "`kernel'"
    opt_criterion = "`opt_criterion'"
    alpha = `alpha'
    beta = `beta'
    se_method = "`se_method'"
    j = `j'
    order = 1 // `order'
	h = `h'
	t0 = `t0'

	
    opt = RDOptionsPrep(m,      ///
	                    kernel, ///
						h, ///
						opt_criterion, ///
                        alpha, ///
						beta, ///
						se_method, ///
						j, ///
						order,
						t0)  

    ret = RDResults()

    if("`kernel'"=="optimal") _error("RDTOpt_fit not implemented yet.")
  
    else ret = RDHonest_fit(df,opt,kernC)
  
  }

  
  restore
  
						/* Output display */
  
  
  //  results taken from Mata to Stata
  Display `y' `runv' `Srd' `alpha' `beta' estimate bias stdErr lEqLimit uEqLimit lLimit uLimit h ///
`opt_criterion' `se_method' effObs Leverage M2 fs M1 `trtmnt' 



                            /* ereturn */
 

     /* esample */
	eret post, esample(`touse')
	//eret post, esample(`touse')
	
	/* names */
	eret local cmd "rdhonestV1"  
    eret local dep_var  "`y'"	
    eret local run_var "`runv'" 
	if (!`Srd') eret local trtmnt_var "`trtmnt'"
	if (`NoCovar') eret local covar "`covar'"
	
	
	/* methods */ 
	eret local se_method = "`se_method'"
	eret local opt_criterion = "`opt_criterion'" 
	if (`Srd') eret local rd = "SharpRD" 
	if (!`Srd') eret local rd = "FuzzyRD" 
	
	/* estimates */
	eret scalar estimate = estimate
	if (!`Srd') eret scalar estimate_fs = fs
	eret scalar se = stdErr
	eret scalar bias = bias 
	eret scalar HLCi = HLCi 
	
	
	/* CI */
	eret scalar TCiL = lEqLimit 
	eret scalar TCiU = uEqLimit
	eret scalar OCiL = lLimit
	eret scalar OCiU = uLimit 
	eret scalar alpha = `alpha'
	if ("`opt_criterion'"=="OCI") eret scalar beta = `beta'
	
	/* other scalars */
	eret scalar bandwidth = h 
	eret scalar M = M1 
	if (!`Srd') eret scalar M_fs = M2
	eret scalar effObs = effObs
	eret scalar Leverage = Leverage
   
	 
  
end 



          *================================================
          *            Sub-routines for Display
          *=================================================



program define Display

  args y runv Srd alpha beta estimate bias stdErr lEqLimit ///
       uEqLimit lLimit uLimit h opt_criterion se_method ///
	   effObs Leverage M2 fs M1 trtmnt
  
  di _newline  
  di "{txt} Dependent variable: ", _continue
  DisplayVar `y'
  di "{txt} Running variable: ", _continue
  DisplayVar `runv' 
  if (!`Srd') {
  	di "{txt} Treatment variable: ", _continue
	DisplayVar `trtmnt' 
  }
  di "{hline 80}"
  di " {txt} Estimate  {col 15} Maximum Bias  {col 30} Std. Error {col 45} [" %2.0f `=(1-`alpha')*100',,"{txt}% Conf. {col 60} intervals]"
  di "{hline 80}" 
  di  " {res}"  `estimate'  "{col 15} {res}" `bias' "{col 30} {res}" `stdErr' "{col 45} {res}" `lEqLimit' "{col 60} {res}"`uEqLimit'
  di "{hline 80}"
  di " " as text %2.0f `=(1-`alpha')*100' "{txt}% One-sided Conf. intervals:"
  di "{txt} (" as result `lLimit' "{txt}, Inf), " , _continue
  di "{txt} (-Inf, " as result `uLimit' "{txt})"
  di "{txt} Bandwidth: " as result `h'
  di "{txt} Optimisation criterion: {res} `opt_criterion'"
  if ("`opt_criterion'"=="OCI") di "{txt} Optimisation criterion: {res} `opt_criterion'" ///
                                "{txt} with beta " as result `beta'
  di "{txt} Std Error method: {res} `se_method'"
  di "{txt} Number of effective observations: " as result  `effObs'
  di "{txt} Maximum leverage for estimated parameter: " as result  `Leverage'
  if (`Leverage'>0.1) {
  	di in r "Warning: Maximal leverage is large."
	di in r "Inference may be inaccurate.Consider using bigger bandwidth."
  }   
  if (!`Srd') {
  	di "{txt} First-stage estimate: " as result  `fs'
	di "{txt} Smootheness const. M [first-stage,reduced-form]: [" as result  `M2' "{txt}," `M1' "{txt}]"
  }
  else{
  	di "{txt} Smootheness const. M: " as result `M1'
  }
  di _newline
  end
  
  
program define DisplayVar

    args vname 

	if (c(linesize) >= 100) local abname = "`vname'"		
	else if (c(linesize) > 80) local abname = abbrev("`vname'", 20+(c(linesize)-80))
	else local abname = abbrev("`vname'", 20)
	
	di "{res} `abname' "
	if `" `:var label `vname' ' "' == "" {
		di "{res} `abname' ", _continue
		di `" {res} (`: var label `vname' ' )"' 
	}
	
end 


          *================================================
          *                       Mata
          *=================================================


*--------------------
*     contents
*--------------------
* 1  RD dataset format 
* 2  
* 3  
* 4  
* 5  
* 6
* 7
* 8  Accessary Functions 

set matastrict on

mata:


//------------------------------------
//-> 1  constructs RD dataset format
//------------------------------------

mata clear

class RDData {
  	
    /* variables */
	real matrix X,Y,sigma2
	real vector cluster, weight, m, p
    real scalar cutoff
	string Class

    /* functions */
    void setup()
  }


  
  void RDData::setup(real matrix X, real matrix Y, real scalar c,  ///
                     real matrix sigma2, real vector weight, real vector cluster, string Class) {

    real matrix s

    /* assign matrices */
	
    this.cutoff = c
	this.X = X :- c ; this.Y = Y 
	this.sigma2 = sigma2 ; this.weight = weight
	this.m = (this.X:<0) ; this.p = (this.X:>=0)
    this.cluster = cluster
	
    /* Supplied variances */
    if(rows(sigma2)!=length(X)) _error("length of supplied variances does not match data");

	
	/* Supplied weights */  
     if(length(weight)!=length(X)) _error("length of supplied weights does not match data")

	/* Supplied weights */ 
	  if(length(cluster)!=length(X)) _error("length of supplied cluster ID does not match data")
	  
	 /* sort data */
    s = sort((this.X,this.m, this.p, this.weight, this.cluster,this.Y,this.sigma2),1)
    this.X = s[.,1];  this.m = s[.,2]; this.p = s[.,3]
	this.weight= s[.,4] ;  this.cluster = s[.,5] 
	this.Y = s[.,(6..(cols(Y)+5))] ; this.sigma2 = s[.,((cols(Y)+6)..(cols(sigma2)+cols(Y)+5))];
	

	/* class */
	this.Class = Class 
    
  }
  


class RDData scalar RDDataPrep(real matrix X, real matrix Y, real scalar c , ///
                               real matrix sigma2, real vector weight, ///
							   real vector cluster, string Class) {
    
	/* constructs df */	
    class RDData scalar df  // create an instance of RDData
    
	df.setup(X,Y,c,sigma2,weight,cluster,Class)
    
	
    return(df)
  }

//------------------------------------
//-> 2  constructs Options
//------------------------------------

// Remarks:
// R.2.1 option c cutoff omiited 
// R.2.2 this class stores all the options input by the user
 
//------------------------
//-2.1 class delcarations
 
  class RDOptions {
    /* variables */
    real scalar  h, alpha, beta, j, order, t0
    string kernel, opt_criterion, se_method
    real vector m
    /* functions */
    void setup()
  }

  
  void RDOptions::setup( real vector u_m , 
    string u_kernel, real scalar u_h,
    string u_opt_criterion, real scalar u_alpha, real scalar u_beta,
    string u_se_method, real scalar u_j,
    real scalar u_order,real scalar u_t0) {

    /* assign values and defaults */
   m = u_m
   kernel = u_kernel
   h = u_h
   opt_criterion = u_opt_criterion
   alpha = u_alpha
   beta = u_beta
   se_method = u_se_method
   j = u_j
   order = u_order
   t0 = u_t0
  }
 
//---------------------
//-2.2 Prep functions
  
   class RDOptions scalar RDOptionsPrep(real vector u_m , 
    string u_kernel, real scalar u_h,
    string u_opt_criterion, real scalar u_alpha, real scalar u_beta,
    string u_se_method, real scalar u_j,
    real scalar u_order, real scalar u_t0) {
		
		
    /* constructs opt */
    class RDOptions scalar opt // create an instance of RDData
	
    opt.setup(u_m , u_kernel, u_h, u_opt_criterion, u_alpha, u_beta,
    u_se_method, u_j,u_order,u_t0)
	
    return(opt)
  }
  
 


//-----------------------------------------
//-> 3  Regression
//-----------------------------------------
  
/*
		==Remarks==
		
		R.3.1. it uses following functions defined seperately:
			(1)
			(2) 
			(3) 
			(4) 
		
		R.3.2. RDLPReg function is based on LPReg function 
*/  
  
  
  
//----------------------------------------  
//-3.1  LPReg--Local Polynomial regression

//-3.1.1 a class for the output
  
  class LPRegOutput {
		/* contains output from LPReg */
		real vector theta /* estimates */
		real matrix sigma2 /* squared residuals */
		real vector var /* sampling variance of theta */
		real vector w /* kern weights plus obs weight */
		real scalar eo /* effective observations */
		real vector p , m /* below and above indicator */
	
		real matrix wgt     /* OLS weight on Y_i that gives beta[1,] aka theta */
		                  
	}

	
//-3.1.2 implement LPReg




  class LPRegOutput scalar LPReg (real matrix X, real matrix Y, real matrix sigma2, 
                                  real matrix weight, real scalar h,
		string kernel, real scalar order, string se_method, real scalar j) {

		real scalar ord, effObs
		real matrix R, Gamma, beta, hsigma2, dsigma2, nsigma2, var
		/* weight: obs weights; */ 
		/* wt: OLS weights give theta ; */
		/* W: kern weights plus obs weights*/
		real vector wgt, w, wgt_unif
		
		real vector p , m
		
		class LPRegOutput scalar output
		
		// ?? 
		//ord = order
		
		/* W=w is a member -- kern weights plus obs weight */
		/* order = 0 ; boundary = 0 */
		w = EqKernWeight(X:/h,kernel,0,0):*weight 
		          // st_matrix("W",W)
                // stata("matlist W")
	
		p = (X:>=0) ; m = (X:<0)
		
		/* form the regression data matrix according to FWL */
		R = J(1,order+1,X) :^ (0..order)
		R = (R:*p,R)
		
	   // printf("The bandwidth for LPRegOutput is %9.0g \n", h) 
		 
		Gamma = quadcross(R,w :* R)
		
		if( h == 0| det(Gamma)==0) {
			printf("{red}Error: bandwidth too small or singular matrix.\n")
			exit(error(3352))
		}
		
		wgt = ((invsym(Gamma) * (w :* R)')[1,.])'
		/* To compute effective observations, rescale against uniform kernel*/
		/*  use weight rather than kernel + weight */
		wgt_unif = ((invsym(quadcross(R,weight:* R)) * (weight :* R)')[1,.])'

		// estimates from using OLS weights
		// squared residuals  allowing for Y being multi-variates
		beta = (invsym(Gamma) * (w :* R)' )* Y
		hsigma2 = (Y - R*beta) 	
		hsigma2=hsigma2[.,(J(1, cols(hsigma2), (1..cols(hsigma2))))]:*
		hsigma2[.,( vec(J(cols(hsigma2),1, (1..cols(hsigma2))))')]
		
		

		/* Robust variance-based formulae */
		
    /* printf("LPReg: Now checking for variance types \n") */
	
      // supplied_var: use use-supplied squared resoduals to compute EHW variance 
		if(strpos(se_method,"supplied_var") > 0) {
      /* printf("Entering LPReg supplied variance conditional \n") */
                        var = wgt:^2:*sigma2
						//printf("sigma2 col is %2.0f \n", cols(sigma2))
                }

		if(strpos(se_method,"EHW") > 0 | strpos(se_method,"ehw") > 0) {
			var = wgt:^2:*hsigma2
		}

		if(strpos(se_method,"nn") > 0 | strpos(se_method,"NN") > 0 ) {
			
			nsigma2 = sigmann(X,Y,j,weight)
			var = wgt:^2:*nsigma2
			//printf("nn col is %2.0f \n",cols(nsigma2))
		}

			
    /* catch error when no variance is provided */
    if( max((var:==.)) ) _error("No variance was computed by LPReg")
  
	

		/* return output */
		output.theta = beta[1,.]
		output.sigma2 = hsigma2
		output.var = colsum(var)
		output.w = w
		output.eo = length(X)*sum(wgt_unif:^2)/sum(wgt:^2) 
		output.wgt = wgt
		output.p = p
		output.m = m
		//printf("var is %9.4f ",output.var)
		return(output)
	}

//-----------------------------------
//-3.2 RD Local Polynomial regression	
		
	
//-3.2.1 a class for the output

  class RDLPregOutput extends LPRegOutput {
  	
			real scalar se /* standard error vector */
			real scalar fs /* 1st stage estimate (for T) */
			real scalar estimate /* treatment effect */
			real vector kw /* kernel weight */
	}

//-3.2.2 impplement LPRegOutput

class RDLPregOutput scalar RDLPreg(class RDData scalar df, real scalar h,
	|string kernel, real scalar order, string se_method,
		real scalar no_warning, real scalar j) {

    /* grep: need to debug for sclass H order 2 */

		/* check for errors */
		if (h <= 0) _error("Non-positive bandwidth h")

	 	/* set defaults */
		if (args() <= 3 & kernel == "") kernel = "tri"
		if (args() <= 4 & order == .) order = 1
	    if (args() <= 5 & se_method == "") se_method = "nn"
		if (args() <= 6 & no_warning == .) no_warning = 0
		if (args() <= 7 & j == .) j = 3

		/* variable declarations */
		real scalar plugin /* to be implemented later */
		real matrix Y, X, w, sigma2, Xm, Xp, weight
		
		class RDLPregOutput scalar output
		class LPRegOutput scalar r1

		
		/* set kernel weights */
		w = (h <= 0 ? (0 :* df.X) : (EqKernWeight(df.X:/h,kernel,0,0)):*df.weight)
		
		/* variance calculations are faster if we only keep data with positive
	 	kernel weights */

		X = select(df.X,w :> 0) ;	Y = select(df.Y,w :> 0)
        sigma2 = select(df.sigma2,w :> 0) ; weight = select(df.weight,w :> 0) 

		Xm = select(X,X :< 0) ; Xp = select(X,X :>= 0)
		if ((length(Xm) < 3*order | length(Xp) < 3*order) & !no_warning) {
			printf("{red}Warning: Too few observations to compute RD estimates.\n")
			printf("{red}Only %g control and %g treated units with positive weights.\n",
			length(Xm),length(Xp))
		}

		if ((length(Xm) <= order | length(Xp) <= order) & !no_warning) {
			printf("{red}Warning: Too few distinct values to compute RD estimates.\n")
			printf("{red}Only %g unique control and %g unique treated values, \n",
			length(uniqrows(Xm)),length(uniqrows(Xp)))
			printf("{red}with positive weights, for the running variable. \n")
		}

    
       /* regression */
	   
	  /* diagnostics */
      /* printf("{red}entering RDLPreg supplied variance conditional \n")
      printf("truncated sigma2 has first entry %9.0g \n", sigma2[1]) */
	  
	  r1 = LPReg(X , Y, sigma2, weight, h, kernel, order, se_method, j)
	  
      // printf("r1.var has value %9.0g \n", r1.var[1]) 
       
	  
	  /* output */

		output.estimate = r1.theta[1]
		output.w = r1.w
		output.kw = w
		output.se = sqrt(r1.var[1])
		output.sigma2 = r1.sigma2
		plugin = .
		output.eo = r1.eo
		output.p = r1.p
		output.m = r1.m
		output.var = r1.var
		output.wgt = r1.wgt
		
		
	    if (df.Class == "frd") {	

		output.fs = r1.theta[2]
		output.estimate = r1.theta[1]/r1.theta[2]
		output.se = sqrt(quadcross(
		           (1,-output.estimate,-output.estimate,output.estimate^2)',r1.var')
				   /output.fs^2)		
	    }

		return(output)
		
	
		
		
	}
	

//-----------------------------------------
//-> 5  Rule of thumb choosing M
//-----------------------------------------
  

//-Remarks:
//-5.(1) Use global quadratic regression to estimate a bound on the second derivative 
//-5.(2) Use only when M is not supplied by users. 

real vector MROT_fit(real matrix X, real matrix Y, | real scalar IsIP) {
		
		if (args()<3) IsIP = 0 

		real vector M, col1, col2, col3, col4, col5, m3coeff
		real matrix  m3mat, Xm ,Xp , Ym, Yp 
		//class RDData scalar dfm ,dfp
		
		
		if (IsIP == 1) {
			
			/* Step 1: Estimate global polynomial regression */
			col1 = X :^ 0; col2 = X; col3 = X :^ 2; col4 = X :^ 3; col5 = X :^ 4;
			m3mat = col1, col2, col3, col4, col5
			m3coeff = invsym(m3mat'm3mat)*m3mat'(Y)
		    
			/* Step 2: maximum occurs either at endpoints, or else at the extremum */ 
			/* the extremum = -r1[4]/(4*r1[5]), if the extremum is in the support */
			M = MROT_aux_endpoint(m3coeff, X)
			
		} 
		else if (cols(Y) == 1) {
			
			/* on either side */
			Xp = select(X,X:>=0); Yp = select(Y,X:>=0) 
			Xm = select(X,X:<0); Ym = select(Y,X:<0) 
			
			M = max((MROT_fit(Xm,Ym,1),MROT_fit(Xp,Yp,1)))
			
		}
		else {
			
			/* for FRD compute M for two reduced form regressions */
			M = (0,0)
	
			/* Reduced form regression */
			Xp = select(X,X:>=0); Yp = select(Y[.,1],X:>=0) 
			Xm = select(X,X:<0); Ym = select(Y[.,1],X:<0) 
			M[1] = max((MROT_fit(Xm,Ym,1),MROT_fit(Xp,Yp,1)))
			/* 1st stage regression */
			Xp = select(X,X:>=0); Yp = select(Y[.,2],X:>=0) 
			Xm = select(X,X:<0); Ym = select(Y[.,2],X:<0) 
			M[2] = max((MROT_fit(Xm,Ym,1),MROT_fit(Xp,Yp,1)))
			
		}
		

		return (M) 
	}
	

       /* Auxilary function for MROT_fit step 2 */
  real scalar MROT_aux_endpoint(real vector coef, real vector X){
		
		real scalar Fmin, Fmax, Fe, M 
		
		if (abs(coef[5])<=1e-10) {
			Fe = 8e307
		}
		else {
			Fe = -coef[4]/(4*coef[5])
		}
		
		Fmin = abs( 2*coef[3] + 6*min(X)*coef[4] + 12*min(X)^2*coef[5])
		Fmax = abs( 2*coef[3] + 6*max(X)*coef[4] + 12*max(X)^2*coef[5])
		
		if ( min(X)<Fe & max(X)>Fe ) {
			M = max((abs( 2*coef[3] + 6*Fe*coef[4] + 12*Fe^2*coef[5]), max((Fmin,Fmax))))
		}
		else {
			M = max((Fmin,Fmax))
		}
		
		return (M) 
  }


//-----------------------------------------
//-> 5  Prelim Variance for Optimal h
//-----------------------------------------



//-------------------------------------------  
//-5.1  IKBW_fit--Local Polynomial regression


  class RDData scalar IKBW_fit(class RDData scalar df, real matrix kernC,| string kernel,
		real scalar verbose) {
		/* Calculate bandwidth for sharp RD based on local linear regression
		using method by Imbens and Kalyanaraman (ReStud 2012) */
		/* note: no option for order != 1 */

		if(args() < 3)  kernel = "triangular"
		if(args() < 4)  verbose = 0

		real scalar order, Nm, Np, N, kernInd, cons, varm, varp, h2m, h2p, m2m, m2p
		real scalar h1, f0, m3, rm, rp
		real vector col1, col2, col3, col4, col5, m3coeff, tempY, Ym, Yp, Xm, Xp
		real matrix s,X, m3mat, tempX, m2mat

		X = df.X
		order = 1
		Nm = sum(df.m)
		Np = sum(df.p)
		N = Nm + Np

		/* kernels must be specified as strings */

		/* translate kernel strings into kernel indices */
		/* uni : 0  tri : 1  epa: 2 */
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
		
		/* run prelim var which asssumes homoskedasticity */
		df = RDPrelimVar(df, kernC, "Silverman")
		h1 = 1.84 * sqrt(variance(X))/(N^(1/5))
		f0 = sum(abs(X) :<= h1)/(2 * N * h1)	
		varm = select(df.sigma2,df.m:==1)[1]; varp = select(df.sigma2,df.p:==1)[1]
	
		/* run a linear regression of y on 1{x >= 0}, x, x^2 and x^3, then extract
		coefficients of x^3 */
		col1 = X :^ 0; col2 = (X :>= 0); col3 = X; col4 = X :^ 2;  col5 = X :^ 3
		m3mat = col1, col2, col3, col4, col5
		m3coeff = invsym(m3mat'm3mat)*m3mat'(df.Y)
		m3 = 6 * m3coeff[5]

		/* left and right bandwidth */
		h2m = 7200^(1/7) * (varm/(f0 * m3^2))^(1/7) * Nm^(-1/7)
		h2p = 7200^(1/7) * (varp/(f0 * m3^2))^(1/7) * Np^(-1/7)

		/* run a linear regression of y_m on x_m and x_m^2, amongst obs that lie
		within the bandwidth */
		Yp = select(df.Y,df.p:==1) ; Ym = select(df.Y,df.m:==1) 
		Xp = select(df.X,df.p:==1) ; Xm = select(df.X,df.m:==1)
		
		tempY = select(Ym,Xm :>= -h2m)
		tempX = select(Xm,Xm :>= -h2m)
		col1 = tempX :^ 0; col2 = tempX; col3 = tempX :^ 2
		m2mat = col1, col2, col3
		m2m = 2 * (invsym(m2mat'm2mat) * m2mat'tempY)[3]

		/* same procedure: y_p on x_p and x_p^2, within bandwidth */
		tempY = select(Yp,Xp :<= h2p)
		tempX = select(Xp,Xp :<= h2p)
		col1 = tempX :^ 0; col2 = tempX; col3 = tempX :^ 2
		m2mat = col1, col2, col3
		m2p = 2 * (invsym(m2mat'm2mat) * m2mat'tempY)[3]

		/* Calculate regularization terms */
		rm = 2160 * varm/(sum(Xm :>= -h2m) * h2m^4)
		rp = 2160 * varp/(sum(Xp :<= h2p) * h2p^4)

		/* output parameters if verbose */
		if (verbose) {
			printf("\n h1: %9.0g \n N_{-}: %9.0g \n  N_{+}: %9.0g", h1, Nm, Np)
			printf("\n f0: %9.0g \n sigma^2_{+}(0): %9.0g^2 \n sigma^2_{-}(0): %9.0g^2",
				f0, sqrt(varp), sqrt(varm))
			printf("\n m3: %9.0g \n h_{2,+}: %9.0g \n h_{2,-}: %9.0g",m3,h2p,h2m)
			printf("\n m^{(2)}_{+}: %9.0g \n m^{(2)}_{-}: %9.0g",m2p,m2m)
			printf("\n r_{+}: %9.0g \n r_{-}: %9.0g \n\n",rp,rm)
		}
		
		// ? a real scalar, then why it's declared to be of RDData class 
		return(cons * ((varp + varm)/(f0 * N * ((m2p - m2m)^2 + rm + rp)))^(1/5))
	}
	
	
//----------------------------------------  
//-5.2  RDPrelimVar--computes prelim var

//- Remarks: 
//- R.5.2.(1) assuming Homoskedasticity
//- R.5.2.(2) difference with Reg functions:
//-           Here: sigma2 stands for estimate for variance under Homoskedasticity
//-           Reg: squared residuals
//- R.5.2.(3) this function updates the matrix sigma2 in the object df. 

  class RDData scalar RDPrelimVar(class RDData scalar df, real matrix kernC,| string se_initial) {
		
        
	
		real matrix X, Xp, Xm
		class RDLPregOutput scalar r1
		real scalar h1, hmin, lm, lp, varm, varp 
		class RDData scalar drf 
		

		/* set defaults: EHW and IK bandwidth*/
		if(args()==2) se_initial = "IKEHW"

		/* concatenate Xm and Xp and compute "rule of thumb" bandwidth */
		X = df.X
		
		Xp = select(df.X,df.X:>=0) ; Xm = select(df.X,df.X:<0)
		
		hmin = max((uniqrows(Xp)[3], uniqrows(Xm)[3],
                    sort(Xp,1)[4], sort(abs(Xm),1)[4]))
		h1 = max((1.84*sqrt(variance(X))/(sum(length(X)))^(1/5), hmin))
		
		/* FRD */
		drf = df
		if (df.Class == "frd") {
			
			drf.Y=df.Y[.,1]
			drf.Class = "srd"
			
		}
					 
        /*  */

		if(strpos(se_initial,"Silverman") > 0) {
					
			if (cols(df.Y) == 1) {
				
			r1 = RDLPreg(df,h1,"uni",0,"EHW")
			
			/* Variance adjustment on either side */
			lp = sum(r1.p)
			lm = sum(r1.m)
			
			varp = sum(r1.sigma2:*r1.p)*1/(lp-1)
			varm = sum(r1.sigma2:*r1.m)*1/(lm-1)
			
			df.sigma2 = (df.X:<0):*varm + (df.X:>=0):*varp 
				
			} 
			else {	
				_error("This method for preliminary variance estimation not supported")
			}
			
		}
		else if(strpos(se_initial,"IKEHW") > 0) {
			
			h1 = IKBW_fit(drf, kernC)	
			r1 = RDLPreg(df,max((h1,hmin)),"tri",1,"EHW")
			
			lp = sum(r1.p)
			lm = sum(r1.m)
			
			varp = colsum(r1.sigma2:*r1.p)/lp
			varm = colsum(r1.sigma2:*r1.m)/lm
			
			
			df.sigma2 = (df.X:<0):*J(rows(df.X),1,varm) + (df.X:>=0):*J(rows(df.X),1,varp) 
			
		}
		else {
			_error("Unknown method for estimating initial variance.")
		}

		return(df)
	}
	 	

//------------------------------------
//-> 6  RDHonest_fit 
//------------------------------------

/*
		==Remarks==
		
		R.6.(1) it uses following functions defined seperately:
			 (a) RDPrelimVar: compute prelim variance 
			 (b) RDOptBW_fit: conpute optimal bandwidth
			 (c) gridSearchOptim: optimasation function for the case of Holder and order 2
			 (d) RDLPreg: RD-local polynomial regression
		R.6.(2)
*/

//-6.1 Declare a class for results
         
struct RDResults {
      real scalar estimate, fs, leverage
      real scalar lff /* not sure what this is */
      real scalar bias, sd, lower, upper, hl, eo, h, naive
  }
  

//-6.2 implement NPRDHonest_fit

//-Remarks
//-6.2.(1) opt.t0: Initial estimate of the treatment effect for calculating the
//-            optimal bandwidth. Only relevant for Fuzzy RD.
//-6.2.(2) T0bias=1: When evaluating the maximum bias of the estimate, use the T0.
//-              =0: use the estimate itself.
//-        Only relevant for Fuzzy RD.   

struct RDResults scalar NPRDHonest_fit (class RDData scalar df,
		class RDOptions scalar opt, real matrix kernC, 
		real scalar T0bias) {
			
//==============================================================================
//-6.2.0 declare variables	

		class RDLPregOutput scalar r1
		real scalar h, z
		real vector wp, wm, M //, hVec
		real matrix XX, XXp, XXm 
		struct RDResults scalar results

		h = opt.h
				
//-6.2.1 deal with missing bandwidths
	
	if(opt.h <= 0) {
        // hacky fix: previously opt.h == .
      printf("Bandwidth (h) missing or invalid. Running RDOptBW_fit \n")
      h = RDOptBW_fit(df, opt, kernC)
		}	
	//printf("Bandwidth (h) computation is completed \n")
			
//-6.2.2 compute the bias 

	/* Local Polinomial regression */
		r1 = RDLPreg(df,h, opt.kernel, opt.order, opt.se_method, 1, opt.j)
		// Suppress warnings about too few observations 
		wp = select(r1.wgt,r1.p:==1)
		wm = select(r1.wgt,r1.m:==1)
		XX = select(df.X,r1.kw:>0)
		XXp = select(XX,r1.p:==1)
		XXm = select(XX,r1.m:==1)
		M = opt.m
		
	/* multiply bias and sd by r1.fs to make if free of first stage */   
		if (df.Class=="frd" & T0bias) {
					 
				 r1.se = r1.se*abs(r1.fs)
                 M = ((M[1]+M[2]*abs(opt.t0)), M)
				 
			}
			else if (df.Class=="frd" & !T0bias) {
				
				  M = ((M[1]+M[2]*abs(r1.estimate)) / abs(r1.fs), M)
				
			}
				

	/* bias and one-sided CI ending porint */
	
		/* If bandwidths are too small, set bias to machine maximum */
		if (sum(wp :>= 0) == 0 | sum(wm :>= 0) == 0) {
			results.bias = results.sd = results.upper = results.hl = sqrt(1.0X+4e/10)
			results.lower = -results.upper
			
		}
		else {
			
			results.bias = M[1]/2 * abs(sum(wm:* XXm:^2)-sum(wp:* XXp:^2))
			results.sd = r1.se
			results.hl = (CVb(results.bias/results.sd, opt.alpha)) * results.sd
			
			
            results.lower = r1.estimate - results.bias - invnormal(1-opt.alpha) * results.sd
			results.upper = r1.estimate + results.bias + invnormal(1-opt.alpha) * results.sd
		}
		
	/* compute naive CIs */			
		
		z = invnormal(1-opt.alpha/2)
		results.naive = normal(z-results.bias/results.sd) - normal(-z-results.bias/results.sd)
	
		results.leverage = max(r1.wgt:^2)/sum(r1.wgt:^2)
		results.h = h
		results.eo = r1.eo
		results.estimate = r1.estimate
		results.fs = r1.fs 
		//printf("estimates is: %4.2f \n",r1.estimate)
		
		return(results)
		}

//==============================================================================

//-6.3 implement RDHonest_fit

struct RDResults scalar RDHonest_fit (class RDData scalar df,
		class RDOptions scalar opt, real matrix kernC ) {
   
	/* declare variables */	
    struct RDResults scalar results



    /* initial se estimate */
		
    if ( max((df.sigma2:==.))
      & (strpos(opt.se_method,"supplied_var") > 0 | opt.h == 0)) {
      printf("Generating initial variance estimates via RDPrelimVar \n")
	  // when no Prevar and no h provided 
      df = RDPrelimVar(df,kernC)
	  // df is updated with PreVar computed 
	  // printf("RDPrelimVar computation completed \n")
	  //printf("cols sigma2 is %2.0f \n",cols(df.sigma2))
    }
	
	/* smoothing class constant */
	if(max(opt.m) < 0) {
	  printf("Smoothing class constant M missing or invalid. Running MROT_fit \n")
      opt.m = MROT_fit(df.X,df.Y)	
	}    
	
	/* optimal bandwidth */
   results =  NPRDHonest_fit(df, opt,kernC, 0)



    /*** numerical scalars ***/ 
      st_numscalar("estimate",results.estimate)
      st_numscalar("bias",results.bias)
      st_numscalar("stdErr",results.sd)
      st_numscalar("h",results.h)
      st_numscalar("effObs",results.eo)
	  st_numscalar("Leverage",results.leverage)
      // one-sided confidence interval limits
      st_numscalar("lLimit",results.lower)
      st_numscalar("uLimit",results.upper)
      // symmetric two-sided confidence interval limits
	  st_numscalar("HLCi",results.hl)
      st_numscalar("lEqLimit",results.estimate - results.hl)
      st_numscalar("uEqLimit",results.estimate + results.hl)
	  // Smoothness Constant 
	  st_numscalar("M1",opt.m[1])
      // *** strings from specified options ***
      //st_strscalar("se_method",opt.se_method)
    
	if (df.Class=="frd"){		
		st_numscalar("fs",results.fs)
		st_numscalar("M1",opt.m[1])
		st_numscalar("M2",opt.m[2])		
	}
	

		return(results)
	}
	
//------------------------
//-> 7  Compute Optimal h
//------------------------

//- 7.1 Compute Optimal h

  real vector RDOptBW_fit (class RDData scalar df, class RDOptions scalar opt,
    real matrix kernC) {

    transmorphic S /* to optimize function */
    real scalar h0, hmin, hmax, h // ROT guess, minimal and maximal h
    real vector h_opt, supp
    real matrix X

    /* First check if sigma2 is supplied */
    //if (! max((df.sigma2:<0))) df = RDPrelimVar(df,kernC)

    /* compute Silverman rule of thumb bandwidth */
    X = df.X
    h0 = 1.06 * sqrt(variance(X)) * length(X)^(-1/5)

    /* initialize bandwidth optimization */
    S = optimize_init()
    optimize_init_argument(S,1,df)
    optimize_init_argument(S,2,opt)
    optimize_init_argument(S,3,kernC)
    optimize_init_which(S, "min")

    /* support of h used in optimisation */
      hmin = max((uniqrows(select(df.X,df.p:==1))[opt.order+1],
	              uniqrows(abs(select(df.X,df.m:==1)))[opt.order+1]))
	
      hmax = max(abs(df.X))
	  
	  
     /* Optimize piecewise-constant function using modified golden section search.
      Careful: criterion may not be unimodal (but appears to be so, even for
      triangular kernel) */
      if (strpos(opt.kernel,"uni") > 0) {
         /* printf("Entering golden section search \n") */
         supp = sort(uniqrows(df.X),1)
         /* changed select(supp,supp :> hmin) */
         h = gss(&RDOptBWGss(), select(supp,supp :>= hmin), df, opt, kernC)
         /* printf("Golden section bandwidth obtained is %9.0g \n",h) */
      }
      else {
        // set up parameters of optimization problem
        optimize_init_argument(S,4,hmin)
        optimize_init_argument(S,5,hmax)
        optimize_init_evaluator(S,&RDOptBWEval())
        optimize_init_params(S,h0)
        optimize_init_technique(S,"bfgs")
        h_opt = optimize(S)
        h = h_opt
      }
    
    /* A REMPLIR */
    return(h)
  }
  
  
 
//-7.2  Objective functions in optimisation 

//-7.2.1 Objective functions
  
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

    optLocal.setup(opt.m, opt.kernel, abs(b), opt.opt_criterion,
      opt.alpha, opt.beta, "supplied_var", opt.j, opt.order, opt.t0)

    r = NPRDHonest_fit(df,optLocal,kernC,1)
	
    if (strpos(optLocal.opt_criterion,"OCI") > 0) val = 2*r.bias + r.sd*(
      invnormal(1-optLocal.alpha) + invnormal(optLocal.beta))
    else if (strpos(optLocal.opt_criterion,"MSE") > 0) val = (r.bias)^2 + (r.sd)^2
    else if (strpos(optLocal.opt_criterion,"FLCI") > 0) val = r.hl
    else {
      printf("optimality criterion %s not implemented yet \n", optLocal.opt_criterion)
      _error("exiting RDOptBWEval()")
         }
  }

//-7.2.1 objective function for golden section search 

  real scalar RDOptBWGss (real scalar h, class RDData scalar df,
    class RDOptions scalar opt, real matrix kernC) {

    struct RDResults scalar r
    class RDOptions scalar optLocal
    real scalar val

    /* options should be that uses PrelimVar stored in sigma2 */
    optLocal.setup(opt.m, opt.kernel, h, opt.opt_criterion,
      opt.alpha, opt.beta, "supplied_var", opt.j, opt.order,opt.t0)

    r = NPRDHonest_fit(df,optLocal,kernC,1)

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
  
 

  
	
//-----------------------------------------
//-> 8  Accessary Fcts
//-----------------------------------------

	
	/* generate KNN standard errors--1*/
	
	real matrix sigmann(real matrix X, real matrix Y, real scalar j, real vector weight) {
		
		real matrix Xm ,Xp, Ym, Yp, sigma2, weightm, weightp
		
		Xm = select(X,X:<0) ; Xp = select(X,X:>=0) 
		Ym = select(Y,X:<0) ; Yp = select(Y,X:>=0)
		weightm = select(weight,X:<0) ; weightp = select(weight,X:>=0)  
		
		/* For RD, compute variance separately on either side of cutoff */
		
		sigma2 = sigmaNN(Xm,Ym,j,weightm)\sigmaNN(Xp,Yp,j,weightp)
			
		
		return(sigma2)	
		
	}
	
	
	/* generate KNN standard errors--2 */

  real matrix sigmaNN(real matrix X, real matrix Y, real scalar j, real vector weight) {
    /* assume X is sorted */
    real scalar n,k,d,Jk
    real matrix sigma2,s,ind,YTemp

    n = length(X)
    sigma2 = J(n,cols(Y)^2,.)
	
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
      Jk = sum(weight:*ind)
      /* compute mean of Y with positive indices */
      YTemp = (sort((ind,weight:*Y),-1))[1..Jk,2..(cols(Y)+1)]:/Jk
	 
      sigma2[k,.] = vec(quadcross(Jk/(Jk+weight[k]) :* (Y[k,.] - colsum(YTemp)),
	                    (Y[k,.] - colsum(YTemp))))'
	}					
		return (sigma2)				
    }
  /* enhance rdrobust_kweight() to EqKern() standards by
		including order and boundary parameters */

  real vector EqKernWeight(real matrix Xh, string kernel,
		 real scalar order, real scalar boundary)
	{
		/* input Xh as X/h */
		real vector u,w, su 
		
		/* support */
		su = (Xh:<=1):*(Xh:>=(-1+boundary))
		u = Xh
		
		if (boundary == 0 & order == 0) {
			if (kernel=="epanechnikov"|kernel=="epa") {
				w = (0.75:*(1:-u:^2):*su)
			}
			else if(kernel=="uniform"|kernel=="uni") {
				w = (0.5:*su)
			}
			else if(kernel=="triangular"|kernel=="tri"){
				w = (1:-abs(u)):*su
			}
			else {
				printf("{red}Error: Unsupported kernel type\n")
				exit(error(198))
			}
		}
		else if (boundary == 0 & order == 1) { 
			if (kernel=="epanechnikov"|kernel=="epa") {
				w = (0.75:*(1:-u:^2):*su)
			}
			else if(kernel=="uniform"|kernel=="uni") {
				w = (0.5:*su)
			}
			else if(kernel=="triangular"|kernel=="tri"){
				w = (1:-abs(u)):*su
			}
			else {
				printf("{red}Error: Unsupported kernel type\n")
				exit(error(198))
			}
		}
		else if (boundary == 0 & order == 2) { 
			if (kernel=="epanechnikov"|kernel=="epa") {
				w = (15/32) :* (3:-7:*u:^2) :* (1:- u:^2):*su
			}
			else if(kernel=="uniform"|kernel=="uni") {
				w = ((9 :- 15 :* u:^2) / 8) :* su
			}
			else if(kernel=="triangular"|kernel=="tri"){
				w = 6/7 :* (2:-5:*u:^2) :* (1:-abs(u)) :* su
			}
			else {
				printf("{red}Error: Unsupported kernel type\n")
				exit(error(198))
			}
		}
		else if (boundary == 1 & order == 0) {
			if (kernel=="epanechnikov"|kernel=="epa") {
				w = 1.5 :* (1 :- u:^2):*su
			}
			else if(kernel=="uniform"|kernel=="uni") {
				w = su
			}
			else if(kernel=="triangular"|kernel=="tri"){
				w = 2 * (1 :- u) :* su
			}
			else {
				printf("{red}Error: Unsupported kernel type \n")
				exit(error(198))
			}
		}
		else if (boundary == 1 & order == 1) {
			if (kernel=="epanechnikov"|kernel=="epa") {
				w = 6/19 * (16 :- 30*u) :* (1 :- u:^2) :* su
			}
			else if(kernel=="uniform"|kernel=="uni") {
				w = (4 :- 6 :* u) :* su
			}
			else if(kernel=="triangular"|kernel=="tri"){
				w = 6 * (1 :- 2 :* u) :* (1 :- u) :* su
			}
			else {
				printf("{red}Error: Unsupported kernel type\n")
				exit(error(198))
			}
		}
		else if (boundary == 1 & order == 2) {
			if (kernel=="epanechnikov"|kernel=="epa") {
				w = 1/8 * (85 :- 400*u :+ 385*u:^2) :* (1 :- u:^2) :* su
			}
			else if(kernel=="uniform"|kernel=="uni") {
				w = (9 :- 36*u :+ 30*u:^2) :* su
			}
			else if(kernel=="triangular"|kernel=="tri"){
				w = 12 * (1 :- 5*u :+ 5*u:^2) :* (1 :- u) :* su
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
 

  /* modified golden section optimizer to use with mata objective functions */ 
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
  

 
 /* implement CVb with noncentral chi-squared distribution */

real scalar CVb(real scalar B,| real scalar alpha) {
    if (args()==1) alpha = 0.05
    
    // declare variables
    real scalar root
    
    // check for null bias/alpha out of range
    if(B == .) return(.)
    if(B < 0 | alpha <= 0 | alpha >= 1) _error("B negative or alpha out of range")

	if (B<10) {
		root = sqrt(invnchi2(1,B^2,1-alpha))
	}
    else{
		root = B + invnormal(1-alpha)
	}
    return(root)
    }

	/* Kernel Constant */ 
	
// load kernel constants into memory. 
// Note: uniform == 0, triangular == 1, epanechnikov == 2.
kernCmu = (

   //kernel		order		boundary		mu0		mu1		mu2		mu3		mu4

	0	,	0	,	1	,	1	,	0.5	,	0.333333343	,	0.25	,	0.200000003	\
	1	,	0	,	1	,	1	,	0.333333343	,	0.166666672	,	0.100000001	,	0.06666667	\
	2	,	0	,	1	,	1	,	0.375	,	0.200000003	,	0.125	,	0.085714288	\
	0	,	1	,	1	,	1	,	0	,	-0.166666672	,	-0.200000003	,	-0.200000003	\
	1	,	1	,	1	,	1	,	0	,	-0.100000001	,	-0.100000001	,	-0.085714288	\
	2	,	1	,	1	,	1	,	0	,	-0.115789473	,	-0.120300755	,	-0.106015034	\
	0	,	2	,	1	,	1	,	0	,	0	,	0.050000001	,	0.085714288	\
	1	,	2	,	1	,	1	,	0	,	0	,	0.028571429	,	0.042857144	\
	2	,	2	,	1	,	1	,	0	,	0	,	0.033482142	,	0.051587302	\
	0	,	0	,	0	,	1	,	0	,	0.333333343	,	0	,	0.200000003	\
	1	,	0	,	0	,	1	,	0	,	0.166666672	,	0	,	0.06666667	\
	2	,	0	,	0	,	1	,	0	,	0.200000003	,	0	,	0.085714288	\
	0	,	1	,	0	,	1	,	0	,	0.333333343	,	0	,	0.200000003	\
	1	,	1	,	0	,	1	,	0	,	0.166666672	,	0	,	0.06666667	\
	2	,	1	,	0	,	1	,	0	,	0.200000003	,	0	,	0.085714288	\
	0	,	2	,	0	,	1	,	0	,	0	,	0	,	-0.085714288	\
	1	,	2	,	0	,	1	,	0	,	0	,	0	,	-0.038775511	\
	2	,	2	,	0	,	1	,	0	,	0	,	0	,	-0.047619049	)

kernCnu = (

	//kernel		order		boundary		nu0		nu1		nu2		nu3		nu4	
	
	0	,	0	,	1	,	1	,	0.5	,	0.333333343	,	0.25	,	0.200000003	\
	1	,	0	,	1	,	1.333333373	,	0.333333343	,	0.13333334	,	0.06666667	,	0.03809524	\
	2	,	0	,	1	,	1.200000048	,	0.375	,	0.171428576	,	0.09375	,	0.057142857	\
	0	,	1	,	1	,	4	,	1	,	0.533333361	,	0.400000006	,	0.342857152	\
	1	,	1	,	1	,	4.800000191	,	0.600000024	,	0.171428576	,	0.085714288	,	0.057142857	\
	2	,	1	,	1	,	4.497982025	,	0.700435281	,	0.235536203	,	0.128215268	,	0.088872902	\
	0	,	2	,	1	,	9	,	1.5	,	0.771428585	,	0.578571439	,	0.485714287	\
	1	,	2	,	1	,	10.28571415	,	0.857142866	,	0.22857143	,	0.114285715	,	0.07272727	\
	2	,	2	,	1	,	9.816468239	,	1.018105149	,	0.322420627	,	0.175161213	,	0.116515428	\
	0	,	0	,	0	,	0.5	,	0	,	0.166666672	,	0	,	0.100000001	\
	1	,	0	,	0	,	0.666666687	,	0	,	0.06666667	,	0	,	0.01904762	\
	2	,	0	,	0	,	0.600000024	,	0	,	0.085714288	,	0	,	0.028571429	\
	0	,	1	,	0	,	0.5	,	0	,	0.166666672	,	0	,	0.100000001	\
	1	,	1	,	0	,	0.666666687	,	0	,	0.06666667	,	0	,	0.01904762	\
	2	,	1	,	0	,	0.600000024	,	0	,	0.085714288	,	0	,	0.028571429	\
	0	,	2	,	0	,	1.125	,	0	,	0.160714284	,	0	,	0.08214286	\
	1	,	2	,	0	,	1.329446077	,	0	,	0.06180758	,	0	,	0.013570104	\
	2	,	2	,	0	,	1.25	,	0	,	0.08116883	,	0	,	0.021228772	)

kernCpi = (

	//kernel		order		boundary		pi0		pi1		pi2		pi3		pi4		
	
	0	,	0	,	1	,	1	,	0.5	,	0.333333343	,	0.25	,	0.200000003	\
	1	,	0	,	1	,	1	,	0.333333343	,	0.166666672	,	0.100000001	,	0.06666667	\
	2	,	0	,	1	,	1	,	0.375	,	0.200000003	,	0.125	,	0.085714288	\
	0	,	1	,	1	,	1.666666627	,	0.592592597	,	0.364197522	,	0.279012352	,	0.235116601	\
	1	,	1	,	1	,	1.5	,	0.375	,	0.1875	,	0.125	,	0.09375	\
	2	,	1	,	1	,	1.566986322	,	0.438184172	,	0.2290048	,	0.155643702	,	0.118335322	\
	0	,	2	,	1	,	2.175755024	,	0.705453038	,	0.43738088	,	0.329359412	,	0.268931001	\
	1	,	2	,	1	,	1.89442718	,	0.429325044	,	0.214662522	,	0.139991507	,	0.102655999	\
	2	,	2	,	1	,	2.005158424	,	0.507928908	,	0.266179353	,	0.177708849	,	0.132114798	\
	0	,	0	,	0	,	1	,	0.5	,	0.333333343	,	0.25	,	0.200000003	\
	1	,	0	,	0	,	1	,	0.333333343	,	0.166666672	,	0.100000001	,	0.06666667	\
	2	,	0	,	0	,	1	,	0.375	,	0.200000003	,	0.125	,	0.085714288	\
	0	,	1	,	0	,	1	,	0.5	,	0.333333343	,	0.25	,	0.200000003	\
	1	,	1	,	0	,	1	,	0.333333343	,	0.166666672	,	0.100000001	,	0.06666667	\
	2	,	1	,	0	,	1	,	0.375	,	0.200000003	,	0.125	,	0.085714288	\
	0	,	2	,	0	,	1.323789954	,	0.487500012	,	0.278854787	,	0.197500005	,	0.157419801	\
	1	,	2	,	0	,	1.205510974	,	0.311559111	,	0.139869452	,	0.084430546	,	0.060140885	\
	2	,	2	,	0	,	1.244526863	,	0.360331625	,	0.171775013	,	0.106710099	,	0.077066191	)


kernCpmse = (	
	
	//kernel		order		boundary		pmse	
	
	0	,	0	,	1	,	1.25992105	\
	1	,	0	,	1	,	1.817120593	\	
	2	,	0	,	1	,	1.621920532	\
	0	,	1	,	1	,	2.701920077	\
	1	,	1	,	1	,	3.437543855	\
	2	,	1	,	1	,	3.199896319	\
	0	,	2	,	1	,	4.161095457	\
	1	,	2	,	1	,	4.976587966	\
	2	,	2	,	1	,	4.724482988	\
	0	,	0	,	0	,	1E+100	\
	1	,	0	,	0	,	1E+100	\
	2	,	0	,	0	,	1E+100	\
	0	,	1	,	0	,	1.350960039	\
	1	,	1	,	0	,	1.888175023	\
	2	,	1	,	0	,	1.718771928	\
	0	,	2	,	0	,	1E+100	\
	1	,	2	,	0	,	1E+100	\	
	2	,	2	,	0	,	1E+100	)

kernC = (kernCmu,  kernCnu[.,4..8] ,  kernCpi[.,4..8] ,  kernCpmse[.,4] )	


  end