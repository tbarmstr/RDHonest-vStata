*! version 1  
*! CYG




program define cyg_rdhonest, eclass

version 14.2

    	syntax [anything(name=0)] [if] [in] [aw fw pw iw/] ///
		   [, m(real 0) ///
		   c(real 0) /// 
		   h(real 0) ///
		   kernel(string) ///
		   opt_criterion(string) ///
		   alpha(real 0.05) beta(real 0.8) ///
		   se_method(string) ///
           j(integer 3) ///
		   order(integer 1)]

 
               /*  find y x treatment and covariates  */
	
	//       DepVar runingVar TreatmentVar CoVar 
	// RD:    y        runv                covar 
    // FRD:   y        runv      trtmnt    covar 
	
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
	

		              /* mark the data sample*/
						   
	      	 
  marksample touse
  preserve
  qui keep if `touse'
  
  // set more off for now
  set more off

  
                        /* set defaults for string options*/
 
  if ("`kernel'" == "") local kernel "tri"
  if ("`opt_criterion'" == "") local opt_criterion "MSE"
  if ("`se_method'" == "") local se_method "nn"

   
                        
	                       /* Mata */
 set matastrict on
  
  
  mata{ 
    Y = st_data(.,("`y'"),0);
	X = st_data(.,("`runv'"),0)	
	if (`Srd' == 0) T = st_data(.,("`trtmnt'"),0);
	if (`NoCovar' == 0) Z = st_data(.,("`covar'"),0);
  
    /* initialize data frame and map Y and X in */
    df = RDDataPrep(X,Y,`c')
	if (`Srd' == 0) dfFS = RDDataPrep(X,T,`c');

    /* printf("`se_method', `j', `sclass', `order', `se_initial'") */

    /* initialize options */
    m = `m'
    kernel = "`kernel'"
    opt_criterion = "`opt_criterion'"
    alpha = `alpha'
    beta = `beta'
    se_method = "`se_method'"
    j = `j'
    order =  `order'
	h = `h'

	
    opt = RDOptionsPrep(m, ///
	                    kernel,h, ///
						opt_criterion, ///
                        alpha,beta, ///
						se_method, ///
						j,///
						order)  

    //ret = RDResults()
	
  }
 

  
end 




*********************************************
//          sub-programs from ivregress ados 
//          extracting endog exog iv and y
*********************************************
program _iv_parse, sclass
	
	local n 0

	gettoken lhs 0 : 0, parse(" ,[") match(paren) bind
	if (strpos("(",`"`lhs'"')) {
		fvunab lhs : `lhs'
		if `:list sizeof lhs' > 1 {
			gettoken lhs rest : lhs
			local 0 `"`rest' `0'"'
		}
	}
	IsStop `lhs'
	if `s(stop)' { 
		error 198 
	}  
	_fv_check_depvar `lhs'
	while `s(stop)'==0 {
		if "`paren'"=="(" {
			local n = `n' + 1
			if `n'>1 {
				capture noi error 198
di as error `"syntax is "(all instrumented variables = instrument variables)""'
				exit 198
			}
			gettoken p lhs : lhs, parse(" =") bind
			while "`p'"!="=" {
				if "`p'"=="" {
					capture noi error 198
di as error `"syntax is "(all instrumented variables = instrument variables)""'
di as error `"the equal sign "=" is required"'
					exit 198
				}
				local end`n' `end`n'' `p'
				gettoken p lhs : lhs, parse(" =") bind
			}
			/* An undocumented feature is that we can specify
			   ( = <insts>) with GMM estimation to impose extra
			   moment conditions 
			*/ 
			if "`end`n''" != "" {
				fvunab end`n' : `end`n''
			}
			fvunab exog`n' : `lhs'
		}
		else {
			local exog `exog' `lhs'
		}
		gettoken lhs 0 : 0, parse(" ,[") match(paren) bind
		IsStop `lhs'
	}
	mata: st_local("0",strtrim(st_local("lhs")+ " " + st_local("0")))

	fvunab exog : `exog'
	fvexpand `exog'
	local exog `r(varlist)'
	tokenize `exog'
	local lhs "`1'"
	local 1 " "
	local exog `*'
	
	// Eliminate vars from `exog1' that are in `exog'
	local inst : list exog1 - exog
	if ("`end1'" != "") {
		fvunab end1 : `end1'
		fvexpand `end1'
		local end1 `r(varlist)'
	}
	
	// `lhs' contains depvar, 
	// `exog' contains RHS exogenous variables, 
	// `end1' contains RHS endogenous variables, and
	// `inst' contains the additional instruments
	// `0' contains whatever is left over (if/in, weights, options)
	
	sret local lhs `lhs'
	sret local exog `exog'
	sret local endog `end1'
	sret local inst `inst'
	sret local zero `"`0'"'

end

// Borrowed from ivreg.ado	
program define IsStop, sclass

	if `"`0'"' == "[" {
		sret local stop 1
		exit
	}
	if `"`0'"' == "," {
		sret local stop 1
		exit
	}
	if `"`0'"' == "if" {
		sret local stop 1
		exit
	}
	if `"`0'"' == "in" {
		sret local stop 1
		exit
	}
	if `"`0'"' == "" {
		sret local stop 1
		exit
	}
	else {
		sret local stop 0
	}

end






**********************************************
//            Mata class declarations     
***********************************************

mata:



//*********************************************
//      constructs RD dataset format     
//*********************************************

class RDData {
  	
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
  


class RDData scalar RDDataPrep(real matrix X, real matrix Y, real scalar c ,| real vector sigma2) {
    /* constructs df */
	
    class RDData scalar df  // create an instance of RDData
    if(args()==4) df.setup(X,Y,c,sigma2)
    else df.setup(X,Y,c)
	
    return(df)
  }


//*********************************************
//      constructs Options     
//      option c cutoff omiited 
//*********************************************


//* this class stores all the options input by the user
 
 
  class RDOptions {
    /* variables */
    real scalar m, h, alpha, beta, j, order
    string kernel, opt_criterion, se_method

    /* functions */
    void setup()
  }

  
  void RDOptions::setup(real scalar u_m ,| 
    string u_kernel, real scalar u_h,
    string u_opt_criterion, real scalar u_alpha, real scalar u_beta,
    string u_se_method, real scalar u_j,
    real scalar u_order) {

    /* assign values and defaults */
		
	//*
    m = u_m
	
	//* 
    if(args() >= 2 & u_kernel != "") kernel = u_kernel
    else kernel = "tri"

    if(args() >= 3 & u_h != .)	h = u_h
    else hp = . /* need to explicitly specify bandwidth is missing */


    if(args() >= 4 & u_opt_criterion != .)	opt_criterion = u_opt_criterion
    else opt_criterion = ""
    /* need to explicitly specify opt_criterion is missing */

    if(args() >= 5)					alpha = u_alpha
    else alpha = 0.05
    if(args() >= 6)					beta = u_beta
    else beta = 0.8
    if(args() >= 7)					se_method = u_se_method
    else se_method = "nn"
    if(args() >= 8)				j = u_j
    else j = 3
    if(args() >= 9)				order = u_order
    else order = 1

  }
  
   class RDOptions scalar RDOptionsPrep(real scalar u_m , 
    string u_kernel, real scalar u_h,
    string u_opt_criterion, real scalar u_alpha, real scalar u_beta,
    string u_se_method, real scalar u_j,
    real scalar u_order) {
    /* constructs opt */
	
    class RDOptions scalar opt // create an instance of RDData
	
    opt.setup(u_m , u_kernel, u_h, u_opt_criterion, u_alpha, u_beta,
    u_se_method, u_j,u_order)
	
    return(opt)
  }
  

//*********************************************
//      constructs Results format     
//*********************************************

struct RDResults {
      real scalar estimate
      real scalar lff /* not sure what this is */
      real scalar bias, sd, lower, upper, hl, eo, hp, hm, naive
  }
  




end 


