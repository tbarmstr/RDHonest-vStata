*! rdhonest v0.1 by KHLee MAY-27-2019
*! rdhonest v0.2 by YGChen and SZhang, Aug 2022

program define rdhonest, eclass byable(onecall)
	if _caller() >= 11 {
		local vv : di "version " string(_caller()) ":"
	}
	version 14

	if _by() {
		local BY `"by `_byvars'`_byrc0':"'
	}

	if replay() {
		if `"`e(cmd)'"' != "rdhonest" { 
			error 301
		}
		else if _by() { 
			error 190 
		}
		else {
			Display `0'
		}
		exit
	}
	`vv' ///
	`BY' Estimate `0'
	ereturn local cmdline `"rdhonest `0'"'
end


*================================================
* main function
*================================================

program Estimate, eclass byable(recall) sortpreserve

	local vv : di "version " string(_caller()) ":"
	version 14

	* parsing ===================================
	rdparse `0'
	
	local depvar `s(depvar)'
	local runvar `s(runvar)'
	local covar `s(covar)'
	local treat `s(treat)'
	local 0 `s(zero)'
	
	syntax [if] [in] /// 
	  [,c(real 0) /// cutoff, default: 0
		m(numlist >0 min=1) /// Bound on second derivative of the conditional mean function, default: ROT calculation (add ">0" to allow only positive inputs)
		kernel(string) /// kernel function used, choose from triangular, uniform, optimal, default: triangular
		opt_criterion(string) /// optimality criterion, choose from MSE, FLCI, OCI, default: MSE
		h(real 0) /// bandwidth, default: optimal level calculated following opt_criterion
		se_method(string) /// choose from "NN", "EHW", if not specified but PVARiance is specified, use user-supplied variance, ) 
		alpha(real 0.05) /// confidence level, default: 0.05
		beta(real 0.8) /// quantile of excess length to optimize, default: 0.8
        j(integer 3) /// number of nearest neighbors if se_method == "NN"
		t0(real 0) /// used for fuzzy FD, initial estimate of the treatment effect for calculating the optimal bandwidth, default: 0
		/// pinf(integer 0) /// =1 if user wants to do inference at cutoff point instead of RD
		PVARiance(varlist) /// 
		CLuster(varlist) /// clustering variables 
		Wgtvar(varname) /// weighting variables
		SAVEWgtest(string) ///
		NOPARAMeter /// not show paramaters used
		ITERLog /// show iteration log
	   ]

	_get_diopts diopts, `options'

	// class of regressions
	local sRD `=("`treat'"=="")' // sharp RD
	local fRD `=("`treat'"!="")' // fuzzy RD
	// fix parsing issues for fuzzy RD
	if `fRD' {
		local treat `s(runvar)'
		local runvar `s(treat)'
	}

	// Indicator: no covariates, nocovar_ind =1
	local nocovar_ind `="`covar'"==""'
    
	// Indicator: user supply preliminary variance, prevar_ind=1
	local prevar_ind `="`pvariance'" != ""'
	
	// Indicator: cluster specified, cluster_ind=1
	local cluster_ind `="`cluster'" != ""'
	
	// Indicator: weights of obs provided, weight_ind=1
	local weight_ind `="`wgtvar'" != ""'	

	marksample touse
	markout `touse' `devpar' `covar' `runvar' `treat' `pvariance' `cluster' `wgtvar', strok

	set more off // set more off

	// populate lists based on fully expanded list
	local totexp `covar' `runvar' `treat'
	fvexpand `totexp' if `touse'
	local totexp `r(varlist)'
	fvexpand `runvar' if `touse'
	local runvar `r(varlist)'
	fvexpand `covar' if `touse'
	local covar `r(varlist)'
	fvexpand `treat' if `touse'
	local treat `r(varlist)'
	foreach vargroup in runvar covar treat {
		local there: list `vargroup' & totexp
		local leftover: list `vargroup' - there
		foreach var of local leftover {
			_ms_put_omit `var'
			local there `there' `s(ospec)'
		}
		local `vargroup' `there'
	}

	// sort by cluster id for clustering
	if (`cluster_ind') {
		sort `cluster'
	}

	* input control =============================

	/* specification errors */
	// If x in both covar and runvar
	local both : list covar & runvar
	foreach x of local both {
		di as err "`x' specified as both a running variable and a covariate"
		exit 498
	}
	// If x in both runvar and treat
	local both : list runvar & treat
	foreach x of local both {
		di as err "`x' specified as both a running variable and a treatment variable"
		exit 498
	}
	// If x on both LHS and (RHS or inst)
	local both : list depvar & runvar
	if "`both'" != "" {
		di as err "`both' specified as both the dependent variable and a running variable"
		exit 498
	}
	// If x in both depvar and covar
	local both : list depvar & covar
	if "`both'" != "" {
		di as err "`both' specified as both the dependent variable and a covariate"
		exit 498
	}
	// If x on both depvar and treat
	local both : list depvar & treat
	if "`both'" != "" {
		di as err "`both' specified as both the dependent variable and a treatment variable"
		exit 498
	}

	/* option errors */ 
	
	* check collinearity ========================
	qui{
		// check collinearity: save for later
	}
	
	* set default values ========================
	if ("`m'" == ""){
		local m -1 // replaced by ROT calculation 
		local m_opted 1 // mark M optimization for display
	}
	if ("`kernel'" == "") local kernel "tri"
	if ("`opt_criterion'" == "") local opt_criterion "MSE"
	if (`h' == 0) local bw_opted 1 // mark bandwidth optimization for display
	
	if (`prevar_ind') {
		local se_method "supplied_var"
		// di "user-supplied variance is used for calculation"
	}
	else if (!`prevar_ind') {
		if ("`se_method'" == "") local se_method "NN"
		else if ("`se_method'" != "NN") & ("`se_method'" != "EHW"){
			noi di in red "se_method() mis-specified: only NN and EHW are allowed"
		}
	}
	
	// cap drop `savewgtest'
	if ("`savewgtest'"!="") & _bylastcall() {
		confirm new var `savewgtest'
	}
	if ("`savewgtest'"=="") local savewgtest "NA"  

	* bylastcall macro 
	local Bylastcall = _bylastcall()
	
	* start mata
	set matastrict on
	mata{
		if (`fRD') Y = st_data(.,("`depvar'","`treat'"));
		if (`sRD') Y = st_data(.,("`depvar'"));
		X = st_data(.,("`runvar'"))	
		if (!`nocovar_ind') Z = st_data(.,("`covar'"));	
		
		sample = st_data(.,"`touse'")
		ID = J(1,1,1..rows(X))'

		rdclass = (`sRD'>0 ? "srd":"frd")
		
		if (`prevar_ind') sigma2 = st_data(.,("`pvariance'"));
		if (!`prevar_ind') sigma2 = J(rows(Y),cols(Y)^2,.);
		if (`cluster_ind') cluster = st_data(.,("`cluster'"));
		if (!`cluster_ind') cluster = J(rows(Y),1,.);
		if (`weight_ind') weight = st_data(.,("`wgtvar'")); 
		if (!`weight_ind') weight = J(rows(Y),1,1);
		rho = J(1,cols(Y)^2,.);

		// select sample
		id = select(ID,sample:==1);
		Y = select(Y,sample:==1);
		X = select(X,sample:==1)
		if (!`nocovar_ind') Z = select(Z,sample:==1)
		sigma2 = select(sigma2,sample:==1)
		cluster = select(cluster,sample:==1)
		weight = select(weight,sample:==1)

		// initialize data frame and map Y and X in
		df = RDDataPrep(id,X,Y,`c',sigma2,weight,cluster,rho,rdclass)

		// initialize options
		m=(`=subinstr("`m'"," ",",",.)')
		kernel = "`kernel'"
		opt_criterion = "`opt_criterion'"
		alpha = `alpha'
		beta = `beta'
		se_method = "`se_method'"
		j = `j'
		order = 1
		h = `h'
		t0 = `t0'
		
		// log display option
		if ("`iterlog'"=="") nolog = 1
		else if ("`iterlog'"!="") nolog = 0

		opt = RDOptionsPrep(m,      ///
							kernel, ///
							opt_criterion, ///
							h, ///
							se_method, ///
							alpha, ///
							beta, ///
							j, ///
							order, ///
							t0,
							nolog)  


		// return
		ret = RDResults()

		if (("`savewgtest'"!="NA") & `Bylastcall'){
			ret = RDHonest_fit(df,opt,kernC,"`savewgtest'",ID,sample)
		}
		else{
			ret = RDHonest_fit(df,opt,kernC)
		}
	}

	* ereturn
	// esample
	eret post, esample(`touse')
	
	// command and variable names
	eret local cmd "rdhonest"  
    eret local depvar  "`depvar'"	
    eret local runvar "`runvar'" 
	if (`fRD') eret local treat "`treat'"
	if (`nocovar_ind') eret local covar "`covar'"
	if (`cluster_ind') eret local cluster "`cluster'"
	if (`weight_ind') eret local weight "`weight'"
	eret local savewgtest "`savewgtest'"

	// RD class
	if (`sRD') eret local rd = "Sharp" 
	if (`fRD') eret local rd = "Fuzzy"
	
	// parameters
	eret scalar M = M1
	if (`fRD') eret scalar M_fs = M2
	
	if ("`kernel'"=="epa"){
		eret local kernel = "{ul:epa}nechnikov"
	}
	else if("`kernel'"=="uni"){
		eret local kernel = "{ul:uni}form"
	}
	else if("`kernel'"=="tri"){
		eret local kernel = "{ul:tri}angular"
	}
	else{
		eret local kernel = "`kernel'"
	}
	
	eret local opt_crit = "`opt_criterion'" 
	eret scalar bandwidth = h
	eret local se_method = "`se_method'" 
	eret scalar alpha = `alpha'
	if ("`opt_criterion'"=="OCI") eret scalar beta = `beta'
	
	// estimates
	eret scalar est = estimate
	if (`fRD') eret scalar est_fs = fs
	eret scalar se = stdErr
	eret scalar bias = bias 
	eret scalar HLCi = HLCi 
	
	// CI
	eret scalar TCiL = lEqLimit 
	eret scalar TCiU = uEqLimit
	eret scalar OCiL = lLimit
	eret scalar OCiU = uLimit
	
	// other scalars
	eret scalar effObs = effObs
	eret scalar Leverage = Leverage
	eret scalar cutoff = `c'

	local title "Honest inference: {res:`=upper("`e(rd)'")'} Regression Discontinuity"
	eret local title `title'

	Display, `noparameter' fRD(`fRD') opt_bw(`bw_opted') opt_m(`m_opted') cluster_ind(`cluster_ind')
end



*==========================================================
* Syntax auxilary functions
*==========================================================

* parsing function ========================================

// Borrrowed from _iv_parse.ado, with modification
/* The syntax for rdhonest is: 
fuzzy RD: Y (X=Z) W, options
RD: Y X W
where
Y: depvar, X: runvar, Z: treat, W: covar
*/

program rdparse, sclass

	//syntax cmd_to_parse [if] [in]
	//noi di "`cmd_to_parse'"
	//marksample touse
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
				di as error `"syntax is "(one treatment variable = one running variable)""'
				exit 198
			}
			gettoken p lhs : lhs, parse(" =") bind
			while "`p'"!="=" {
				if "`p'"=="" {
					capture noi error 198
					di as error `"syntax is "(one treatment variable = one running variable)""'
					di as error `"the equal sign "=" is required"'
					exit 198
				}
				local runvar`n' `runvar`n'' `p'
				gettoken p lhs : lhs, parse(" =") bind
			}
			if "`runvar`n''" != "" {
				fvunab runvar`n' : `runvar`n''
				if `:list sizeof runvar`n'' > 1 {
				capture noi error 198
				di as error `"syntax is "(one treatment variable = one running variable)""'
				di as error `"only one treatment variable allowed"'
				}
			}
			fvunab covar`n' : `lhs'
			if `:list sizeof covar`n'' > 1 {
			capture noi error 198
			di as error `"syntax is "(one treatment variable = one running variable)""'
			di as error `"only one running variable allowed"'
			}
		}
		else {
			local covar `covar' `lhs'
		}
		gettoken lhs 0 : 0, parse(" ,[") match(paren) bind
		IsStop `lhs'
	}
	mata: st_local("0",strtrim(st_local("lhs")+ " " + st_local("0")))

	fvunab covar : `covar'
	fvexpand `covar'
	local covar `r(varlist)'
	tokenize `covar'
	local lhs "`1'"
	local 1 " "
	local covar `*'
	
	// Eliminate vars from `covar1' that are in `covar'
	local treat : list covar1 - covar
	if ("`runvar1'" != "") {
		fvunab runvar1 : `runvar1'
		fvexpand `runvar1'
		local runvar `r(varlist)'
	}
	else {
		tokenize `covar'
		local runvar "`1'"
		local 1 " "
		local covar `*'
	}
	
	// `lhs' contains the dependent variable
	// `runvar' contains the running variable
	// `treat' contains the treatment variable
	// `covar' contains RHS covariates
	// `0' contains whatever is left over (if/in, options)
	
	sret local depvar `lhs'
	sret local covar `covar'
	sret local runvar `runvar'
	sret local treat `treat'
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

* display function ========================================

program Display

	syntax , [NOPARAMeter fRD(real 0) opt_bw(real 0) opt_m(real 0) cluster_ind(real 0)*]
	
	local C1 "_col(2)"
	local C2 "_col(16)"
	local C3 "_col(34)"
	local C4 "_col(50)"
	local C5 "_col(67)"

	di
	di _col(1) as text "`e(title)'"
	/* main table */
	di "{hline 78}"
	di `C1' as text %-14s "Estimate" ///
	   `C2' as text %-14s "Maximum Bias" ///
	   `C3' as text %-14s "Std. Error" ///
	   `C4' as text %-2s "[" %2.0f `=(1-`e(alpha)')*100' "% Conf." ///
	   `C5' as text %9s "intervals ]"
	di "{hline 78}"
	di `C1' as res %-11.0g e(est) ///
	   `C2' as res %-11.0g e(bias) ///
	   `C3' as res %-11.0g e(se) ///
	   `C4' as res %-11.0g e(TCiL) ///
	   `C5' as res %11.0g e(TCiU)
	di "{hline 78}"
	/* other estimation results: one-sided CIs, optimal h, effective observations */
	di `C1' as text %2.0f `=(1-`e(alpha)')*100' "% One-sided Conf. intervals:" " (" ///
			as res %-11.0g e(OCiL) as text ", Inf), " ///
			as text "(-Inf, " ///
			as res %11.0g e(OCiU) as text ")"	
	if (`fRD') {
		di `C1' as text "First-stage estimate: " as res %-11.0g e(est_fs)
	}
	di `C1' as text "Bandwidth" _continue
	if (`opt_bw') {
		di as text " (optimized): " as res %-11.0g e(bandwidth)
	}
	else {
		di as text ": " as res %-11.0g e(bandwidth)
	}

	di `C1' as text "Number of effective observations: " as res %-11.0g e(effObs)
	
	/* if asked, display other paramters */
	if "`noparameter'" == ""  {
		di
		di `C1' as text "{it:Parameters}:"
		di `C1' as text "Cutoff: " as res %-11.0g e(cutoff)
		di `C1' as text "Kernel: " as res "`e(kernel)'"
		if ("`e(opt_crit)'"=="OCI"){
			di `C1' as text "Optimization criterion: " as res "`e(opt_crit)'" ///
					as text ", with beta " as res %-11.0g e(beta)
		}
		else {
			di `C1' as text "Optimization criterion: " as res "`e(opt_crit)'" 
		}
		
		if ("`e(se_method)'"=="supplied_var"){
			di `C1' as text "Standard error estimation method: " as res "user-supplied variance"
		}
		else{
			di `C1' as text "Standard error estimation method: " as res "`e(se_method)'"
		}
		
		di `C1' as text "Maximum leverage for estimated parameter: " as res %-11.0g e(Leverage)
		if (e(Leverage)>0.1) {
			di `C1' in r "Warning: Maximal leverage is large."
			di `C1' in r "Inference may be inaccurate. Consider using bigger bandwidth."
		}
		
		if (`fRD') {
			di `C1' as text "Smoothness constant M (first-stage" _continue
			if (`opt_m'){
				di as text ", rule of thumb): " as res %-11.0g e(M_fs)
			}
			else {
				di as text "): " as res %-11.0g e(M_fs)
			}
			di `C1' as text "Smoothness constant M (reduced-form" _continue
			if (`opt_m'){
				di as text ", rule of thumb): " as res %-11.0g e(M)
			}
			else {
				di as text "): " as res %-11.0g e(M)
			}
		}
		else{
			di `C1' as text "Smoothness constant M" _continue
			if (`opt_m'){
				di as text " (rule of thumb): " as res %-11.0g e(M)
			}
			else {
				di as text ": "  as res %-11.0g e(M)
			}
		}
	}
	
	/* variables */
	di "{hline 78}"
	di `C1' as text %21s "Dependent variable: " _continue
	DisplayVar `e(depvar)'
	di `C1' as text %21s "Running variable: " _continue
	DisplayVar `e(runvar)'

	if (`fRD') {
		di `C1' as text %21s "Treatment variable: " _continue
		DisplayVar `e(treat)'
	}

	if (`cluster_ind'){
		di `C1' as text %21s "Clustered by: " _continue
		DisplayVar `e(cluster)'
	}
	
	if ("`e(savewgtest)'"!="NA"){
		di
		di `C1' as text "{it:Generated variables}:"
		di `C1' as text %21s "Estimation weight" _continue
		if _by() {
			di as text " (only for the last by-group): " as res "`e(savewgtest)'"
		}
		else{
			di as text ": " as res "`e(savewgtest)'"
		}
	}
end

* auxilary display functions: for single variable
program define DisplayVar
	
	local Cvar "_col(22)"

    args vname 
	local abname = abbrev("`vname'", 20)
	di `Cvar' as res %-40s "`abname'"
	
end 

*==========================================================
* Mata codes
*==========================================================
set matastrict on

mata:
	mata clear

	// 1. constructs RD dataset format==========================
	class RDData {
		/* variables */
		real matrix X,Y,sigma2
		real vector cluster, weight, m, p, ID
		real vector rho
		real scalar cutoff
		string rdclass

		/* functions */
		void setup()
	}

	void RDData::setup(real vector ID, real matrix X, real matrix Y, real scalar c, ///
		real matrix sigma2, real vector weight, real vector cluster, real vector rho, string rdclass) {

		real matrix s

		/* assign matrices */
		this.cutoff = c
		this.X = X :- c ; this.Y = Y 
		this.sigma2 = sigma2 ; this.weight = weight
		this.m = (this.X:<0) ; this.p = (this.X:>=0)
		this.cluster = cluster
		this.ID = ID
		
		/* supplied variances */
		if(rows(sigma2)!=length(X)) _error("length of supplied variances does not match data");

		/* supplied weights */  
		if(length(weight)!=length(X)) _error("length of supplied weights does not match data")

		/* supplied weights */ 
		if(length(cluster)!=length(X)) _error("length of supplied cluster ID does not match data")
		  
		/* sort data: by cluster id and X */
		s = sort((this.X, this.m, this.p, this.weight, this.cluster, this.Y, this.sigma2, this.ID),(5,1))
		this.X = s[.,1];  this.m = s[.,2]; this.p = s[.,3]
		this.weight= s[.,4] ;  this.cluster = s[.,5] 
		this.Y = s[.,(6..(cols(Y)+5))] ; this.sigma2 = s[.,((cols(Y)+6)..(cols(sigma2)+cols(Y)+5))];
		this.ID = s[.,cols(s)]

		/* class */
		this.rdclass = rdclass 

		this.rho = rho
	}

	class RDData scalar RDDataPrep(real vector ID, real matrix X, real matrix Y, real scalar c , ///
								   real matrix sigma2, real vector weight, ///
								   real vector cluster, real vector rho, string rdclass) {
		
		/* constructs df */	
		class RDData scalar df  // create an instance of RDData
		
		df.setup(ID,X,Y,c,sigma2,weight,cluster,rho,rdclass)
	  
		return(df)
	}

	// 2. constructs RD options=================================

	class RDOptions {
		/* variables */
		real scalar  h, alpha, beta, j, order, t0, nolog
		string kernel, opt_criterion, se_method
		real vector m
		/* functions */
		void setup()
	}

	void RDOptions::setup(real vector u_m , 
		string u_kernel, string u_opt_criterion, real scalar u_h,
		string u_se_method, real scalar u_alpha, real scalar u_beta,
		real scalar u_j, real scalar u_order,real scalar u_t0, 
		real scalar u_nolog) {

		/* assign values and defaults */
		m = u_m
		kernel = u_kernel
		opt_criterion = u_opt_criterion
		h = u_h
		se_method = u_se_method
		alpha = u_alpha
		beta = u_beta
		j = u_j
		order = u_order
		t0 = u_t0
		nolog = u_nolog
	}

	class RDOptions scalar RDOptionsPrep(real vector u_m , 
		string u_kernel, string u_opt_criterion, real scalar u_h,
		string u_se_method, real scalar u_alpha, real scalar u_beta,
		real scalar u_j, real scalar u_order, real scalar u_t0,
		real scalar u_nolog) {
		
		/* constructs opt */
		class RDOptions scalar opt // create an instance of RDData
		
		opt.setup(u_m, u_kernel, u_opt_criterion, u_h, u_se_method, 
				u_alpha, u_beta, u_j, u_order, u_t0, u_nolog)
		
		return(opt)
	}

	// 3. Regression============================================

	// 3.1 LPReg: local polynomial regression

	class LPRegOutput {
		/* contains output from LPReg */
		real vector theta /* estimates */
		real matrix sigma2 /* squared residuals */
		real vector var /* sampling variance of theta */
		real vector w /* kern weights plus obs weight */
		real scalar eo /* effective observations */
		real vector p , m /* below and above indicator */

		real matrix wgt  /* OLS weight on Y_i that gives theta */

		real matrix res /* residuals */		
		real matrix clu_setup /*clustering set up*/		  
	}

	class LPRegOutput scalar LPReg (real matrix X, real matrix Y, real matrix sigma2, 
									real matrix weight, real matrix cluster, real scalar h, 
									string kernel, real scalar order, string se_method, 
									real scalar j, real vector rho) {

		real scalar effObs
		real matrix R, Gamma, beta, res, hsigma2, dsigma2, nsigma2, var, clu_setup, Vaug
		/* weight: obs weights */ 
		/* wt: OLS weights give theta */
		/* W: kern weights plus obs weights */
		/* cluster_res: cluster selected by non-missing residual */
		real vector wgt, w, wgt_unif
		real vector p, m
		
		class LPRegOutput scalar output
		
		/* w is a vector -- kern weights plus obs weight */
		/* order = 0 ; boundary = 0 */
		w = EqKernWeight(X:/h,kernel,0,0):*weight 
			// st_matrix("w",w)
			// stata("matlist w")

		p = (X:>=0) ; m = (X:<0)
		
		/* form the regression data matrix according to FWL */
		R = J(1,order+1,X) :^ (0..order)
		R = (R:*p,R)
		
		Gamma = quadcross(R,w :* R)
		
		if(h == 0| det(Gamma)==0) {
			//printf("{red}Error: bandwidth too small or singular matrix.\n")
			output.theta = J(1,cols(Y),0)
			output.sigma2 = J(rows(Y),cols(Y)^2,0)
			output.var = J(1,cols(Y)^2,0)
			output.w = w
			output.eo = 0
			output.wgt = J(rows(X),1,0)
			output.res = J(rows(Y),1,0)
			output.p = p
			output.m = m
			output.clu_setup = panelsetup(cluster,1)
		}
		else{
			wgt = ((invsym(Gamma) * (w :* R)')[1,.])'
			/* To compute effective observations, rescale against uniform kernel*/
			/* use weight rather than kernel + weight */
			wgt_unif = ((invsym(quadcross(R,weight:* R)) * (weight :* R)')[1,.])'

			/* estimates from using OLS weights */
			/* squared residuals allowing for Y being multi-variates */
			beta = (invsym(Gamma) * (w :* R)' )* Y
			res = (Y - R*beta)

			/* Robust variance-based formula */
			// supplied_var: use use-supplied squared residuals to compute EHW variance 
			if ( max((cluster:==.)) ) {
				if(strpos(se_method,"supplied_var") > 0) {
					var = colsum(wgt:^2:*sigma2)
					output.sigma2 = sigma2
				}

				if(strpos(se_method,"EHW") > 0 | strpos(se_method,"ehw") > 0) {
					hsigma2=res[.,(J(1, cols(res), (1..cols(res))))]:*
					res[.,( vec(J(cols(res),1, (1..cols(res))))')]
					var = colsum(wgt:^2:*hsigma2)
					output.sigma2 = hsigma2
				}

				if(strpos(se_method,"NN") > 0 | strpos(se_method,"nn") > 0 ) {
					nsigma2 = sigmann(X,Y,j,weight)
					var = colsum(wgt:^2:*nsigma2)
					output.sigma2 = nsigma2
				}
			}
			else {
				clu_setup = panelsetup(cluster,1)
				if(strpos(se_method,"supplied_var") > 0) {
					output.sigma2 = sigma2
					var = colsum(wgt:^2:*sigma2) :+ rho*( sum(panelsum(wgt,clu_setup):^2) - sum(wgt:^2) )
				}
				
				else {
					hsigma2=res[.,(J(1, cols(res), (1..cols(res))))]:*
					res[.,( vec(J(cols(res),1, (1..cols(res))))')]
					output.sigma2 = hsigma2
					Vaug = panelsum(wgt:*res, clu_setup)
					var = vec( cross(Vaug,Vaug) )'
				}
			}
			
			/* catch error when no variance is provided */
			if( max((var:==.)) ) _error("No variance was computed by LPReg")

			/* return output */
			output.theta = beta[1,.]
			output.var = var
			output.w = w
			output.eo = length(X)*sum(wgt_unif:^2)/sum(wgt:^2) 
			output.wgt = wgt
			output.res = res
			output.p = p
			output.m = m
			output.clu_setup = clu_setup
		}
		
		//printf("var is %9.4f ",output.var)
		return(output)
	}


	// 3.2 RDLRReg: RD local polynomial regression

	class RDLPregOutput extends LPRegOutput {
		real scalar se /* standard error vector */
		real scalar fs /* 1st stage estimate (for T) */
		real scalar estimate /* treatment effect */
		real vector kw /* kernel weight */
	}

	class RDLPregOutput scalar RDLPreg(class RDData scalar df, real scalar h,
		| string kernel, real scalar order, string se_method,
			real scalar no_warning, real scalar j) {

		/* check for errors */
		if (h <= 0) _error("Non-positive bandwidth h")

		/* set defaults */
		if (args() <= 3 & kernel == "") kernel = "tri"
		if (args() <= 4 & order == .) order = 1
		if (args() <= 5 & se_method == "") se_method = "NN"
		if (args() <= 6 & no_warning == .) no_warning = 0
		if (args() <= 7 & j == .) j = 3

		/* variable declarations */
		real scalar plugin /* to be implemented later */
		real matrix Y, X, w, sigma2, Xm, Xp, weight, cluster
		real vector rho
		
		class RDLPregOutput scalar output
		class LPRegOutput scalar r1

		/* set rho */
		rho = df.rho

		/* set kernel weights */
		w = (h <= 0 ? (0 :* df.X) : (EqKernWeight(df.X:/h,kernel,0,0)):*df.weight)
		
		/* variance calculations are faster if we only keep data with positive
		kernel weights */
		X = select(df.X,w :> 0) ; Y = select(df.Y,w :> 0)
		sigma2 = select(df.sigma2,w :> 0) ; weight = select(df.weight,w :> 0) 
		cluster = select(df.cluster,w :> 0)

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
		r1 = LPReg(X, Y, sigma2, weight, cluster, h, kernel, order, se_method, j, rho)
	  
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
		output.res = r1.res
		output.clu_setup = r1.clu_setup

		if (df.rdclass == "frd") {	
			output.fs = r1.theta[2]
			output.estimate = r1.theta[1]/r1.theta[2]
			output.se = sqrt(quadcross(
					   (1,-output.estimate,-output.estimate,output.estimate^2)',r1.var')
					   /output.fs^2)		
		}

		return(output)
		
	}


	// 4. Rule of thumb choosing M =============================

	real vector MROT_fit(real matrix X, real matrix Y, | real scalar IsIP) {
			
		if (args()<3) IsIP = 0 

		real vector M, col1, col2, col3, col4, col5, m3coeff
		real matrix  m3mat, Xm ,Xp , Ym, Yp 
		
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
			/* sharp RD: on either side */
			Xp = select(X,X:>=0); Yp = select(Y,X:>=0) 
			Xm = select(X,X:<0); Ym = select(Y,X:<0) 
			
			M = max((MROT_fit(Xm,Ym,1),MROT_fit(Xp,Yp,1)))
		}
		else {
			/* fuzzy RD: compute M for two reduced form regressions */
			M = (0,0)

			/* reduced form regression */
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
	
	// Auxilary function for MROT_fit step 2
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

	// 5. Preliminary variance used for h ======================

	class RDData scalar IKBW_fit(class RDData scalar df, real matrix kernC,| string kernel,
		real scalar verbose) {
		/* Calculate bandwidth for sharp RD based on local linear regression
		using method by Imbens and Kalyanaraman (ReStud 2012) */

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

	// 5.1 build preliminary variance estimation

	class RDPrelimVarOutput{
		real vector p , m, moul /* below and above indicator */
		real matrix res, sigma2 /* residuals */
	}

	class RDPrelimVarOutput scalar RDPrelimEst(class RDData scalar df, real matrix kernC, string se_initial) {
		
		real matrix X, Xp, Xm
		class RDLPregOutput scalar r1
		class RDPrelimVarOutput scalar output
		real scalar h1, hmin 
		class RDData scalar drf 
		real matrix clu_setup
		real vector moul

		/* concatenate Xm and Xp and compute "rule of thumb" bandwidth */
		X = df.X
		
		Xp = select(df.X,df.X:>=0) ; Xm = select(df.X,df.X:<0)
		
		hmin = max((uniqrows(Xp)[3], uniqrows(Xm)[3],
					sort(Xp,1)[4], sort(abs(Xm),1)[4]))
		h1 = max((1.84*sqrt(variance(X))/(sum(length(X)))^(1/5), hmin))
		
		/* fuzzy RD */
		drf = df
		if (df.rdclass == "frd") {
			drf.Y=df.Y[.,1]
			drf.rdclass = "srd"
		}
					 
		if (strpos(se_initial,"Silverman") > 0) {	
			if (cols(df.Y) == 1) {
				r1 = RDLPreg(df,h1,"uni",0,"EHW")
			} 
			else {	
				_error("This method for preliminary variance estimation is not supported.")
			}
		}
		else if (strpos(se_initial,"IKEHW") > 0) {
			h1 = IKBW_fit(drf, kernC)
			r1 = RDLPreg(df,max((h1,hmin)),"tri",1,"EHW")
		}
		else {
			_error("Unknown method for estimating initial variance.")
		}

		if ( max((df.cluster:!=.)) ){
			moul = moulton_est(r1.res, r1.clu_setup)
		}

		/* output */
		output.p = r1.p
		output.m = r1.m
		output.res = r1.res
		output.sigma2 = r1.sigma2
		output.moul = moul

		return(output)
	}

	class RDData scalar RDPrelimVar(class RDData scalar df, real matrix kernC,| string se_initial) {
		
		real matrix X, Xp, Xm
		real vector moul
		real scalar lm, lp, varm, varp
		class RDPrelimVarOutput scalar r1
		
		/* set defaults: EHW and IK bandwidth*/
		if(args()==2) se_initial = "IKEHW"

		/* concatenate Xm and Xp and compute "rule of thumb" bandwidth */
		X = df.X
		
		Xp = select(df.X,df.X:>=0) ; Xm = select(df.X,df.X:<0)
		r1 = RDPrelimEst(df, kernC, se_initial)
					 
		if (strpos(se_initial,"Silverman") > 0) {	
			if (cols(df.Y) == 1) {
				/* variance adjustment on either side */
				lp = sum(r1.p)
				lm = sum(r1.m)

				varp = sum(r1.sigma2:*r1.p)*1/(lp-1)
				varm = sum(r1.sigma2:*r1.m)*1/(lm-1)
				
				df.sigma2 = (df.X:<0):*varm + (df.X:>=0):*varp
			}
		}
		else if (strpos(se_initial,"IKEHW") > 0) {
			lp = sum(r1.p)
			lm = sum(r1.m)

			varp = colsum(r1.sigma2:*r1.p)/lp
			varm = colsum(r1.sigma2:*r1.m)/lm
			
			df.sigma2 = (df.X:<0):*J(rows(df.X),1,varm) + (df.X:>=0):*J(rows(df.X),1,varp)	
		}
		df.rho = r1.moul
		
		return(df)
	}

	// 5.2 Moulton estimate of rho for clustering ==============
	real vector moulton_est(real matrix res, real matrix clu_setup) {

		real scalar den
		real matrix res_aug
		real vector moul

		den = sum( panelsum(J(rows(res),1,1),clu_setup):^2 ) - rows(res)
		
		if (den > 0){
			res_aug = panelsum(res,clu_setup)
			moul = vec(cross(res_aug,res_aug)-cross(res,res)):/den
		}
		else{
			moul = J(cols(res)^2,1,0)
		}

		return (moul')
	}

	// 6. compute optimal h=====================================

	real vector RDOptBW_fit (class RDData scalar df, class RDOptions scalar opt,
		real matrix kernC) {

		transmorphic S /* to optimize function */
		real scalar h0, hmin, hmax, h /* ROT guess, minimal and maximal h */
		real scalar errorcode
		real vector h_opt, supp
		real matrix X

		/* compute Silverman rule of thumb bandwidth */
		X = df.X
		h0 = 1.06 * sqrt(variance(X)) * length(X)^(-1/5)

		/* initialize bandwidth optimization */
		S = optimize_init()
		optimize_init_argument(S,1,df)
		optimize_init_argument(S,2,opt)
		optimize_init_argument(S,3,kernC)
		optimize_init_which(S, "min")
		if (opt.nolog == 1) {
			optimize_init_tracelevel(S, "none")
		}

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
			 h = gss(&RDOptBWGss(), select(supp,supp :>= hmin), df, opt, kernC)
			 /* printf("Golden section bandwidth obtained is %9.0g \n",h) */
		}
		else {
			// set up parameters of optimization problem
			optimize_init_argument(S,4,hmin)
			optimize_init_argument(S,5,hmax)
			optimize_init_evaluator(S,&RDOptBWEval())
			optimize_init_params(S,h0)
			optimize_init_technique(S,"bfgs 10 nr 30") 
			optimize_init_verbose(S, 0)
			errorcode = _optimize(S)
			
			if (errorcode == 0){
			   // no error occurs
			   h_opt = optimize_result_params(S) 
			   h = h_opt
			}
			else {
			  // error occurs then call grid search the optimal bandwidth 
			   printf("{red}Grid search optimizer is called \n")
			   supp = sort(uniqrows(df.X),1)
			   h = gs(&RDOptBWGss(), select(supp,supp :>= hmin), df, opt, kernC)	
			}
		}
		return(h)
	}

	// 7. main function: RDHonest_fit===========================

	/* declare a class for results */
	struct RDResults {
		real scalar estimate, fs, leverage
		real scalar bias, sd, lower, upper, hl, eo, h, naive
	}

	// 7.1 NPRDHonest_fit

	struct RDResults scalar NPRDHonest_fit (class RDData scalar df,
		class RDOptions scalar opt, real matrix kernC, real scalar T0bias) {

		class RDLPregOutput scalar r1
		real scalar h, z
		real vector wp, wm, M
		real matrix XX, XXp, XXm 
		struct RDResults scalar results

		h = opt.h
				
		if(opt.h <= 0) {
			//printf("Bandwidth (h) missing or invalid, running RDOptBW_fit \n")
			h = RDOptBW_fit(df, opt, kernC)
		}

		/* run RD local polinomial regression */
		// Suppress warnings about too few observations 

		r1 = RDLPreg(df, h, opt.kernel, opt.order, opt.se_method, 1, opt.j)

		wp = select(r1.wgt,r1.p:==1)
		wm = select(r1.wgt,r1.m:==1)
		XX = select(df.X,r1.kw:>0)
		XXp = select(XX,r1.p:==1)
		XXm = select(XX,r1.m:==1)
		M = opt.m
		
		/* multiply bias and sd by r1.fs to make if free of first stage */   
		if (df.rdclass=="frd" & T0bias) { 
			r1.se = r1.se*abs(r1.fs)
			M = ((M[1]+M[2]*abs(opt.t0)), M) 
		}
		else if (df.rdclass=="frd" & !T0bias) {
			M = ((M[1]+M[2]*abs(r1.estimate)) / abs(r1.fs), M)
		}
				

		/* bias and one-sided CI ending point */
		// If bandwidths are too small, set bias to machine maximum
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

		return(results)
	}

	// 7.2 RDHonest_fit

	struct RDResults scalar RDHonest_fit (class RDData scalar df,
		class RDOptions scalar opt, real matrix kernC,
		| string optionVname, real vector ID, real vector sample) {

		/* declare variables */	
		struct RDResults scalar results
		class RDLPregOutput scalar regoutput 
		real matrix kw, wgt, tempID, Sample, est_w

		/* initial se estimate */
		if ( max((df.sigma2:==.))
		& (strpos(opt.se_method,"supplied_var") > 0 | opt.h == 0)) {
		//printf("Generating initial variance estimates via RDPrelimVar \n")
			// when no Prevar and no h provided 
			df = RDPrelimVar(df,kernC)
		}

		/* smoothing class constant */
		if(max(opt.m) < 0) {
		printf("Using Armstrong and Kolesar (2020) rule of thumb for smoothness constant M \n")
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
		// smoothness Constant 
		st_numscalar("M1",opt.m[1])

		if (df.rdclass=="frd"){		
		st_numscalar("fs",results.fs)
		st_numscalar("M1",opt.m[1])
		st_numscalar("M2",opt.m[2])		
		}

		/* generate a variable from the option savewgtest */
		if (args()== 6) {
			regoutput = RDLPreg(df, results.h, opt.kernel, opt.order, "ehw", 1, opt.j)
								
			wgt = regoutput.wgt
			kw = regoutput.kw

			/* put wgt into kw */
			tempID = J(1,1,1..rows(kw))'
			kw = ((kw:>0),tempID,kw)
			kw = sort(kw,(-1,2))
			kw[1..rows(wgt),3] = wgt
			kw = sort(kw,2)
			kw = kw[.,3]

			/* put kw into orginal sample to get vector savewgtest */
			kw = sort((df.ID,kw),1)
			Sample = sort((ID,sample),(-2,1))
			Sample[1..rows(kw),2] = kw[.,2]
			est_w = sort(Sample,1)
			est_w = est_w[.,2]

			/* back to Stata */	
			real scalar Vid 
			Vid = st_addvar("double",optionVname)
			st_store(.,optionVname,est_w) 
		}
	return(results)
	}
	
	// 7.3 	Optimalization function
	
	void RDOptBWEval (real scalar todo, real scalar b, class RDData scalar df,
		class RDOptions scalar opt, real matrix kernC, real scalar hmin, real scalar hmax,
		val, grad, hess) {

		struct RDResults scalar r
		class RDOptions scalar optLocal

		/* check to make sure bandwidth is feasible, else kicks the function back into range */
		if (b < hmin) b = hmin + 0.2*(hmax-hmin)/2
		else if (b > hmax) b = hmax - 0.2*(hmax-hmin)/2

		optLocal.setup(opt.m, opt.kernel, opt.opt_criterion, abs(b), "supplied_var",
			opt.alpha, opt.beta, opt.j, opt.order, opt.t0, opt.nolog)

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

	real scalar RDOptBWGss (real scalar h, class RDData scalar df,
		class RDOptions scalar opt, real matrix kernC) {

		struct RDResults scalar r
		class RDOptions scalar optLocal
		real scalar val

		/* options should be that uses PrelimVar stored in sigma2 */
		 optLocal.setup(opt.m, opt.kernel, opt.opt_criterion, h, 
			"supplied_var", opt.alpha, opt.beta, opt.j, opt.order,opt.t0, opt.nolog)

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


	// 8. accessary functions===================================

	real matrix sigmann(real matrix X, real matrix Y, real scalar j, real vector weight) {
			
		real matrix Xm ,Xp, Ym, Yp, sigma2, weightm, weightp
		
		Xm = select(X,X:<0) ; Xp = select(X,X:>=0) 
		Ym = select(Y,X:<0) ; Yp = select(Y,X:>=0)
		weightm = select(weight,X:<0) ; weightp = select(weight,X:>=0)  
		
		/* For RD, compute variance separately on either side of cutoff */
		sigma2 = sigmaNN(Xm,Ym,j,weightm)\sigmaNN(Xp,Yp,j,weightp)
			
		return(sigma2)	
		
	}


	real matrix sigmaNN(real matrix X, real matrix Y, real scalar j, real vector weight) {
		
		/* assume X is sorted */
		real scalar n,k,d,Jk
		real matrix sigma2,s,ind,YTemp

		n = length(X)
		sigma2 = J(n,cols(Y)^2,.)
		
		for(k = 1; k <= n; k++) {
		s = max((k-j,1))..k

		if (length(s)==1) {
			/* avoid error when trying to join a nonexistent vector with one that exists */
			d = (sort(abs(X[(k+1)..min((k+j,n))]:-X[k]),1))[j]
		}
		else if (k == n) {
			/* avoid error when trying to join a nonexistent vector with one that exists */
			d = (sort(abs((X[s[1..length(s)-1]]):-X[k]),1))[j]
		}
		else {
			d = (sort(abs((X[s[1..length(s)-1]]\X[(k+1)..min((k+j,n))]):-X[k]),1))[j]
		}
		
		ind = (abs(X :- X[k]) :<= d)
		ind[k] = 0
		Jk = sum(weight:*ind)
		/* compute mean of Y with positive indices */
		YTemp = (sort((ind,weight:*Y),-1))[1..sum(ind),2..(cols(Y)+1)]:/Jk
		sigma2[k,.] = vec(quadcross(Jk/(Jk+weight[k]) :* (Y[k,.] - colsum(YTemp)), (Y[k,.] - colsum(YTemp))))'
		}					
		return (sigma2)				
	}


	real vector EqKernWeight(real matrix Xh, string kernel,
			 real scalar order, real scalar boundary) {
		/* input Xh as X/h */
		real vector u, w, su 
		
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

		// iteratively apply f to support points to find optimal bandwidth
		h_opt = 1.0X+4e
		for (i=1; i<=length(supp); i++) {
			val = (*f)(supp[i],df,opt,kernC)
			if (val <= h_opt) {
				h_opt = val
				i_opt = i
			}
		}
		return(supp[i_opt])
	}
	  
	 /* grid search optimizer to use with mata objective functions */ 
	real scalar gs(pointer scalar f, real vector xs, class RDData scalar df,
		class RDOptions scalar opt, real matrix kernC) {

		real scalar  i, h_opt, i_opt, val
		real vector supp
		
		h_opt = 1.0X+4e

		// iteratively apply f to support points to find optimal bandwidth
		for (i=1; i<=length(xs); i++) {
			val = (*f)(xs[i],df,opt,kernC)
			if (val <= h_opt) {
				h_opt = val
				i_opt = i
			}
		}
		return(xs[i_opt])
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
