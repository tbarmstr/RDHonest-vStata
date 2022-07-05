*! v0.2.0 KHLee JUN-29-2018
/* This file contains the definitions of RDData, RDOptions and RDResults, and
forces Stata to load them into memory. */

version 14.2
set matastrict on
set more off

/*
clear all

use lee08, replace
*/

mata mata clear

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
    if(args() >= 2) {
      bw_equal = (u_bw_equal > 0 ? 1 : 0)
    }

    if(args() >= 3 & u_kernel != "") kernel = u_kernel
    else kernel = "tri"

    if(args() >= 4 & u_hp != .)	hp = u_hp
    else hp = . /* need to explicitly specify bandwidth is missing */

    hm = ((bw_equal == 1 | (args() >= 5 & u_hm == .)) ? hp : u_hm)

    if(args() >= 6 & u_opt_criterion != .)	opt_criterion = u_opt_criterion
    else opt_criterion = ""
    /* need to explicitly specify opt_criterion is missing */

    if(args() >= 7)					alpha = u_alpha
    else alpha = 0.05
    if(args() >= 8)					beta = u_beta
    else beta = 0.8
    if(args() >= 9)					se_method = u_se_method
    else se_method = "nn"
    if(args() >= 10)				j = u_j
    else j = 3
    if(args() >= 11)				sclass = u_sclass
    else sclass = "H"
    if(args() >= 12)				order = u_order
    else order = 1
    if(args() >= 13)				se_initial = u_se_initial
    else se_initial = "IKEHW"

  }
  struct RDResults {
      real scalar estimate
      real scalar lff /* not sure what this is */
      real scalar bias, sd, lower, upper, hl, eo, hp, hm, naive
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
  
  struct RDOptBWOutput {
    real scalar hp, hm
    real vector sigma2p, sigma2m
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
      kernWeights = EqKernWeight(x/h,kernel,order,0)
      w = ((invsym(u_Gamma) * (kernWeights :* R)')[1,])'
      return(w)
  }

  class RDLPregOutput extends LPRegOutput {

      public:
          real vector wm /* weights for negative observations */
          real vector wp /* weights for positive observations */
          real vector sigma2m, sigma2p /* squared residuals */
          real vector se /* standard error vector */
          real matrix GammaM, GammaP

      private:
          real vector sigma2
          real vector var
          real vector w
          real matrix Gamma
  }
end

/*
RDHonest voteshare margin, m(0.1) hp(10) kernel("uni") ///
se_method("nn") j(3)

supplied_var test:
RDHonest voteshare margin, m(0.1) hp(10) kernel("uni") ///
se_method("supplied_var") j(3)

*/

// command to run this do file:
// do RDHonest_build.do
