*! testRDHonest.do v2.0 KHLee MAY-27-2019

clear all
set more off

use lee08, replace

/*
local y voteshare
local x margin
local c = 0

mata {
  Y = st_data(.,("`y'"),0);	X = st_data(.,("`x'"),0)
  /* initialize data frame and map Y and X in */
  df = RDDataPrep(X,Y,`c')

  /* initialize options */
  m = 0.1
  bw_equal = 1
  kernel = "uniform"
  hp = 10;  hm = 10
  opt_criterion = "MSE"
  alpha = 0.05
  beta = 0.8
  se_method = "supplied_var"
  j = 3
  sclass = "H"
  order = 1
  se_initial = "IKEHW"

  opt = RDOptionsPrep(m,bw_equal,kernel,hp,hm,opt_criterion,
      alpha,beta,se_method,j,sclass,order,se_initial)

  ret = RDResults()
}
*/

// standard test
rdhonest voteshare margin, m(0.1) hp(10) ///
kernel("uni") se_method("nn") j(3)

// what if bandwidths are not equal?
rdhonest voteshare margin, m(0.1) bw_equal(0) hp(10) hm(8) ///
kernel("uni") se_method("nn") j(3)

// supplied variance test
rdhonest voteshare margin, m(0.1) hp(10) ///
kernel("uni") se_method("supplied_var") j(3)

// test for sclass H, order 2
rdhonest voteshare margin, m(0.1) hp(10) ///
kernel("uni") se_method("nn") j(3) order(2)

// continue with Oreopoulos data
use cghs, replace
gen logEarn = log(earnings)

// uniform kernel
rdhonest logEarn yearat14, m(0.04) c(1947) ///
kernel("uni") opt_criterion("FLCI") sclass("H")

// triangular kernel
rdhonest logEarn yearat14, m(0.04) c(1947) ///
kernel("tri") opt_criterion("FLCI") sclass("H")
