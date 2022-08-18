// Run from main directory

clear all

global datadir = "data"
global logdir = "tests"

set seed 2022

cap log close 
log using `"${logdir}/rdhonest_test.log"', replace

********************************************************************************
//  1. Lee08
qui: use `"${datadir}/lee08.dta"', clear

// 1.1 uni kernel
rdhonest voteshare margin, m(0.1) kernel("uni") h(10) 

// 1.1.1 uni kernel + optimal h 
rdhonest voteshare margin, m(0.1) kernel("uni") 
rdhonest voteshare margin, m(0.1) kernel("uni") opt_criterion("OCI")
rdhonest voteshare margin, m(0.1) kernel("uni") opt_criterion("FLCI") 

// 1.1.2 uni kernel + optimal h + without M + est_w saved as wgt 
rdhonest voteshare margin, kernel("uni") savew(wgt)

// 1.2 tri kernel 
rdhonest voteshare margin, m(0.1) kernel("tri") h(10) 

// 1.2.1 tri kernel + optimal h 
rdhonest voteshare margin, m(0.1) kernel("tri") 
rdhonest voteshare margin, m(0.1) kernel("tri") opt_criterion("OCI")
rdhonest voteshare margin, m(0.1) kernel("tri") opt_criterion("FLCI") 

cap drop wgt /*for the dofile to run*/
// 1.2.2 tri kernel + optimal h + without M + est_w saved as wgt 
rdhonest voteshare margin, kernel("tri") savew(wgt)

// 1.3 display option tests and by option tests
// 1.3.1 by option (psuedo categories)
qui{
	gen random = runiform()
	gen bygrp = 1 if random<0.5
	replace bygrp = 0 if mi(bygrp)
	cap drop random wgt
}
bys bygrp: rdhonest voteshare margin, m(0.1) kernel("uni") h(10)
bys bygrp: rdhonest voteshare margin, m(0.1) kernel("uni") h(10) savew(wgt)

// 1.3.2 hide parameters
rdhonest voteshare margin, m(0.1) kernel("uni") h(10) noparam

// 1.3.3 show iteration log
rdhonest voteshare margin, m(0.1) kernel("tri") iterl

********************************************************************************
// 2. rcp
qui: use `"${datadir}/rcp.dta"', clear 

// 2.1 uni kernel
rdhonest cn (retired=elig_year), m(4 0.4) kernel("uni") h(3) 

// 2.1.1 uni kernel + optimal h 
rdhonest cn (retired=elig_year), m(4 0.4) kernel("uni") t0(0)
rdhonest cn (retired=elig_year), m(4 0.4) kernel("uni") t0(0) opt_criterion("OCI")
rdhonest cn (retired=elig_year), m(4 0.4) kernel("uni") t0(0) opt_criterion("FLCI") 

// 2.1.2 uni kernel + optimal h + without M + est_w saved as wgt 
rdhonest cn (retired=elig_year), kernel("uni") t0(0) savew(wgt)

// 2.2 tri kernel 
rdhonest cn (retired=elig_year), m(4 0.4) kernel("tri") h(3) 

// 2.2.1 tri kernel + optimal h 
rdhonest cn (retired=elig_year), m(4 0.4) kernel("tri") t0(0)
rdhonest cn (retired=elig_year), m(4 0.4) kernel("tri") opt_criterion("OCI") t0(0)
rdhonest cn (retired=elig_year), m(4 0.4) kernel("tri") opt_criterion("FLCI") t0(0)

cap drop wgt /*for the dofile to run*/
// 2.2.1 tri kernel + optimal h + without M + est_w saved as wgt 
rdhonest cn (retired=elig_year), kernel("tri") t0(0) savew(wgt)

// 2.3 display option tests and by option tests
// 2.3.1 by option (psuedo categories)
qui{
	gen random = runiform()
	gen bygrp = 1 if random<0.5
	replace bygrp = 0 if mi(bygrp)
	cap drop random wgt
}
bys bygrp: rdhonest cn (retired=elig_year), m(4 0.4) kernel("uni") h(3)
bys bygrp: rdhonest cn (retired=elig_year), m(4 0.4) kernel("uni") h(3) savew(wgt)

// 2.3.2 hide parameters
rdhonest cn (retired=elig_year), m(4 0.4) kernel("uni") h(3) noparam

// 2.3.3 show iteration log
rdhonest cn (retired=elig_year), m(4 0.4) kernel("tri") t0(0) iterl

log close
clear
