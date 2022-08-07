clear all

adopath + "."

cap log close 
log using "rdhonest.txt",replace 

********************************************************************************
//  1. Lee08
qui: use "./data/lee08.dta",clear 


// 1.1 uni kernel

rdhonestV1 voteshare margin, m(0.1) kernel("uni") h(10) 
rdhonestV1 voteshare margin, m(0.1) kernel("uni") 
rdhonestV1 voteshare margin, kernel("uni") 

// 1.2 tri kernel 

rdhonestV1 voteshare margin, m(0.1) kernel("tri") h(10) 

// 1.2.1 tri kernel + optimal h 

rdhonestV1 voteshare margin, m(0.1) kernel("tri") 
rdhonestV1 voteshare margin, m(0.1) kernel("tri")   opt_criterion("OCI")
rdhonestV1 voteshare margin, m(0.1) kernel("tri")   opt_criterion("FLCI") 

// 1.2.1 tri kernel + optimal h + without M
rdhonestV1 voteshare margin, kernel("tri") 
rdhonestV1 voteshare margin, kernel("tri")   opt_criterion("OCI")
rdhonestV1 voteshare margin, kernel("tri")   opt_criterion("FLCI") 

********************************************************************************
// 2. rcp
qui: use "./data/rcp.dta",clear 

// 2.1 uni kernel

rdhonestV1  cn (retired=elig_year), m(4 0.4) kernel("uni") h(3) 
rdhonestV1  cn (retired=elig_year), m(4 0.4) kernel("uni") t0(0)
rdhonestV1  cn (retired=elig_year), kernel("uni") t0(0)

// 2.2 tri kernel 

rdhonestV1  cn (retired=elig_year), m(4 0.4) kernel("tri") h(3) 

// 2.2.1 tri kernel + optimal h 

rdhonestV1  cn (retired=elig_year), m(4 0.4) kernel("tri") t0(0)
rdhonestV1  cn (retired=elig_year), m(4 0.4) kernel("tri")   opt_criterion("OCI") t0(0)
rdhonestV1  cn (retired=elig_year), m(4 0.4) kernel("tri")   opt_criterion("FLCI")  t0(0)

// 2.2.1 tri kernel + optimal h + without M
rdhonestV1  cn (retired=elig_year), kernel("tri") t0(0)
rdhonestV1  cn (retired=elig_year), kernel("tri")   opt_criterion("OCI") t0(0)
rdhonestV1  cn (retired=elig_year), kernel("tri")   opt_criterion("FLCI") t0(0)
 
log close 
