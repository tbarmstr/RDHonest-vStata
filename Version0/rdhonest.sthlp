{smcl}
{* *! version 0.2.2  KHLee 31jan2019}{...}

{viewerjumpto "Syntax" "examplehelpfile##syntax"}{...}
{viewerjumpto "Description" "examplehelpfile##description"}{...}
{viewerjumpto "Options" "examplehelpfile##options"}{...}
{viewerjumpto "Remarks" "examplehelpfile##remarks"}{...}
{viewerjumpto "Examples" "examplehelpfile##examples"}{...}
{title:Title}

{phang}
{bf:rdhonest} {hline 2} Calculate RDHonest statistic, maximum bias, standard errors and confidence intervals, as in Kolesar and Rothe (2018)

{marker syntax}{...}
{title:Syntax}

{p 8 17 2}
{cmdab:rdhonest}
{depvar}
[{indepvars}]
{ifin}
m(real)
[{cmd:,} {it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{synopt:{opt c(#)}}set regression discontinuity cutoff; default is {cmd:c(0)}{p_end}
{synopt:{opt hp(#)}}set the bandwidth above the cutoff; if set to zero, and an optimality condition is specified, RDHonest will find the optimal bandwidth {p_end}
{synopt:{opt kernel(string)}}set the kernel for the local regression: default is {it:{ul:tri}angular}; other options are {it:{ul:epa}nechnikov} and {it:{ul:uni}form} {p_end}
{synopt:{opt hm(#)}} set the bandwidth below the cutoff; if not specified, will use the value for {cmd:hp(#)}{p_end}
{synopt:{opt opt_criterion(string)}}set the optimality criterion for the optimal bandwidth: default is {it:MSE}; other options are {it:OCI} and {it:FLCI}{p_end}
{synopt:{opt bw_equal(#)}}integer indicating if bandwidths above and below the cutoff are equal; default is {cmd:bw_equal(1)}{p_end}
{synopt:{opt alpha(#)}}sets the worst-case coverage (1-alpha) of the confidence interval; default alpha is 0.05{p_end}
{synopt:{opt beta(#)}}sets the quantile of excess length when optimizing the objective function in the OCI option; default is 0.8{p_end}
{synopt:{opt se_method(string)}}sets method used to compute standard errors; default is {it:nn}, other options are {it:EHW}, {it:demeaned} and {it:supplied_var} {p_end}
{synopt:{opt j(#)}}sets number of nearest neighbors used when {it:nn} option is selected for {cmd:se_method()}; default is 3 {p_end}
{synopt:{opt sclass(string)}}sets the smoothness class used to compute the worst-case bias; default is {it:H}, other option is {it:T} {p_end}
{synopt:{opt order(#)}}sets the order of the local regression; default is 1, other option is 2 {p_end}
{synopt:{opt se_initial(string)}}sets method used to compute preliminary estimates of variances (for optimal bandwidth calculations); default is {it:IKEHW}, other options are {it:Silverman}, {it:SilvermanNN}, {it:IKdemeaned} and {it:NN} {p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
{it:m(#)} is the bound for the Taylor/Holder approximation error; see the RDHonest PDF manual for details.{p_end}
{p 4 6 2}
{*: {cmd:fweight}s are NOT allowed; DO NOT see {help weight}.}

{marker description}{...}
{title:Description}

{pstd}
{cmd:rdhonest} calculates estimators and one- and two-sided CIs based on local polynomial estimators in RD under the second-order Taylor or Holder smoothness class.

{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt c(#)} sets the regression discontinuity cutoff. The default value for the cutoff is zero.

{phang}
{opt hp(#)} sets the bandwidth above the cutoff. If {opt hm(#)} is not specified, then {opt hp(#)} sets the bandwidth on both sides of the cutoff.
If neither {opt hp(#)} nor {opt hm(#)} is specified, then an optimality condition needs to be provided; in this case, rdhonest will find the optimal bandwidth.

{phang}
{opt kernel(string)} sets the kernel for the local regression (in the context of standard sharp regression discontinuity).
The default option is the triangular ({it:tri}) kernel, k(u) = (1-|u|)_{+}.
Other options include the epanechnikov ({it:epa}, 0.75(1-u^2)_{+}) and uniform ({it:uni}, (|u|<1)/2) kernels.

{phang}
{opt hm(#)} sets the bandwidth below the cutoff. If {opt hm(#)} is not specified, rdhonest will use the bandwidth in {opt hp(#)}.

{phang}
{opt opt_criterion(string)} sets the optimality criterion used to compute the optimal bandwidth.
The default is the finite-sample maximum mean-squared error ({it:MSE}).
Other options are the given quantile of the excess length of one-sided confidence intervals ({it:OCI}) and the length of two-sided fixed-length confidence intervals ({it:FLCI}).

{phang}
{opt bw_equal(#)} specifies if the bandwidths above and below the cutoff are equal, when optimal bandwidths have to be calculated. If bw_equal is 1, then the computed bandwidths are equal; else if bw_equal is 0, then the bandwidths are not.

{phang}
{opt alpha(#)} sets the worst-case coverage of the confidence interval. rdhonest computes confidence intervals that guarantee that the true parameter is covered with probability (1-alpha). The default alpha is 0.05.

{phang}
{opt beta(#)} sets the quantile of excess length.
This is only used when the optimal bandwidth is chosen using the one-sided confidence interval criterion ({it:OCI}).
The default value of beta is 0.8.

{phang}
{opt se_method(#)} sets the method used to compute standard errors when running local regressions above and below the cutoff.
The default method is that of k-nearest neighbors ({it:nn}).
Others are Eicker-Huber-White standard errors ({it:EHW}), the demeaned analog of EHW standard errors obtained by subtracting the estimated slope ({it:demeaned}), and if the standard errors are supplied by the user ({it:supplied_var}).

{phang}
{opt j(#)} sets the number of nearest neighbors used to compute standard errors when {opt se_method("nn")} is specified.
The default j is 3.

{phang}
{opt sclass(string)} sets the smoothness class used to compute the worst-case bias.
The default is the Holder class ({it:H}), which bounds the second derivative of the estimated function globally.
The other option is the Taylor class ({it:T}), which bounds the second derivative of the estimated function at the cutoff by the value {it:m(#)} specified above.

{phang}
{opt order(#)} sets the order of the local regression.
The default is 1 (local linear regression), and the other option is 2 (local quadratic regression).

{phang}
{opt se_initial(string)} sets the method used to compute preliminary estimates of variances (for optimal bandwidth calculations).
The default option is {it:IKEHW}, based on residuals from a local linear regression using a triangular kernel and Imbens and Kalyanaraman bandwidth.
Other options include:
{it:IKdemeaned}, Based on sum of squared deviations of outcome from estimate of intercept in local linear regression with triangular kernel and Imbens and Kalyanaraman bandwidth;
{it:Silverman}, using residuals from a local constant regression with uniform kernel and bandwidth selected using Silverman's rule of thumb, as in Equation (14) in Imbens and Kalyanaraman (2012);
{it:SilvermanNN}, using nearest neighbor estimates, rather than residuals;
{it:NN}, using nearest neighbor estimates, without assuming homoscedasticity.

{marker remarks}{...}
{title:Remarks}

{pstd}
For detailed information on the rdhonest package, see the associated PDF manual, {bf: Honest Inference in Sharp Regression Discontinuity in STATA} (RDHonest.pdf).


{marker examples}{...}
{title:Examples}

{phang}
{cmd:. rdhonest voteshare margin, m(0.1) hp(10) kernel("uni") se_method("nn") sclass("H") j(3)}{p_end}

{phang}
{cmd:. rdhonest voteshare margin, m(0.1) kernel("tri") opt_criterion("MSE") sclass("H")}{p_end}

{phang}
{cmd:. rdhonest logEarn yearat14, m(0.04) hp(2) c(1947) kernel("uni") sclass("H")}{p_end}


{marker references}{...}
{title:References}

{marker AK2018a}{...}
{phang}
Armstrong, Timothy B., and Michal Kolesar. 2018. Optimal inference in a class of regression models. {it:Econometrica}, 86, 655-683. 

{marker AK2018b}{...}
{phang}
Armstrong, Timothy B., and Michal Kolesar. 2018. Simple and Honest Confidence Intervals in Nonparametric Regression.

{marker KR2018}{...}
{phang}
Kolesar, Michal, and Christoph Rothe. 2018. Inference in regression discontinuity designs with a discrete running variable, {it:American Economic Review}, 108, 2277–2304.
{p_end}

{title:End}